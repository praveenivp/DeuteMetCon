%% input files
MeasPath='/ptmp/pvalsala/deuterium/20241021_phantomtest';
sn=fullfile(MeasPath,'TWIX');

dirst_csi=dir(fullfile(sn,"*rpcsi_fid*.dat"));
dirst_csi=dirst_csi(2);
dirst_csi_ssfp=dir(fullfile(sn,"*rpcsi_ssfp*.dat"));
% dirst_csi_ssfp(2:3);
dirst_me=dir(fullfile(sn,"*pvrh_trufi_5E_*.dat"));
dirst_me=dirst_me(2);

pn=fullfile(MeasPath,sprintf('proc/csi_GRE_%s',datetime('today','Format','yyyyMMMdd')));


addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))


metabolites=getMetaboliteStruct('phantom',0);

%% process noise anf get fieldmap
fn_noise=dir(fullfile(sn,"*oise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true,'fm','IDEAL','Solver','IDEAL-modes'};

CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL','fm',[]};

CSI_setting_ssfp={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL-modes','fm','IDEAL'};


%% process all CSI
mcobj_csi=cell(length(dirst_csi),1);
for cf=1:length(dirst_csi)
       fn=fullfile(sn,dirst_csi(cf).name);
    mcobj_csi{cf}=MetCon_CSI(fn,CSI_setting{:});
end
mcobj_csi_ssfp=cell(length(dirst_csi_ssfp),1);
for cf=1:length(dirst_csi_ssfp)
       fn=fullfile(sn,dirst_csi_ssfp(cf).name);
    mcobj_csi_ssfp{cf}=MetCon_CSI(fn,CSI_setting_ssfp{:});
end
% process all ME data
mcobj_me=cell(length(dirst_me),1);
for cf=1:length(dirst_me)
       fn=fullfile(sn,dirst_me(cf).name);
    mcobj_me{cf}=MetCon_ME(fn,ME_setting{:});
end
%  mcobj_me{1}.PlotResults
%% reslicing and registration
mkdir(pn)
cd(pn)
mcobj_all=[mcobj_csi,mcobj_csi_ssfp,mcobj_me];
cellfun(@(x) x.WriteImages(),mcobj_all,'UniformOutput',false)
dirst_nii=dir(fullfile(pwd,"Metcon_SNR*.nii"));
resliced_vol=myspm_reslice(dirst_nii(1).name,dirst_nii,'nearest','r');
[realigned_vol]=realign_vol(resliced_vol);
data_label={'CSI','CSI-PC-SSFP','ME-PC-bSSFP'};


im_sf=3;
realigned_vol3=NDimscale(realigned_vol,im_sf);
vox_vol=cellfun(@(x)prod(x.DMIPara.resolution_PSF(1:3)*1e2),mcobj_all,'UniformOutput',true);
%%
AllDescription=cellfun(@(obj) getDescription(obj),mcobj_all,'UniformOutput',false);
%mask 
slcSel=28*im_sf;
  % 
% figure,imMask=CreateMask(squeeze(realigned_vol3(:,:,slcSel,1,1)));

%  save('processced_data.mat','dirst*','CSI_setting','ME_setting','mcobj_me','mcobj_csi','slcSel','tubes','realigned_vol3','vox_vol','imMask','AllDescription','-v7.3')

%% get ROIs and plot
slcSel=28*im_sf;%+(-20:-15); % more slices for better stats
getRotm= @(theta)[cosd(theta) sind(theta); -sind(theta) cosd(theta)];
tubes={[46,115],[110,42]}; %diagonally opposite 9th tube and 
center=0.5*(tubes{1}+tubes{2});
%figure,
% imagesc(squeeze(realigned_vol3(:,:,slcSel,1,1))),hold on, axis image
for i=1:9
mask=zeros(size(realigned_vol3,1),size(realigned_vol3,3),'logical');

ROIs{i}=round(getRotm(36*(i-1))*(tubes{1}-center)'+center');
mask(ROIs{i}(1),ROIs{i}(2))=true;
mask=imdilate(mask,strel('sphere',round(3*im_sf)));
%  contour(mask,1,'color','red','linewidth',2)
end

% plotting
figure(41),clf
tt=tiledlayout(4,2,'TileSpacing','compact','Padding','compact');

picked=[9,5,7,2]; % selected tubes
ROIs{9}=ROIs{9}-[0 ;2];
ROIs{7}=ROIs{7}-[-1 ;2];
ROIs{2}=ROIs{2}-[1 ;1];

PlotSel=1:length(mcobj_all);    
stat_mean=zeros(size(realigned_vol,4),size(realigned_vol,5));
for ii=1:size(realigned_vol,4)
    nexttile()
    realigned_vol4=realigned_vol3./reshape(vox_vol,1,1,1,1,[]);
im_plot=reshape(mean(realigned_vol4(:,:,slcSel,ii,:),3),size(realigned_vol4,1),[]);

%black sorounding
mask=reshape(repmat(imMask(:,:),[1 1 1 size(realigned_vol4,5) ]),size(realigned_vol4,1),[]);
im_plot(~mask)=nan;

imagesc(im_plot,'AlphaData',mask),colormap("turbo"),colorbar,axis image,hold on
title(metabolites(ii).name)
set(gca,'color','black')

% plot ROI
cMask=zeros(size(realigned_vol4,1),size(realigned_vol4,3),'logical');
cMask(ROIs{picked(ii)}(1),ROIs{picked(ii)}(2))=true;

% picked2=[8,6,8,1];
% cMask(ROIs{picked2(ii)}(1),ROIs{picked2(ii)}(2))=true;
% picked2=[7,4,6,1];
% cMask(ROIs{picked2(ii)}(1),ROIs{picked2(ii)}(2))=true;


cMask=imdilate(cMask,strel('disk',round(3.5*im_sf)),1);
cMet=reshape(mean(realigned_vol4(:,:,slcSel,ii,:),3),numel(cMask),[]);
cMet=cMet(cMask,:);
% stat_mean(ii,:)=mean(cMet,1);
  stat_mean(ii,:)=median(cMet,1);
 stat_std(ii,:)=std(cMet,[],1);

cMask=repmat(cMask,[1 size(realigned_vol4,5)]);
  contour(cMask,1,'EdgeColor',[0,0,0],'EdgeAlpha',1,'linewidth',1.5)

yticks([]),xticks(round((0.5:1:(0.5+size(realigned_vol3,5)))*size(realigned_vol3,3)))
 xticklabels(data_label)
end

% colors_met=[[0,0,255];[0, 128,0]; [255,212,42];[255,102,0]]; 
% colors_datalabels=([[127,184,221];[235,168,139];[245,215,143];]-10)./250;
nexttile(5,[2 2])
hb_handle=bar(stat_mean,'LineWidth',0.75);
% set(gca,"ColorOrder",colors_datalabels);
arrayfun(@(x)set(x,'FaceAlpha',0.6),hb_handle);
%
hold on;
for cMethod = 1:size(stat_mean,2)
    % get x positions per group
    xpos = hb_handle(cMethod).XEndPoints% + hb_handle(cMethod).XOffset;
    % draw errorbar
    eb_h=errorbar(xpos, stat_mean(:,cMethod), stat_std(:,cMethod), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1.5);
end
ylabel('SNR'),grid on,grid minor
 set(gca,'ylim',get(gca,'ylim')+[0 10])
% xticklabels(AllDescription(PlotSel))
 legend(data_label),xticklabels({metabolites.name})
fontsize(gcf,'scale',1.9);

    set(gcf,'Color','w','Position',[397 64 1154 943],'InvertHardcopy','off')
disp(stat_mean./(stat_mean(:,1))) %disp mean ratios

%% plot signal levels
p1=squeeze(mcobj_me {1}.img(1,24,32,27,1,:))';
p2=squeeze(mcobj_csi_ssfp{1}.img(1,25,10,26,1,:))';

vol_fac=(vox_vol(3)/vox_vol(2));
dc_fac=sqrt(mcobj_me{1}.DMIPara.DutyCycle./mcobj_csi_ssfp{1}.DMIPara.DutyCycle);
acq_fac=sqrt((10/(5*18))/(10/(64*4)));

me_amp_fac= (1/vol_fac).*dc_fac.*acq_fac;

figure(9),clf,hold on
plot(mcobj_me{1}.DMIPara.PC_deg,abs(p1(:)),'-*')
plot(mcobj_csi_ssfp{1}.DMIPara.PC_deg,me_amp_fac.*abs(p2(:)),'-*')
legend('ME','CSI')
%%


function [realigned_vol]=realign_vol(vols)
[optimizer,metric] = imregconfig("multimodal");
realigned_vol=zeros(size(vols));
for i=1:size(vols,5)
    %calcualte translation from water
[tform]=imregtform(vols(:,:,:,1,i),vols(:,:,:,1,1),"translation",optimizer,metric);
%apply for all metabolites
for j=1:size(vols,4)
realigned_vol(:,:,:,j,i)=imwarp(vols(:,:,:,j,i),tform,"OutputView",imref3d(size(vols(:,:,:,1,1))));
end

end

end

function imout=NDimscale(im_nd,scale)
sz=size(im_nd);
im_nd=reshape(im_nd,sz(1),sz(2),sz(3),[]);
imout=zeros([sz(1:3)*scale,sz(4:end)]);
for ii=1:size(im_nd,4)
    imout(:,:,:,ii)=imresize3(im_nd(:,:,:,ii),scale,'nearest');
end
imout=reshape(imout,[size(imout,1:3),sz(4:end)] );
end

function DesStr=getDescription(mcobj)

  if(mcobj.DMIPara.isCSI)
    seqType='CSI-';
  else
    seqType='ME-';
  end
  if(mcobj.DMIPara.TR*1e3<20)
    seqType=[seqType,'bSSFP'];
    isbSSFP=true;
  else
    seqType=[seqType,'FISP'];
    isbSSFP=false;
  end

DesStr=sprintf('%s|TR %.0fms',seqType,mcobj.DMIPara.TR*1e3);

if(isbSSFP)
    DesStr=[DesStr,sprintf('|%d PC',length(mcobj.DMIPara.PC_deg))];
end
end
