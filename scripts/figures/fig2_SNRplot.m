%% input files
sn='/ptmp/pvalsala/deuterium/20240813_spectral';
pn='/ptmp/pvalsala/deuterium/20240813_spectral/proc/SNR';
dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));

dirst_me=dir(fullfile(sn,"*pvrh_trufi_5E_*.dat"));


addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))


metabolites=getMetaboliteStruct('phantom',-4);

%% process noise anf get fieldmap
fn_noise=dir(fullfile(sn,"*Noise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

fn_fm=(fullfile(sn,dirst_me(1).name));
ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true};
mcobj_ideal=MetCon_bSSFP(fn_fm,'Solver','IDEAL',ME_setting{:},'mask',80);
fm_meas_Hz=mcobj_ideal.Experimental.fm_est*(-2*pi)/(6.536 /42.567);
%% final setttings
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','AMARES'};
ME_setting=[ME_setting,{'fm',fm_meas_Hz,'Solver','pinv'}];


%% process all CSI
mcobj_csi=cell(length(dirst_csi),1);
for cf=1:length(dirst_csi)
       fn=fullfile(sn,dirst_csi(cf).name);
    mcobj_csi{cf}=MetCon_CSI(fn,CSI_setting{:});
end
%% process all ME data
mcobj_me=cell(length(dirst_me),1);
for cf=1:length(dirst_me)
       fn=fullfile(sn,dirst_me(cf).name);
    mcobj_me{cf}=MetCon_bSSFP(fn,ME_setting{:});
end
%  mcobj_me{1}.PlotResults
%% reslicing and registration
%  cellfun(@(x) x.WriteImages,mcobj_csi,'UniformOutput',false)
%  cellfun(@(x) x.WriteImages,mcobj_me,'UniformOutput',false)
dirst_nii=dir(fullfile(pwd,"Metcon*.nii"));
resliced_vol=myspm_reslice(dirst_nii(1).name,dirst_nii,'linear','r');
[realigned_vol]=realign_vol(resliced_vol);

im_sf=3;
realigned_vol3=NDimscale(realigned_vol,im_sf);
vox_vol=[cellfun(@(x)prod(x.DMIPara.resolution_PSF(1:3)*1e2),mcobj_me(1:2),'UniformOutput',true);...
cellfun(@(x)prod(x.DMIPara.resolution_PSF(1:3)*1e2),mcobj_csi(1:5),'UniformOutput',true)];

AllDescription=cellfun(@(obj) getDescription(obj),[mcobj_me(1:2);mcobj_csi(1:5)],'UniformOutput',false);
%mask 
% imMask=CreateMask(squeeze(realigned_vol3(:,slcSel,:,1,1)));

% save('processced_data.mat','dirst*','CSI_setting','ME_setting','mcobj_me','mcobj_csi','slcSel','tubes','realigned_vol3','vox_vol','imMask','AllDescription','-v7.3')

%% get ROIs
getRotm= @(theta)[cosd(theta) sind(theta); -sind(theta) cosd(theta)];
tubes={[49,63],[105,33]}; %diagonally opposite 9 thtube and 
center=0.5*(tubes{1}+tubes{2});
figure,
imagesc(squeeze(realigned_vol3(:,slcSel,:,1,1))),hold on, axis image
for i=1:9
mask=zeros(size(realigned_vol3,1),size(realigned_vol3,3),'logical');

ROIs{i}=round(getRotm(36*(i-1))*(tubes{1}-center)'+center');
mask(ROIs{i}(1),ROIs{i}(2))=true;
mask=imdilate(mask,strel('sphere',round(1.3*im_sf)));
contour(mask,1,'color','red','linewidth',2)
end

%% plotting
figure(34),clf
tt=tiledlayout(4,2);



slcSel=32*im_sf;
stat_mean=zeros(size(realigned_vol,4),size(realigned_vol,5));
for ii=1:size(realigned_vol,4)
    nexttile()
    realigned_vol4=realigned_vol3./reshape(vox_vol./vox_vol(1),1,1,1,1,[]);
im_plot=reshape(realigned_vol4(:,slcSel,:,ii,:),size(realigned_vol4,1),[]);

%black sorounding
mask=reshape(repmat(imMask(:,:),[1 1 1 size(realigned_vol4,5) ]),size(realigned_vol4,1),[]);
im_plot(~mask)=nan;

imagesc(im_plot,'AlphaData',mask),colormap("jet"),colorbar,axis image,hold on
title(metabolites(ii).name)
set(gca,'color','black')

% plot ROI
picked=[9,5,7,2]; % selected tubes

cMask=zeros(size(realigned_vol4,1),size(realigned_vol4,3),'logical');
cMask(ROIs{picked(ii)}(1),ROIs{picked(ii)}(2))=true;

% picked2=[8,6,8,1];
% cMask(ROIs{picked2(ii)}(1),ROIs{picked2(ii)}(2))=true;
% picked2=[7,4,6,1];
% cMask(ROIs{picked2(ii)}(1),ROIs{picked2(ii)}(2))=true;


cMask=imdilate(cMask,strel('sphere',round(2*im_sf)));
cMet=reshape(realigned_vol4(:,slcSel,:,ii,:),numel(cMask),[]);
cMet=cMet(cMask,:);
stat_mean(ii,:)=mean(cMet,1);

cMask=repmat(cMask,[1 size(realigned_vol3,5)]);
  contour(cMask,1,'color','red','linewidth',0.1)

yticks([]),xticks(round((0.5:1:(0.5+size(realigned_vol3,5)))*size(realigned_vol3,3)))
xticklabels([1:7])

end
nexttile(5,[2 2])
bar(stat_mean')
legend({metabolites.name})
grid on
xticklabels(AllDescription)



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
    imout(:,:,:,ii)=imresize3(im_nd(:,:,:,ii),scale);
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
    seqType=[seqType,'GRE'];
    isbSSFP=false;
  end

DesStr=sprintf('%s|%s|%s|TR %.0fms',seqType,mcobj.flags.doCoilCombine,mcobj.flags.Solver,mcobj.DMIPara.TR*1e3);

if(isbSSFP)
    DesStr=[DesStr,sprintf('|%d PC',length(mcobj.DMIPara.PC_deg))];
end
end
