%% run both fig6_CSI_bSSFP script and Fig8_ME_bSSFP script

% read all mM vols replace water with SNR vol.
% read anat with skull for showing
% read anat_brain for masking
% vioin plots are still in native resolution


%set up paths and data specifics

addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))

MeasPath='/ptmp/pvalsala/deuterium/dataForPublication/sub-02-DMI';
pn=fullfile(MeasPath,sprintf('proc_all_%s',datetime('today','Format','yyyyMMMdd')));
mkdir(pn);cd(pn);
dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));
dirst_me=dir(fullfile(sn,"*rh_trufi_5E_18PC*.dat"));

metabolites=getMetaboliteStruct('invivo');
type_label_s2=[1     3     2     3     1     2     3];
intake_time_all=[43,55,67,87,122,132,143]; % minutes
data_label={'CSI-FISP','CSI-bSSFP','ME-bSSFP'};

%% process noise and init settings
fn_noise=dir(fullfile(sn,"*rh_trufi*noise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true};
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL-modes','fm','IDEAL'};
% fm_meas_Hz=mcobj_ideal.Experimental.fm_est*(-2*pi)/(6.536 /42.567);
%% process all data and export nifti
mcobj_me=cell(length(dirst_me),1);
for cf=1:length(dirst_me)

       fn=fullfile(sn,dirst_me(cf).name);   
    mcobj_me{cf}=MetCon_ME(fn,ME_setting{:},'fm','IDEAL','Solver','IDEAL-modes');
end
mcobj_me{3}.ShiftMetcon([0,-10,0])% shift ~10 mm in read for up in IS direction: motion corr
% CSI data

mcobj_csi=cell(length(dirst_csi),1);
for cf=1:length(dirst_csi)
       fn=fullfile(sn,dirst_csi(cf).name);
    mcobj_csi{cf}=MetCon_CSI(fn,CSI_setting{:});
end
%write images

[intake_time_s2,order_idx_s2]=sort(cellfun(@(x) x.getMinutesAfterIntake('08:36'),[mcobj_csi;mcobj_me],'UniformOutput',true));


%% reslice
cd(pn)
cellfun(@(x) x.WriteImages('',{'snr','mM'}),[mcobj_csi;mcobj_me],'UniformOutput',false);
reslice_all_sag=myspm_reslice('../anat/*anat*_sag*.nii',dir('Metcon_*.nii'), 'nearest','rs');
reslice_all_tra=myspm_reslice('../anat/*anat*_tra*.nii',dir('Metcon_*.nii'), 'nearest','rt');

%% read mifti and try imcompose
cd(pn)
anat_nii_sag=double(MyNiftiRead("../anat/*anat_sag.nii",'IPR'));
brain_nii_sag=double(MyNiftiRead("../anat/*brain_sag.nii",'IPR'));
brain_nii_sag=brain_nii_sag./max(brain_nii_sag);

[met_snr_me_sag]=MyNiftiRead('rsMetcon_SNR*.nii','IPR');
[met_mm_sag]=MyNiftiRead('rsMetcon_mM*.nii','IPR');
 met_mm_sag(:,:,:,1,:)=met_snr_me_sag(:,:,:,1,:);

anat_nii_tra=double(MyNiftiRead("../anat/*anat_tra.nii",'PRS'));
brain_nii_tra=double(MyNiftiRead("../anat/*brain_tra.nii",'PRS'));
brain_nii_tra=brain_nii_tra./max(brain_nii_tra);

[met_snr_me_tra]=MyNiftiRead('rtMetcon_SNR*.nii','PRS');
[met_mm_tra]=MyNiftiRead('rtMetcon_mM*.nii','PRS');
met_mm_tra(:,:,:,1,:)=met_snr_me_tra(:,:,:,1,:);
%%
slct=14;
slcs=18;
% crop_val  sag             tra : [idx1:end-idx2,idx3:end-idx4] 
crop_val={[20 40 45 25],[20 40 45 50]};

data_label={'CSI','CSI-PC-bSSFP','ME-PC-bSSFP'};


[im_anat]=ComposeImage({anat_nii_sag,anat_nii_tra},{slcs,slct},crop_val);
[im_brain]=ComposeImage({brain_nii_sag,brain_nii_tra},{slcs,slct},crop_val);

midlabel= @(Sz,N) Sz/2:Sz:Sz*N*1.5;
title_str={'Water [SNR]','Glc [mM]','Glx [mM]'};
%  as(im_anat)

%
[imComposed_csi]=ComposeImage({met_mm_sag(:,:,:,:,(type_label_s2==1)),met_mm_tra(:,:,:,:,(type_label_s2==1))}, ...
    {slcs,slct},crop_val,[0 5 0]);
% as(imComposed2)

%
SNR_scl=ones(size(met_mm_tra,5),1);
SNR_scl(type_label_s2==2)=4.0425/4.2831;
SNR_scl(type_label_s2==3)=4.0425/2.6670;


load('/ptmp/pvalsala/deuterium/H4DK-64GH/proc/combined/masks.mat')
 % mask_me=(im_brain./max(im_brain(:)))>0.05;
 mask_me=m_tra|m_sag;
% mask_me=imComposed_csi(:,:,1,1,1)>8;
 mask_me=imerode(mask_me,strel('sphere',10));
 mask_me=imdilate(mask_me,strel('sphere',8));

 common_settting={'prctile',99,'SlcSel',1,'transform',@(x) x,'cax_im',[0 1],'Mask',mask_me};
cax_met={[0 70],[0 2.2],[0 2.2]};
figure(45),clf
    set(gcf,"Position",[50 50 1080 1080])
tt=tiledlayout(9,3,'TileSpacing','compact','Padding','compact');
sgtitle('CSI vs CSI-bSSFP','fontsize',24)
cb_all={};cax_all={};cnt=1;
for i=1:3
    nexttile(tt,[2 1]);
[cb_all{cnt},cax_all{cnt}]=overlayplot(im_anat,imComposed_csi,'MetIdx',i,'cax',cax_met{i},common_settting{:});
if(i==1),yticks(midlabel(size(im_anat,1),3)),yticklabels(intake_time_all(type_label_s2==1)),else,yticks([]),end
title(title_str{i})
 cb_all{cnt}=colorbar;
cnt=cnt+1;

end



% cax_met{1}=[0 40]
[imComposed_ssfp]=ComposeImage({met_mm_sag(:,:,:,:,(type_label_s2==2)),met_mm_tra(:,:,:,:,(type_label_s2==2))}, ...
    {slcs,slct},crop_val,[0,0,0]);
imComposed_ssfp(:,:,:,1,:)=imComposed_ssfp(:,:,:,1,:)*SNR_scl(3);
for i=1:3
    nexttile(tt,[2 1]);
[cb_all{cnt},cax_all{cnt}]=overlayplot(im_anat,imComposed_ssfp,'MetIdx',i,'cax',cax_met{i},common_settting{:});
if(i==1),yticks(midlabel(size(im_anat,1),4)),yticklabels(intake_time_all(type_label_s2==2)),else,yticks([]),end
title(title_str{i})
cb_all{cnt}=colorbar;
cnt=cnt+1;

end
cax_met{1}=[0 50];
[imComposed_me]=ComposeImage({met_mm_sag(:,:,:,:,(type_label_s2==3)),met_mm_tra(:,:,:,:,(type_label_s2==3))}, ...
    {slcs,slct},crop_val,[0,0,0]);
imComposed_me(:,:,:,1,:)=imComposed_me(:,:,:,1,:)*SNR_scl(2);
for i=1:3
    nexttile(tt,[3 1]);
[cb_all{cnt},cax_all{cnt}]=overlayplot(im_anat,imComposed_me,'MetIdx',i,'cax',cax_met{i},common_settting{:});
if(i==1),yticks(midlabel(size(im_anat,1),4)),yticklabels(intake_time_all(type_label_s2==3)),else,yticks([]),end
title(title_str{i})
cb_all{cnt}=colorbar;
cnt=cnt+1;

end



for i=1:length(cb_all)
makeColorbar(cb_all{i},cax_all{i});
end

annotation(gcf,'textbox',...
    [0.03 0.5 0.222177412973777 0.0359869131020137],...
    'String',{'Time after intake [min]'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',14,'EdgeColor','none','FitBoxToText','on');

set(gcf,'Color','w','InvertHardcopy','off')


% stats from registered maps

mask=brain_nii_tra(:,:,:,1,1)>1e-3;
 mask=imerode(mask,strel('disk',10,0));
mask=imdilate(mask,strel('disk',10,0));
% as(mask.*anat_nii_tra)

col_tra=reshape(met_mm_tra,[],4,size(met_mm_tra,5));

colors_label=lines(3);

mask=mask(:) & all(col_tra(:,2,:)<3,3);
mask=mask(:) & all(col_tra(:,3,:)<3,3);

% figure,
% tt=tiledlayout(1,3);
violin_format={'mc',[],'medc',[]};
nexttile(tt,[2 1])

violin(squeeze(col_tra(mask(:),1,:)).*SNR_scl',violin_format{:},'facecolor',colors_label(type_label_s2,:));
grid on,grid minor, title(title_str{1}),xticks(1:7),xticklabels(intake_time_all),set(gca,'FontWeight','bold')

mean_val=mean(squeeze(col_tra(mask(:),1,:)).*SNR_scl');
hold on,plot(find(type_label_s2==1),mean_val(type_label_s2==1),'Marker','x','Color',colors_label(1,:),'LineWidth',1.5)
hold on,plot(find(type_label_s2==2),mean_val(type_label_s2==2),'Marker','x','Color',colors_label(2,:),'LineWidth',1.5)
hold on,plot(find(type_label_s2==3),mean_val(type_label_s2==3),'Marker','x','Color',colors_label(3,:),'LineWidth',1.5)


nexttile(tt,[2 1])
violin(squeeze(col_tra(mask(:),2,:)),violin_format{:},'facecolor',colors_label(type_label_s2,:));
grid on,grid minor, title(title_str{2}),xticks(1:7),xticklabels(intake_time_all),set(gca,'FontWeight','bold')

mean_val=mean(squeeze(col_tra(mask(:),2,:)));
hold on,plot(find(type_label_s2==1),mean_val(type_label_s2==1),'Marker','x','Color',colors_label(1,:),'LineWidth',1.5)
hold on,plot(find(type_label_s2==2),mean_val(type_label_s2==2),'Marker','x','Color',colors_label(2,:),'LineWidth',1.5)
hold on,plot(find(type_label_s2==3),mean_val(type_label_s2==3),'Marker','x','Color',colors_label(3,:),'LineWidth',1.5)

xlabel('Time after intake [min]')


nexttile(tt,[2 1])
vh=violin(squeeze(col_tra(mask(:),3,:)),violin_format{:},'facecolor',colors_label(type_label_s2,:));
grid on,grid minor, title(title_str{3}),xticks(1:7),xticklabels(intake_time_all),set(gca,'FontWeight','bold')

mean_val=mean(squeeze(col_tra(mask(:),3,:)));
hold on,plot(find(type_label_s2==1),mean_val(type_label_s2==1),'Marker','x','Color',colors_label(1,:),'LineWidth',1.5)
hold on,plot(find(type_label_s2==2),mean_val(type_label_s2==2),'Marker','x','Color',colors_label(2,:),'LineWidth',1.5)
hold on,plot(find(type_label_s2==3),mean_val(type_label_s2==3),'Marker','x','Color',colors_label(3,:),'LineWidth',1.5)



legend(gca,vh([1 3 2]),data_label([1 2 3]),'Position',[0.620349794317177 0.00128943966208529 0.370370363554469 0.020370369873665],...
    'NumColumns',3)



%%

function makeColorbar(cb_handle,cax)



cb_handle.Visible='off';
ax2 = axes('Position',cb_handle.Position);

imagesc(linspace(0,1,100)'),colormap(ax2,'turbo')
set(ax2,'YAxisLocation','right','FontSize',10,'FontWeight','bold','YDir','normal')
yticks(linspace(0,100,4))
yticklabels(round(linspace(0,1,4)*cax(2),1))
xticks([])
end



function [imComposed]=ComposeImage(Vol,slcSel,crop_val,shift)


if(~exist('shift','var'))
    shift=[0,0,0];
end
for i=1:length(Vol)
    Vol{i}=ndCircShift(Vol{i},shift,1:3);
    im{i}=Vol{i}(crop_val{i}(1):end-crop_val{i}(2),crop_val{i}(3):end-crop_val{i}(4),slcSel{i},:,:);
end

imsz1=size(im{1});
imsz2=size(im{2});

diff=abs(imsz1(2)-imsz2(2));

if(imsz1(1)>imsz2(1))   
    im{2}=padarray(im{2},[0 floor(diff/2)],0,'pre');
    im{2}=padarray(im{2},[0 ceil(diff/2)],0,'post');
else
    im{1}=padarray(im{1},[0 floor(diff/2)],0,'pre');
    im{1}=padarray(im{1},[0 ceil(diff/2)],0,'post');


end
imComposed=cat(2,im{:});
end
