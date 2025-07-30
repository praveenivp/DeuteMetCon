%% run both fig6_CSI_bSSFP script and Fig8_ME_bSSFP script

% read all mM vols replace water with SNR vol.
% read anat with skull for showing
% read anat_brain for masking
% vioin plots are still in native resolution


% type_label_s2=[1     3     2     3     1     2     3];
% intake_time_s2=[43,55,67,87,122,132,143];



addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))

metabolites=getMetaboliteStruct('invivo');
type_label_s2=[1     3     2     3     1     2     3];
intake_time_s2=[43,55,67,87,122,132,143]; % minutes
data_label={'CSI-FISP','CSI-bSSFP','ME-bSSFP'};

%% process sub-01
MeasPath='/ptmp/pvalsala/deuterium/dataForPublication/sub-01-DMI';
pn_s1=fullfile(MeasPath,sprintf('proc_all_%s',datetime('today','Format','yyyyMMMdd')));
mkdir(pn_s1);cd(pn_s1);
dirst_csi=dir(fullfile(MeasPath,"*rpcsi*.dat"));
dirst_me=dir(fullfile(MeasPath,"*rh_trufi_5E_18PC*.dat"));

% process noise and init settings
fn_noise=dir(fullfile(MeasPath,"*rh_trufi*noise*.dat"));
twix_noise=mapVBVD(fullfile(MeasPath,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true};
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL-modes','fm','IDEAL'};
% fm_meas_Hz=mcobj_ideal.Experimental.fm_est*(-2*pi)/(6.536 /42.567);
%% process all data and export nifti
mcobj_me_s1=cell(length(dirst_me),1);
for cf=1:length(dirst_me)
       fn=fullfile(MeasPath,dirst_me(cf).name);   
    mcobj_me_s1{cf}=MetCon_ME(fn,ME_setting{:},'fm','IDEAL','Solver','IDEAL-modes');
end
% mcobj_me_s2{3}.ShiftMetcon([0,-10,0])% shift ~10 mm in read for up in IS direction: motion corr
% CSI data

mcobj_csi_s1=cell(length(dirst_csi),1);
for cf=1:length(dirst_csi)
       fn=fullfile(MeasPath,dirst_csi(cf).name);
    mcobj_csi_s1{cf}=MetCon_CSI(fn,CSI_setting{:});
end

%write images
[intake_time_s1,order_idx_s1]=sort(cellfun(@(x) x.getMinutesAfterIntake('08:13'),[mcobj_csi_s1;mcobj_me_s1],'UniformOutput',true));

% reslice
cellfun(@(x) x.WriteImages('',{'snr','mM'}),[mcobj_csi_s1;mcobj_me_s1],'UniformOutput',false);
reslice_all_sag=myspm_reslice('../anat/*anat_sag*.nii',dir('Metcon_*.nii'), 'nearest','rs');
reslice_all_tra=myspm_reslice('../anat/*anat_tra*.nii',dir('Metcon_*.nii'), 'nearest','rt');
type_label_s1 =[1     2     3     1     3     2     3];
%% process sub-03
MeasPath='/ptmp/pvalsala/deuterium/dataForPublication/sub-03-DMI';
pn_s3=fullfile(MeasPath,sprintf('proc_all_%s',datetime('today','Format','yyyyMMMdd')));
mkdir(pn_s3);cd(pn_s3);
dirst_csi=dir(fullfile(MeasPath,"*rpcsi*.dat"));
dirst_me=dir(fullfile(MeasPath,"*rh_trufi_5E_18PC*.dat"));

% process noise and init settings
fn_noise=dir(fullfile(MeasPath,"*rh_trufi*noise*.dat"));
twix_noise=mapVBVD(fullfile(MeasPath,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true};
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL-modes','fm','IDEAL'};
% fm_meas_Hz=mcobj_ideal.Experimental.fm_est*(-2*pi)/(6.536 /42.567);
%% process all data and export nifti
mcobj_me_s3=cell(length(dirst_me),1);
for cf=1:length(dirst_me)

       fn=fullfile(MeasPath,dirst_me(cf).name);   
    mcobj_me_s3{cf}=MetCon_ME(fn,ME_setting{:},'fm','IDEAL','Solver','IDEAL-modes');
end
% CSI data

mcobj_csi_s3=cell(length(dirst_csi),1);
for cf=1:length(dirst_csi)
       fn=fullfile(MeasPath,dirst_csi(cf).name);
    mcobj_csi_s3{cf}=MetCon_CSI(fn,CSI_setting{:});
end
%write images
[intake_time_s3,order_idx_s3]=sort(cellfun(@(x) x.getMinutesAfterIntake('10:10'),[mcobj_csi_s3;mcobj_me_s3],'UniformOutput',true));

% reslice
cellfun(@(x) x.WriteImages('',{'snr','mM'}),[mcobj_csi_s3;mcobj_me_s3],'UniformOutput',false);
reslice_all_sag=myspm_reslice('../anat/*anat*_sag*.nii',dir('Metcon_*.nii'), 'nearest','rs');
reslice_all_tra=myspm_reslice('../anat/*anat*_tra*.nii',dir('Metcon_*.nii'), 'nearest','rt');

type_label_s3 =[1     2     3     1     2     3 1];
%%
% load('/ptmp/pvalsala/deuterium/paper/sub2_H4_modes/data_S2_modes.mat')
% load('/ptmp/pvalsala/deuterium/paper/sub3_modes/data_S3_modes.mat','intake_time_s3','vox_vol_s3','type_label_s3')

% read mifti and try imcompose
for cSub=[1 3]
    switch(cSub)
        case 1
             cd(pn_s1);
        case 3
            cd(pn_s3);
    end

    anat_nii_sag{cSub}=double(MyNiftiRead("../anat/*anat_sag.nii",'IPR'));
    brain_nii_sag{cSub}=double(MyNiftiRead("../anat/*brain_sag.nii",'IPR'));
    brain_nii_sag{cSub}=brain_nii_sag{cSub}./max(brain_nii_sag{cSub},[],'all');
     anat_nii_sag{cSub}=anat_nii_sag{cSub}./max(anat_nii_sag{cSub},[],'all');


    anat_nii_tra{cSub}=double(MyNiftiRead("../anat/*anat_tra.nii",'PRS'));
    brain_nii_tra{cSub}=double(MyNiftiRead("../anat/*brain_tra.nii",'PRS'));
    brain_nii_tra{cSub}=brain_nii_tra{cSub}./max(brain_nii_tra{cSub},[],'all');
     anat_nii_tra{cSub}=anat_nii_tra{cSub}./max(anat_nii_tra{cSub},[],'all');

    met_mm_me_sag{cSub}=MyNiftiRead('rsMetcon_mM*.nii','IPR');
    met_mm_me_tra{cSub}=MyNiftiRead('rtMetcon_mM*.nii','PRS');
    %show two of each method
   met_mm_me_sag{cSub}(:,:,:,:,7)=[];
   met_mm_me_tra{cSub}(:,:,:,:,7)=[];
end


%There is some misregistration
met_mm_me_tra{1}=(circshift(met_mm_me_tra{1},-7,2));
met_mm_me_tra{3}=(circshift(met_mm_me_tra{3},-5,2));

type_label_s1 =[1     2     3     1     3     2     3];
type_label_s3 =[1     2     3     1     2     3 1];
   type_label_s1(7)=[];
   type_label_s3(7)=[];

% met_mm_me_tra(:,:,:,1,:)=met_snr_me_tra(:,:,:,1,:);
%%

cSub=1;
slct_s1=15;
slcs_s1=20;
% crop_val  sag             tra : [idx1:end-idx2,idx3:end-idx4]
crop_val_s1={[45 70 40 30],[25 35 20 55]};

[im_anat{cSub}]=ComposeImage2({anat_nii_sag{cSub},anat_nii_tra{cSub}},{slcs_s1,slct_s1},crop_val_s1);

[im_brain{cSub}]=ComposeImage2({brain_nii_sag{cSub},brain_nii_tra{cSub}},{slcs_s1,slct_s1},crop_val_s1);
[imComposed_csi{cSub}]=ComposeImage2({met_mm_me_sag{cSub},met_mm_me_tra{cSub}}, ...
    {slcs_s1,slct_s1},crop_val_s1,[0 5 0]);




cSub=3;
slct_s3=14;
slcs_s3=18;
% crop_val  sag             tra : [idx1:end-idx2,idx3:end-idx4]
crop_val_s3={[45 60 30 20],[35 25 35 55]};

[im_anat{cSub}]=ComposeImage2({anat_nii_sag{cSub},anat_nii_tra{cSub}},{slcs_s3,slct_s3},crop_val_s3);

[im_brain{cSub}]=ComposeImage2({brain_nii_sag{cSub},brain_nii_tra{cSub}},{slcs_s3,slct_s3},crop_val_s3);
[imComposed_csi{cSub}]=ComposeImage2({met_mm_me_sag{cSub},met_mm_me_tra{cSub}}, ...
    {slcs_s3,slct_s3},crop_val_s3,[0 5 0]);


data_label={'CSI','CSI-PC-bSSFP','ME-PC-bSSFP'};
midlabel= @(Sz,N) Sz/2:Sz:Sz*N*1.5;
title_str={'Water [SNR]','Glucose [mM]','Glx [mM]'};




% plot figure
cSub=1;
load('/ptmp/pvalsala/deuterium/fig8_sub1and3/sub-01/MC/masks_s1.mat')
% mask_me=(im_brain{cSub}./max(im_brain{cSub}(:)))>0.01;
mask_me=m_tra|m_sag;
% mask_me=imComposed_csi(:,:,1,1,1)>8;
mask_me=imerode(mask_me,strel('sphere',10));
mask_me=imdilate(mask_me,strel('sphere',8));
mask_me=repmat(mask_me,[1 size(met_mm_me_sag{cSub},5)]);
cax_met={[0 70],[0 2.2],[0 2.2]}
figure(45),clf
set(gcf,"Position",[50 35 1425 1080])
tt=tiledlayout(4,2,'TileSpacing','compact','Padding','loose');
sgtitle('CSI vs CSI-bSSFP','fontsize',24)
cb_all={};cax_all={};cnt=1;

common_setttings={'cax_im',[0 0.5],'prctile',98,'SlcSel',1,'transform',@(x) x,'Mask',mask_me};

i=2;
nexttile(tt,1,[1 1]);
[cb_all{cnt},cax_all{cnt}]=overlayplot(repmat(im_anat{cSub},[1 size(met_mm_me_sag{cSub},5)]), ...
    createImMontage(squeeze(imComposed_csi{cSub}(:,:,1,2,:)),size(met_mm_me_sag{cSub},5)),'MetIdx',1,...
    'cax',cax_met{i},common_setttings{:});
% xticks(midlabel(size(im_anat{cSub},2),4)),xticklabels(intake_time_s1)
yticks([])
title(title_str{i})
cb_all{cnt}=colorbar;
cnt=cnt+1;

i=3;
nexttile(tt,5,[1 1]);
[cb_all{cnt},cax_all{cnt}]=overlayplot(repmat(im_anat{cSub},[1 size(met_mm_me_sag{cSub},5)]),createImMontage(squeeze(imComposed_csi{cSub}(:,:,1,3,:)),size(met_mm_me_sag{cSub},5)),'MetIdx',1,'prctile',98,'SlcSel',1,'transform',@(x) x,...
    'cax',cax_met{i},common_setttings{:});
yticks([])
% if(i==1),yticks(midlabel(size(im_anat{cSub},1),3)),yticklabels(intake_time_s1),else,yticks([]),end
title(title_str{i})
cb_all{cnt}=colorbar;
cnt=cnt+1;


cSub=3;
load('/ptmp/pvalsala/deuterium/fig8_sub1and3/sub-03/MC/masks_s3.mat')
mask_me=m_tra|m_sag;

% mask_me=imComposed_csi(:,:,1,1,1)>8;
mask_me=imerode(mask_me,strel('sphere',10));
mask_me=imdilate(mask_me,strel('sphere',8));
mask_me=repmat(mask_me,[1 size(met_mm_me_sag{cSub},5)]);
common_setttings{end}=mask_me   ;

i=2;
nexttile(tt,2,[1 1]);
[cb_all{cnt},cax_all{cnt}]=overlayplot(repmat(im_anat{cSub},[1 size(met_mm_me_sag{cSub},5)]),createImMontage(squeeze(imComposed_csi{cSub}(:,:,1,2,:)),size(met_mm_me_sag{cSub},5)),'MetIdx',1,...
    'cax',cax_met{i},common_setttings{:});
if(i==1),yticks(midlabel(size(im_anat{cSub},1),3)),yticklabels(intake_time_fisp),else,yticks([]),end
title(title_str{i}),yticks([])
cb_all{cnt}=colorbar;
cnt=cnt+1;

i=3;
nexttile(tt,6,[1 1]);
[cb_all{cnt},cax_all{cnt}]=overlayplot(repmat(im_anat{cSub},[1 size(met_mm_me_sag{cSub},5)]),createImMontage(squeeze(imComposed_csi{cSub}(:,:,1,3,:)),size(met_mm_me_sag{cSub},5)),'MetIdx',1,...
    'cax',cax_met{i},common_setttings{:});
title(title_str{i}),yticks([])
cb_all{cnt}=colorbar;
cnt=cnt+1;

for i=1:length(cb_all)
    makeColorbar(cb_all{i},cax_all{i});
end

set(gcf,'Color','w','InvertHardcopy','off')
%

% stats from registered maps
cSub=1;
mask=brain_nii_tra{cSub}(:,:,:,1,1)>1e-3;
mask=imerode(mask,strel('disk',10,0));
mask=imdilate(mask,strel('disk',10,0));
% as(mask.*anat_nii_tra{cSub})

col_tra=reshape(met_mm_me_tra{cSub},[],4,size(met_mm_me_tra{cSub},5));

colors_label=lines(3);

mask=mask(:) & all(col_tra(:,2,:)<3,3);
mask=mask(:) & all(col_tra(:,3,:)<3,3);

% figure,
% tt=tiledlayout(1,3);
violin_format={'mc',[],'medc',[]};

nexttile(tt,3,[1 1])
violin(squeeze(col_tra(mask(:),2,:)),violin_format{:},'facecolor',colors_label(type_label_s1,:));
grid on,grid minor,xticks(1:7),xticklabels(intake_time_s1),set(gca,'FontWeight','bold')

mean_val=mean(squeeze(col_tra(mask(:),2,:)));
hold on,plot(find(type_label_s1==1),mean_val(type_label_s1==1),'Marker','x','Color',colors_label(1,:),'LineWidth',1.5)
hold on,plot(find(type_label_s1==2),mean_val(type_label_s1==2),'Marker','x','Color',colors_label(2,:),'LineWidth',1.5)
hold on,plot(find(type_label_s1==3),mean_val(type_label_s1==3),'Marker','x','Color',colors_label(3,:),'LineWidth',1.5)




nexttile(tt,7,[1 1])
vh=violin(squeeze(col_tra(mask(:),3,:)),violin_format{:},'facecolor',colors_label(type_label_s1,:));
grid on,grid minor, xticks(1:7),xticklabels(intake_time_s1),set(gca,'FontWeight','bold')

mean_val=mean(squeeze(col_tra(mask(:),3,:)));
hold on,plot(find(type_label_s1==1),mean_val(type_label_s1==1),'Marker','x','Color',colors_label(1,:),'LineWidth',1.5)
hold on,plot(find(type_label_s1==2),mean_val(type_label_s1==2),'Marker','x','Color',colors_label(2,:),'LineWidth',1.5)
hold on,plot(find(type_label_s1==3),mean_val(type_label_s1==3),'Marker','x','Color',colors_label(3,:),'LineWidth',1.5)
xlabel('Time after intake [min]')

%%%
cSub=3;
mask=brain_nii_tra{cSub}(:,:,:,1,1)>1e-3;
mask=imerode(mask,strel('disk',10,0));
mask=imdilate(mask,strel('disk',10,0));
% as(mask.*anat_nii_tra{cSub})

col_tra=reshape(met_mm_me_tra{cSub},[],4,size(met_mm_me_tra{cSub},5));

colors_label=lines(3);

mask=mask(:) & all(col_tra(:,2,:)<3,3);
mask=mask(:) & all(col_tra(:,3,:)<3,3);

% figure,
% tt=tiledlayout(1,3);
violin_format={'mc',[],'medc',[]};

nexttile(tt,4,[1 1])
violin(squeeze(col_tra(mask(:),2,:)),violin_format{:},'facecolor',colors_label(type_label_s3,:));
grid on,grid minor,xticks(1:7),xticklabels(intake_time_s3),set(gca,'FontWeight','bold')

mean_val=mean(squeeze(col_tra(mask(:),2,:)));
hold on,plot(find(type_label_s3==1),mean_val(type_label_s3==1),'Marker','x','Color',colors_label(1,:),'LineWidth',1.5)
hold on,plot(find(type_label_s3==2),mean_val(type_label_s3==2),'Marker','x','Color',colors_label(2,:),'LineWidth',1.5)
hold on,plot(find(type_label_s3==3),mean_val(type_label_s3==3),'Marker','x','Color',colors_label(3,:),'LineWidth',1.5)



nexttile(tt,8,[1 1])
vh=violin(squeeze(col_tra(mask(:),3,:)),violin_format{:},'facecolor',colors_label(type_label_s3,:));
grid on,grid minor,xticks(1:7),xticklabels(intake_time_s3),set(gca,'FontWeight','bold')

mean_val=mean(squeeze(col_tra(mask(:),3,:)));
hold on,plot(find(type_label_s3==1),mean_val(type_label_s3==1),'Marker','x','Color',colors_label(1,:),'LineWidth',1.5)
hold on,plot(find(type_label_s3==2),mean_val(type_label_s3==2),'Marker','x','Color',colors_label(2,:),'LineWidth',1.5)
hold on,plot(find(type_label_s3==3),mean_val(type_label_s3==3),'Marker','x','Color',colors_label(3,:),'LineWidth',1.5)
xlabel('Time after intake [min]')
legend(gca,vh([1 2 3]),data_label([1 2 3]),'Position',[0.620349794317177 0.00128943966208529 0.370370363554469 0.020370369873665],...
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



function [imComposed]=ComposeImage2(Vol,slcSel,crop_val,shift)


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

if(imsz1(2)>imsz2(2))
    im{2}=padarray(im{2},[0 floor(diff/2)],0,'pre');
    im{2}=padarray(im{2},[0 ceil(diff/2)],0,'post');
else
    im{1}=padarray(im{1},[0 floor(diff/2)],0,'pre');
    im{1}=padarray(im{1},[0 ceil(diff/2)],0,'post');


end
imComposed=cat(1,im{:});
end
