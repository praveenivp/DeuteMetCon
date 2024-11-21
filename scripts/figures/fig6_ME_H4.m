%% process and write images
%% input files


MeasPath='/ptmp/pvalsala/deuterium/H4DK-64GH';
sn=fullfile(MeasPath,'TWIX');

dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));
% dirst_csi=dirst_csi([1 3 4]);
pn=fullfile(MeasPath,sprintf('proc/ME_bSSFP_%s',datetime('today','Format','yyyyMMMdd')));
mkdir(pn);
dirst_me=dir(fullfile(sn,"*rh_trufi*.dat"));
dirst_me=dirst_me([1 3 4]);
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))

metabolites=getMetaboliteStruct('invivo');


%% 
% dirst_fm=dir(fullfile(sn,'*piv_gre_B0mapping*'));
%  fmobj=B0map(fullfile(sn,dirst_fm(1).name),'UnwrapMode','SpatialUnwrap','doRegularization',false);
% fmobj.exportNIFTI()
%% process noise anf get fieldmap
fn_noise=dir(fullfile(sn,"*rh_trufi*noise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

fn_fm=(fullfile(sn,dirst_me(3).name));
ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',true,'doZeropad',[1 1 1]*0.5,'parfor',true};
% mcobj_ideal=MetCon_bSSFP(fn_fm,'Solver','IDEAL',ME_setting{:},'mask',80);
% fm_meas_Hz=mcobj_ideal.Experimental.fm_est*(-2*pi)/(6.536 /42.567);
%% process all ME- data

% ME_setting=[ME_setting,{'fm',fm_meas_Hz,'Solver','pinv'}];
mcobj_me=cell(length(dirst_me),1);
for cf=1:length(dirst_me)

       fn=fullfile(sn,dirst_me(cf).name);   
    mcobj_me{cf}=MetCon_bSSFP(fn,ME_setting{:},'fm','IDEAL','Solver','IDEAL-modes');
%     mcobj_me{cf-1}=MetCon_bSSFP(fn,ME_setting{:},'fm','M01169_piv_gre_B0mapping_5Echoes_fmap.nii','Solver','pinv');
end

%% Process CSI
dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));
dirst_csi=dirst_csi([1 3]);

CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL','fm','IDEAL'};


mcobj_csi=cell(length(dirst_csi),1);
for cf=1:length(dirst_csi)
       fn=fullfile(sn,dirst_csi(cf).name);
    mcobj_csi{cf}=MetCon_CSI(fn,CSI_setting{:});
end


%%
cd(pn)
intake_time_me=cellfun(@(x) x.getMinutesAfterIntake('08:36'),mcobj_me,'UniformOutput',true);
ylabel_str_me=strsplit(sprintf('%.0f min\n',intake_time_me(:)),'\n');

intake_time_csi=cellfun(@(x) x.getMinutesAfterIntake('08:36'),mcobj_csi,'UniformOutput',true);
ylabel_str_csi=strsplit(sprintf('%.0f min\n',intake_time_csi(:)),'\n');
[intake_time_all,intake_order]=sort([intake_time_csi;intake_time_me]);

norm_mat_csi=abs(mcobj_csi{1}.getNormalized());
norm_mat_csi=imgaussfilt3(norm_mat_csi(:,:,:,1),2);
cellfun(@(x) x.WriteImages([],{'snr','mM'}, 1./norm_mat_csi),mcobj_csi,'UniformOutput',false);
norm_mat_me=abs(mcobj_me{1}.getNormalized());
 norm_mat_me=imgaussfilt3(norm_mat_me(:,:,:,1),2);
%  norm_mat_me=medfilt3(norm_mat_me(:,:,:,1),[5 5 5]);
% norm_mat_me(norm_mat_me<0.5)=0.5;
  cellfun(@(x) x.WriteImages([],{'snr','mM'},1./norm_mat_me),mcobj_me,'UniformOutput',false);
%   cellfun(@(x) x.WriteImages([],{'snr','mM'}),mcobj_me,'UniformOutput',false);
SNR_scl=prod(mcobj_csi{1}.DMIPara.resolution_PSF(1:3))./prod(mcobj_me{1}.DMIPara.resolution_PSF(1:3));
% reslice images to anatomy
dirst_nii=dir('Metcon_*.nii');
reslice_all_sag=myspm_reslice('../*anat*_sag*.nii',dirst_nii, 'nearest','rs');
reslice_all_tra=myspm_reslice('../*anat*_tra*.nii',dirst_nii, 'nearest','rt');


%% try imcompose
anat_nii_sag=double(MyNiftiRead("../*anat*_sag*.nii",'IPR'));

[met_snr_me_sag]=MyNiftiRead('rsMetcon_SNR*.nii','IPR');
[met_mm_me_sag]=MyNiftiRead('rsMetcon_mM*.nii','IPR');
 met_mm_me_sag(:,:,:,1,:)=met_snr_me_sag(:,:,:,1,:);

anat_nii_tra=double(MyNiftiRead("../*anat*_tra*.nii",'PRS'));


[met_snr_me_tra]=MyNiftiRead('rtMetcon_SNR*.nii','PRS');
[met_mm_me_tra]=MyNiftiRead('rtMetcon_mM*.nii','PRS');
  met_mm_me_tra(:,:,:,1,:)=met_snr_me_tra(:,:,:,1,:);
%%
slct=15;
slcs=18;
% crop_val  sag             tra : [idx1:end-idx2,idx3:end-idx4] 
crop_val={[20 40 45 25],[20 40 45 50]};
[im_anat]=ComposeImage({anat_nii_sag,anat_nii_tra},{slcs,slct},crop_val);

midlabel= @(Sz,N) Sz/2:Sz:Sz*N*1.5;
title_str={'Water [SNR]','Glc [mM]','Glx [mM]'};
%  as(im_anat)

%
[imComposed_csi]=ComposeImage({met_mm_me_sag(:,:,:,:,[1 4]),met_mm_me_tra(:,:,:,:,[1 4])}, ...
    {slcs,slct},crop_val,[0 5 0]);
% as(imComposed2)
% as(imComposed)

%  mask_me=(im_anat./max(im_anat(:)))>0.05;
mask_me=imComposed_csi(:,:,1,1,1)>8;
mask_me=imerode(mask_me,strel('sphere',10));
mask_me=imdilate(mask_me,strel('sphere',10));

cax_met={[0 70],[0 3],[0 3]}
figure(45),clf
    set(gcf,"Position",[50 50 1080 1080])
tt=tiledlayout(7,3,'TileSpacing','compact','Padding','compact');
sgtitle('CSI-FISP vs ME-bSSFP','fontsize',24)
cb_all={};cax_all={};cnt=1;
for i=1:3
    nexttile(tt,[2 1]);
[cb_all{cnt},cax_all{cnt}]=overlayplot(im_anat,imComposed_csi,'MetIdx',i,'prctile',98,'SlcSel',1,'transform',@(x) x,...
  'cax',cax_met{i},'Mask',mask_me,'alpha_overlay',0.6,'cax_im',[0 0.35]);
if(i==1),yticks(midlabel(size(im_anat,1),3)),yticklabels(intake_time_csi),else,yticks([]),end
title(title_str{i})
 cb_all{cnt}=colorbar;
cnt=cnt+1;

end

cax_met{1}=[0 50]
[imComposed_me]=ComposeImage({met_mm_me_sag(:,:,:,:,[2 3 5]),met_mm_me_tra(:,:,:,:,[2 3 5])}, ...
    {slcs,slct},crop_val,[0,0,0]);
imComposed_me(:,:,:,1,:)=imComposed_me(:,:,:,1,:)*SNR_scl;
for i=1:3
    nexttile(tt,[3 1]);
[cb_all{cnt},cax_all{cnt}]=overlayplot(im_anat,imComposed_me,'MetIdx',i,'prctile',98,'SlcSel',1,'transform',@(x) x,...
  'cax',cax_met{i},'Mask',mask_me,'alpha_overlay',0.6,'cax_im',[0 0.35]);
if(i==1),yticks(midlabel(size(im_anat,1),4)),yticklabels(intake_time_me),else,yticks([]),end
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

mask=met_mm_me_tra(:,:,:,1,1)>8;
mask=imerode(mask,strel('sphere',10));
mask=imdilate(mask,strel('sphere',10));
% as(mask.*anat_nii_tra)

col_tra=reshape(met_mm_me_tra,[],4,5);
data_label={'CSI-FISP','CSI-bSSFP','ME-bSSFP'};
type_label_s2=[1 3 3 1 3];
colors_label=lines(3);

mask=mask(:) & all(col_tra(:,2,:)<4,3);
mask=mask(:) & all(col_tra(:,3,:)<4,3);

% figure,
% tt=tiledlayout(1,3);
violin_format={'mc',[],'medc',[]};
nexttile(tt,[2 1])

violin(squeeze(col_tra(mask(:),1,:)).*[1;SNR_scl;SNR_scl;1;SNR_scl  ]',violin_format{:},'facecolor',colors_label(type_label_s2,:));
grid on, title(title_str{1}),xticklabels(intake_time_all),set(gca,'FontWeight','bold')

mean_val=mean(squeeze(col_tra(mask(:),1,:)).*[1;SNR_scl;SNR_scl;1;SNR_scl  ]');
hold on,plot(find(type_label_s2==1),mean_val(type_label_s2==1),'Marker','x','Color',colors_label(1,:),'LineWidth',1.5)
hold on,plot(find(type_label_s2==3),mean_val(type_label_s2==3),'Marker','x','Color',colors_label(3,:),'LineWidth',1.5)

nexttile(tt,[2 1])
violin(squeeze(col_tra(mask(:),2,:)),violin_format{:},'facecolor',colors_label(type_label_s2,:));
grid on, title(title_str{2}),xticklabels(intake_time_all),set(gca,'FontWeight','bold')

mean_val=mean(squeeze(col_tra(mask(:),2,:)).*[1;SNR_scl;SNR_scl;1;SNR_scl  ]');
hold on,plot(find(type_label_s2==1),mean_val(type_label_s2==1),'Marker','x','Color',colors_label(1,:),'LineWidth',1.5)
hold on,plot(find(type_label_s2==3),mean_val(type_label_s2==3),'Marker','x','Color',colors_label(3,:),'LineWidth',1.5)
xlabel('Time after intake [min]')


nexttile(tt,[2 1])
vh=violin(squeeze(col_tra(mask(:),3,:)),violin_format{:},'facecolor',colors_label(type_label_s2,:));
grid on, title(title_str{3}),xticklabels(intake_time_all),set(gca,'FontWeight','bold')

mean_val=mean(squeeze(col_tra(mask(:),3,:)).*[1;SNR_scl;SNR_scl;1;SNR_scl  ]');
hold on,plot(find(type_label_s2==1),mean_val(type_label_s2==1),'Marker','x','Color',colors_label(1,:),'LineWidth',1.5)
hold on,plot(find(type_label_s2==3),mean_val(type_label_s2==3),'Marker','x','Color',colors_label(3,:),'LineWidth',1.5)


legend(gca,vh([1 2]),data_label([1 3]),'Location','northwest')


%%

% 
% norm_mat_csi=abs(mcobj_csi{1}.getNormalized());
% norm_mat_csi=norm_mat_csi(:,:,:,1);
%  norm_mat_csi=imgaussfilt3(norm_mat_csi(:,:,:,1),1);
% cellfun(@(x) x.WriteImages([],{'snr','mM'}, 1./norm_mat_csi),mcobj_csi,'UniformOutput',false);
% norm_mat_me=abs(mcobj_me{1}.getNormalized());
%  norm_mat_me=norm_mat_me(:,:,:,1);
%    norm_mat_me=imgaussfilt3(norm_mat_me(:,:,:,1),1);
% cellfun(@(x) x.WriteImages([],{'snr','mM'},1./norm_mat_me),mcobj_me,'UniformOutput',false);
% 
% figure(2),clf
% tt=tiledlayout(2,2)
% mask_me=mcobj_me{1}.getMask(95);
% mask_me=imerode(mask_me,strel('sphere',2));
% 
% 
% mask_csi=mcobj_csi{1}.getMask(97);
% mask_csi=imerode(mask_csi,strel('sphere',2));
% % 
% % col_snr=cellfun(@(x)x.getNormalized(),mcobj_me,'UniformOutput',false);
% % col_snr=cat(5,col_snr{:});
% %  as(col_snr.*mask)
% % col_snr=reshape(col_snr,[],size(col_snr,4),size(col_snr,5));
% % col_snr=abs(col_snr(mask_me(:),:,:));
% 
% col_mM_csi=cellfun(@(x)x.getmM(1./norm_mat_csi),mcobj_csi,'UniformOutput',false);
% col_mM_csi=cellfun(@(x)reshape(x,[],size(x,4),size(x,5)),col_mM_csi,'UniformOutput',false);
% % mask(col_mM(:,:,:,2,1)>3)=false;
% col_mM_csi=cellfun(@(x)x(mask_csi(:),:,:),col_mM_csi,'UniformOutput',false);
% 
% 
% col_mM_me=cellfun(@(x)x.getmM(1./norm_mat_me),mcobj_me,'UniformOutput',false);
% col_mM_me=cellfun(@(x)reshape(x,[],size(x,4),size(x,5)),col_mM_me,'UniformOutput',false);
% % mask(col_mM(:,:,:,2,1)>3)=false;
% col_mM_me=cellfun(@(x)x(mask_me(:),:,:),col_mM_me,'UniformOutput',false);
% 
% 
% intake_time_all=[intake_time_csi;intake_time_me;];
% [intake_time_all,idx_order]=sort(intake_time_all);
% % 
% % col_mM=cat(5,col_mM{:});
% % mask(col_mM(:,:,:,2,1)>3)=false;
% %  as(col_mM.*mask)
% % col_mM=reshape(col_mM,[],size(col_mM,4),size(col_mM,5));
% % col_mM=abs(col_mM(mask(:),:,:));
% 
% col_mM=[col_mM_csi;col_mM_me];
% col_mM=col_mM(idx_order);
% 
% % for cMet=1:3
% % nexttile(tt,[2 1]),cla
% % violin(squeeze(col_snr(:,cMet,:)));
% % grid on
% % title([metabolites(cMet).name ,' [SNR]' ])
% % end
% acceptable_range ={[5 10],[0 5],[0 4]};
% for cMet=2:3
% nexttile(tt,[2 1]),cla
% % col_m<col_mM(:,cMet,:)
% 
% col_mM_cmet=cellfun(@(x)x(  x(:,cMet)>acceptable_range{cMet}(1) & x(:,cMet)<acceptable_range{cMet}(2) ,cMet),col_mM,'UniformOutput',false);
% 
% 
% violin2(col_mM_cmet')
% xticklabels(intake_time_me),
% 
% % ylim([0,4]),grid on
% title([metabolites(cMet).name ,' [mM]' ])
% end

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


%old code
% %% plot overlay
% anat_nii_sag=double(MyNiftiRead("../*anat*_sag*.nii",'IPR'));
% [met_snr_me_sag]=MyNiftiRead('rsMetcon_SNR*.nii','IPR');
% [met_mm_me_sag]=MyNiftiRead('rsMetcon_mM*.nii','IPR');
% met_mm_me_sag(:,:,:,1,:)=met_snr_me_sag(:,:,:,1,:);
% % close all
% figure(16),clf
% tt=tiledlayout(1,6,'TileSpacing','tight','Padding','compact');
% transform=@(x) ndflip(permute(x(1:end,1:end,:,:,:),[1:5]),[]);
% % title_str={'Glx','Glc','D20'};
% title_str={'Water [SNR]','Glc [mM]','Glx [mM]'};
% cax_met={[0 40],[0 3],[0 3]}
% 
% cb_all={};cax_all={};cnt=1;
% slct=16;
% slcs=30;
% for i=1:3
%     nexttile(1+(i-1)*2);
% [cb_all{cnt},cax_all{cnt}]=overlayplot(anat_nii_sag,met_mm_me_sag,'MetIdx',i,'prctile',95,'SlcSel',slcs,'transform',transform,...
%     'cax',cax_met{i});
% if(i==1),yticklabels(ylabel_str_me(1:end-1)); else, yticks([]); end
% title(title_str{i})
% cnt=cnt+1;
% end
% %
% 
% 
% anat_nii_tra=double(MyNiftiRead("../*anat*_tra*.nii",'PRS'));
% [met_snr_me_tra]=MyNiftiRead('rtMetcon_SNR*.nii','PRS');
% [met_mm_me_tra]=MyNiftiRead('rtMetcon_mM*.nii','PRS');
%  met_mm_me_tra(:,:,:,1,:)=met_snr_me_tra(:,:,:,1,:);
% % transform=@(x) ndflip(permute(x(30:end-20,1:end,:,:,:),[1 2 3 4 5]),[]);
% % title_str={'Glx','Glc','D20'};
% % title_str={'Glx/Water_{14 min}','Glc/Water_{14 min}','Water/Water_{-9 min}'};
% for i=1:3
%     nexttile(2+(i-1)*2);
% [cb_all{cnt},cax_all{cnt}]=overlayplot(anat_nii_tra,met_mm_me_tra,'MetIdx',i,'prctile',98,'SlcSel',slct,'transform',transform,...
%   'cax',cax_all{i});
% yticks([]);
% title(title_str{i})
% cb_all{cnt}=colorbar;
% cnt=cnt+1;
% 
% end
% 
% 
% 
% for i=4:length(cb_all)
% makeColorbar(cb_all{i},cax_all{i});
% end
% 
% %% stats
% % figure(13),clf
% mask_me=mcobj_me{1}.getMask(90);
% mask_me=imerode(mask_me,strel('sphere',2));
% col_snr=cellfun(@(x)x.getNormalized,mcobj_me,'UniformOutput',false);
% col_snr=cat(5,col_snr{:});
% % as(col_snr.*mask)
% col_snr=reshape(col_snr,[],size(col_snr,4),size(col_snr,5));
% col_snr=abs(col_snr(mask_me(:),:,:));
% 
% col_mM=cellfun(@(x)x.getmM,mcobj_me,'UniformOutput',false);
% col_mM=cat(5,col_mM{:});
% mask_me(col_mM(:,:,:,2,1)>3)=false;
%  as(col_mM.*mask_me)
% col_mM=reshape(col_mM,[],size(col_mM,4),size(col_mM,5));
% col_mM=abs(col_mM(mask_me(:),:,:));
% 
% 
% for cMet=1:3
% nexttile(tt,6+cMet),cla
% violin(squeeze(col_snr(:,cMet,:)));
% grid on
% title([metabolites(cMet).name ,' [SNR]' ])
% end
% for cMet=1:3
% nexttile(tt,9+cMet),cla
% % col_m<col_mM(:,cMet,:)
% violin(squeeze(col_mM(:,cMet,:)))
% xticklabels(intake_time_me),
% 
% ylim([0,4]),grid on
% title([metabolites(cMet).name ,' [mM]' ])
% end

