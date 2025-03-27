%%
load('/ptmp/pvalsala/deuterium/paper/sub1_HOSJ_modes/data_S1_modes.mat')
load('/ptmp/pvalsala/deuterium/paper/sub2_H4_modes/data_S2_modes.mat')
load('/ptmp/pvalsala/deuterium/paper/sub3_modes/data_S3_modes.mat')
%%
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'));


%%
figure(2),clf
tt=tiledlayout(3,1,'Padding','compact','TileSpacing','none');
slc_csl_me=20:2:42;

nexttile()
mcg=flip(permute(mcobj_me{1}.getNormalized().*mcobj_me{2}.coilNormMat,[1 3 2 4]),3);
imT=@(xin)flip(permute(xin,[1 3 2 4 5]),3);
%  mcg=flip(permute(mcobj_me{2}.getNormalized(),[1 3 2 4]),3);
mcg=imT(resliced_s1(:,:,:,:,type_label_s1==3));
cIntakeTime=intake_time_s1(type_label_s1==3);
imagesc(createImMontage(mcg(:,:,slc_csl_me,2,3),length(slc_csl_me)))
axis image
xticks([]),yticks(48/2),yticklabels(sprintf('S1 (%d min)',cIntakeTime(3)))
title('ME-bSSFP glucose maps [SNR]')
% colorbar,
clim([0 30])

nexttile()
mcg=imT(resliced_s2(:,:,:,:,type_label_s2==3));
cIntakeTime=intake_time_s2(type_label_s2==3);
imagesc(createImMontage(mcg(:,:,slc_csl_me,2,2),length(slc_csl_me)))
axis image
xticks([]),yticks(48/2),yticklabels(sprintf('S2 (%d min)',cIntakeTime(2)))
colorbar,
clim([0 30])

nexttile()
mcg=imT(resliced_s3(:,:,:,:,type_label_s3==3));
cIntakeTime=intake_time_s3(type_label_s3==3);
imagesc(createImMontage(mcg(:,:,slc_csl_me,2,1),length(slc_csl_me)))
axis image
xticks([]),yticks(48/2),yticklabels(sprintf('S3 (%d min)',cIntakeTime(2)))
% colorbar,
clim([0 30])

fontsize(gcf,"scale",1.5)
set(gcf,'color','w','Position',[73 459 1612 552])

%% Better figure with Anatomical overlay

%sub1
s1=fullfile('/ptmp/pvalsala/deuterium/paper/sub1_HOSJ_modes', ...
    'Metcon_SNR_m01011_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_IDEAL-modes.nii');
s1_anat=fullfile('/ptmp/pvalsala/deuterium/dataForPublication/anat/sub-01','anat_tra.nii');
resliced_s1=myspm_reslice(dir(s1_anat),dir(s1), 'nearest','rt_');


%sub2
s1=fullfile('/ptmp/pvalsala/deuterium/paper/sub2_H4_modes', ...
    'Metcon_SNR_m01182_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_IDEAL-modes.nii');
s2_anat=fullfile('/ptmp/pvalsala/deuterium/dataForPublication/anat/sub-02','anat_tra.nii');
resliced_s2=myspm_reslice(dir(s2_anat),dir(s1), 'nearest','rt_');

%sub3
s3=fullfile('/ptmp/pvalsala/deuterium/paper/sub3_modes', ...
    'Metcon_SNR_m00118_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_IDEAL-modes.nii');
s3_anat=fullfile('/ptmp/pvalsala/deuterium/dataForPublication/anat/sub-03','anat_tra.nii');
resliced_s3=myspm_reslice(dir(s3_anat),dir(s3), 'nearest','rt_');


%% reslice
s1=MyNiftiRead(fullfile('/ptmp/pvalsala/deuterium/paper/sub1_HOSJ_modes', ...
    'rt_Metcon_SNR_m01011_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_IDEAL-modes.nii'),'PRI');
s1_anat=MyNiftiRead(fullfile('/ptmp/pvalsala/deuterium/dataForPublication/anat/sub-01','anat_tra.nii'),'PRI');

%sub2
s2=MyNiftiRead(fullfile('/ptmp/pvalsala/deuterium/paper/sub2_H4_modes', ...
    'rt_Metcon_SNR_m01182_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_IDEAL-modes.nii'),'PRI');
s2_anat=MyNiftiRead(fullfile('/ptmp/pvalsala/deuterium/dataForPublication/anat/sub-02','anat_tra.nii'),'PRI');

%sub3
s3=MyNiftiRead(fullfile('/ptmp/pvalsala/deuterium/paper/sub3_modes', ...
    'rt_Metcon_SNR_m00118_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_IDEAL-modes.nii'),'PRI');
s3_anat=MyNiftiRead(fullfile('/ptmp/pvalsala/deuterium/dataForPublication/anat/sub-03','anat_tra.nii'),'PRI');

%%
figure(2)
tt=tiledlayout(3,1,'TileSpacing','none','Padding','compact');
cb_all={};cax_all={};cnt=1;
% for i=1:3
    nexttile(tt,[1 1]);

    sel_s1=@(x)createImMontage( x(:,25:end-60,4:3:end,2),10);
    sel_s1_anat=@(x)createImMontage( x(:,25:end-60,4:3:end),10);
[cb_all{cnt},cax_all{cnt}]=overlayplot(sel_s1_anat(s1_anat),sel_s1(s1),'MetIdx',1,'prctile',98,'SlcSel',1,'transform',@(x) x,...
  'cax',[0 20],'Mask',sel_s1_anat(s1_anat)>0.1,'alpha_overlay',0.4,'cax_im',[0 0.6]);

cIntakeTime=intake_time_s1(type_label_s1==3);
axis image
xticks([]),yticklabels(sprintf('S1 (%d min)',cIntakeTime(end)))
title('ME-bSSFP Glucose maps [SNR]','fontsize',14)

    nexttile(tt,[1 1]);
cnt=cnt+1;
    sel_s2=@(x)createImMontage( x(15:end-30,45:end-50,2:3:end-3,2),10);
    sel_s2_anat=@(x)createImMontage( x(15:end-30,45:end-50,2:3:end-3),10);
[cb_all{cnt},cax_all{cnt}]=overlayplot(sel_s2_anat(s2_anat),sel_s2(s2),'MetIdx',1,'prctile',98,'SlcSel',1,'transform',@(x) x,...
  'cax',[0 20],'Mask',sel_s2_anat(s2_anat)>0.1,'alpha_overlay',0.4,'cax_im',[0 0.6]);
cb_all{cnt}=colorbar();
cIntakeTime=intake_time_s2(type_label_s2==3)-23; % some reason there is 23 min offset;
axis image
xticks([]),yticklabels(sprintf('S2 (%d min)',cIntakeTime(end)))


    nexttile(tt,[1 1]);
cnt=cnt+1;
    sel_s3=@(x)createImMontage( x(:,35:end-55,4:3:end,2),10);
    sel_s3_anat=@(x)double(createImMontage( x(:,35:end-55,4:3:end),10));
[cb_all{cnt},cax_all{cnt}]=overlayplot(sel_s3_anat(s3_anat),sel_s3(s3),'MetIdx',1,'prctile',98,'SlcSel',1,'transform',@(x) x,...
  'cax',[0 20],'Mask',sel_s3_anat(s3_anat)>0.1,'alpha_overlay',0.4,'cax_im',[0 0.6]);


cIntakeTime=intake_time_s3(type_label_s3==3);
axis image
xticks([]),yticklabels(sprintf('S3 (%d min)',cIntakeTime(end)))

makeColorbar(cb_all{2},cax_all{2})
set(gcf,'color','w')

%% %%%%%%%%%%%%%%%%%   supporting CSI data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('/ptmp/pvalsala/deuterium/DDRM-T4EU/proc/skull')
metabolites=getMetaboliteStruct('invivo');
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','AMARES','fm',[]};
fn=fullfile('/ptmp/pvalsala/deuterium/DDRM-T4EU/TWIX','allData#S94Tuebingen#F52148#M991#D220523#T102734#rpcsi_fid_Stan_Res19_20mm.dat');
mcobj_csi_DDRM=MetCon_CSI(fn,CSI_setting{:});
mc=(mcobj_csi_DDRM.getNormalized);
inTime=mcobj_csi_DDRM.getMinutesAfterIntake('09:15');  %78 min

%%
anat_tra=myspm_reslice(dir("M00991_rpcsi_fid_Stan_Res19_20mm.nii"),dir("anat_tra*.nii"),'nearest','rt_');
anat_tra=double(anat_tra)./double(max(anat_tra(:)));
%%
spec=myfft( mcobj_csi_DDRM.img,5);
skull_mask=(anat_tra(:,:,:,1)>0.05&anat_tra(:,:,:,2)==0);
brain_mask=(anat_tra(:,:,:,2)>0.1);

slcSel=30;
skull_mask_sel=zeros(size(skull_mask),'logical');
skull_mask_sel(:,:,slcSel)=skull_mask(:,:,slcSel);

brain_mask_sel=zeros(size(brain_mask),'logical');
brain_mask_sel(:,:,slcSel)=brain_mask(:,:,slcSel);



%%
[expParams,pk]=getAMARES_structs(mcobj_csi_DDRM);
skull_ind=find(skull_mask_sel);
for j=1:length(skull_ind)
    [pxl_idx(1),pxl_idx(2),pxl_idx(3)]=ind2sub(size(mcobj_csi_DDRM.img,2:4),skull_ind(j));
    fid=squeeze(mcobj_csi_DDRM.img(1,pxl_idx(1),pxl_idx(2),pxl_idx(3),:));
    fm_est=0;

    expParams.offset=fm_est;
    [fitResults, fitStatus,~,CRBResults] = AMARES.amaresFit(double(fid(:)), expParams, pk,33,'quiet',false);

    % get data from axes
    fh=gcf();
    ax=fh.Children(4);
    lines = findobj(ax, 'type', 'line');
    ydata_skull=[];
    for i = 1:length(lines)
        xdata = get(lines(i), 'XData');
        ydata_skull(:,i,j) = get(lines(i), 'YData');

    end
end
brain_ind=find(brain_mask_sel);
for j=1:5:length(brain_ind)
    [pxl_idx(1),pxl_idx(2),pxl_idx(3)]=ind2sub(size(mcobj_csi_DDRM.img,2:4),brain_ind(j));
    fid=squeeze(mcobj_csi_DDRM.img(1,pxl_idx(1),pxl_idx(2),pxl_idx(3),:));
    fm_est=0;

    expParams.offset=fm_est;
    [fitResults, fitStatus,~,CRBResults] = AMARES.amaresFit(double(fid(:)), expParams, pk,33,'quiet',false);

    % get data from axes
    fh=gcf();
    ax=fh.Children(4);
    lines = findobj(ax, 'type', 'line');
    ydata_brain=[];
    for i = 1:length(lines)
        xdata = get(lines(i), 'XData');
        ydata_brain(:,i,j) = get(lines(i), 'YData');

    end
end

%


%%

img = double(anat_tra(:,:,slcSel,1));
mask_1 = double(skull_mask(:,:,slcSel));
mask_2 = double(brain_mask(:,:,slcSel));


% Create a red transparent mask
red_mask = cat(3, ones(size(mask_1)), zeros(size(mask_1)), zeros(size(mask_1)));
blue_mask = cat(3, zeros(size(mask_1)), zeros(size(mask_1)), ones(size(mask_1)));
overlay1 = uint8(red_mask .* 255);
overlay2 = uint8(blue_mask .* 255);

% Overlay the red transparent mask on the original image
figure(9),clf;
tt=tiledlayout(2,4);
% nexttile()
% imagesc(img);
% title('Structural image')
% xticks([]),yticks([])

nexttile(1,[2,2])
imagesc(img);
xticks([]),yticks([]),colormap(gca,'gray')
hold on;
h1 = imshow(overlay1);
set(h1, 'AlphaData', mask_1*0.2);
h2 = imshow(overlay2);
set(h2, 'AlphaData', mask_2*0.2);
title('Segmentation')
% for idx=1:4
% nexttile()
% imagesc((squeeze(mc(:,:,slcSel,idx))))
% xticks([]),yticks([])
% title([mcobj_csi_DDRM.metabolites(idx).name,' [SNR]'])
% end


nexttile(3,[1 2])
plot(xdata+4.8,mean(ydata_brain,3),'LineWidth',1.5)
ax=gca;
ax.XAxis.Direction='reverse';
title('Brain voxels (blue)')
xlabel('frequency [ppm]')
xlim([0.5 7])
legend('fit','data')
nexttile(7,[1 2])
plot(xdata+4.8,mean(ydata_skull,3),'LineWidth',1.5)
title('skull voxels (red)')
fontsize(gcf,'scale',1.3)
ax=gca;
ax.XAxis.Direction='reverse';
set(gcf,'color','w')
xlabel('frequency [ppm]')
xlim([0.5 7])
legend('fit','data')
%%
function makeColorbar(cb_handle,cax)



cb_handle.Visible='off';
ax2 = axes('Position',cb_handle.Position);

imagesc(linspace(0,1,100)'),colormap(ax2,'turbo')
set(ax2,'YAxisLocation','right','FontSize',10,'FontWeight','bold','YDir','normal')
yticks(linspace(0,100,4))
yticklabels(round(linspace(0,1,4)*cax(2),0))
xticks([])
end

