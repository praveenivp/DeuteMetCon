%%
load('/ptmp/pvalsala/deuterium/paper/sub1_HOSJ_modes/data_S1_modes.mat')
load('/ptmp/pvalsala/deuterium/paper/sub2_H4_modes/data_S2_modes.mat')
load('/ptmp/pvalsala/deuterium/paper/sub3_modes/data_S3_modes.mat')

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