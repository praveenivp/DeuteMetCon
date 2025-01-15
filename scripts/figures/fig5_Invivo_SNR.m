
%subject 2
% MeasPath='/ptmp/pvalsala/deuterium/H4DK-64GH';
% sn=fullfile(MeasPath,'TWIX');
% dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));
% dirst_me=dir(fullfile(sn,"*rh_trufi*.dat"));
% dirst_me=dirst_me([1 3 4]);
% pn='/ptmp/pvalsala/deuterium/paper/sub2_H4_modes';
% type_label_s2=[1     3     2     3     1     2     3];
% subject 3
MeasPath='/ptmp/pvalsala/deuterium/I3BL-CJ5O';
sn=fullfile(MeasPath,'TWIX');
dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));
dirst_csi=dirst_csi(1:4);
dirst_me=dir(fullfile(sn,"*rh_trufi*.dat"));
dirst_me=dirst_me([1 3]);
pn='/ptmp/pvalsala/deuterium/paper/sub3_modes';
type_label_s3=[1     2     3     1     2     3    ];

mkdir(pn)
cd(pn)

addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))

metabolites=getMetaboliteStruct('invivo');


%% process ME and CSI
fn_noise=dir(fullfile(sn,"*rh_trufi*noise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true};
% process all ME- data
mcobj_me=cell(length(dirst_me),1);
for cf=1:length(dirst_me)

       fn=fullfile(sn,dirst_me(cf).name);   
    mcobj_me{cf}=MetCon_ME(fn,ME_setting{:},'fm','IDEAL','Solver','IDEAL-modes');
end


% CSI data
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL-modes','fm','IDEAL'};
mcobj_csi=cell(length(dirst_csi),1);
for cf=1:length(dirst_csi)
       fn=fullfile(sn,dirst_csi(cf).name);
    mcobj_csi{cf}=MetCon_CSI(fn,CSI_setting{:});
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBJECT 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mask_s2=mcobj_me{end}.getMask(90);
% % as(squeeze(mask_s2(:,slcSel_s2,:)))
% [intake_time_s2,order_idx_s2]=sort(cellfun(@(x) x.getMinutesAfterIntake('08:13'),[mcobj_csi;mcobj_me],'UniformOutput',true));
% cellfun(@(x) x.WriteImages('',{'snr'}),[mcobj_csi;mcobj_me],'UniformOutput',false);
% 
% %
% vox_vol_s2=cellfun(@(x) prod(x.DMIPara.resolution_PSF(1:3)*1e2), [mcobj_csi;mcobj_me],'UniformOutput',true);
% vox_vol_s2=vox_vol_s2(order_idx_s2);
% resliced_s2=myspm_reslice(dir('Metcon_SNR_m01156_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_*.nii'),dir('Metcon*.nii'),'linear','r');
% 
% save('data_S2_modes.mat','vox_vol_s2','resliced_s2','mask_s2','intake_time_s2','type_label_s2')


%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBJECT 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_s3=mcobj_me{end}.getMask(90);
% as(squeeze(mask_s2(:,slcSel_s2,:)))
[intake_time_s3,order_idx_s3]=sort(cellfun(@(x) x.getMinutesAfterIntake('10:10'),[mcobj_csi;mcobj_me],'UniformOutput',true));
cellfun(@(x) x.WriteImages('',{'snr'}),[mcobj_csi;mcobj_me],'UniformOutput',false);

%
vox_vol_s3=cellfun(@(x) prod(x.DMIPara.resolution_PSF(1:3)*1e2), [mcobj_csi;mcobj_me],'UniformOutput',true);
vox_vol_s3=vox_vol_s3(order_idx_s3);
resliced_s3=myspm_reslice(dir('Metcon_SNR_m00102_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_IDEAL-modes.nii'),dir('Metcon*.nii'),'linear','r');

save('data_S3_modes.mat','vox_vol_s3','resliced_s3','mask_s3','intake_time_s3','type_label_s3')


%%  load and plot results

load('/ptmp/pvalsala/deuterium/paper/sub1_HOSJ_modes/data_S1_modes.mat')
load('/ptmp/pvalsala/deuterium/paper/sub2_H4_modes/data_S2_modes.mat')
load('/ptmp/pvalsala/deuterium/paper/sub3_modes/data_S3_modes.mat')

 resliced_s1_av= AverageByType(resliced_s1,type_label_s1);
  resliced_s2_av= AverageByType(resliced_s2,type_label_s2);
   resliced_s3_av= AverageByType(resliced_s3,type_label_s3);

% mask_s1=imerode(mask_s1,strel('sphere',2));
% %  mask_s2=mcobj_me{1}.getMask(80);
%  mask_s2=resliced_s2_av(:,:,:,1,1)>12;
%  mask_s2=imdilate(mask_s2,strel('sphere',1));
slcSel_s1=36;
slcSel_s2=34;
slcSel_s3=38;

mask_s3=resliced_s3_av(:,:,:,1,1)>10;
mask_s2=resliced_s2_av(:,:,:,1,1)>11;
mask_s1=resliced_s1_av(:,:,:,1,1)>10;
 % mask_s2(:,1:slcSel_s2-1,:)=false;
% mask_s2(:,slcSel_s2+1:end,:)=false;
% 
% mask_s1(:,1:slcSel_s1-1,:)=false;
% mask_s1(:,slcSel_s1+1:end,:)=false;


imPlot=cat(6,resliced_s1_av(9:42,slcSel_s1,:,:,:), ...
    resliced_s2_av(9:42,slcSel_s2,:,:,:),resliced_s3_av(9:42,slcSel_s3,:,:,:));
% dim1xdim2x seqType x subject x metabolite
imPlot=permute(imPlot ,[1 3 5 6 4 2]);
%scale with voxel volume (type x subjects )
vox_scal=ones(1,1,3);
vox_scal(1,1,:,1)=vox_vol_s1([4,6,7])./vox_vol_s1(4);
% vox_scal(1,1,:,2)=vox_scal(1,1,:,1);
imPlot=imPlot./vox_scal;


mask_all=squeeze(cat(4,mask_s1(9:42,slcSel_s1,:),mask_s2(9:42,slcSel_s2,:),mask_s3(9:42,slcSel_s3,:)));
mask_all=permute(repmat(mask_all,[1 1 1 3 4]),[1 2 4 3 5]);

imPlot=reshape(imPlot,size(imPlot,1),size(imPlot,2),[],size(imPlot,5));
mask_all=reshape(mask_all,size(mask_all,1),size(mask_all,2),[],size(mask_all,5));
%
% mask_s2=mcobj_me{1}.getMask(90);

figure(67),clf

tt=tiledlayout(5,12,"TileSpacing","tight","Padding","compact");

data_label={'CSI-FISP','CSI-bSSFP','ME-bSSFP'};


clim_all={[0 75],[0 35],[0 25],[0 15]};

for cMet=1:4
nexttile([2 3])
imagesc(createImMontage(imPlot(:,:,:,cMet),3))
hold on
if(cMet<2)
contour(createImMontage(mask_all(:,:,:,cMet),3),'LineWidth',0.5,'EdgeColor',[0.98 0.45 0.73],'EdgeAlpha',0.7)
end
xticks([size(imPlot,1)*0.5:size(imPlot,1):size(imPlot,1)*4])
xticklabels(data_label)
yticks([size(imPlot,2)*0.5:size(imPlot,2)+2:size(imPlot,2)*4])
yticklabels({'S1','S2','S3'})
title([metabolites(cMet).name, ' [SNR]'])
axis image
fontsize(gca,"scale",1.2)
colorbar
clim(clim_all{cMet})
colormap(gca,'turbo')
end



colors_label=lines(3);



yax_all={[-5 100],[-5 50],[-5 40]};




violin_format={'bw',1,'mc','rx','medc',[]};



data_label_s1=data_label(type_label_s1);
data_s1=reshape(resliced_s1_av,[],4,size(resliced_s1_av,5))./reshape(vox_scal,1,1,[]);
for cMet=1:3
nexttile([1 4])
violin2( squeeze(data_s1(mask_s1(:),cMet,:)),'facecolor',colors_label,violin_format{:})
xticks([])%  xticks(1:3),xticklabels(data_label),
grid on,grid minor
ylim(yax_all{cMet})
title([metabolites(cMet).name,' [SNR]'])
legend off
 if(cMet==1),ylabel('S1'),end
%create percentage
mv=mean( squeeze(data_s1(mask_s1(:),cMet,:)));
P40=prctile( squeeze(data_s1(mask_s1(:),cMet,:)),99,1)*1.2;
text(2.02,P40(2),sprintf('%+.0f%%',mv(2)/mv(1)*100-100),'HorizontalAlignment','left','FontWeight','bold')
text(3.02,P40(3),sprintf('%+.0f%%',mv(3)/mv(1)*100-100),'HorizontalAlignment','left','FontWeight','bold')

end



data_label_s2=data_label(type_label_s2);
data_s2=reshape(resliced_s2_av,[],4,size(resliced_s2_av,5))./reshape(vox_scal,1,1,[]);
for cMet=1:3
nexttile([1 4])
vh=violin2( squeeze(data_s2(mask_s2(:),cMet,:)),'facecolor',colors_label,violin_format{:});
%  xticklabels(data_label)
grid on,grid minor
ylim(yax_all{cMet})
legend off
xticks([])% xticks(1:3),xticklabels(data_label)
 if(cMet==1),ylabel('S2'),end
 %create percentage
mv=mean( squeeze(data_s2(mask_s2(:),cMet,:)));
P40=prctile( squeeze(data_s2(mask_s2(:),cMet,:)),99,1)*1.2;
text(2.02,P40(2),sprintf('%+.0f%%',mv(2)/mv(1)*100-100),'HorizontalAlignment','left','FontWeight','bold')
text(3.02,P40(3),sprintf('%+.0f%%',mv(3)/mv(1)*100-100),'HorizontalAlignment','left','FontWeight','bold')
end


data_label_s3=data_label(type_label_s3);
data_s3=reshape(resliced_s3_av,[],4,size(resliced_s3_av,5))./reshape(vox_scal,1,1,[]);
for cMet=1:3
nexttile([1 4])
vh=violin2( squeeze(data_s3(mask_s3(:),cMet,:)),'facecolor',colors_label,violin_format{:});
%  xticklabels(data_label)
grid on,grid minor
ylim(yax_all{cMet})
legend off
 xticks(1:3),xticklabels(data_label)
 if(cMet==1),ylabel('S3'),end
 %create percentage
mv=mean( squeeze(data_s3(mask_s3(:),cMet,:)));
P40=prctile( squeeze(data_s3(mask_s3(:),cMet,:)),99,1)*1.2;
text(2.02,P40(2),sprintf('%+.0f%%',mv(2)/mv(1)*100-100),'HorizontalAlignment','left','FontWeight','bold')
text(3.02,P40(3),sprintf('%+.0f%%',mv(3)/mv(1)*100-100),'HorizontalAlignment','left','FontWeight','bold')
end

% legend(gca,vh([1 3 2]),data_label,'Location','northwest')

fontsize(gcf,"scale",1.3)
set(gcf,'color','w','Position', [148 76 1570 1232])



%% 
function Vol_av= AverageByType(Vol,type_idx)

Vol_av=cat(5,mean(Vol(:,:,:,:,type_idx==1),5),mean(Vol(:,:,:,:,type_idx==2),5),mean(Vol(:,:,:,:,type_idx==3),5));
end