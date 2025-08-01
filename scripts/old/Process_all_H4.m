MeasPath='/ptmp/pvalsala/deuterium/H4DK-64GH';
sn=fullfile(MeasPath,'TWIX');

dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));
pn=fullfile(MeasPath,'proc',sprintf('proc_all_%s',datetime('today','Format','yyyyMMMdd')));
mkdir(pn);
dirst_me=dir(fullfile(sn,"*rh_trufi*.dat"));
dirst_me=dirst_me([1 3 4]);
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
metabolites=getMetaboliteStruct('invivo');


%% process noise anf get fieldmap
fn_noise=dir(fullfile(sn,"*rh_trufi*noise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true};
% fm_meas_Hz=mcobj_ideal.Experimental.fm_est*(-2*pi)/(6.536 /42.567);
%% process all ME- data
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
%write images
cd(pn)
[intake_time_s2,order_idx_s2]=sort(cellfun(@(x) x.getMinutesAfterIntake('08:13'),[mcobj_csi;mcobj_me],'UniformOutput',true));
cellfun(@(x) x.WriteImages('',{'snr','mM'}),[mcobj_csi;mcobj_me],'UniformOutput',false);

%% for SNR
vox_vol_s2=cellfun(@(x) prod(x.DMIPara.resolution_PSF(1:3)*1e2), [mcobj_csi;mcobj_me],'UniformOutput',true);
vox_vol_s2=vox_vol_s2(order_idx_s2);
resliced_s2=myspm_reslice(dir('Metcon_SNR_m01156_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_*.nii'),dir('Metcon_SNR*.nii'),'linear','r');
mask_s2=mcobj_me{end}.getMask(90);
% as(squeeze(mask_s2(:,slcSel_s2,:)))
type_label_s2=[1     3     2     3     1     2     3];
save('data_S2_modes.mat','vox_vol_s2','resliced_s2','mask_s2','intake_time_s2','type_label_s2')

%% 


% load('/ptmp/pvalsala/deuterium/paper/sub1_HOSJ/data_S1.mat')
% load('/ptmp/pvalsala/deuterium/paper/sub2_H4/data_S2.mat')


load('/ptmp/pvalsala/deuterium/paper/sub1_HOSJ_modes/data_S1_modes.mat')
load('/ptmp/pvalsala/deuterium/paper/sub2_H4_modes/data_S2_modes.mat')


mask_s1=imerode(mask_s1,strel('sphere',2));
% mask_s2=mcobj_me{end}.getMask(93);
mask_s2=imerode(mask_s2,strel('sphere',1));
slcSel_s1=36;
slcSel_s2=32;

% mask_s2(:,1:slcSel_s2-1,:)=false;
% mask_s2(:,slcSel_s2+1:end,:)=false;
% 
% mask_s1(:,1:slcSel_s1-1,:)=false;
% mask_s1(:,slcSel_s1+1:end,:)=false;


imPlot=(cat(6,resliced_s1(9:40,slcSel_s1,:,:,[4 6 7]),resliced_s2(9:40,slcSel_s2,:,:,[5,6,7])));
% dim1xdim2x seqType x subject x metabolite
imPlot=permute(imPlot ,[1 3 5 6 4 2]);
%scale with voxel volume (type x subjects )
vox_scal=ones(1,1,3,2);
vox_scal(1,1,:,1)=vox_vol_s1([4,6,7])./vox_vol_s1(4);
vox_scal(1,1,:,2)=vox_vol_s2([5,6,7])./vox_vol_s2(5);
imPlot=imPlot./vox_scal;


mask_all=squeeze(cat(4,mask_s1(9:40,slcSel_s1,:),mask_s2(9:40,slcSel_s2,:)));
mask_all=permute(repmat(mask_all,[1 1 1 3 4]),[1 2 4 3 5]);

imPlot=reshape(imPlot,size(imPlot,1),size(imPlot,2),[],size(imPlot,5));
mask_all=reshape(mask_all,size(mask_all,1),size(mask_all,2),[],size(mask_all,5));
%
% mask_s2=mcobj_me{1}.getMask(90);

figure(67),clf

tt=tiledlayout(3,12,"TileSpacing","compact","Padding","compact");

data_label={'CSI-FISP','CSI-bSSFP','ME-bSSFP'};


clim_all={[0 75],[0 30],[0 30],[0 15]};

for cMet=1:4
nexttile([1 3])
imagesc(createImMontage(imPlot(:,:,:,cMet),3))
hold on
contour(createImMontage(mask_all(:,:,:,cMet),3),'LineWidth',0.5,'EdgeColor','r','EdgeAlpha',0.1)
xticks([size(imPlot,1)*0.5:size(imPlot,1):size(imPlot,1)*4])
xticklabels(data_label)
yticks([size(imPlot,2)*0.5:size(imPlot,2):size(imPlot,2)*4])
yticklabels({'S1','S2'})
title(metabolites(cMet).name)
axis image
fontsize(gca,"scale",1.2)
colorbar
clim(clim_all{cMet})
colormap(gca,'parula')
end












colors_label=lines(3);



yax_all={[-5 100],[-5 50],[-5 40]};




violin_format={'bw',1,'mc','rx','medc',[]};



data_label_s1=data_label(type_label_s1);
data_s1=reshape(resliced_s1,[],4,size(resliced_s1,5))./reshape(vox_vol_s1./vox_vol_s1(1),1,1,[]);
for cMet=1:3
nexttile([1 4])
violin2( squeeze(data_s1(mask_s1(:),cMet,:)),'facecolor',colors_label(type_label_s1,:),violin_format{:})
xticklabels(sort(intake_time_s1))
grid on,grid minor
ylim(yax_all{cMet})
title(metabolites(cMet).name)
legend off
end



data_label_s2=data_label(type_label_s2);
data_s2=reshape(resliced_s2,[],4,size(resliced_s2,5))./reshape(vox_vol_s2./vox_vol_s2(1),1,1,[]);
for cMet=1:3
nexttile([1 4])
vh=violin2( squeeze(data_s2(mask_s2(:),cMet,:)),'facecolor',colors_label(type_label_s2,:),violin_format{:});
xticklabels(sort(intake_time_s2))
grid on,grid minor
ylim(yax_all{cMet})
legend off
end

legend(gca,vh([1 3 2]),data_label,'Location','northwest')

fontsize(gcf,"scale",1.3)
set(gcf,'color','w','Position',[198 56 1500 900])

% [hleg,hico]=legend(data_label{1},'del',data_label{2},'del',data_label{3},'mean')
% % Delete objects associated with del
% istxt = strcmp(get(hico, 'type'), 'text');
% hicot = hico(istxt);
% hicol = hico(~istxt);
% delete(hicot(ismember(get(hicot, 'String'), {'del'})));
% delete(hicol(ismember(get(hicol, 'Tag'),    {'del'})));

% save('data_S2.mat','vox_vol_s2','resliced_s2','mask_s2','intake_time_s2','type_label_s2')