MeasPath='/ptmp/pvalsala/deuterium/H4DK-64GH';
sn=fullfile(MeasPath,'TWIX');

dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));
% dirst_csi=dirst_csi([1 3 4]);
% pn=fullfile(MeasPath,'proc',sprintf('proc_all_%s',datetime('today','Format','yyyyMMMdd')));
% mkdir(pn);
dirst_me=dir(fullfile(sn,"*rh_trufi*.dat"));
dirst_me=dirst_me([1 3 4]);
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))

metabolites=getMetaboliteStruct('invivo3');


%% process noise anf get fieldmap
fn_noise=dir(fullfile(sn,"*rh_trufi*noise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

% fn_fm=(fullfile(sn,dirst_me(3).name));
ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true};
% mcobj_ideal=MetCon_bSSFP(fn_fm,'Solver','IDEAL',ME_setting{:},'mask',80);
% fm_meas_Hz=mcobj_ideal.Experimental.fm_est*(-2*pi)/(6.536 /42.567);
%% process all ME- data
% ME_setting=[ME_setting,{'fm',fm_meas_Hz,'Solver','pinv'}];
mcobj_me=cell(length(dirst_me),1);
for cf=1:length(dirst_me)

       fn=fullfile(sn,dirst_me(cf).name);   
    mcobj_me{cf}=MetCon_bSSFP(fn,ME_setting{:},'fm','IDEAL','Solver','pinv');
end


%% CSI data
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','pinv','fm','IDEAL'};


mcobj_csi=cell(length(dirst_csi),1);
for cf=1:length(dirst_csi)
       fn=fullfile(sn,dirst_csi(cf).name);
    mcobj_csi{cf}=MetCon_CSI(fn,CSI_setting{:});
end

%%
pn='/ptmp/pvalsala/deuterium/paper/sub2_H4';
cd(pn)
[intake_time_s2,order_idx_s2]=sort(cellfun(@(x) x.getMinutesAfterIntake('08:13'),[mcobj_csi;mcobj_me],'UniformOutput',true));
cellfun(@(x) x.WriteImages('',{'snr'}),[mcobj_csi;mcobj_me],'UniformOutput',false);

%%
vox_vol_s2=cellfun(@(x) prod(x.DMIPara.resolution_PSF(1:3)*1e2), [mcobj_csi;mcobj_me],'UniformOutput',true);
vox_vol_s2=vox_vol_s2(order_idx_s2);

%%
resliced_s2=myspm_reslice(dir('Metcon_SNR_m01156_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_pinv.nii'),dir('Metcon*.nii'),'linear','r');

%% 
slcSel_s1=32;
slcSel_s2=40;

imPlot=(cat(6,resliced_s1(9:40,slcSel_s2,:,:,[4 6 5]),resliced_s2(9:40,slcSel_s2,:,:,[6,5,4])));
% dim1xdim2x seqType x subject x metabolite
imPlot=permute(imPlot ,[1 3 5 6 4 2]);
imPlot=reshape(imPlot,size(imPlot,1),size(imPlot,2),[],size(imPlot,5));
%
mask_s2=mcobj_me{1}.getMask(90);

figure(67),clf

tt=tiledlayout(3,12,"TileSpacing","compact","Padding","compact");

data_label={'CSI-FISP','CSI-bSSFP','ME-bSSFP'};




for cMet=1:4
nexttile([1 3])
imagesc(createImMontage(imPlot(:,:,:,cMet),3))
xticks([size(imPlot,1)*0.5:size(imPlot,1):size(imPlot,1)*4])
xticklabels(data_label)
yticks([size(imPlot,2)*0.5:size(imPlot,2):size(imPlot,2)*4])
yticklabels({'S1','S2'})
title(metabolites(cMet).name)
axis image
end











colors_label=lines(3);



yax_all={[-5 100],[-5 50],[-5 40]};




violin_format={'bw',2,'mc','rx','medc',[]};



load('/ptmp/pvalsala/deuterium/paper/sub1_HOSJ/data_S1.mat')
data_label_s1=data_label(type_label_s1);
data_s1=reshape(resliced_s1,[],4,size(resliced_s1,5))./reshape(vox_vol_s1./vox_vol_s1(1),1,1,[]);
for cMet=1:3
nexttile([1 4])
violin2( squeeze(data_s1(mask_s1(:),cMet,:)),'facecolor',colors_label(type_label_s1,:),violin_format{:})
xticklabels(sort(intake_time_s1))
grid on,grid minor
ylim(yax_all{cMet})
title(metabolites(cMet).name)
end


load('/ptmp/pvalsala/deuterium/paper/sub2_H4/data_S2.mat')
data_label_s2=data_label(type_label_s2);
data_s2=reshape(resliced_s2,[],4,size(resliced_s2,5))./reshape(vox_vol_s2./vox_vol_s2(1),1,1,[]);
for cMet=1:3
nexttile([1 4])
violin2( squeeze(data_s2(mask_s2(:),cMet,:)),'facecolor',colors_label(type_label_s2,:),violin_format{:})
xticklabels(sort(intake_time_s2))
grid on,grid minor
ylim(yax_all{cMet})
end

% save('data_S2.mat','vox_vol_s2','resliced_s2','mask_s2','intake_time_s2','type_label_s2')