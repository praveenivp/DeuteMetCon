MeasPath='/ptmp/pvalsala/deuterium/HOSJ-D6P2';
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
    mcobj_me{cf}=MetCon_bSSFP(fn,ME_setting{:},'fm','IDEAL','Solver','IDEAL-modes');
end


%% CSI data
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL-modes','fm','IDEAL'};


mcobj_csi=cell(length(dirst_csi),1);
for cf=1:length(dirst_csi)
       fn=fullfile(sn,dirst_csi(cf).name);
    mcobj_csi{cf}=MetCon_CSI(fn,CSI_setting{:});
end

%%
pn='/ptmp/pvalsala/deuterium/paper/sub1_HOSJ_modes';
mkdir(pn)
cd(pn)
[intake_time_s1,order_idx]=sort(cellfun(@(x) x.getMinutesAfterIntake('08:13'),[mcobj_csi;mcobj_me],'UniformOutput',true));
cellfun(@(x) x.WriteImages('',{'snr'}),[mcobj_csi;mcobj_me],'UniformOutput',false);

%%
vox_vol_s1=cellfun(@(x) prod(x.DMIPara.resolution_PSF(1:3)*1e2), [mcobj_csi;mcobj_me],'UniformOutput',true);
vox_vol_s1=vox_vol_s1(order_idx);

%%
resliced_s1=myspm_reslice(dir('Metcon_SNR_m01011_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180_*.nii'),dir('Metcon*.nii'),'linear','r');

%% 
mask_s1=mcobj_me{1}.getMask(92);

figure(63),clf


tt=tiledlayout(3,3,"TileSpacing","compact","Padding","compact");

data_label={'CSI-FISP','CSI-bSSFP','ME-bSSFP'};
colors_label=lines(3);
type_label_s1=[1 2 3 1 3 2 3];
data_label_s1=data_label(type_label_s1);

data_s1=reshape(resliced_s1,[],4,size(resliced_s1,5))./reshape(vox_vol_s1./vox_vol_s1(1),1,1,[]);
for cMet=1:3
nexttile()
violin( squeeze(data_s1(mask_s1(:),cMet,:)),'facecolor',colors_label(type_label_s1,:))
xticklabels(sort(intake_time_s1))
grid on
end

save('data_S1_modes.mat','vox_vol_s1','resliced_s1','mask_s1','intake_time_s1','type_label_s1')