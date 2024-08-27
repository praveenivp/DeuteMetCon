%% input files
sn='/ptmp/pvalsala/deuterium/20240813_spectral';
pn='/ptmp/pvalsala/deuterium/20240813_spectral/proc/SNR';
dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));

dirst_me=dir(fullfile(sn,"*pvrh_trufi_5E_*.dat"));


addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))


metabolites=getMetaboliteStruct('phantom');

%% process noise anf get fieldmap
fn_noise=dir(fullfile(sn,"*Noise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

fn_fm=(fullfile(sn,dirst_me(1).name));
ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true};
mcobj_ideal=MetCon_bSSFP(fn_fm,'Solver','IDEAL',ME_setting{:});
fm_meas_Hz=mcobj_ideal.Experimental.fm_est*(-2*pi)/(6.536 /42.567);
%% final setttings
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','wsvd','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','AMARES'};
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

%% reslicing and registration
% cellfun(@(x) x.WriteImages,mcobj_csi,'UniformOutput',false)
% cellfun(@(x) x.WriteImages,mcobj_me,'UniformOutput',false)
dirst_nii=dir(fullfile(pn,"Metcon*.nii"));
resliced_vol=myspm_reslice(dirst_nii(1).name,dirst_nii,'linear','r');
[realigned_vol]=realign_vol(resliced_vol);

%% plotting
figure,
tt=tiledlayout(4,1);

slcSel=32;
for ii=1:size(realigned_vol,4)
    nexttile()
    realigned_vol2=realigned_vol./reshape(vox_vol,1,1,1,1,[]);
im_plot=reshape(realigned_vol2(:,slcSel,:,ii,:),size(realigned_vol,1),[]);

imagesc(im_plot),colormap("jet"),colorbar,axis image
title(metabolites(ii).name)
end

vox_vol=[cellfun(@(x)prod(x.DMIPara.resolution_PSF(1:3)*1e2),mcobj_me(1:2),'UniformOutput',true),...
cellfun(@(x)prod(x.DMIPara.resolution_PSF(1:3)*1e2),mcobj_csi(1:5),'UniformOutput',true)];

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
