%%
%% input files
MeasPath='/ptmp/pvalsala/deuterium/20241021_phantomtest';
sn=fullfile(MeasPath,'TWIX');

dirst_csi=dir(fullfile(sn,"*rpcsi_fid*.dat"));
dirst_csi=dirst_csi(2);
dirst_csi_ssfp=dir(fullfile(sn,"*rpcsi_ssfp*.dat"));
% dirst_csi_ssfp(2:3);
dirst_me=dir(fullfile(sn,"*pvrh_trufi_5E_*.dat"));
dirst_me=dirst_me(2);

pn=fullfile(MeasPath,sprintf('proc/csi_GRE_%s',datetime('today','Format','yyyyMMMdd')));


addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))


metabolites=getMetaboliteStruct('phantom',0);

%% process noise anf get fieldmap
fn_noise=dir(fullfile(sn,"*oise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true,'fm','IDEAL','Solver','IDEAL-modes'};

CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL','fm',[]};

CSI_setting_ssfp={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL-modes','fm','IDEAL'};


%% process all CSI
mcobj_csi=cell(length(dirst_csi),1);
for cf=1:length(dirst_csi)
       fn=fullfile(sn,dirst_csi(cf).name);
    mcobj_csi{cf}=MetCon_CSI(fn,CSI_setting{:});
end
mcobj_csi_ssfp=cell(length(dirst_csi_ssfp),1);
for cf=1:length(dirst_csi_ssfp)
       fn=fullfile(sn,dirst_csi_ssfp(cf).name);
    mcobj_csi_ssfp{cf}=MetCon_CSI(fn,CSI_setting_ssfp{:});
end
% process all ME data
mcobj_me=cell(length(dirst_me),1);
for cf=1:length(dirst_me)
       fn=fullfile(sn,dirst_me(cf).name);
    mcobj_me{cf}=MetCon_ME(fn,ME_setting{:});
end

%%
[fw64_fispcsi,PSF_fispcsi,W_fispcsi]=getPSF_CSI(mcobj_csi{1}.twix);
[fw64_csibssfp,PSF_csibssfp,W_csibssfp]=getPSF_CSI(mcobj_csi_ssfp{1}.twix);
[fw64_mebssfp,PSF_mebssfp,W_mebssfp]=getPSF_CSI(mcobj_me{1}.twix);

%%
figure,
hold on
plot(PSF_fispcsi.FOV_E1,abs(PSF_fispcsi.PSF_E1)./max(abs(PSF_fispcsi.PSF_E1)),'LineWidth',1.5)
plot(PSF_csibssfp.FOV_E1,abs(PSF_csibssfp.PSF_E1)./max(abs(PSF_csibssfp.PSF_E1)),'LineWidth',1.5)
plot(PSF_mebssfp.FOV_E1,abs(PSF_mebssfp.PSF_E1)./max(abs(PSF_mebssfp.PSF_E1)),'LineWidth',1.5)
% plot(PSF_mebssfp.FOV_E2,abs(PSF_mebssfp.PSF_E2)./max(abs(PSF_mebssfp.PSF_E2)))
% plot(PSF_mebssfp.FOV_E3,abs(PSF_mebssfp.PSF_E3)./sum(abs(PSF_mebssfp.PSF_E3)))
grid minor,grid on
legend('CSI-FISP','CSI-bSSFP','ME-bSSFP')
xlim([-1 1]*50)

title('PSF')
xlabel('distance [mm]'),ylabel('amplitude [a.u]')

set(gcf,'Color','w')
fontsize(gcf,'scale',1.3)
box on
