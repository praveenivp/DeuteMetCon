%%
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/MATLAB'))

addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))

metabolites=getMetaboliteStruct('invivo1');
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','wsvd','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','AMARES'};

%% invivo 1
CSIdataset='/ptmp/pvalsala/deuterium/EAZW-GUMK/TWIX/allData#S94Tuebingen#F5070#M412#D220124#T093615#rpcsi_fid_Stan_Res25_156mm.dat';
mcobj=MetCon_CSI(CSIdataset,CSI_setting{:});
mask=mcobj.getMask(60);
mask=imerode(mask,strel("sphere",5));
% T2* [ms]  24.3812 17.0160 25.2516 23.8431
% cs [Hz]  -0.6647 -55.7530 -144.9229 -205.0421
% cs [ppm]  4.7000 3.8022 2.3490 1.3692

%% invivo 2

CSIdataset='/ptmp/pvalsala/deuterium/DA77-F3UY/TWIX/allData#S94Tuebingen#F10215#M696#D220524#T092924#rpcsi_fid_Stan_res156_moreopti.dat';
mcobj=MetCon_CSI(CSIdataset,CSI_setting{:});
mask=mcobj.getMask(60);
mask=imerode(mask,strel("sphere",5));
% T2* [ms]  20.0680 12.7326 22.9249 25.2875
% cs [Hz]  9.1671 -46.2185 -134.4634 -195.6384
% cs [ppm]  4.7000 3.7974 2.3592 1.3622


%% invivo 2.1(longer readout)

CSIdataset='/ptmp/pvalsala/deuterium/DA77-F3UY/TWIX/allData#S94Tuebingen#F10214#M695#D220524#T091915#rpcsi_fid_Stan_res156_optimal.dat';
mcobj=MetCon_CSI(CSIdataset,CSI_setting{:});
mask=mcobj.getMask(60);
mask=imerode(mask,strel("sphere",5));

% T2* [ms]  20.5969 14.2529 22.2796 25.2529
% cs [Hz]  9.0229 -46.9185 -135.4832 -198.4703
% cs [ppm]  4.7000 3.7883 2.3449 1.3184
%% phantom
CSIdataset='/ptmp/pvalsala/deuterium/20240813_spectral/meas_MID00857_FID14528_rpcsi_fid_Stan_res15_6_optimal.dat';

% cs [Hz]  -3.4080 -63.5005 -150.2146 -205.1276
% cs [ppm]  4.7000 3.7206 2.3074 1.4125

mcobj=MetCon_CSI(CSIdataset,CSI_setting{:});
mask=mcobj.getMask(60);
mask=imerode(mask,strel("sphere",5));

% as metabolites are highly localized, we use different masking
mask4=(mcobj.Metcon)>5;
[median(T2Star1(mask4(:,1),1)) median(T2Star1(mask4(:,2),2)) median(T2Star1(mask4(:,3),3)) median(T2Star1(mask4(:,4),4))]*1e3
% T2* [ms]   77.9295   18.8752   14.0095   15.9714 % lot of variablility
% depening on the vial. Pure vials have double theses numbers except
% water!0

%% plot and get statistics 
FWHM_hz=mcobj.Experimental.linewidth;
 T2Star=1./(pi*FWHM_hz);  %s
 as(T2Star.*mask*1e3) %ms
T2Star1=reshape(T2Star,[],size(T2Star,4));
T2Star1=T2Star1(mask(:),:);
cs=reshape(mcobj.Experimental.chemicalshift,[],length(metabolites));
cs=cs(mask(:),:);
 cs=cs*mcobj.DMIPara.ImagingFrequency_MHz;
  freq_shift=median(cs);
  freq_shift_ppm=(freq_shift-freq_shift(1))./(mcobj.DMIPara.ImagingFrequency_MHz)+4.7;

clc
  fprintf("T2* [ms]  %.4f %.4f %.4f %.4f\n",median(T2Star1*1e3)) 
  fprintf("cs [Hz]  %.4f %.4f %.4f %.4f\n",freq_shift) 
  fprintf("cs [ppm]  %.4f %.4f %.4f %.4f\n",freq_shift_ppm) 
  
  figure,
  subplot(121),violin(T2Star1*1e3);ylabel('T2* [ms]'),ylim([0 40]),grid on
  subplot(122),violin(cs);ylabel('chemical shift [Hz]'),grid on