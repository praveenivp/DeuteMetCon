%%
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/MATLAB'))

addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))

metabolites=getMetaboliteStruct('invivo');
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

%% invivo 3
CSIdataset='/ptmp/pvalsala/deuterium/HOSJ-D6P2/TWIX/allData#S94Tuebingen#F16603#M997#D230924#T093336#rpcsi_fid_Stan_res156_moreopti.dat';
mcobj=MetCon_CSI(CSIdataset,CSI_setting{:});
mask=mcobj.getMask(60);
mask=imerode(mask,strel("sphere",5));
mask=mcobj.Experimental.relativeNorm<0.7;
mask=imdilate(mask,strel("sphere",1));
% T2Star=[21.0052 13.9556 22.3202 31.8583]*1e-3; %s

%% invivo 4
CSIdataset='/ptmp/pvalsala/deuterium/H4DK-64GH/TWIX/allData#S94Tuebingen#F17843#M1178#D101024#T103246#rpcsi_fid_Stan_res156_moreopti.dat';
mcobj=MetCon_CSI(CSIdataset,CSI_setting{:});
mask=mcobj.getMask(60);
mask=imerode(mask,strel("sphere",5));
mask=mcobj.Experimental.relativeNorm<0.7;
mask=imdilate(mask,strel("sphere",1));

%% invivo 5
CSIdataset='/ptmp/pvalsala/deuterium/I3BL-CJ5O/TWIX/allData#S94Tuebingen#F3446#M114#D281124#T120217#rpcsi_fid_Stan_res156_moreopti.dat';
mcobj=MetCon_CSI(CSIdataset,CSI_setting{:});
mask=mcobj.getMask(60);
mask=imerode(mask,strel("sphere",5));
mask=mcobj.Experimental.relativeNorm<0.7;
mask=imdilate(mask,strel("sphere",1));
% T2* [ms]  19.4120 13.9282 22.3839 20.7921
% T2* std  [ms]  9.6127 11.3821 48.6879 67.1582
% cs [Hz]  3.7134 -51.1761 -139.8969 -203.9604
% cs [ppm]  4.7000 3.8054 2.3595 1.3154

%% phantom
metabolites=getMetaboliteStruct('phantom');
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','wsvd','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','AMARES'};
CSIdataset='/ptmp/pvalsala/deuterium/20240813_spectral/meas_MID00857_FID14528_rpcsi_fid_Stan_res15_6_optimal.dat';
mcobj=MetCon_CSI(CSIdataset,CSI_setting{:});
% as metabolites are highly localized, we use different masking
FWHM_hz=mcobj.Experimental.linewidth;
  T2Star=1./(pi*FWHM_hz);  %s
  cs=reshape(mcobj.Experimental.chemicalshift,[],length(metabolites));

mask4=(mcobj.getNormalized)>20;
 as(T2Star.*mask4*1e3)
mask4=reshape(mask4,[],4);
T2Star1=reshape(T2Star,[],size(T2Star,4));
[median(T2Star1(mask4(:,1),1)) median(T2Star1(mask4(:,2),2)) median(T2Star1(mask4(:,3),3)) median(T2Star1(mask4(:,4),4))]*1e3
[std(T2Star1(mask4(:,1),1)) std(T2Star1(mask4(:,2),2)) std(T2Star1(mask4(:,3),3)) std(T2Star1(mask4(:,4),4))]*1e3

 freq_shift=[median(cs(mask4(:,1),1)) median(cs(mask4(:,2),2)) median(cs(mask4(:,3),3)) median(cs(mask4(:,4),4))];
  freq_shift_std=[std(cs(mask4(:,1),1)) std(cs(mask4(:,2),2)) std(cs(mask4(:,3),3)) std(cs(mask4(:,4),4))];
  freq_shift_ppm=(freq_shift-freq_shift(1))./(mcobj.DMIPara.ImagingFrequency_MHz)+4.7;
freq_shift_std=freq_shift_std./(mcobj.DMIPara.ImagingFrequency_MHz)
% T2Star=[74.7374   22.0817   53.4715   75.7936]*1e-3; % s
% T2Star_std=[38.9049    2.8619   14.3097   28.7368]*1e-3; % s
% depening on the vial. Pure vials have double theses numbers except
% water
% cs [Hz]  -3.7858  -68.1933 -157.3928 -216.7620
% cs [ppm]   4.7000    3.6503    2.1966    1.2290
%cs[ppm] std 0.0700    0.0531    0.0654    0.0796

%% plot and get statistics 
FWHM_hz=mcobj.Experimental.linewidth;
 T2Star=1./(pi*FWHM_hz);  %s
 as(T2Star.*mask*1e3) %ms
T2Star1=reshape(T2Star,[],size(T2Star,4));
T2Star1=T2Star1(mask(:),:);
cs=reshape(mcobj.Experimental.chemicalshift,[],length(metabolites));
cs=cs(mask(:),:);
 cs=cs;%*mcobj.DMIPara.ImagingFrequency_MHz;
  freq_shift=median(cs);
  freq_shift_ppm=(freq_shift-freq_shift(1))./(mcobj.DMIPara.ImagingFrequency_MHz)+4.7;

clc
  fprintf("T2* [ms]  %.4f %.4f %.4f %.4f\n",median(T2Star1*1e3)) 
  fprintf("T2* std  [ms]  %.4f %.4f %.4f %.4f\n",std(T2Star1*1e3)) 
  fprintf("cs [Hz]  %.4f %.4f %.4f %.4f\n",freq_shift) 
  fprintf("cs [ppm]  %.4f %.4f %.4f %.4f\n",freq_shift_ppm) 
  
  figure,
  subplot(121),violin(T2Star1*1e3);ylabel('T2* [ms]'),ylim([0 40]),grid on
  subplot(122),violin(cs);ylabel('chemical shift [Hz]'),grid on