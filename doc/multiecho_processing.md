# Multi-echo data processing

# Calculate noise decorrelation matrix
As our trufi does noise scans only with parallel imaging. We acquire noise scans (imaging with 0 flip angle) with same BW and resolution. 
```
twix_noise=mapVBVD(fullfile(sn,'meas_MID00280_FID08611_pvrh_trufi_5E_noise_12P5mm.dat'),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);
```
As noise decorrelation is a crucial step for reconstruction in SNR units, take care of the following:
- remove oversampling to remove anti-alias filter transition band (factor 0.79?)
- `Metcon_bSSFP` work with oversampling removed data (factor sqrt(2))
- readout bandwidth is same between imaging and nosie (factor sqrt(dwell time ratio))

# calculate fieldmaps

Reconstruct proton ME-GRE data amd calcualte field maps. For this, we need `recoVBVD` for cartesian PI reconstrucion. 

```
 fmobj=B0map(fullfile(sn,dirst(i).name)...
         ,'UnwrapMode','SpatialUnwrap','doRegularization',true);
 save(fullfile(pn2,sprintf('fmobj_measuid%d.mat',fmobj.reco_obj.twix.hdr.Config.MeasUID)),'fmobj');
 fmobj.reco_obj.WriteNIFTI();
 V=double(fmobj.mask.*fmobj.Fmap(:,:,:,1));
V=flip(permute(V,[2 1 3 4 ]),2); %9.4 T transverter
 MyNIFTIWrite(V,fmobj.reco_obj.twix,sprintf('M%d_B0map_rad_s.nii',fmobj.reco_obj.twix.hdr.Config.MeasUID),'1H fielmap rad/s');
```

# registering fieldmaps


# check resonances

# set T1/T2 times

## new phantom relaxometry values
```
%new phantom: /ptmp/pvalsala/deuterium/20240102_new2Hphantom
met_name={'water','glucose','Glx','lactate'};
freq_shift=[ -65 -154 -214]; %Hz
T1=[444.30, 69.53, 157.65, 244.38]*1e-3; %s
T2=[272.75,52.09,114.31,245.17]*1e-3; %s
```
for invivo we have all values except lactate

```
%invivo: /ptmp/pvalsala/deuterium/EAZW-GUMK/proc
met_name={'water','glucose','Glx','lactate'};
freq_shift=[3.9 -58.7  -152 -200]; % Hz
T1=[342.04 56.62  151.55 165.02]*1e-3;%s
T2=[53.93 45.96 80.10 96.47]*1e-3;%s
```
```
%invivo: /ptmp/pvalsala/deuterium/DA77-F3UY
met_name={'water','glucose','Glx','lactate'};
freq_shift=[3.9 -58.7  -152 -200]; %Hz
T1=[343.58 94.21 169.32 190.64]*1e-3;%s
T2=[287 65.88 124.5 180]*1e-3;%s
```

```
clear metabolites;
for i=1:length(freq_shift)
    metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift(i),'name',met_name{i});
end
```

