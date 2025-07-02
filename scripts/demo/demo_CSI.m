%%
CSI_setting1={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','AMARES'};
CSI_setting2={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL'};

fn='/ptmp/pvalsala/deuterium/20240813_spectral/meas_MID00855_FID14526_rpcsi_ssfp_Stan25_15_6mm.dat';

mcobj1=MetCon_CSI(fn,CSI_setting1{:});
mcobj2=MetCon_CSI(fn,CSI_setting2{:});