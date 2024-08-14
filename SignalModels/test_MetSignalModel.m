%% test case : bSSFP T1==T2 and FA=90, Steady state signal =0.5
TR=10e-3;
FA=90;
TE=TR/2;
metabolites=struct("T1_s",0.4,"T2_s",0.4,"freq_shift_Hz",0,"name","water","T2star_s",0.45);

 [Msig_all]=MetSignalModel   (metabolites,TE,pi,TR,0,deg2rad(FA),'bSSFP',1);
 assert(abs(squeeze(Msig_all))-0.5 < 1e-3,"steady state signal should be 0.5 for T1/T2=1 and FA=90")
  [Msig_all]=MetSignalModel   (metabolites,TE,pi,TR,0,deg2rad(FA),'bSSFP-peters',1);
 assert(abs(squeeze(Msig_all))-0.5 < 1e-3,"steady state signal should be 0.5 for T1/T2=1 and FA=90")

 %% for FA=0 , signal should be zero
 FA=0;
 [Msig_all]=MetSignalModel   (metabolites,TE,pi,TR,0,deg2rad(FA),'bSSFP',1);
 assert(abs(squeeze(Msig_all)) < eps,"steady state signal should be 0 for  FA=0")
  [Msig_all]=MetSignalModel   (metabolites,TE,pi,TR,0,deg2rad(FA),'bSSFP-peters',1);
 assert(abs(squeeze(Msig_all)) < eps,"steady state signal should be 0 for  FA=0")


 %% test case: GRE : TR >>>T1 and FA=90, Steady state signal =1
TR=20;
FA=90;
TE=0e-3;
metabolites=struct("T1_s",.4,"T2_s",.4,"freq_shift_Hz",0,"name","water","T2star_s",20);

 [Msig_all]=MetSignalModel   (metabolites,TE,pi,TR,0,deg2rad(FA),'GRE',1e-4);
 assert(abs(abs(squeeze(Msig_all))-1) <1e-4,"steady state signal should be 1 for TR>>T1 and FA=90")
  [Msig_all]=MetSignalModel   (metabolites,TE,pi,TR,0,deg2rad(FA),'GRE-peters',1e-5);
 assert(abs(abs(squeeze(Msig_all))-1) < 1e-4,"steady state signal should be 1 for TR>>T1 and FA=90")

 %% ernst angle should give maximum signal
 clc
TR=100e-3;
FA=0:5:90;
TE=1e-3;
metabolites=struct("T1_s",.4,"T2_s",.4,"freq_shift_Hz",0,"name","water","T2star_s",0.9);

FA_ernst=rad2deg(acos(exp(-TR/metabolites(1).T1_s)));
Msig_ersnst=MetSignalModel   (metabolites,TE,pi,TR,0,deg2rad(FA_ernst),'GRE',0.9);
 [Msig_all]=squeeze(MetSignalModel   (metabolites,TE,pi,TR,0,deg2rad(FA),'GRE',0.9));
 assert(all(abs((Msig_all))<Msig_ersnst),"Ernst angle gives maximum flip angle")

Msig_ersnst=MetSignalModel   (metabolites,TE,pi,TR,0,deg2rad(FA_ernst),'GRE-peters',0.9);
  [Msig_all]=squeeze(MetSignalModel(metabolites,TE,pi,TR,0,deg2rad(FA),'GRE-peters',0.9));
 assert(all(abs((Msig_all))<Msig_ersnst),"Ernst angle gives maximum flip angle")



fprintf('all (8/8) assertion passed\n')
