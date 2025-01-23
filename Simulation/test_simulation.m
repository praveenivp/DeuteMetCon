%% test
clearvars
metabolites=getMetaboliteStruct('invivo',0);
%T2* messes up spectral seperation.
T2_star={Inf,Inf,Inf,Inf};
[metabolites.T2star_s]=T2_star{:};


TR=19e-3;
N=5;

TE=2.4e-3+(0:N-1)*3.6e-3;
assert(max(TE)<(TR-2e-3))
FA=deg2rad(50);

MetSimObj_bssfp=MetSim(metabolites,TE,pi,TR,FA,'GRE','MatSize',32,'Noisefac',1e-3,'B0range',50);
sig=MetSimObj_bssfp.sig;
fm=MetSimObj_bssfp.FieldMap_Hz;
fm_fac=((2*pi)*(42.567/6.536)); % 2H [Hz] to 1H [rad/s]
sig2=mean(sig,4);
% as(sig2);
% as(myfft(MetSimObj.sig,3))

cMask=sum(MetSimObj_bssfp.experimental.Phantom,3)>0;
 IDEALobj_pinv=IDEAL(metabolites,TE,'solver','phaseonly','fm',0*fm*fm_fac,'mask',cMask,'PhaseCorr',false,'mask',cMask);
 metabol_con_pinv=IDEALobj_pinv'*permute(sig2,[1 2 4 3]);
as(metabol_con_pinv,'select',':,:,:','title','pinv','complexSelect','m','colormap','jet')
%%
IDEALOptions={'SmoothFM',2,'maxit',10,'tol',1};
 IDEALobj_ideal=IDEAL(metabolites,TE,'solver','IDEAL','fm',[],'mask',cMask,'PhaseCorr',false,IDEALOptions{:});
 metabol_con_ideal=IDEALobj_ideal'*permute(sig2,[1 2 4 3]);
as(metabol_con_ideal,'select',':,:,:','title','IDEAL','complexSelect','m','colormap','jet')

%%

TR=20e-3;
TE=linspace(2e-3,TR-2e-3,5);
dfreq=0;
PhaseCyles=linspace(0,2*pi-(2*pi/20),1);
FA=deg2rad(50);
Msig_all=MetSignalModel(metabolites,TE,0,TR,dfreq,FA,'GRE');

figure,plot(TE*1e3,squeeze(abs(Msig_all(:,:)))')
yyaxis('right')
plot(TE*1e3,rad2deg(angle(squeeze(Msig_all(:,:,1)))'))
xlabel('echo time [ms]')
legend(sp_name{1:3},sp_name{1:3})

%% plot conditioning/NSA
 figure, 
deltaT_={(0.1:0.01:6)*1e-3,(0.1:0.01:0.5)*1e-3};
N_=[5, 64];

for i=1:2
    deltaT=deltaT_{i};
    N=N_(i)
NSA=zeros(length(metabolites),length(deltaT));
for cDT=1:length(deltaT)
    TE=2e-3+(0:N-1)*deltaT(cDT);
 IDEALobj_pinv=IDEAL(metabolites,TE,'solver','phaseonly','PhaseCorr',true);
 A=(IDEALobj_pinv.getA);
NSA(:,cDT)=abs(diag(inv(A'*A)));
end

subplot(1,2,i),plot(deltaT*1e3,1./(NSA),'LineWidth',2)
xlabel('echo spacing [ms]'),ylabel('NSA')
title('multi-echo | N=5'),grid on
end
title('CSI | N=64')
legend({metabolites.name},'Location','southeast')
fontsize(gcf,'scale',1.5)
% set(gcf,'color','w')
%% test sim class
clc
TR=19e-3;
N=5;

TE=2e-3+(0:N-1)*3.6e-3;
assert(max(TE)<(TR-2e-3))



PhaseCyles=linspace(0,2*pi-(2*pi/20),1);
FA=deg2rad(50);
MetSimObj=MetSim(metabolites,TE,PhaseCyles,TR,FA,'GRE','NoiseFac',0.001,'B0Range',10,'MatSize',32);
sig=MetSimObj.sig;
% sig=sig./abs(sig);

cMask=sum(MetSimObj.experimental.Phantom,3)>0;
fm=1*MetSimObj.FieldMap_Hz;
% as(sig);
% as(myfft(MetSimObj.sig,3))

% test ideal solver



 IDEALobj_pinv=IDEAL(metabolites,TE,'mask',cMask,'solver','pinv','fm',fm,'PhaseCorr',true);
  metabol_con_pinv=IDEALobj_pinv'*permute(sig,[1 2 4 3]);
as(metabol_con_pinv,'select',':,:,:','title','pinv','complexSelect','m','colormap','jet')

A=(IDEALobj_pinv.getA);
 as(inv(A'*A))

