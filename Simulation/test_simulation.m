%% test
clearvars
    freq_shift_WGX=[3.9 -58.7  -152 -215];
    T1=[432.88 69.65 147.6 190.64]*1e-3;%s
    T2=[287 65.88 124.5 180]*1e-3;%s

    sp_name={'D20','Glucose','Glutamate','Lactate'};
    clear metabolites;
    for i=1:4
        metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift_WGX(i),'name',sp_name{i});
    end

TR=19e-3;
N=5;

TE=2e-3+(0:N-1)*3.6e-3;
assert(max(TE)<(TR-2e-3))

dfreq=0;
PhaseCyles=linspace(0,2*pi-(2*pi/20),40);
FA=deg2rad(50);
Msig_all=MetSignalModel(metabolites,TE,PhaseCyles,TR,dfreq,FA,'bSSFP');

figure,
tiledlayout("flow")
for i=1:length(TE)
nexttile()
plot(squeeze(abs(Msig_all(:,1,:)))')
yyaxis('right')
plot(TE*1e3,rad2deg(angle(squeeze(Msig_all(:,:,1)))'))
xlabel('echo time [ms]')
legend(sp_name{:},sp_name{:})
end

%% test sim class

MetSimObj_bssfp=MetSim(metabolites,TE,PhaseCyles,TR,FA,'bSSFP','MatSize',32,'Noisefac',0.01,'B0range',50);
sig=MetSimObj_bssfp.sig;
fm=MetSimObj_bssfp.FieldMap_Hz;
sig2=mean(sig,4);
as(sig2);
% as(myfft(MetSimObj.sig,3))


%% test ideal solver


cMask=sum(MetSimObj_bssfp.experimental.Phantom,3)>0;
 IDEALobj_pinv=IDEAL(metabolites,TE,'solver','phaseonly','fm',0*fm,'mask',cMask,'PhaseCorr',true,'mask',cMask);
 metabol_con_pinv=IDEALobj_pinv'*permute(sig2,[1 2 4 3]);
as(metabol_con_pinv,'select',':,:,:','title','pinv','complexSelect','m','colormap','jet')
%%
IDEALOptions={'FilterSize',5,'maxit',100,'tol',0.1};
 IDEALobj_ideal=IDEAL(metabolites,TE,'solver','IDEAL','fm',0*fm,'mask',cMask,'PhaseCorr',true,IDEALOptions{:});
 metabol_con_ideal=IDEALobj_ideal'*permute(sig2,[1 2 4 3]);
as(metabol_con_ideal,'select',':,:,:','title','IDEAL','complexSelect','m','colormap','jet')

%%
clearvars
    cs_ppm=[0 1.0256 ]*1e-6; %D20 Glu Glx Lac/fat
    freq_shift_WGX=[3.9 -58.7  -152 -215]+5;
    T1=[432.88 69.65 147.6 190.64]*1e-3;%s
    T2=[287 65.88 124.5 180]*1e-3;%s

%      T2star=[287 65.88 124.5 180]*1e-3*(1/10);%s
     T2star=[20 13 15 15]*1e-3;
clear metabolites
    sp_name={'D20','Glucose','Glutamate','Lactate'};
    for i=1:4
        metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift_WGX(i),'T2star_s',T2star(i),'name',sp_name{i});
    end

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
clearvars
    cs_ppm=[0 1.0256 ]*1e-6; %D20 Glu Glx Lac/fat
%     freq_shift_WGX=[3.9 -58.7  -152 -215];
freq_shift_WGX=[3.9 -56.1  -145.2 -201.4];
    T1=[432.88 69.65 147.6 190.64]*1e-3;%s
    T2=[287 65.88 124.5 180]*1e-3;%s

%      T2star=[287 65.88 124.5 180]*1e-3*(1/10);%s
     T2star=[20 13 15 15]*1e-3;
clear metabolites
    sp_name={'D20','Glucose','Glutamate','Lactate'};
    for i=1:4
        metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift_WGX(i),'T2star_s',T2star(i),'name',sp_name{i});
    end
% figure, 
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


%%
IDEALOptions={'FilterSize',7,'maxit',100,'tol',0.1};
 IDEALobj_ideal=IDEAL(metabolites,TE,'solver','IDEAL','fm',0*fm,'PhaseCorr',true,'mask',cMask,IDEALOptions{:});%,'mask',sum(MetSimObj.experimental.Phantom,3)>0);
 metabol_con_ideal=IDEALobj_ideal'*permute(sig,[1 2 4 3]);
as(metabol_con_ideal,'select',':,:,:','complexSelect','m','colormap','jet')
