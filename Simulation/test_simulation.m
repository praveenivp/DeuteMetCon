%% test
    cs_ppm=[0 1.0256 ]*1e-6; %D20 Glu Glx Lac/fat
    freq_shift_WGX=[3.9 -58.7  -152 -215];
    T1=[432.88 69.65 147.6 190.64]*1e-3;%s
    T2=[287 65.88 124.5 180]*1e-3;%s

    sp_name={'D20','Glucose','Glutamate','Lactate'};
    clear metabolites;
    for i=1:3
        metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift_WGX(i),'name',sp_name{i});
    end

TR=23e-3;
TE=linspace(2e-3,TR-1e-3,5);
dfreq=0;
PhaseCyles=linspace(0,2*pi-(2*pi/20),20);
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

MetSimObj=MetSim(metabolites,TE,PhaseCyles,TR,FA,'bSSFP');

MetSimObj.getSig();
as(MetSimObj.sig);
as(myfft(MetSimObj.sig,3))


%%

%%
clearvars
    cs_ppm=[0 1.0256 ]*1e-6; %D20 Glu Glx Lac/fat
    freq_shift_WGX=[3.9 -58.7  -152 -215];
    T1=[432.88 69.65 147.6 190.64]*1e-3;%s
    T2=[287 65.88 124.5 180]*1e-3;%s

     T2star=[287 65.88 124.5 180]*1e-3*(1/15);%s
clear metabolites
    sp_name={'D20','Glucose','Glutamate','Lactate'};
    for i=1:4
        metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift_WGX(i),'T2star_s',T2star(i),'name',sp_name{i});
    end

TR=23e-3;
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


%% test sim class
clc
TR=23e-3;
TE=linspace(2e-3,TR-2e-3,10);
dfreq=0;
PhaseCyles=linspace(0,2*pi-(2*pi/20),1);
FA=deg2rad(50);
MetSimObj=MetSim(metabolites,TE,PhaseCyles,TR,FA,'GRE','NoiseFac',0.02,'MatSize',32);

MetSimObj.getSig();
as(MetSimObj.sig);
% as(myfft(MetSimObj.sig,3))

%% test ideal solver
clc
 IDEALobj_pinv=IDEAL(metabolites,TE,'mask',sum(MetSimObj.experimental.Phantom,3)>0,'solver','pinv');
 sig=MetSimObj.sig;
%  sig=sig./abs(sig);
 metabol_con_pinv=IDEALobj_pinv'*permute(sig,[1 2 4 3]);

as(metabol_con_pinv,'select',':,:,:','title','pinv')
%%

 IDEALobj_ideal=IDEAL(metabolites,TE,'solver','IDEAL','maxit',100);%,'mask',sum(MetSimObj.experimental.Phantom,3)>0);
 metabol_con_ideal=IDEALobj_ideal'*permute(sig,[1 2 4 3]);
as(metabol_con_ideal,'select',':,:,:')
