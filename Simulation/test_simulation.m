%% test
    cs_ppm=[0 1.0256 ]*1e-6; %D20 Glu Glx Lac/fat
    freq_shift_WGX=[3.9 -58.7  -152 -215];
    T1=[432.88 69.65 147.6 190.64]*1e-3;%s
    T2=[287 65.88 124.5 180]*1e-3;%s

    sp_name={'D20','Glucose','Glutamate','Lactate'};
    clear metabolites;
    for i=1:4
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

TR=30e-3;
TE=linspace(2e-3,TR-1e-3,20);
dfreq=0;
PhaseCyles=linspace(0,2*pi-(2*pi/20),1);
FA=deg2rad(50);
Msig_all=MetSignalModel(metabolites,TE,0,TR,dfreq,FA,'GRE');

figure,plot(TE*1e3,squeeze(abs(Msig_all(:,:)))')
yyaxis('right')
plot(TE*1e3,rad2deg(angle(squeeze(Msig_all(:,:,1)))'))
xlabel('echo time [ms]')
legend(sp_name{:},sp_name{:})


%% test sim class

MetSimObj=MetSim(metabolites,TE,PhaseCyles,TR,FA,'GRE');

MetSimObj.getSig();
as(MetSimObj.sig);
% as(myfft(MetSimObj.sig,3))

%% test ideal solver

 IDEALobj_ideal=IDEAL(metabolites,TE,'mask',sum(MetSimObj.experimental.Phantom,3)>0,'solver','IDEAL');
 sig=MetSimObj.sig;
%  sig=sig./abs(sig);
 metabol_con_ideal=IDEALobj_ideal'*permute(sig,[1 2 4 3]);

as(metabol_con_ideal,'select',':,:,:')
%%

 IDEALobj_ideal=IDEALsolver(metabolites,TE_s,0*fm_Hz,'IDEAL',300);
 metabol_con_ideal=IDEALobj_ideal'*sig;

