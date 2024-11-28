%% echo spacing
addpath(genpath('/ptmp/pvalsala/Packages/pulseq'));
res=12.5e-3; %m
 GammaH2=6.536e6 ; %Hz/T
 Nx=32;
 BW=480; %Hz/pixel

SR=200; % mT/m/ms 

 DW=1/(Nx*BW); %s
 kmax=pi/res; % rad/m
 Gread=1e3/(GammaH2*res*(1/BW)); %mT/m


    %% using pulseq
BW_All=370;
clear data
for i=1:length(BW_All)
Nx=32; Ny=Nx; Nz=28;            % Define FOV and resolution
res=12.5e-3;
fov=[Nx Ny Nz]*res;     % Define FOV
Tread=1/BW_All(i);
% define system properties
sys=mr.opts('maxGrad',42,'gradUnit','mT/m','maxSlew',200,'slewUnit','mT/m/ms',...
            'rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6,'gamma',6.536e6,'B0',9.38);
seq=mr.Sequence(sys);           % Create a new sequence object


% Create non-selective pulse
[rf, rfDelay] = mr.makeBlockPulse(8*pi/180,sys,'Duration',1.4e-3);

% Define other gradients and ADC events
deltak=1./fov;
gx = mr.makeTrapezoid('x',sys,'FlatArea',Nx*deltak(1),'FlatTime',Tread);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2);
gxfly = mr.makeTrapezoid('x',sys,'Area',-gx.area);

Necho=5;
ES=1*((mr.calcDuration(gx)+mr.calcDuration(gxfly)))
TE=1*(mr.calcDuration(rf)/2+mr.calcDuration(gxPre)+mr.calcDuration(gx)/2);
TE2=TE+ES*(0:4);

TR=mr.calcDuration(rf)/2+mr.calcDuration(gxPre) + ES*Necho+mr.calcDuration(gxPre);

DC=Necho*mr.calcDuration(adc)/TR;
data.TR(i)=TR;
data.ES(i)=ES;
data.TE(i)=TE;
data.DC(i)=DC;
data.BW(i)=BW_All(i);
end

%
figure,tiledlayout("flow")
nexttile()
plot(data.BW,data.DC*1e2),title('Duty [%]')
nexttile()
plot(data.BW,data.TR*1e3),title('TR [ms]')
nexttile()
plot(data.BW,data.ES*1e3),title('ES [ms]')
nexttile()
plot(data.BW,data.TE*1e3),title('TE1 [ms]')
%%
seq.addBlock(rf)
seq.addBlock(gxPre)
seq.addBlock(gx,adc)
seq.addBlock(gxfly)
seq.addBlock(gx,adc)
seq.addBlock(gxfly)
seq.addBlock(gx,adc)
seq.addBlock(gxfly)
seq.addBlock(gx,adc)
seq.addBlock(gxfly)
seq.addBlock(gxPre)
seq.plot
%% plot conditioning/NSA
metabolites=getMetaboliteStruct('invivo');

deltaT=(0.01:0.01:0.5)*1e-3;
N=64;
NSA_CSI=zeros(length(metabolites),length(deltaT));
for cDT=1:length(deltaT)
    TE=2e-3+(0:N-1)*deltaT(cDT);
    IDEALobj_pinv=IDEAL(metabolites,TE,'solver','phaseonly','PhaseCorr',true);
    A=(IDEALobj_pinv.getA);
    NSA_CSI(:,cDT)=abs(diag(inv(A'*A)));
end
deltaTE=(0.1:0.01:6)*1e-3;
NE=5;
NSA_ME=zeros(length(metabolites),length(deltaTE));
for cDT=1:length(deltaTE)
    TE=2e-3+(0:NE-1)*deltaTE(cDT);
    IDEALobj_pinv=IDEAL(metabolites,TE,'solver','phaseonly','PhaseCorr',true);
    A=(IDEALobj_pinv.getA);
    NSA_ME(:,cDT)=abs(diag(inv(A'*A)));
end


figure,
tt=tiledlayout(1,2,"TileSpacing","compact",'Padding','compact');

nexttile()
plot(deltaTE*1e3,1./(NE*NSA_ME),'LineWidth',2)
hold on
line([1 1]*3.7,[0 1],'linewidth',2,'color','magenta')
xlabel('echo spacing [ms]'),ylabel('Normalized NSA')
title('multi-echo | N=5'),grid on
legend([{metabolites.name},'ME-bSSFP'],'Location','northwest')

nexttile()
plot(deltaT*1e3*N,1./(N*NSA_CSI),'LineWidth',2)
title('CSI '),grid on
xlabel('readout length [ms]'),ylabel('Normalized NSA')
hold on
line([1 1]*14.9,[0 1],'linewidth',2,'color','magenta' ),line([1 1]*28.44,[0 1],'linewidth',2,'color','cyan')
legend([{metabolites.name},'CSI-bSSFP','CSI-FISP'],'Location','northwest')
fontsize(gcf,'scale',1.5)
set(gcf,'color','w','Position', [335 381 1289 556])