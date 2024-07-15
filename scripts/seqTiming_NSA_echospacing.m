%% echo spacing
addpath(genpath('C:\Users\pvalsala\Documents\Packages2\pulseq\matlab'));
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
%invivo: /ptmp/pvalsala/deuterium/EAZW-GUMK/proc
met_name={'water','glucose','Glx','lactate'};
freq_shift_WGX=[3.9 -58.7  -152 -210];
T1=[432.88 69.65 147.6 190.64]*1e-3;%s
T2=[287 65.88 124.5 180]*1e-3;%s

clear metabolites;
for i=1:length(freq_shift_WGX)
    metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift_WGX(i),'name',met_name{i});
end

figure, 
deltaT_={(0.1:0.01:6)*1e-3,(0.01:0.01:0.5)*1e-3};
N_=[5, 64*3];

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

if(i==1)
subplot(1,2,i),plot(deltaT*1e3,1./(N*NSA),'LineWidth',2)
hold on
line([4.4 4.4],[0 1],'linewidth',2),line([1 1]*3.7,[0 1],'linewidth',2)
xlabel('echo spacing [ms]'),ylabel('NSA')
title('multi-echo | N=5'),grid on
else
    subplot(1,2,i),plot(deltaT*1e3*N,1./(N*NSA),'LineWidth',2)
title('CSI '),grid on
xlabel('readout length [ms]'),ylabel('NSA')
legend({metabolites.name},'Location','southeast')
hold on
line([1 1]*14.9,[0 1],'linewidth',2),line([1 1]*94.8,[0 1],'linewidth',2)
end

end



fontsize(gcf,'scale',1.5)
set(gcf,'color','w')