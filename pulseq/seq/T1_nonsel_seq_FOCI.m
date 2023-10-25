%% Load FOCI pulse
load('S:\Deuterium\20230920_xpulseq\foci.mat')
taxis2=0:1e-3:5.2-1e-3;
 data2=interp1(taxis,data,taxis2,"linear",0);
addpath(genpath('C:\Users\pvalsala\Documents\Packages2\pulseq\matlab'));
%% input parameters
Nx=512;

FA=pi/2;
pulse_dur=2e-3;
dwell_s=250e-6; %[s]
TE=pulse_dur; %[s]
TR= 2500e-3; %[s]

% 10 mins with 8 averages TR 2.5 s
%  TI_array=[20:15:280 300:50:600 700:250:1500 2000]*1e-3; %[s]
%  5 mins with 8 averages TR 2.5s
 TI_array=[ 25:20:280 300:100:600 700:300:1000]*1e-3; %[s]
Nav=8;
wait_s=0;
dummy_scans=0;
Nrep=length(TI_array);




system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 200e-6,'maxSlew',140,'slewUnit','T/m/s','maxGrad',30,'GradUnit','mT/m','gamma',6.536e6,'B0',9.4);

seq=mr.Sequence(system);              % Create a new sequence object

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration',pulse_dur, 'system', system);
%
crusherx=mr.makeTrapezoid('x','amplitude',10e-3*system.gamma,'duration',10e-3);
crushery=mr.makeTrapezoid('y','amplitude',10e-3*system.gamma,'duration',10e-3);
crusherz=mr.makeTrapezoid('z','amplitude',10e-3*system.gamma,'duration',10e-3);
% rf_foci = mr.makeAdiabaticPulse('hypsec');

rf_foci=mr.makeArbitraryRf(data2(:,1).*exp(1i*2*pi*data2(:,2)),pi,...
    'FreqOffset', 0,'PhaseOffset',0,'system',system);%,'dwell',1e-6

% Define delays and ADC events
adc = mr.makeAdc(Nx,'Duration',dwell_s*Nx, 'system', system,'delay',system.adcDeadTime);

%
spoilerx=mr.makeTrapezoid('x','amplitude',20e-3*system.gamma,'duration',10e-3);
spoilery=mr.makeTrapezoid('y','amplitude',20e-3*system.gamma,'duration',10e-3);
spoilerz=mr.makeTrapezoid('z','amplitude',20e-3*system.gamma,'duration',10e-3);
%
assert(TE>mr.calcDuration(rf)/2);
assert(TR >= (mr.calcDuration(adc)+mr.calcDuration(spoilerx)));

time_until_exc= TI_array(1)-mr.calcDuration(rf_foci)/2 -mr.calcDuration(rf)/2-mr.calcDuration(crusherx);
assert(time_until_exc>1e-3);
% Loop over repetitions and define sequence blocks

        

for rep=1:Nrep
    for i=1:(Nav+dummy_scans)
        seq.addBlock(rf_foci);

        % inversion block
        seq.addBlock(crusherx,crushery,crusherz);
        time_until_exc= TI_array(rep)-mr.calcDuration(rf_foci)/2 -mr.calcDuration(rf)/2-mr.calcDuration(crusherx);
        seq.addBlock(mr.makeDelay(time_until_exc));


        rand_phase = mod(117*(i^2 + i + 2), 360)*pi/180;
        rf = mr.makeBlockPulse(FA,'Duration',pulse_dur, 'system', system,'phaseOffset', 0);
        seq.addBlock(rf);
        time_until_adc= TE-mr.calcDuration(rf)/2;
        seq.addBlock(mr.makeDelay(time_until_adc));


        if(i>dummy_scans)
            seq.addBlock(adc);
        else
            seq.addBlock(mr.makeDelay(mr.calcDuration(adc)))
        end
        % fill the rest of TR
        TR_fill= TR- (mr.calcDuration(rf_foci)/2+TI_array(rep)+ ...
            TE+mr.calcDuration(adc)+mr.calcDuration(spoilerx));
        seq.addBlock(spoilerx,spoilery,spoilerz);
        seq.addBlock(mr.makeDelay(TR_fill));
    end
    % inter rep delay
    if(wait_s>0)
        seq.addBlock(mr.makeDelay(wait_s));
    end
end
% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
seq.setDefinition('Name', 'T1_nonsel');
seq.setDefinition('dwell', dwell_s);
seq.setDefinition('TR',TR);
seq.setDefinition('TI_array',TI_array);
seq.setDefinition('rf_dur',2e-3);
seq.setDefinition('averages',Nav);
seq.setDefinition('repetitions',Nrep);
seq.setDefinition('MeasurementTime',seq.duration)

% display helpfull pararmeters
fprintf('readout: %.2f ms\n',dwell_s*Nx*1e3);
fprintf('Spectral resolution %.2f Hz\n',1/(Nx*dwell_s))
fprintf('Nrep: %d , Nav: %d, %.1f min \n',Nrep,Nav,seq.duration/60)


seq.plot()
pn='\\mrz10\upload9t\USERS\Praveen\20230920_xpulseq';
pulseqFileName=sprintf('%s_%dmin.seq','T1_nonsel_FOCI_2H',round(seq.duration/60));
 delete(pulseqFileName);
 seq.write(fullfile(pn,pulseqFileName))       % Write to pulseq file