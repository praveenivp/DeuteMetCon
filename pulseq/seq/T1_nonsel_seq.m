%input parameters
Nx=512;
Nav=16;
FA=pi/2;
dwell_s=200e-6; %[s]
delayTE=2.5e-3; %[s]
TR_array=[110:10:190  200:20:400 450:50:600 700:200:1500 2000]*1e-3-delayTE; %[s]
wait_s=5;
dummy_scans=8;
Nrep=length(TR_array);


fprintf('readout: %.2f ms\n',dwell_s*Nx*1e3);
fprintf('Spectral resolution %.2f Hz\n',1/(Nx*dwell_s))
system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration',2e-3, 'system', system);

% Define delays and ADC events
adc = mr.makeAdc(Nx,'Duration',dwell_s*Nx, 'system', system,'delay',system.adcDeadTime);

%
spoilerx=mr.makeTrapezoid('x','amplitude',20,'duration',5e-3);
spoilery=mr.makeTrapezoid('y','amplitude',20,'duration',5e-3);
spoilerz=mr.makeTrapezoid('z','amplitude',20,'duration',5e-3);
%
assert(delayTE>=mr.calcDuration(rf));
assert(TR_array(1)>= (mr.calcDuration(adc)+mr.calcDuration(spoilerx)));
% Loop over repetitions and define sequence blocks

        

for rep=1:Nrep

for i=1:(Nav+dummy_scans)
    rf = mr.makeBlockPulse(FA,'Duration',2e-3, 'system', system,'phaseOffset', 0);
    adc = mr.makeAdc(Nx,'Duration',dwell_s*Nx, 'system', system,'delay',system.adcDeadTime,'phaseOffset', 0);
    seq.addBlock(rf,mr.makeDelay(delayTE));
    if(i>dummy_scans)
         seq.addBlock(adc);
         seq.addBlock(spoilerx,spoilery,spoilerz);
         seq.addBlock(mr.makeDelay(TR_array(rep)- ...
             (mr.calcDuration(adc)+mr.calcDuration(spoiler))));
    else
        seq.addBlock(mr.makeDelay(mr.calcDuration(adc)))
         seq.addBlock(spoilerx,spoilery,spoilerz);
         seq.addBlock(mr.makeDelay(TR_array(rep)- ...
             (mr.calcDuration(adc)+mr.calcDuration(spoiler))));
    end
end
seq.addBlock(mr.makeDelay(wait_s));
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
seq.setDefinition('TR_array',TR_array);
seq.setDefinition('rf_dur',2e-3);
seq.setDefinition('averages',Nav);
seq.setDefinition('repetitions',Nrep);


pn='\\mrz10\upload9t\USERS\Praveen\20230920_xpulseq';
delete('T1_nonsel_2H.seq');
seq.write(fullfile(pn,'T1_nonsel_2H.seq'))       % Write to pulseq file