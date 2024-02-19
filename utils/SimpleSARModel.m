function [SAR,maxFA]=SimpleSARModel(FA,PulseDur,TR,RefAmp,kFactor)
% simple SAR model to calculate average SAR and maximum possible flip angle
% for rectangular pulse.  
%
% [SAR,maxFA]=SimpleSARModel(FA,PulseDur,TR,RefAmp=447,kFactor=1)
%
% INPUTS
% FA - flip angle in degrees
% PulseDur - RF pulse duration in seconds
% TR - repetition time in seconds
% RefAmp - Referenece amplitude in V from transmitter adjust
% kFactor - from coil file [in (W/kg) / W]
%
%OUTPUT
% SAR [%]
% maxFA [deg]


if(~exist('RefAmp','var'))
    RefAmp=447; %v
end

if(~exist('kFactor','var'))
    kFactor=1; %[1/kg]
end

% A 0.5 ms block pulse with reference ampltitude yields 90 deg flipangle.
pulseVoltage= RefAmp*(FA/90)*(0.5e-3/PulseDur); %V
if(pulseVoltage>387)
    warning('Input Pulse voltage is probably too high : %.2f V',pulseVoltage)
end

inputPower= pulseVoltage.*pulseVoltage/50; %W

RF_dutyCycle= PulseDur./TR;

SAR=kFactor*inputPower*RF_dutyCycle;

% 6 min SAR limit. for some reason I can do only 7.5 W/kg not 10 W/kg
maxSAR_6min=7.5; %W/kg
% maxSAR_10s=20; %W/kg % not relevant mostly in our case

%assuming RF duty cycle is the same
maxInputPower=inputPower.*(maxSAR_6min./SAR);

maxPulseVoltage= sqrt(maxInputPower*50);
if(maxPulseVoltage>387)
    warning('maximum Pulse voltage is probably too high : %.2f V : clipping',maxPulseVoltage)
    maxPulseVoltage=387;
end
maxFA= (maxPulseVoltage/RefAmp)*(PulseDur/0.5e-3)*90;


%% verbose output if no output argumetns
if(nargout==0)
fprintf('Input pulse Voltage = %.2f V\n',pulseVoltage);
fprintf('Estimated SAR = %.2f W/kg\n',SAR);
fprintf('maximum possible Pulse Voltage = %.2f V\n',maxPulseVoltage);
fprintf('maximum possible flipangle = %.2f deg\n',maxFA);
end

end
