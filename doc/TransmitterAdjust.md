# Calibration
## transmitter adjust
We try to determine the RF pulse amplitude required for 0.5 ms rectangular pulse to get 90 deg flip angle or (1 ms/180 deg flip angle).  

The reference voltage and 
```matlab 
para.RefVoltage=twix.hdr.Spice.TransmitterReferenceAmplitude; % [volt]
para.FlipAngle=deg2rad(twix.hdr.Dicom.adFlipAngleDegree); % [rad]
para.pulseVoltage=twix.hdr.Phoenix.sTXSPEC.aRFPULSE{1}.flAmplitude; % [volt]
```

## Dico Adjust
In addition to measuring forward and reflected power for safety,  we also measure the relationship between the voltage set in software and RFPA output measurement.  The maximum pulse voltage is determined from here (only if it less than maximum pulse voltage from the coil file?). 

## what is flip angle

          Flip angle ->         pulse voltage                  ->           RFPA                         -> time varying field.
problems:                 Reference voltage limit                  max pulse voltage -> pulse clipping

