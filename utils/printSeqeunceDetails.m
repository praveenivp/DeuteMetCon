function [s]=printSeqeunceDetails(twix_obj)
c=1;
p=twix_obj.hdr.Phoenix;
s{c}=sprintf('Meas ID: %d  | %dD sequence   ',twix_obj.hdr.Meas.MeasUID); c=c+1;
try
if(twix_obj.hdr.Meas.ulEnableRFSpoiling)
    s{c}=sprintf('Without RFSpoiling');  c=c+1;
else 
    s{c}=sprintf('With RF Spoiling');  c=c+1;
end
catch
    s{c}=sprintf('RFSpoiling unknown'); c=c+1;
end
s{c}=sprintf(' |  %s |  %s.dll \n',twix_obj.hdr.Config.SequenceDescription,twix_obj.hdr.Config.SequenceFileName); c=c+1;
s{c}=sprintf('TR (ms): %.2f \n',1e-3*p.alTR{1}); c=c+1;
s{c}=[strcat(sprintf('TE (ms): '),sprintf('%.2f ', 1e-3*[p.alTE{1:p.lContrasts}])) newline ]; c=c+1;
s{c}=sprintf('Encoding matrix (RxPxS) : %d x %d x %d\n',p.sKSpace.lBaseResolution,p.sKSpace.lPhaseEncodingLines,p.sKSpace.lPartitions); c=c+1;
if isfield(twix_obj.hdr.Dicom,'adFlipAngleDegree')
    s{c}=sprintf('Flip angle(deg): %d\n',twix_obj.hdr.Dicom.adFlipAngleDegree);  c=c+1;
end
s{c}=sprintf('%d Repetitions and %d Averages\n', twix_obj.image.NRep,twix_obj.image.NAve);  c=c+1;
s{c}=sprintf('Bandwidth: %d Hz/pxl and dwell time(with Read OS): %d ns \n', round(1/(2e-9*p.sKSpace.lBaseResolution*twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1})),twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1});  c=c+1;

p=p.sSliceArray;
for i=1:size(p.asSlice,2)
s{c}=sprintf('\n**********Slice #%d ************\n',i);     c=c+1;
s{c}=sprintf('FOV(RxPxS mm):%d X %d X %d \n',p.asSlice{i}.dReadoutFOV,p.asSlice{i}.dPhaseFOV,p.asSlice{i}.dThickness); c=c+1;
s{c}=sprintf('Orientation: %s \n',string(fieldnames(p.asSlice{i}.sNormal))); c=c+1;

if(isfield(p.asSlice{i},'sPosition'))
    temp=fieldnames(p.asSlice{i}.sPosition); 
    for j=1:size(temp,1)
        s{c}=sprintf('off-center slice: %s = %d mm \n',temp{j},eval(strcat('p.asSlice{i}.sPosition.',temp{j})));  c=c+1;
    end
else
         s{c}=sprintf('Position: Isocenter\n');  c=c+1;
end
end
if(nargout==0)
fprintf('%s',s{:})
end


% pixel shift
p=twix_obj.hdr.Phoenix;
BW_pxl=round(1/(2e-9*p.sKSpace.lBaseResolution*twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1})); % Hz/pxl
BR=p.sKSpace.lBaseResolution;
DW=twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*2e-9; %s no readout OS
Tacq=1/BW_pxl; %s
p=p.sSliceArray;
FOV_read=p.asSlice{i}.dReadoutFOV*1e-3; %m no readout OS
res_read=FOV_read/BR;
% gamma=42.567e6; %Hz/T
gamma=6.536e6 ; %2H Hz/T
Gread=(1/(gamma*FOV_read*DW)); % T/m
offresonace= 150;% Hz
pxl_shift=offresonace/(gamma*Gread); %m
pxl_shift=pxl_shift/res_read; %should be same as offresonace/BW_pxl

end

