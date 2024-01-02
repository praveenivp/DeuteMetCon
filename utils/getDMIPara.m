function para=getDMIPara(twix)

if(iscell(twix))
twix=twix{1};
end



% try to keep all times in s and all angular stuff in rad
para.isCSI= contains(twix.hdr.Phoenix.tSequenceFileName,'csi');

if(~para.isCSI)
% echoes and phase cycles (seq specific)
EchoSel=1:twix.hdr.Phoenix.lContrasts;
% para.PCSel=1:twix.image.NRep;
para.TE=[twix.hdr.Phoenix.alTE{EchoSel}]*1e-6; %s
para.TR=twix.hdr.Phoenix.alTR{1}*1e-6; %s
Range=twix.hdr.Phoenix.sWipMemBlock.alFree{5};
NRep=twix.image.NRep;
Shift=twix.hdr.Phoenix.sWipMemBlock.alFree{6};
if(isempty(Shift));Shift=0;end
para.PC_deg = (Range/(2*NRep):Range/NRep:Range-Range/(2*NRep))+Shift;
para.PhaseCycles=deg2rad(para.PC_deg);

%transmitter
para.FlipAngle=deg2rad(twix.hdr.Dicom.adFlipAngleDegree); %rad
para.pulseVoltage=twix.hdr.Phoenix.sTXSPEC.aRFPULSE{1}.flAmplitude;
     %only work for rect pulse
    para.RefVoltage=twix.hdr.Spice.TransmitterReferenceAmplitude;
    para.PulseDur=(para.RefVoltage/para.pulseVoltage)*(0.5*45/90);
    para.FrequencySystem=twix.hdr.Dicom.lFrequency;
% get timings probably only work for VE12U
 para.acq_duration_s=(twix.image.timestamp(end)-twix.image.timestamp(1))*2.5e-3; %s 
   tok=regexp(twix.hdr.Phoenix.tReferenceImage0,'.(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})\d+$','tokens');
   tok=cellfun(@(x)(str2double(x)),tok,'UniformOutput',false);
   para.StudyDateTime= datetime(tok{1});
   para.SeriesDateTime_start=datetime(tok{1}(1:3)) +seconds(twix.image.timestamp(1)*2.5e-3);
   para.SeriesDateTime_end=datetime(tok{1}(1:3))+seconds(twix.image.timestamp(end)*2.5e-3);

   para.seq_details=printSeqeunceDetails(twix);

else

CSI_VectorSize   = twix.image.NCol;
CSI_SamplingTime = twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*CSI_VectorSize*1e-9;
CSI_TimePoints   = [0:(CSI_VectorSize-1)]'*CSI_SamplingTime/CSI_VectorSize;
CSI_DwellTime    = CSI_SamplingTime/CSI_VectorSize;
para.Bandwidth    = 1/CSI_DwellTime;
% CSI_ppmRange     = CSI_Bandwidth/CSI_sMrprot.Dicom.lFrequency*1e6;
para.FreqAxis=linspace(-0.5/CSI_DwellTime,0.5/CSI_DwellTime,(CSI_VectorSize));

    % treat all time axis as echo dimension
EchoSel=1:twix.hdr.Phoenix.sSpecPara.lVectorSize;
% take rmos flag into account
if(twix.image.NCol == size(twix.image(:,:,1),1)) 
para.dwell=twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9; %s
else %'rmos' flag active
para.dwell=twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9*2; %s
end


% para.PCSel=1:twix.image.NRep;
para.TE=linspace(0,para.dwell*EchoSel(end-1),EchoSel(end)); %s
para.TE=para.TE+twix.hdr.Phoenix.alTE{1}*1e-6; %s
para.TR=twix.hdr.Phoenix.alTR{1}*1e-6; %s
Range=0; %twix.hdr.Phoenix.sWipMemBlock.alFree{5};
NRep=twix.image.NRep;
Shift=[]; %twix.hdr.Phoenix.sWipMemBlock.alFree{6};
if(isempty(Shift));Shift=0;end
para.PC_deg = 0; %(Range/(2*NRep):Range/NRep:Range-Range/(2*NRep))+Shift;
para.PhaseCycles=1; %deg2rad(para.PC_deg);

%transmitter
para.FlipAngle=deg2rad(twix.hdr.Dicom.adFlipAngleDegree); %rad
para.PulseVoltage=twix.hdr.Phoenix.sTXSPEC.aRFPULSE{1}.flAmplitude;
     %only work for rect pulse
    para.RefVoltage=twix.hdr.Spice.TransmitterReferenceAmplitude;
    para.PulseDur=(para.RefVoltage/para.PulseVoltage)*(0.5*45/90);

% get timings probably only work for VE12U
 para.AcquitionDuration=(twix.image.timestamp(end)-twix.image.timestamp(1))*2.5e-3; %s 
   tok=regexp(twix.hdr.Phoenix.tReferenceImage0,'.(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})\d+$','tokens');
   tok=cellfun(@(x)(str2double(x)),tok,'UniformOutput',false);
   para.StudyDateTime= datetime(tok{1});
   para.SeriesDateTime_start=datetime(tok{1}(1:3)) +seconds(twix.image.timestamp(1)*2.5e-3);
   para.SeriesDateTime_end=datetime(tok{1}(1:3))+seconds(twix.image.timestamp(end)*2.5e-3);

   para.SeqDetails=printSeqeunceDetails(twix);

end

end