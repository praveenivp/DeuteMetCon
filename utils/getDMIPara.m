function para=getDMIPara(twix)
%para=getDMIPara(twix)
% fucntion to colllect and calculate some sequence parameters for
% convinience. works for both CSI and ME-bSSFP twix.


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
    try
    Range=twix.hdr.Phoenix.sWipMemBlock.alFree{5};
    NRep=twix.image.NRep;
    Shift=twix.hdr.Phoenix.sWipMemBlock.alFree{6};
    if(isempty(Shift));Shift=0;end
    para.PC_deg = (Range/(2*NRep):Range/NRep:Range-Range/(2*NRep))+Shift;
    catch %for fisp
      para.PC_deg =ones(twix.image.NRep,1)*180;  
    end
    para.PhaseCycles=deg2rad(para.PC_deg);

    %transmitter
    para.FlipAngle=deg2rad(twix.hdr.Dicom.adFlipAngleDegree)/1.07; %rad
    try
        para.pulseVoltage=twix.hdr.Phoenix.sTXSPEC.aRFPULSE{1}.flAmplitude;
        %only work for rect pulse
        para.RefVoltage=twix.hdr.Spice.TransmitterReferenceAmplitude;
        para.PulseDur=(para.RefVoltage/para.pulseVoltage)*(0.5*45/90);
    end
    para.FrequencySystem=twix.hdr.Dicom.lFrequency;

    % get timings probably only work for VE12U
    para.acq_duration_s=(twix.image.timestamp(end)-twix.image.timestamp(1))*2.5e-3; %s
    tok=regexp(twix.hdr.Phoenix.tReferenceImage0,'.(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})\d+$','tokens');
    tok=cellfun(@(x)(str2double(x)),tok,'UniformOutput',false);
    para.StudyDateTime= datetime(tok{1});
    para.SeriesDateTime_start=datetime(tok{1}(1:3)) +seconds(twix.image.timestamp(1)*2.5e-3);
    para.SeriesDateTime_end=datetime(tok{1}(1:3))+seconds(twix.image.timestamp(end)*2.5e-3);

    para.seq_details=printSeqeunceDetails(twix);

    % take rmos flag into account
    para.rmos_flag=twix.image.NCol == size(twix.image(:,:,1),1);
    if(para.rmos_flag)
        para.dwell=twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9; %s
    else %'rmos' flag active
        para.dwell=twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9*2; %s
    end

    %resolution
    sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};
    kp=twix.hdr.Phoenix.sKSpace;

    try %if slice oversampling!
        FOV=[sa.dReadoutFOV sa.dPhaseFOV*kp.dPhaseResolution sa.dThickness + sa.dThickness*kp.dSliceOversamplingForDialog]*1e-3; %m
    catch
        FOV=[sa.dReadoutFOV sa.dPhaseFOV*kp.dPhaseResolution  sa.dThickness]*1e-3; %m
    end

    para.MatrixSize=[kp.lBaseResolution  kp.lPhaseEncodingLines kp.lPartitions ];
    para.resolution=FOV./para.MatrixSize; %m
    para.ShortDescription=sprintf('M%d|TR %.0f ms| %.0f deg | %.2f mm | %d rep | %d echoes',twix.hdr.Config.MeasUID,para.TR*1e3,rad2deg(para.FlipAngle),para.resolution(1)*1e3,length(para.PhaseCycles),length(para.TE));

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
    para.rmos_flag=twix.image.NCol == size(twix.image(:,:,1),1);
    if(para.rmos_flag)
        para.dwell=twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9; %s
    else %'rmos' flag active
        para.dwell=twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9*2; %s
    end


    % para.PCSel=1:twix.image.NRep;
    para.TE=linspace(0,para.dwell*EchoSel(end-1),EchoSel(end)); %s
%     para.TE=para.TE+twix.hdr.Phoenix.alTE{1}*1e-6; %s
    para.TE=para.TE+twix.hdr.Phoenix.sSpecPara.lAcquisitionDelay*1e-6;
    para.TR=twix.hdr.Phoenix.alTR{1}*1e-6; %s
    Range=0; %twix.hdr.Phoenix.sWipMemBlock.alFree{5};
    NRep=twix.image.NRep;
    Shift=[]; %twix.hdr.Phoenix.sWipMemBlock.alFree{6};
    if(isempty(Shift));Shift=0;end
    para.PC_deg = 0; %(Range/(2*NRep):Range/NRep:Range-Range/(2*NRep))+Shift;
    para.PhaseCycles=1; %deg2rad(para.PC_deg);

    %transmitter
    para.FlipAngle=deg2rad(twix.hdr.Dicom.adFlipAngleDegree); %rad

    %only work for rect pulse
    para.RefVoltage=twix.hdr.Spice.TransmitterReferenceAmplitude;
    try
        para.PulseVoltage=twix.hdr.Phoenix.sTXSPEC.aRFPULSE{1}.flAmplitude;
        para.PulseDur=(para.RefVoltage/para.PulseVoltage)*(0.5*para.FlipAngle/90);
    catch

        para.PulseDur=twix.hdr.Phoenix.sWipMemBlock.alFree{4}*1e-6; %s
        para.PulseVoltage=para.RefVoltage*(0.5/para.PulseDur)*(para.FlipAngle/90);
    end

    % get timings probably only work for VE12U
    para.AcquitionDuration=(twix.image.timestamp(end)-twix.image.timestamp(1))*2.5e-3; %s
    tok=regexp(twix.hdr.Phoenix.tReferenceImage0,'.(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})\d+$','tokens');
    tok=cellfun(@(x)(str2double(x)),tok,'UniformOutput',false);
    para.StudyDateTime= datetime(tok{1});
    para.SeriesDateTime_start=datetime(tok{1}(1:3)) +seconds(twix.image.timestamp(1)*2.5e-3);
    para.SeriesDateTime_end=datetime(tok{1}(1:3))+seconds(twix.image.timestamp(end)*2.5e-3);

    para.SeqDetails=printSeqeunceDetails(twix);

        %resolution
    sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};
    kp=twix.hdr.Phoenix.sKSpace;

    try %if slice oversampling!
        FOV=[sa.dReadoutFOV sa.dPhaseFOV*kp.dPhaseResolution sa.dThickness + sa.dThickness*kp.dSliceOversamplingForDialog]*1e-3; %m
    catch
        FOV=[sa.dReadoutFOV sa.dPhaseFOV*kp.dPhaseResolution  sa.dThickness]*1e-3; %m
    end

    para.MatrixSize=[kp.lBaseResolution  kp.lPhaseEncodingLines kp.lPartitions ];
    para.resolution=FOV./para.MatrixSize; %m
    para.ShortDescription=sprintf('M%d|TR %.0f ms| %.0f deg | %.2f mm | %d rep | %d echoes',twix.hdr.Config.MeasUID,para.TR*1e3,rad2deg(para.FlipAngle),para.resolution(1)*1e3,length(para.PhaseCycles),length(para.TE));

  
end

end