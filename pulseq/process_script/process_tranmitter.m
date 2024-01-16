%% load sequence file
system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
seq.read('transmitter.seq')

st.dwell_s=seq.getDefinition('dwell');
st.fa_array=seq.getDefinition('fa');
st.rf_dur=seq.getDefinition('rf_dur');
st.averages=seq.getDefinition('averages');
st.repetitions=seq.getDefinition('repetitions');
%%
dirst=dir('*transmitter*.dat');
st.filename=dirst(end).name %load last file
twix=mapVBVD(st.filename);
twix=twix{end};
data=twix.image{''};
  data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions);
Spectrum=squeeze(sos(myfft(data,1),[2 3]));
st.Cfreq=twix.hdr.Dicom.lFrequency; %Hz
st.hdr=twix.hdr;
st.RefVoltage= twix.hdr.Spice.TransmitterReferenceAmplitude;
st.pulseVoltage= (st.fa_array./st.rf_dur).*((st.RefVoltage*1e-3)/pi);


faxis=linspace(-0.5/st.dwell_s,0.5/st.dwell_s,length(Spectrum));


figure
    subplot(121),plot(faxis,Spectrum),xlabel('Freq [Hz]'),title(st.filename,'Interpreter','none');
       xlim([-200 200])
    subplot(122),plot(st.pulseVoltage,max(Spectrum,[],1)),xlabel('pulse Voltage [V]')
% subplot(122),plot(rad2deg(st.fa_array),max(Spectrum,[],1))

% pulse voltage with maximum amplitude is 90 deg 
[~,max_idx]=max(max(Spectrum,[],1));
optimal_referencevoltage=st.pulseVoltage(max_idx)*(st.rf_dur/0.5e-3)
title(sprintf('current/optimal RefVolatge %.1f V/%.1f V',st.RefVoltage,optimal_referencevoltage))
