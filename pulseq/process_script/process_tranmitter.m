%% setup variables
sn='/ptmp/pvalsala/deuterium/EAZW-GUMK/TWIX/';

addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'));
addpath(genpath('/ptmp/pvalsala/Packages/pulseq'));

%% load sequence file
system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
seq.read(fullfile(sn,'../EXPDATA/pulseq/trans_adjust_2H_1min.seq'))

st.dwell_s=seq.getDefinition('dwell');
st.fa_array=seq.getDefinition('fa_array');
st.rf_dur=seq.getDefinition('rf_dur');
st.averages=seq.getDefinition('averages');
st.repetitions=seq.getDefinition('repetitions');
%%
dirst=dir(fullfile(sn,'allData#S94Tuebingen#F5081#M423#D220124#T103144#pulseq_transmitter.dat'));
st.filename=dirst(end).name %load last file
twix=mapVBVD(st.filename);
if(iscell(twix)),twix=twix{end}; end
data=twix.image{''};
  data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions);
Spectrum=squeeze(sos(myfft(data,1),[2 3]));
st.Cfreq=twix.hdr.Dicom.lFrequency; %Hz
st.hdr=twix.hdr;
st.RefVoltage= twix.hdr.Spice.TransmitterReferenceAmplitude;
st.pulseVoltage= (st.fa_array./st.rf_dur).*((st.RefVoltage*1e-3)/pi);


faxis=linspace(-0.5/st.dwell_s,0.5/st.dwell_s,length(Spectrum));


figure,tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
   nexttile() ,plot(faxis,Spectrum),xlabel('frequency[Hz]'),title(st.filename(27:end),'Interpreter','none');
     set(gca,'ColorOrder',parula(size(Spectrum,2)));
       xlim([-200 200])
    nexttile(),plot(st.pulseVoltage,max(Spectrum,[],1),'LineWidth',1.5),xlabel('pulse voltage [V]')
grid on
text(1.1*min(xlim),1.1*min(ylim),sprintf('Pulse duration= %.1f ms\n\n  ',st.rf_dur*1e3),'FontSize',12)
% pulse voltage with maximum amplitude is 90 deg 
[~,max_idx]=max(max(Spectrum,[],1));
optimal_referencevoltage=st.pulseVoltage(max_idx)*(st.rf_dur/0.5e-3)
title(sprintf('current/optimal RefVolatge %.1f V/%.1f V',st.RefVoltage,optimal_referencevoltage))
 fontsize(gcf,'scale',1.3)

 set(gcf,"Position",[323 325 1516 782])