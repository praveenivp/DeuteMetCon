function [expParam,pk]=getAMARES_structs_T1inv(twix,metabolites,pulseq_st)
%[expParam,pk]=getAMARES_structs_T1inv(twix,metabolites)
% function to get experimental parameters and initital/bounds/prior
% knowledge struct for AMARES fit: allows negative amplitude for inversion
% recovery data

% DMIPara=getDMIPara(twix);

nMet=length(metabolites);

DW= pulseq_st.dwell_s;
BW=1/DW; %Hz
samples=pulseq_st.VectorSize;
timeAxis=0:DW:(samples-1)*DW;
imagingFrequency=twix.hdr.Dicom.lFrequency*1e-6; %MHz
ppmAxis=calcFreqAxis(DW,samples)/imagingFrequency;

expParam=struct('imagingFrequency', imagingFrequency,...
    'BW', BW,...
    'timeAxis', timeAxis(:), ...
    'dwellTime', DW,...
    'ppmAxis',ppmAxis(:), ...
    'beginTime',pulseq_st.AcqDelay_s, ...
    'offset',0,...
    'samples',samples);

phase=zeros(1,nMet);
amplitude=ones(1,nMet);%*max(abs(mobj.img),[],'all');
linewidth=ones(1,nMet).*10;
freq=[metabolites.freq_shift_Hz];
chemShift=freq./imagingFrequency;


initialValues=[];
bounds=[];
priorKnowledge=[];
for i=1:nMet
    initialValues=[initialValues,...
        struct("peakName",metabolites(i).name,"chemShift",chemShift(i),"linewidth",linewidth(i),"amplitude",amplitude(i),"phase",phase(i))];
    bounds=[bounds,...
        struct("peakName",metabolites(i).name,"chemShift",chemShift(i)+[-1 1]*0.1,"linewidth",[6 30],"amplitude",amplitude(i)*[-Inf Inf],"phase",[0 360],"chemShiftDelta",[],"amplitudeRatio","[]")];
    priorKnowledge=[priorKnowledge,...
        struct("peakName",metabolites(i).name,"multiplet",[],"chemShiftDelta",[],"amplitudeRatio",[],"G_linewidth",[],"G_amplitude",[],"G_phase",1,"RelPhase",[],"G_chemShiftDelta",[],"refPeak",0)];


end
priorKnowledge(1).refPeak=1; %water
pk=struct("bounds",bounds,"initialValues",initialValues,"priorKnowledge",priorKnowledge);


end