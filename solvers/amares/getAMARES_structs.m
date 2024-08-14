function [expParam,pk]=getAMARES_structs(mcobj)
%[expParam,pk]=getAMARES_structs(mcobj)
% function to get experimental parameters and initital/bounds/prior
% knowledge struct for AMARES fit

metabolites=mcobj.metabolites;
DMIPara=mcobj.DMIPara;

nMet=length(metabolites);

DW= DMIPara.dwell;
BW=1/DW; %Hz
samples=DMIPara.VectorSize;
timeAxis=0:DW:(samples-1)*DW;
imagingFrequency=DMIPara.ImagingFrequency_MHz;
ppmAxis=linspace(-BW/2,BW/2-BW/samples,samples)/imagingFrequency;

expParam=struct('imagingFrequency', imagingFrequency,...
    'BW', BW,...
    'timeAxis', timeAxis(:), ...
    'dwellTime', DW,...
    'ppmAxis',ppmAxis(:), ...
    'beginTime',DMIPara.AcqDelay_s, ...
    'offset',0,...
    'samples',samples);

phase=zeros(1,nMet);
amplitude=ones(1,nMet);%*max(abs(mobj.img),[],'all');
linewidth=ones(1,nMet).*12;
freq=[metabolites.freq_shift_Hz];
chemShift=freq./imagingFrequency;


initialValues=[];
bounds=[];
priorKnowledge=[];
for i=1:nMet
    initialValues=[initialValues,...
        struct("peakName",metabolites(i).name,"chemShift",chemShift(i),"linewidth",linewidth(i),"amplitude",amplitude(i),"phase",phase(i))];
    bounds=[bounds,...
        struct("peakName",metabolites(i).name,"chemShift",chemShift(i)+[-1 1]*0.2,"linewidth",[2 25],"amplitude",amplitude(i)*[0 Inf],"phase",[0 360],"chemShiftDelta",[],"amplitudeRatio","[]")];
    priorKnowledge=[priorKnowledge,...
        struct("peakName",metabolites(i).name,"multiplet",[],"chemShiftDelta",[],"amplitudeRatio",[],"G_linewidth",[],"G_amplitude",[],"G_phase",1,"RelPhase",[],"G_chemShiftDelta",[],"refPeak",0)];


end
priorKnowledge(1).refPeak=1; %water
pk=struct("bounds",bounds,"initialValues",initialValues,"priorKnowledge",priorKnowledge);


end