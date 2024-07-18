function pk=PriorKnowledge_DMI(amares_struct)
   



%% initial values
if(nargin==0)
initialValues=[struct("peakName",{'D20'},"chemShift",0,"linewidth",12,"amplitude",1,"phase",23),...
    struct("peakName",{'Glc'},"chemShift",-1.01,"linewidth",24,"amplitude",1,"phase",23),...
    struct("peakName",{'Glx'},"chemShift",-2.5,"linewidth",24,"amplitude",1,"phase",23+90)];


% set bounds
bounds=[struct("peakName",{"D20"},"chemShift",[-Inf Inf],"linewidth",[0 Inf],"amplitude",[0 Inf],"phase",[0 360],"chemShiftDelta",[],"amplitudeRatio","[]"),...
    struct("peakName",{'Glc'},"chemShift",[-2 0],"linewidth",[0 Inf],"amplitude",[0 Inf],"phase",[0 360],"chemShiftDelta",[],"amplitudeRatio","[]"),...
    struct("peakName",{'Glx'},"chemShift",[-Inf Inf],"linewidth",[0 Inf],"amplitude",[0 Inf],"phase",[0 360],"chemShiftDelta",[],"amplitudeRatio","[]"),...
    ];

priorKnowledge=[struct("peakName","D20","multiplet",[],"chemShiftDelta",[],"amplitudeRatio",[],"G_linewidth",[],"G_amplitude",[],"G_phase",1,"RelPhase",[],"G_chemShiftDelta",[],"refPeak",1),...
    struct("peakName","Glc","multiplet",[],"chemShiftDelta",[],"amplitudeRatio",[],"G_linewidth",[],"G_amplitude",[],"G_phase",1,"RelPhase",[],"G_chemShiftDelta",[],"refPeak",1),...
    struct("peakName","Glx","multiplet",[],"chemShiftDelta",[],"amplitudeRatio",[],"G_linewidth",[],"G_amplitude",[],"G_phase",1,"RelPhase",[],"G_chemShiftDelta",[],"refPeak",1),...
    ];

else
% initialValues=[struct("peakName",{'D20'},"chemShift",amares_struct.chemShift(1),"linewidth",amares_struct.linewidth(1),"amplitude",amares_struct.amplitude(1),"phase",amares_struct.phase(1)),...
%     struct("peakName",{'Glc'},"chemShift",amares_struct.chemShift(2),"linewidth",amares_struct.linewidth(2),"amplitude",amares_struct.amplitude(2),"phase",amares_struct.phase(2)),...
%     struct("peakName",{'Glx'},"chemShift",amares_struct.chemShift(3),"linewidth",amares_struct.linewidth(3),"amplitude",amares_struct.amplitude(3),"phase",amares_struct.phase(3))];

initialValues=[];
bounds=[];
priorKnowledge=[];
for i=1:length(amares_struct.chemShift)
    initialValues=[initialValues,...
        struct("peakName",amares_struct.peakName(i),"chemShift",amares_struct.chemShift(i),"linewidth",amares_struct.linewidth(i),"amplitude",amares_struct.amplitude(i),"phase",amares_struct.phase(i))];
bounds=[bounds,...
    struct("peakName",amares_struct.peakName(i),"chemShift",amares_struct.chemShift(i)+[-1 1]*2,"linewidth",[0 20],"amplitude",amares_struct.amplitude(i)*[-1 1]*100,"phase",[-180 180],"chemShiftDelta",[],"amplitudeRatio","[]")];
priorKnowledge=[priorKnowledge,...
    struct("peakName",amares_struct.peakName{i},"multiplet",[],"chemShiftDelta",[],"amplitudeRatio",[],"G_linewidth",[],"G_amplitude",[],"G_phase",1,"RelPhase",[],"G_chemShiftDelta",[],"refPeak",1)];


end
end


pk=struct("bounds",bounds,"initialValues",initialValues,"priorKnowledge",priorKnowledge);


end