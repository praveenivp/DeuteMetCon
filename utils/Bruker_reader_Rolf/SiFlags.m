function flaginds = SiFlags( flags )
% Determines what flags are set in a MdhBitField

%% Flaglabels:
flaglabel{1} = 'AcqEnd';
flaglabel{2} = 'RTFeedback';
flaglabel{3} = 'HPFeedback';
flaglabel{4} = 'Online';
flaglabel{5} = 'Offline';
flaglabel{6} = 'Syncdata';
flaglabel{7} = '';
flaglabel{8} = '';
flaglabel{9} = 'LastScanInConcat';
flaglabel{10} = '';
flaglabel{11} = 'RawDataCorrection';
flaglabel{12} = 'LastScanInMeas';
flaglabel{13} = 'ScanScaleFactor';
flaglabel{14} = '2ndHadamardPulse';
flaglabel{15} = 'RefPhaseStabScan';
flaglabel{16} = 'PhaseStabScan';
flaglabel{17} = '3DFFT';
flaglabel{18} = 'SignRev';
flaglabel{19} = 'PhaseFFT';
flaglabel{20} = 'Swapped';
flaglabel{21} = 'PostSharedLine';
flaglabel{22} = 'PhasCor';
flaglabel{23} = 'PatRefScan';
flaglabel{24} = 'PatRefAndImaScan';
flaglabel{25} = 'Reflect';
flaglabel{26} = 'NoiseAdjScan';
flaglabel{27} = 'ShareNow';
flaglabel{28} = 'LastMeasuredLine';
flaglabel{29} = 'FirstScanInSlice';
flaglabel{30} = 'LastScanInSlice';
flaglabel{31} = 'TREffectiveBegin';
flaglabel{32} = 'TREffectiveEnd';
flaglabel{32} = 'MDSRefPosition';
flaglabel{34} = 'SlcAveraged';
flaglabel{35} = '';
flaglabel{36} = '';
flaglabel{37} = '';
flaglabel{38} = '';
flaglabel{39} = '';
flaglabel{40} = '';
flaglabel{41} = '';
flaglabel{42} = '';
flaglabel{41} = 'FirstScanInBlade';
flaglabel{42} = 'LastScanInBlade';
flaglabel{43} = 'LastBladeInTR';
flaglabel{44} = '';
flaglabel{45} = '';
flaglabel{46} = 'RetroLastphase';
flaglabel{47} = 'RetroEndOfMeas';
flaglabel{48} = 'RepeatThisHeartbeat';
flaglabel{49} = 'RepeatPrevHeartbeat';
flaglabel{50} = 'AbortScanNow';
flaglabel{51} = 'LastHeartBeat';
flaglabel{52} = 'DummyScan';
flaglabel{53} = 'ArrDetDisabled';



lastbit = 53;
if isa(flags, 'uint32')
    lastbit = 32;
end
flaginds = [];
ind = 1;
for cnt = 1:lastbit
    if bitget(flags, cnt)
        display( ['Flag ', num2str(cnt-1),': ',flaglabel{cnt}]);
        flaginds(ind) = cnt-1;
        ind = ind+1;
    end
end

end

