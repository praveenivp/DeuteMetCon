function [Mdhs, Scans, Samples] = SiReadallMdhs(file,FilePos)
%%Open file
if nargin == 0
    file = -1;
end
if file <0
    path = FilePos;
    FilePos = -1;
    file = fopen(path);
end
Scans.data = 0;
Scans.PhaseCorr = 0;
Scans.PhaseStab = 0;
Scans.NoiseAdj = 0;
Scans.RTFeedback = 0;
Scans.PhaseStabRef = -1;
Samples.data = 0;
Samples.PhaseCorr = 0;
Samples.NoiseAdj = 0;
Samples.RTFeedback = 0;
Samples.PhaseStab = 0;
%% Read header length
if FilePos<0
    fseek(file, 0, 'bof');
    FilePos = fread(file, 1,'uint32');
end
Currpos = ftell(file);
fseek(file,FilePos,'bof');
%Calculate approximate number of Mdhs
%Read the last Mdh
c = ftell(file);
fseek(file, -(128+16*2*4),'eof');
LastMdh = SiReadMdh(file,0,0);
fseek(file,c,'bof');
nMdhs = LastMdh.Scan-1;
cnt = 1;
Mdh = SiReadMdh(file,0,0);
Mdhs(1) = Mdh;
Mdhs(nMdhs) = Mdh;
Channel0 = Mdh.IChannel;
display(['Reading ',int2str(nMdhs),' MDHs'])
fseek(file,c,'bof');
NSamples = Mdh.NSamples;
%for cnt=2:nMdhs
progress_str = '';
while NSamples > 16 
    if (round(cnt/1000) == cnt/1000)
        prevLength      = numel(progress_str);
        progress_str    = sprintf('%d %% done. \n', round(cnt/nMdhs*100));
        fprintf([repmat('\b',1,prevLength) '%s'],progress_str);
        %display([int2str(round(cnt/nMdhs*100)),'% done']);
    end
    if cnt == 427
    cnt
    end
    %Read next Mdh
    %Skip first 8 bytes
    fseek(file, 8, 'cof');
    %Now read
    Scan = fread(file,1,'uint32');
    Tim = fread(file,1,'uint32');
    PMUTime = fread(file,1,'uint32');
    flags = uint64(fread(file,1,'uint64'));
    %Flags:
    % Bit 4: Acquired scan
    % Bit 9: Last Scan in ??
    % Bit 12: Last Scan in ??
    % Bit 19: Last Scan in ??
    % Bit 22: Phase correction scan (EPI)
    % Bit 23: Calibration region for GRAPPA
    % Bit 24: No GRAPPA calibration scan 
    % Bit 25: Line reflection (EPI)
    % Bit 26: First Scan in ??
    % Bit 29: First Scan in Slice
    % Bit 30: Last Scan in Slice
    %Skip next 4 bytes
    %fseek(file, 4, 'cof');
    %Go on reading
    NSamples = fread(file,1,'uint16');
    if NSamples > 16
        Mdhs(cnt).Scan = Scan;
        Mdhs(cnt).Time = Tim;
        Mdhs(cnt).PMUTime = PMUTime;
        Mdhs(cnt).flags = flags;
        Mdhs(cnt).NSamples = NSamples;
        Mdhs(cnt).NChannels = fread(file,1,'uint16');
        Mdhs(cnt).ILine = fread(file,1,'uint16');
        Mdhs(cnt).IAcq = fread(file,1,'uint16');
        Mdhs(cnt).ISlice = fread(file,1,'uint16');
        Mdhs(cnt).IPartition = fread(file,1,'uint16');
        Mdhs(cnt).IEcho = fread(file,1,'uint16');
        Mdhs(cnt).IPhase = fread(file,1,'uint16');
        Mdhs(cnt).IMeas = fread(file,1,'uint16');
        Mdhs(cnt).ISet = fread(file,1,'uint16');
        Mdhs(cnt).ISegment = fread(file,1,'uint16');
        %Usually not needed
        Mdhs(cnt).IceDim=fread(file,5,'uint16');
        Mdhs(cnt).CutOffData=fread(file,2,'uint16');
        %--------------------------------------------------
        fseek(file, 8, 'cof');
        %fseek(file, 22, 'cof');
        Mdhs(cnt).lTimeSinceRF = fread(file,1,'uint32');
        Mdhs(cnt).ICenterLine = fread(file,1,'uint16');
        Mdhs(cnt).ICenterPartition = fread(file,1,'uint16');
        fseek(file, 8, 'cof');
        Mdhs(cnt).Free1 = fread(file,1,'uint16');
        Mdhs(cnt).Free2 = fread(file,1,'uint16');
        Mdhs(cnt).Free3 = fread(file,1,'uint16');
        Mdhs(cnt).Free4 = fread(file,1,'uint16');
        fseek(file, 28, 'cof');
        Mdhs(cnt).IChannel = fread(file,1,'uint32')-Channel0;
    if bitget(Mdhs(cnt).flags,22)
        Scans.PhaseCorr = Scans.PhaseCorr+1;
        Samples.PhaseCorr = Mdhs(cnt).NSamples;
    elseif bitget(Mdhs(cnt).flags,26)    %Noise adjust scans (flag 25)
        Scans.NoiseAdj = Scans.NoiseAdj+1;
        Samples.NoiseAdj = Mdhs(cnt).NSamples;
        Mdhs(cnt).Data = 0;                                                 % Data = 0 means that this scan is not part of the final image (e.g. noise adjustment scan)
    elseif bitget(Mdhs(cnt).flags,2)    %Realtime feedback scans (flag 1)
        Scans.RTFeedback = Scans.RTFeedback+1;
        Samples.RTFeedback = Mdhs(cnt).NSamples;
    elseif bitget(Mdhs(cnt).flags,15)   %Phase Stabilisation Scan (flag 14)
        Scans.PhaseStab = Scans.PhaseStab + 1
        Samples.PhaseStab = Mdhs(cnt).NSamples;
        if bitget(Mdhs(cnt).flags, 16)  % Phase Stabilisation Reference Scan (Flag 15)
            Scans.PhaseStabRef = cnt;
        end
    else
        Scans.data = Scans.data+1;
        Samples.data = Mdhs(cnt).NSamples;
    end
    %Skip the data
    %fseek(file, Mdhs(cnt-1).NSamples*2*4, 'cof');
    %Skip NChannels times the data plus NChannels-1 times the mdh
    fseek(file, Mdhs(cnt).NSamples*2*4*Mdhs(cnt).NChannels + 128*(Mdhs(cnt).NChannels-1), 'cof');  
    NSamples = Mdhs(cnt).NSamples;
    cnt = cnt + 1;
    end
end
fseek(file, Currpos, 'bof');



return
end

