function [raw, corr, flags, Mdhs,params] = SiReadRawDataBroken(path,channel,average,removeemptyedges)
% Reads data from interrupted and thus incomplete measurements
% Reads Siemens raw data files in either VB or VE versions
% Called: e.g. [raw, corr, flags, Mdhs] = SiReadRawData(filename,[],1,0);
% removeemptyedges: Removes all lines/partitions/  that are empty at the edges of the k space. Mainly for weighted CSI, where there are usually empty lines at the edges
VD = false;
corr = struct;
flags = 0;
if nargin < 4
    removeemptyedges = 0;
end
if nargin < 3
    average = 1;
end

%%Open file
file = fopen(path);
%% Read header length
tic
fseek(file,0,'eof');
e = ftell(file);
fseek(file,0,'bof');
length = fread(file, 1,'uint32');

%% VD version
if length < 1000
    VD = true;
    % Read number of scans in the file
    % At the beginning of the data file, it states the numer of scans in
    % this file. After that, each scan is described by Scan and Fid
    % Numbers, the offset of the beginning of the scan and the length of
    % this scan, followed by a 128 Byte header. We now only want to scan
    % the last of the scans in a file
    NumScans = fread(file, 1,'uint32');
    for scanNum = 1:NumScans
        ScanNum = fread(file, 1,'uint32');
        FidNum =  fread(file, 1,'uint32');
        Offset = fread(file, 1,'uint64');
        DataLength = fread(file, 1,'uint64');
        if scanNum<NumScans
            fseek(file,128,'cof');
        end
    end
    fseek(file, Offset,'bof');
    
    %Now start reading the header of this scan
    %length = fread(file, 1,'uint32')
    %disp(['Header file length: ',num2str(length)]);
    %fseek(file,length-4,'cof');
    %c = ftell(file)
    if NumScans == 1
        [Mdhs, Scans, Samples] = SiReadallMdhsVEUBroken(file,DataLength);
    else
        [Mdhs, Scans, Samples] = SiReadallMdhsVD(file,DataLength);
    end
    %end  %!!!!
else
    disp(['Header file length: ',num2str(length)]);
    fseek(file,length,'bof');
    c = ftell(file);
    nbytes = e-c;
    cnt = 1;
    [Mdhs, Scans, Samples] = SiReadallMdhs(file,c);
end
toc
%Mdh = SiReadMdh(file,0,0);
%Mdh.Time
%nMdhs = floor((nbytes - 256)/(Mdh.NSamples*2*4+128));
nMdhs = numel(Mdhs);
%% Read header
%linecount = 1;
%while ftell(file) < length
%    line = fgetl(file);
    %disp([num2str(linecount), line]);
    %linecount = linecount+1;
%end
%fseek(file,length,'bof');
%% Read all MDHs to get info about file

%%Read MDH
%fseek(file,128,'cof');
NotAllChannels = 0;
ChannelsToRead = Mdhs(1).NChannels;
SlicesToRead = -1;
if nargin > 1
    if numel(channel)>0
        ChannelsToRead = numel(channel);
    end
end
if ChannelsToRead < Mdhs(1).NChannels
    nMdhs = nMdhs / Mdhs(1).NChannels * ChannelsToRead;
    NotAllChannels = 1;
end
channelcount = 0;
scancount = 0;


% Generate data array
if removeemptyedges == 1
    minLine = min([Mdhs([Mdhs.Data] == 1).ILine]);
    minPar = min([Mdhs([Mdhs.Data] == 1).IPartition]);
    minSeg = min([Mdhs([Mdhs.Data] == 1).ISegment]);
else
    minLine = 0;
    minPar = 0;
    minSeg = 0;
end
NSlices = max([Mdhs.ISlice])+1;
NPhases = max([Mdhs.IPhase])+1;
NLines = max([Mdhs.ILine])+1 - minLine;
NPart = max([Mdhs.IPartition])+1 - minPar;
NMeas = max([Mdhs.IMeas])+1;
NEchos = max([Mdhs.IEcho])+1;
NSets = max([Mdhs.ISet])+1;
NAcqs = max([Mdhs.IAcq])+1;
NSegments = max([Mdhs.ISegment])+1 - minSeg;
NChannels = Mdhs(1).NChannels;
NSamples = Samples.data;
if Scans.PhaseCorr > 0
    corr.Phasecorr.data = complex(zeros(Scans.Phasecorr, Samples.Phasecorr),zeros(Scans.Phasecorr, Samples.Phasecorr));
    corr.Phasecorr.n = 0;
end
if Scans.NoiseAdj > 0
    corr.Noiseadj.data = complex(zeros(Samples.NoiseAdj, Scans.NoiseAdj, NChannels),zeros(Samples.NoiseAdj, Scans.NoiseAdj,NChannels ));
    corr.Noiseadj.n = 0;
end
if Scans.RTFeedback > 0
    corr.RTFeedback.data = complex(zeros(Scans.RTFeedback, Samples.RTFeedback, NChannels),zeros(Scans.RTFeedback, Samples.RTFeedback, NChannels));
    corr.RTFeedback.n = 0;
end
if isfield(Scans,'PhaseStab') == true && Scans.PhaseStab > 0
    corr.PhaseStab.data = complex(zeros(Scans.PhaseStab, Samples.PhaseStab, NChannels),zeros(Scans.PhaseStab, Samples.PhaseStab, NChannels));
    corr.PhaseStab.n = 0;
    corr.PhaseStab.ref = complex(zeros(Samples.PhaseStab,NChannels),zeros(Samples.PhaseStab,NChannels));
end
raw = 0;
flags = 0;
 if average == 0
     raw = complex(zeros(Samples.data, NLines,NPart, NSlices, NEchos, NPhases, NMeas, NSegments, NAcqs, NSets, ChannelsToRead),zeros(Samples.data, NLines,NPart, NSlices, NEchos, NPhases, NMeas, NSegments, NAcqs, NSets, ChannelsToRead));
     flags = zeros(NLines,NPart, NSlices, NEchos, NPhases, NMeas, NSegments, NAcqs, NSets, NChannels,'uint64');
 else
     raw = complex(zeros(NSamples, NLines,NPart, NSlices, NEchos, NPhases, NMeas, NSegments, ChannelsToRead),zeros(Samples.data, NLines,NPart, NSlices, NEchos, NPhases, NMeas, NSegments, ChannelsToRead));
     flags = zeros(NLines,NPart, NSlices, NEchos, NPhases, NMeas, NSegments, NAcqs, NSets, NChannels,'uint64');
 end

% Generate params structure
params.Filename = path;
params.StartTime = Mdhs(1).Time*2.5/1000/60/60;
params.EndTime = Mdhs(end).Time*2.5/1000/60/60;
params.ScanDuration = (Mdhs(end).Time-Mdhs(1).Time)*2.5/1000;
if (NPart == 1)
    if NLines == 1
        params.Dim = 1;
    else
        params.Dim = 2;
    end
else
    params.Dim = 3;
end
params.Matrix = [NSamples,NLines,NPart];
params.Channels = NChannels;



ncorrs = 0;
nnoiseadj = 0;
currchannel=0;
Firstscan = true;
%% Read data
display('Reading data.');
progress_str = '';
if ~VD
for cnt = 1:nMdhs*NChannels
    MdhInd = ceil(cnt/NChannels);
    if (round(cnt/1000) == cnt/1000)
        prevLength      = numel(progress_str);
        progress_str    = sprintf('%d %% done. \n', round(cnt/(nMdhs*NChannels)*100));
        fprintf([repmat('\b',1,prevLength) '%s'],progress_str);
        %display([int2str(round(cnt/(nMdhs*NChannels)*100)),'% done']);
    end
    currchannel = currchannel+1;
    if currchannel > NChannels
        currchannel = 1;
    end
    if NotAllChannels == 1
        channelind = find(channel==currchannel);
        if numel(channelind) == 0
            continue;
        end
    else
        channelind = currchannel;
    end
    %display([num2str(Mdhs(MdhInd).ILine+1),', ', num2str(Mdhs(MdhInd).IPartition+1),', ', num2str(Mdhs(MdhInd).ISlice+1),', ', num2str(Mdhs(MdhInd).IEcho+1),', ', num2str(Mdhs(MdhInd).IPhase+1),', ', num2str(Mdhs(MdhInd).IMeas+1),', ', num2str(Mdhs(MdhInd).ISegment+1),', ', num2str(Mdhs(MdhInd).IAcq+1),', ', num2str(Mdhs(MdhInd).ISet+1),', ', num2str(currchannel)]);
    fseek(file,128,'cof');
    fid = fread(file, Mdhs(MdhInd).NSamples*2, 'float32');
    %% Convert to complex
    fid = reshape(fid,2,[]);
    dat = complex(fid(1,:)',fid(2,:)');
    
    % For EPI-scans, the first three acquisitions are phase correction
    % data, marked by bit 22 in the Flags
    if bitget(Mdhs(MdhInd).flags,22)
        if isfield(corr, 'Phasecorr')
            corr.Phasecorr.n = corr.Phasecorr.n+1;
        else
            corr.Phasecorr.n = 1;
        end
        corr.Phasecorr.data(corr.Phasecorr.n,1:Mdhs(MdhInd).NSamples) = dat;
        corr.Phasecorr.Mdh(corr.Phasecorr.n) = Mdhs(MdhInd);
        
    elseif bitget(Mdhs(MdhInd).flags,26)    %Noise adjust scans (flag 25)
        if isfield(corr, 'Noiseadj')
            corr.Noiseadj.n = corr.Noiseadj.n+1;
        else
            corr.Noiseadj.n = 1;
        end
        corr.Noiseadj.data(1:Mdhs(MdhInd).NSamples,corr.Noiseadj.n,currchannel) = dat;
        corr.Noiseadj.Mdh(corr.Noiseadj.n) = Mdhs(MdhInd);
        
    elseif bitget(Mdhs(MdhInd).flags,2)    %Realtime feedback scans (flag 1)
        if isfield(corr, 'RTFeedback')
            corr.RTFeedback.n = corr.RTFeedback.n+1;
        else
            corr.RTFeedback.n = 1;
        end
        corr.RTFeedback.data(ceil(corr.RTFeedback.n/NChannels),1:Mdhs(MdhInd).NSamples, channelind) = dat;
        corr.RTFeedback.Mdh(ceil(corr.RTFeedback.n/NChannels),channelind) = Mdhs(MdhInd);
    elseif bitget(Mdhs(MdhInd).flags,17)    %PhaseStab scan (flag 16)
        if isfield(corr, 'PhaseStab')
            corr.PhaseStab.n = corr.PhaseStab.n+1;
        else
            corr.PhaseStab.n = 1;
        end
        corr.PhaseStab.data(ceil(corr.PhaseStab.n/NChannels),1:Mdhs(MdhInd).NSamples, channelind) = dat;
        corr.PhaseStab.Mdh(ceil(corr.PhaseStab.n/NChannels),channelind) = Mdhs(MdhInd);
        if bitget(Mdhs(MdhInd).flags,16)    %PhaseStab reference scan (flag 15)
            corr.PhaseStab.ref(:,channelind) = dat;
        end
    elseif bitget(Mdhs(MdhInd).flags,16)    %PhaseStab scan (flag 15)
        if isfield(corr, 'PhaseStabRef')
            corr.PhaseStabRef.n = corr.PhaseStabRef.n+1;
        else
            corr.PhaseStabRef.n = 1;
        end
        corr.PhaseStab.ref(:,channelind) = dat;

     else
          if average == 0   %line, partition, slice, echo, phase, meas, segment, acq, set
                 raw(:,Mdhs(MdhInd).ILine+1-minLine, Mdhs(MdhInd).IPartition+1-minPar,Mdhs(MdhInd).ISlice+1,Mdhs(MdhInd).IEcho+1,Mdhs(MdhInd).IPhase+1,Mdhs(MdhInd).IMeas+1,Mdhs(MdhInd).ISegment+1-minSeg,Mdhs(MdhInd).IAcq+1,Mdhs(MdhInd).ISet+1,currchannel)=dat;
                 flags(Mdhs(MdhInd).ILine+1-minLine, Mdhs(MdhInd).IPartition+1-minPar,Mdhs(MdhInd).ISlice+1,Mdhs(MdhInd).IEcho+1,Mdhs(MdhInd).IPhase+1,Mdhs(MdhInd).IMeas+1,Mdhs(MdhInd).ISegment+1-minSeg,Mdhs(MdhInd).IAcq+1,Mdhs(MdhInd).ISet+1,currchannel)=Mdhs(MdhInd).flags;
          else
                 if Mdhs(MdhInd).ISet == 0 && Mdhs(MdhInd).IAcq == 0
                     raw(:,Mdhs(MdhInd).ILine+1-minLine, Mdhs(MdhInd).IPartition+1-minPar,Mdhs(MdhInd).ISlice+1,Mdhs(MdhInd).IEcho+1,Mdhs(MdhInd).IPhase+1,Mdhs(MdhInd).IMeas+1,Mdhs(MdhInd).ISegment+1-minSeg,currchannel)=dat;
                     flags(Mdhs(MdhInd).ILine+1-minLine, Mdhs(MdhInd).IPartition+1-minPar,Mdhs(MdhInd).ISlice+1,Mdhs(MdhInd).IEcho+1,Mdhs(MdhInd).IPhase+1,Mdhs(MdhInd).IMeas+1,Mdhs(MdhInd).ISegment+1-minSeg,Mdhs(MdhInd).IAcq+1,Mdhs(MdhInd).ISet+1,currchannel)=Mdhs(MdhInd).flags;
                 else
                     raw(:,Mdhs(MdhInd).ILine+1-minLine, Mdhs(MdhInd).IPartition+1-minPar,Mdhs(MdhInd).ISlice+1,Mdhs(MdhInd).IEcho+1,Mdhs(MdhInd).IPhase+1,Mdhs(MdhInd).IMeas+1,Mdhs(MdhInd).ISegment+1-minSeg,currchannel)=raw(:,Mdhs(MdhInd).ILine+1, Mdhs(MdhInd).IPartition+1,Mdhs(MdhInd).ISlice+1,Mdhs(MdhInd).IEcho+1,Mdhs(MdhInd).IPhase+1,Mdhs(MdhInd).IMeas+1,Mdhs(MdhInd).ISegment+1-minSeg,currchannel)+dat;
                     flags(Mdhs(MdhInd).ILine+1-minLine, Mdhs(MdhInd).IPartition+1-minPar,Mdhs(MdhInd).ISlice+1,Mdhs(MdhInd).IEcho+1,Mdhs(MdhInd).IPhase+1,Mdhs(MdhInd).IMeas+1,Mdhs(MdhInd).ISegment+1-minSeg,Mdhs(MdhInd).IAcq+1,Mdhs(MdhInd).ISet+1,currchannel)=Mdhs(MdhInd).flags;
                 end
          end
    end
    
end
else  % if VD
    for cnt = 1:nMdhs
        MdhInd = cnt; %ceil(cnt/NChannels);
        if (round(cnt/1000) == cnt/1000)
            display([int2str(round(cnt/nMdhs*100)),'% done']);
        end
        fseek(file,192,'cof');
        %fseek(file,pos(cnt),'bof');
        for currchannel = 1:NChannels
            type = fread(file,1,'uint8');
            fseek(file, 3, 'cof');
            UID = fread(file,1,'uint32');
            Scan  = fread(file,1,'uint32');
            fseek(file, 12, 'cof');
            ChanID = fread(file,1,'uint16');
            fseek(file, 6, 'cof');  
            fid = fread(file, Mdhs(MdhInd).NSamples*2, 'float32');
            %% Convert to complex
            fid = reshape(fid,2,[]);
            dat = complex(fid(1,:),fid(2,:)).';
            % For EPI-scans, the first three acquisitions are phase correction
            % data, marked by bit 22 in the Flags
            if bitget(Mdhs(MdhInd).flags,22)
                ncorrs = ncorrs+1;
                corr.data(ncorrs,1:Mdhs(MdhInd).NSamples) = dat;
                corr.Mdh(ncorrs) = Mdhs(cnt);
            elseif bitget(Mdhs(MdhInd).flags,26)
                if currchannel == 1
                    if isfield(corr, 'Noiseadj')
                        corr.Noiseadj.n = corr.Noiseadj.n+1;
                    else
                        corr.Noiseadj.n = 1;
                    end
                    corr.Noiseadj.Mdh(corr.Noiseadj.n) = Mdhs(MdhInd);
                end
                corr.Noiseadj.data(1:Mdhs(MdhInd).NSamples,corr.Noiseadj.n,currchannel) = dat;
            else
                if average == 0
                        raw(:,Mdhs(cnt).ILine+1-minLine, Mdhs(cnt).IPartition+1-minPar,Mdhs(cnt).ISlice+1,Mdhs(cnt).IEcho+1,Mdhs(cnt).IPhase+1,Mdhs(cnt).IMeas+1,Mdhs(cnt).ISegment+1-minSeg,Mdhs(cnt).IAcq+1,Mdhs(cnt).ISet+1,currchannel)=dat;
                        flags(Mdhs(cnt).ILine+1-minLine, Mdhs(cnt).IPartition+1-minPar,Mdhs(cnt).ISlice+1,Mdhs(cnt).IEcho+1,Mdhs(cnt).IPhase+1,Mdhs(cnt).IMeas+1,Mdhs(cnt).ISegment+1-minSeg,Mdhs(cnt).IAcq+1,Mdhs(cnt).ISet+1,currchannel)=Mdhs(cnt).flags;
                else
                    if Mdhs(cnt).ISet == 0 && Mdhs(cnt).IAcq == 0
                        raw(:,Mdhs(cnt).ILine+1-minLine, Mdhs(cnt).IPartition+1-minPar,Mdhs(cnt).ISlice+1,Mdhs(cnt).IEcho+1,Mdhs(cnt).IPhase+1,Mdhs(cnt).IMeas+1,Mdhs(cnt).ISegment+1-minSeg,currchannel)=dat;
                        flags(Mdhs(cnt).ILine+1-minLine, Mdhs(cnt).IPartition+1-minPar,Mdhs(cnt).ISlice+1,Mdhs(cnt).IEcho+1,Mdhs(cnt).IPhase+1,Mdhs(cnt).IMeas+1,Mdhs(cnt).ISegment+1-minSeg,Mdhs(cnt).IAcq+1,Mdhs(cnt).ISet+1,currchannel)=Mdhs(cnt).flags;
                    else
                        raw(:,Mdhs(cnt).ILine+1-minLine, Mdhs(cnt).IPartition+1-minPar,Mdhs(cnt).ISlice+1,Mdhs(cnt).IEcho+1,Mdhs(cnt).IPhase+1,Mdhs(cnt).IMeas+1,Mdhs(cnt).ISegment+1-minSeg,currchannel)=raw(:,Mdhs(cnt).ILine+1-minLine, Mdhs(cnt).IPartition+1-minPar,Mdhs(cnt).ISlice+1,Mdhs(cnt).IEcho+1,Mdhs(cnt).IPhase+1,Mdhs(cnt).IMeas+1,Mdhs(cnt).ISegment+1-minSeg,currchannel)+dat;
                        flags(Mdhs(cnt).ILine+1-minLine, Mdhs(cnt).IPartition+1-minPar,Mdhs(cnt).ISlice+1,Mdhs(cnt).IEcho+1,Mdhs(cnt).IPhase+1,Mdhs(cnt).IMeas+1,Mdhs(cnt).ISegment+1-minSeg,Mdhs(cnt).IAcq+1,Mdhs(cnt).ISet+1,currchannel)=Mdhs(cnt).flags;
                    end
                end
            end
        end
    end
end
    toc

return
end
