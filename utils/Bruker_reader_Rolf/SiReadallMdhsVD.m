function [Mdhs, Scans, Samples] = SiReadallMdhsVD(file, DataLength)
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
Scans.NoiseAdj = 0;
Scans.RTFeedback = 0;
Samples.data = 0;
Samples.PhaseCorr = 0;
Samples.NoiseAdj = 0;
Samples.RTFeedback = 0;
%% Read header length
Currpos = ftell(file);
headerlength = fread(file, 1,'uint32');

%Find number of Mdhs
%fseek(file, FilePos+DataLength-128-32-192, 'bof');
LastMdh = SiReadMdhVD(file, Currpos+DataLength-128-32-192, 0);
%LastMdh = SiReadMdhVD(file, Currpos+DataLength-192, 0);

fseek(file, Currpos, 'bof');
nMdhs = LastMdh.Scan-1;
cnt = 1;
Mdh = SiReadMdhVD(file,Currpos+headerlength,0);
Mdhs(1) = Mdh;
display(['Reading ',int2str(nMdhs),' MDHs'])
datpos = zeros(floor(nMdhs),1,'uint64');
datpos(1) = ftell(file);
for cnt=2:nMdhs
    if (round(cnt/1000) == cnt/1000)
        display([int2str(round(cnt/nMdhs*100)),'% done']);
    end
    %Skip the data and the Channel headers
if Mdhs(cnt-1).NSamples < 4
    display(Mdhs(cnt-1).NSamples)
end
    fseek(file, (32+Mdhs(cnt-1).NSamples*2*4)*Mdh.NChannels, 'cof');
    Mdhs(cnt) = SiReadMdhVD(file, ftell(file),0);
    %Read next Mdh
    %Skip first 8 bytes
    %fseek(file, 8, 'cof');
    %Now read
    %Mdhs(cnt).Scan = fread(file,1,'uint32');
    %Mdhs(cnt).Time = fread(file,1,'uint32');
    %Mdhs(cnt).PMUTime = fread(file,1,'uint32');
    %Mdhs(cnt).flags = fread(file,1,'uint32');
    %Flags:
    % Bit 22: Phase correction scan (EPI)
    % Bit 25: Line reflection (EPI)
    % Bit 29: First Scan in Slice
    % Bit 30: Last Scan in Slice
    %Skip next 4 bytes
    %fseek(file, 4, 'cof');
    %Go on reading
    %Mdhs(cnt).NSamples = fread(file,1,'uint16');
    %Mdhs(cnt).NChannels = fread(file,1,'uint16');
    %Mdhs(cnt).ILine = fread(file,1,'uint16');
    %Mdhs(cnt).IAcq = fread(file,1,'uint16');
    %Mdhs(cnt).ISlice = fread(file,1,'uint16');
    %Mdhs(cnt).IPartition = fread(file,1,'uint16');
    %Mdhs(cnt).IEcho = fread(file,1,'uint16');
    %Mdhs(cnt).IPhase = fread(file,1,'uint16');
    %Mdhs(cnt).IMeas = fread(file,1,'uint16');
    %Mdhs(cnt).ISet = fread(file,1,'uint16');
    %Mdhs(cnt).ISegment = fread(file,1,'uint16');
    %fseek(file, 38, 'cof');
    %Mdhs(cnt).Free1 = fread(file,1,'uint16');
    %Mdhs(cnt).Free2 = fread(file,1,'uint16');
    %Mdhs(cnt).Free3 = fread(file,1,'uint16');
    %Mdhs(cnt).Free4 = fread(file,1,'uint16');
    %fseek(file, 28, 'cof');
    %Mdhs(cnt).IChannel = fread(file,1,'uint32')-Channel0;
    datpos(cnt) = ftell(file);
end
fseek(file, Currpos+headerlength, 'bof');
Scans.data = nMdhs;
Samples.data = Mdhs(1).NSamples;
return
end

