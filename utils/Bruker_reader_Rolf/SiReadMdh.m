function Mdh = SiReadMdh(file,FilePos,MdhIndex)
%%Open file
if nargin == 0
    file = -1;
end
if file <0
    path = FilePos;
    FilePos = -1;
    file = fopen(path);
end
%% Read header length
if FilePos<0
    fseek(file, 0, 'bof');
    FilePos = fread(file, 1,'uint32');
end
fseek(file,FilePos,'cof');
%%Define Mdh structure

Mdh.Scan= uint32(0);
Mdh.Time= uint32(0);
Mdh.PMUTime= uint32(0);
Mdh.flags = uint64(0);
Mdh.NSamples=uint16(0);
Mdh.NChannels=uint16(0);
Mdh.ILine=uint16(0);
Mdh.IAcq=uint16(0);
Mdh.ISlice=uint16(0);
Mdh.IPartition=uint16(0);
Mdh.IEcho=uint16(0);
Mdh.IPhase=uint16(0);
Mdh.IMeas=uint16(0);
Mdh.ISet=uint16(0);
Mdh.ISegment=uint16(0);

%Usually not needed
Mdh.IceDim=zeros(5,1,'uint16');
Mdh.CutOffData = zeros(2,1,'uint16');
%--------------------------------------------------

Mdh.lTimeSinceRF = uint32(0);
Mdh.ICenterLine= uint16(0);
Mdh.ICenterPartition= uint16(0);
Mdh.Free1 = uint16(0);
Mdh.Free2 = uint16(0);
Mdh.Free3 = uint16(0);
Mdh.Free4 = uint16(0);
Mdh.IChannel = uint32(0);
Mdh.Data = 1;
%%Read first Mdh
MdhCount = 0;
%Skip first 8 bytes
fseek(file, 8, 'cof');
%Now read
Mdh.Scan = fread(file,1,'uint32');
Mdh.Time = fread(file,1,'uint32');
Mdh.PMUTime = fread(file,1,'uint32');
Mdh.flags = uint64( fread(file,1,'uint64'));
%Skip next 4 bytes
%fseek(file, 4, 'cof');
%Go on reading
Mdh.NSamples = fread(file,1,'uint16');
Mdh.NChannels = fread(file,1,'uint16');
Mdh.ILine = fread(file,1,'uint16');
Mdh.IAcq = fread(file,1,'uint16');
Mdh.ISlice = fread(file,1,'uint16');
Mdh.IPartition = fread(file,1,'uint16');
Mdh.IEcho = fread(file,1,'uint16');
Mdh.IPhase = fread(file,1,'uint16');
Mdh.IMeas = fread(file,1,'uint16');
Mdh.ISet = fread(file,1,'uint16');
Mdh.ISegment = fread(file,1,'uint16');
%Usually not needed
Mdh.IceDim=fread(file,5,'uint16');
Mdh.CutOffData=fread(file,2,'uint16');
%--------------------------------------------------
%fseek(file, 22, 'cof');
fseek(file, 8, 'cof');
Mdh.lTimeSinceRF = fread(file,1,'uint32');
Mdh.ICenterLine = fread(file,1,'uint16');
Mdh.ICenterPartition = fread(file,1,'uint16');
fseek(file, 8, 'cof');
Mdh.Free1 = fread(file,1,'uint16');
Mdh.Free2 = fread(file,1,'uint16');
Mdh.Free3 = fread(file,1,'uint16');
Mdh.Free4 = fread(file,1,'uint16');
fseek(file, 28, 'cof');
Mdh.IChannel = fread(file,1,'uint32');


return
end