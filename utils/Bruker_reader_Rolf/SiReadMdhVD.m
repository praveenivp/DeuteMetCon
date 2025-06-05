function Mdh = SiReadMdhVD(file,FilePos,MdhIndex)
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
fseek(file,FilePos,'bof');
%%Define Mdh structure


Mdh.Scan= uint32(0);
Mdh.Time= uint32(0);
Mdh.PMUTime= uint32(0);
Mdh.flags = uint32(0);
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
Mdh.Free1 = uint16(0);
Mdh.Free2 = uint16(0);
Mdh.Free3 = uint16(0);
Mdh.Free4 = uint16(0);
Mdh.IChannel = uint32(0);
Mdh.Data = 1;
%%Read first Mdh
MdhCount = 0;
%Skip first 4 bytes

fseek(file, 4, 'cof');%Should be 4???
%Now read
UID = fread(file,1,'uint32');
Mdh.Scan = fread(file,1,'uint32');
Mdh.Time = fread(file,1,'uint32');
Mdh.PMUTime = fread(file,1,'uint32');
%Mdh.flags = fread(file,1,'uint32');
%Skip next 24 bytes
fseek(file, 20, 'cof');
Mdh.flags = fread(file,1,'uint32');
Mdh.Test = fread(file,1,'uint32');
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
%fseek(file, 74, 'cof');
seg = ftell(file);
fseek(file, 66, 'cof');
Mdh.Free1 = fread(file,1,'uint16');
Mdh.Free2 = fread(file,1,'uint16');
Mdh.Free3 = fread(file,1,'uint16');
Mdh.Free4 = fread(file,1,'uint16');
fseek(file, 40, 'cof');
%fseek(file, 114, 'cof');  %Should be 122!!!
%Mdh.IChannel = fread(file,1,'uint32');
Mdh.Test = fread(file,1,'uint32');
fseek(file, 4, 'cof');  %Should be 122!!!
%Don't read channel header
%type = fread(file,1,'uint8')
%fseek(file, 3, 'cof');
%UID = fread(file,1,'uint32')
%Scan  = fread(file,1,'uint32')
%fseek(file, 12, 'cof');
%ChanID = fread(file,1,'uint16')
%fseek(file, 6, 'cof');
return
end