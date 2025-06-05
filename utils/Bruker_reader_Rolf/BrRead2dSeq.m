function im = BrRead2dSeq(path, bits)
% Reads Bruker processed data file.
% Reads the 2dSeq-file path without reading or expecting any additional
% information. The size of the resulting row vector depends only on the
% size of the file.
if nargin < 2
    bits = 16;
end

if nargin < 1
    [f,p] = uigetfile([{'2dseq*', ...
            'Bruker image data files';'*.*','All Files'}], ... 
            'Select image file for reading')
    if f == 0
        arr = -1;
        return;
    end
    path = strcat(p,f);
end
path
%% Open file and start reading
if exist(path,'file')== 2
    [file,errormsg] = fopen(path,'r');
    if file == -1
        disp(errormsg);
        im = -2;
        return
    end
    if bits == 16
        im = fread(file,'int16');
    elseif bits == 32
        im = fread(file,'int32');
    end
    fclose(file);
else 
    im = -2;
end

return;
end
