function fid = BrReadFid(path)

% Reads Bruker raw data file.
% Reads the fid-file path without reading or expecting any additional
% information. The size of the resulting row vector depends only on the
% size of the file.

% Written by Rolf Pohmann
% Max Planck Institute for Biological Cybernetics, Tübingen, Germany
path
if nargin < 1
    [f,p] = uigetfile([{'fid*;ser*;rawdata*', ...
            'Bruker raw data files';'*.*','All Files'}], ... 
            'Select fid-file for reading')
    if f == 0
        arr = -1;
        return;
    end
    path = strcat(p,f);
end
%% Open file and start reading
file = fopen(path,'r');
if file == -1
    arr = -2;
    return
end
fid = fread(file,'int32');
fclose(file);
fid = reshape(fid,2,[]);
fid = complex(fid(1,:),fid(2,:));

return;
end
