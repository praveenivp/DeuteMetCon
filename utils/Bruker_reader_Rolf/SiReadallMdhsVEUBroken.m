function [Mdhs, Scans, Samples] = SiReadallMdhsVEUBroken(file, DataLength)
% Reads Mdhs for the case that a measurement was interrupted and thus data is incomplete
%%Open file
if nargin == 0
    file = -1;
end

if ischar(file) 
    file = fopen(file);
end
Scans.data = 0;
Scans.PhaseCorr = 0;
Scans.NoiseAdj = 0;
Scans.RTFeedback = 0;
Scans.PhaseStab = 0;
Samples.data = 0;
Samples.PhaseCorr = 0;
Samples.NoiseAdj = 0;
Samples.RTFeedback = 0;
%% Read header length
Currpos = ftell(file);
headerlength = fread(file, 1,'uint32');

%Find number of Mdhs
%fseek(file, FilePos+DataLength-128-32-192, 'bof');
%LastMdh = SiReadMdhVD(file, Currpos+DataLength-128-32-192, 0);
%LastMdh = SiReadMdhVD(file, Currpos+DataLength-192, 0);
fseek(file,0,'eof');
fend = ftell(file);
fseek(file,Currpos+headerlength, 'bof');
flength = fend - ftell(file);
display(['Reading  MDHs'])
cnt = 1;
Mdhs(cnt) = SiReadMdhVD(file, ftell(file),0);
defMdhs = floor(1.5*flength/(192+(32+Mdhs(1).NSamples*2*4)*Mdhs(1).NChannels));  % first guess for the number of Mdhs
Mdhs(defMdhs) = Mdhs(1);
fileend = 0;
while ~feof(file) && ~fileend   
    if (round(cnt/1000) == cnt/1000)
        fprintf('%d Mdhs read (approximately %d %%)\n',cnt,defMdhs);
    end
    %Skip the data and the Channel headers
    fseek(file, (32+Mdhs(cnt).NSamples*2*4)*Mdhs(cnt).NChannels, 'cof');
    cnt = cnt+1;
    %fprintf('Scan %d of %d\n',cnt,defMdhs);

    Mdhs(cnt) = SiReadMdhVD(file, ftell(file),0);
    if numel(Mdhs(cnt).Free4)== 0
        fileend = 1;
        cnt = cnt-1;
    else
        datpos(cnt) = ftell(file);
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
            Scans.PhaseStab = Scans.PhaseStab + 1;
            Samples.PhaseStab = Mdhs(cnt).NSamples;
            if bitget(Mdhs(cnt).flags, 16)  % Phase Stabilisation Reference Scan (Flag 15)
                Scans.PhaseStabRef = cnt;
            end
        else
            Scans.data = Scans.data+1;
            Samples.data = Mdhs(cnt).NSamples;
        end
    end
end
Mdhs = Mdhs(1:cnt);
fseek(file, Currpos+headerlength, 'bof');
return
end

