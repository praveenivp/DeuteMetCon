function [Mdhs, Scans, Samples] = SiReadallMdhsVEU(file, DataLength)
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
LastMdh = SiReadMdhVD(file, Currpos+DataLength-192, 0);

fseek(file, Currpos, 'bof');
nMdhs = LastMdh.Scan-1;
fseek(file, Currpos+headerlength, 'bof');
datpos = zeros(floor(nMdhs),1,'uint64');
fprintf('Reading %d MDHs: %3d %% done.',nMdhs, 0);
for cnt=1:nMdhs
    if (round(cnt/1000) == cnt/1000)
        fprintf('\b\b\b\b\b\b\b\b\b\b\b%3d %% done.', round(cnt/nMdhs*100));
    end
    %Skip the data and the Channel headers
    if cnt > 1
        fseek(file, (32+Mdhs(cnt-1).NSamples*2*4)*Mdhs(cnt-1).NChannels, 'cof');
    end
    Mdhs(cnt) = SiReadMdhVD(file, ftell(file),0);
    if cnt == 1
        Mdhs(nMdhs) = Mdhs(1);
    end
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
        Scans.PhaseStab = Scans.PhaseStab + 1
        Samples.PhaseStab = Mdhs(cnt).NSamples;
        if bitget(Mdhs(cnt).flags, 16)  % Phase Stabilisation Reference Scan (Flag 15)
            Scans.PhaseStabRef = cnt;
        end
    else
        Scans.data = Scans.data+1;
        Samples.data = Mdhs(cnt).NSamples;
    end
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b%3d %% done.\n', 100);

fseek(file, Currpos+headerlength, 'bof');
return
end

