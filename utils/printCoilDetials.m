function CoilSelectMap=printCoilDetials(filename)
% CoilSelectMap=printCoilDetials(filename)
% function to extract raw Data correction and FFT scale factors from twix


twix=mapVBVD(filename);
if(iscell(twix)), twix=twix{end};end
sRx=twix.hdr.Phoenix.sCoilSelectMeas.aRxCoilSelectData{1};
for ii=1:length(sRx.aFFT_SCALE)
st.aFFT_SCALE(ii)=sRx.aFFT_SCALE{ii}.flFactor;
st.lADCChannel(ii)=sRx.aFFT_SCALE{ii}.lADCChannel;
st.tElement(ii)=sRx.asList{ii}.sCoilElementID.tElement;
st.lADCChannel2(ii)=sRx.asList{ii}.lADCChannelConnected;
end

textHeader=gettextHeader(filename);
nCha=twix.image.NCha;        
findThis = '{\s*{\s*{\s*"[^"]+"[^\n]+';
pStart = regexp(textHeader, findThis, 'start');
if length(pStart) > 0
   pEnd = regexp(textHeader(pStart(1):end), '}[\n\s]*}[\n\s]*}', 'end');
   if length(pEnd) > 0
      allCoilsInHeaderCell = textHeader(pStart(1):pStart(1)+pEnd(1)+1);
      findThis = '{\s*{\s*"(?<name>[^"]+)"\s*}\s*{\s*(?<fft>[\d\.]+)\s*}\s*{\s*(?<re>[\d\.-]+)\s*}\s*{\s*(?<im>[\d\.-]+)\s*}\s*}';
      CoilStructArray = regexp(allCoilsInHeaderCell, findThis,'names');
      if length(CoilStructArray) == nCha
         for c=1:nCha
            CoilSelectMap(c).tName = CoilStructArray(c).name;
            CoilSelectMap(c).fftScale = sscanf(CoilStructArray(c).fft,'%f');
            CoilSelectMap(c).rawDataCorrectionFactor = complex(sscanf(CoilStructArray(c).re,'%f'),sscanf(CoilStructArray(c).im,'%f'));
            CoilSelectMap(c).lADCChannel=st.lADCChannel(c);
            CoilSelectMap(c).tName2=st.tElement(c);
            CoilSelectMap(c).fftScale2=st.aFFT_SCALE(c);

         end
      end
   end
end

end

function textHeader=gettextHeader(filename)

fid = fopen(filename, 'r', 'ieee-le');
% start of actual measurement data (sans header)
fseek(fid,0,'bof');

firstInt  = fread(fid, 1, 'uint32');
secondInt = fread(fid, 1, 'uint32');

% lazy software version check (VB or VD?)
if and(firstInt < 10000, secondInt <= 64)
    version = 'vd';
    disp('Software version: VD (!?)');

    % number of different scans in file stored in 2nd in
    NScans = secondInt;
    measID = fread(fid,1,'uint32');
    fileID = fread(fid,1,'uint32');
    measOffset = cell(1, NScans);
    measLength = cell(1, NScans);
    for k=1:NScans
        measOffset{k} = fread(fid,1,'uint64');
        measLength{k} = fread(fid,1,'uint64'); 
        fseek(fid, 152 - 16, 'cof');
    end
else
    % in VB versions, the first 4 bytes indicate the beginning of the
    % raw data part of the file
    version  = 'vb';
    disp('Software version: VB (!?)');
    measOffset{1} = 0;
    measLength{1} = fileSize;
    NScans     = 1; % VB does not support multiple scans in one file
end

s=1;
   cPos = measOffset{s};
    fseek(fid,cPos,'bof');
    hdr_len = fread(fid, 1,'uint32');

textHeader         = fread(fid, hdr_len, 'uint8=>char').';

fclose(fid)
end