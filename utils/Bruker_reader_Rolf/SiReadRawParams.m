function parout = SiReadRawParams( filename, pars, labels, silent )
if nargin < 4
    silent = 0;
end
if nargin < 3 || numel(labels) < numel(pars)
    labels = pars;
end


labels = replace(labels,'.','_');
labels = replace(labels,'[','_');
labels = replace(labels,']','_');


parout = [];
% Reads the header of a raw data file and returns values
OutputAll = false;
if strcmpi(pars,'def') == 1
    params.strings = {'alTR[0]', 'alTE[0]', 'sSliceArray.asSlice[0].dThickness', 'sSliceArray.asSlice[0].dPhaseFOV','sSliceArray.asSlice[0].dReadoutFOV', 'sSliceArray.asSlice[0].dInPlaneRot', 'lContrasts','sKSpace.lBaseResolution','sKSpace.lPhaseEncodingLines','sKSpace.lPartitions','sRXSPEC.alDwellTime[0]','adFlipAngleDegree[0]','sPat.lAccelFactPE','sPat.lAccelFact3D'};
    params.labels = {'TR','TE','SliceThickness','FOVPhase','FOVRead', 'SliceRotation', 'Contrasts','BaseResolution','PhaseEncodingLines','Partitions','DwellTime','FlipAngle','Acceleration2D', 'Acceleration3D'};
elseif strcmpi(pars,'wip') == 1
    params.strings = {'sWiPMemBlock.alFree[0]','sWiPMemBlock.alFree[1]','sWiPMemBlock.alFree[2]','sWiPMemBlock.alFree[3]','sWiPMemBlock.alFree[4]','sWiPMemBlock.alFree[5]','sWiPMemBlock.alFree[6]','sWiPMemBlock.alFree[7]','sWiPMemBlock.alFree[8]','sWiPMemBlock.alFree[9]','sWiPMemBlock.alFree[10]','sWiPMemBlock.alFree[11]','sWiPMemBlock.alFree[12]','sWiPMemBlock.alFree[13]',
        'sWiPMemBlock.adFree[0]','sWiPMemBlock.adFree[1]','sWiPMemBlock.adFree[2]','sWiPMemBlock.adFree[3]','sWiPMemBlock.adFree[4]','sWiPMemBlock.adFree[5]','sWiPMemBlock.adFree[6]','sWiPMemBlock.adFree[7]','sWiPMemBlock.adFree[8]','sWiPMemBlock.adFree[9]','sWiPMemBlock.adFree[10]','sWiPMemBlock.adFree[11]','sWiPMemBlock.adFree[12]','sWiPMemBlock.adFree[13]'};
    params.labels = {'WIP long 0','WIP long 1', 'WIP long 2','WIP long 3','WIP long 4','WIP long 5','WIP long 6','WIP long 7','WIP long 8','WIP long 9','WIP long 10','WIP long 11','WIP long 12','WIP long 13','WIP double 0','WIP double 1','WIP double 2','WIP double 3','WIP double 4','WIP double 5','WIP double 6','WIP double 7','WIP double 8','WIP double 9','WIP double 10','WIP double 11','WIP double 12','WIP double 13','WIP double 14'};
elseif strcmpi(pars,'pulse') == 1
    params.strings = {'sTXSPEC.asNucleusInfo[0].tNucleus','sTXSPEC.asNucleusInfo[0].flReferenceAmplitude', 'sTXSPEC.asNucleusInfo[0].lFrequency','sTXSPEC.aRFPULSE[0].flAmplitude','alTD[0]','adFlipAngleDegree[0]'};
    params.labels = {'Nucleus','ReferenceVoltage','Frequency','PulseAmplitude','PulseDuration','FlipAngle'};
elseif strcmpi(pars,'all') == 1
    OutputAll = true;
    parout = 0;
    params.strings = {};
    params.labels = {};
else
    if isstruct(pars)
        params = pars;
    else
        params.strings = pars;        
        params.labels = labels;
    end
    for cnt = 1:numel(pars)
        if strcmpi(params.strings{cnt}, 'Contrasts')
            params.strings{cnt} = 'lContrasts';
            params.labels{cnt} = 'Contrasts';
        elseif strcmpi(params.strings{cnt}, 'TE')
            params.string{cnt} = 'alTE[0]';
            params.labels{cnt} = 'TE';
        elseif  strcmpi(params.strings{cnt}, 'TR')
            params.string{cnt} = 'alTR[0]';
            params.labels{cnt} = 'TR';
        else
            %params.strings{cnt} = [params.strings{cnt},char(9)];    % char(9) is a tab: entries have a tab between parameter name and '='
        end   
    end
end
nparams = numel(params.strings);
if nparams >1
    for cnt = 1:nparams
        parout.(params.labels{cnt}) = [];
    end
end


parstofind = nparams;
file = fopen(filename);
%% Read header length
fseek(file,0,'eof');
e = ftell(file);
fseek(file,0,'bof');
length = fread(file, 1,'uint32');
if length < 1000
    VD = true;
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
end

u = strtrim(fgetl(file));
ind = 1;
while strncmp(u,'ulVersion',9) == 0
    u = strtrim(fgetl(file));
end
blank = u(10);        % In old .dat files, there are blanks separating the parameter name and the '='. In new ones it's a tab. To separate parameters with similar names, the separator has to be added to the parameter name
if iscell(params.strings)
for cnt = 1:numel(params.strings)
    params.strings{cnt} = [params.strings{cnt},blank];
end
end
while strncmp(u,'### ASCCONV END ###',20) == 0
    u = strtrim(fgetl(file));
    if OutputAll
        disp(u);
    else
    for cnt = 1:nparams
        if strncmpi(u,params.strings{cnt},numel(params.strings{cnt})) == 1
            val=sscanf(u,'%*s = %f');
            if numel(val)==0
                val = sscanf(u,'%*s = %s');
                val = strip(val,'both','"');
            end
            parstofind = parstofind - 1;
            if silent == 0
                display([params.labels{cnt}, ' = ', num2str(val)]);
            end
            if nparams == 1
                parout = val;
                fclose(file);
                return;
            end
            parout.(params.labels{cnt}) = val;
%             parout(ind).value = val;
%             parout(ind).string = params.strings{cnt};
%             parout(ind).label = params.labels{cnt};
            ind = ind+1;
            break;
        end
        if parstofind == 0
            break;
        end
    end  
    end
end
fclose(file);

end

