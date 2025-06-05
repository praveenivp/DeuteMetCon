classdef DataClass 
    % Read and Process MR data
    % Methods:
    %       DataClass(path)                Constructor. If path given, reads data
    %       SiReadDICOM(filename)   Reads Siemens DICOM file
    %       SiReadRaw(filename)        Read raw data from Siemens raw data file
    %       BrReadRaw(filename)       Reads raw data from Bruker raw data file  
    %       BrReadPro(path, rnum)    Reads Bruker 2dseq file
    %       BrReadRawParams(path)  Reads parameters from Bruker method and acqp files
    %       SiReadRawParams(file)     Reads parameters from Siemens .dat file
    %       DoOrientations                  Handles orientations of different file types
    %       CoilCombine                     Performs coil combination
    %       ReadWeighting           Read Weighting from weighting.txt file
    %       CSIRearrange                    Reorders data for 3D-CSI
    %       CSIApplyWeighting           Applies retrospective weighting correction
    %       CSISpatialFT(zf)                CSI: fourier-transform spatial directions
    %       CSICalcPSF                      CSI: Calculates PSF
    %       CSIRotate                         CSI: Rotates data in right orientation
    %       ChooseSets
    %       SpecShift(index)            CSI: Modifies fid so spectra are centered
    %       PhaseCor(index)             CSI:
    %       CropFID(cropPoints)
    %       fidExtrp (numpoints, ProOrIndex)
    %       ShowRest
    %       ShowStim
    %       Show

    properties
        AcqParameters
        RawData
        ProData
        ProDataStim
        ProDataRest
        ProDataFull
        RecoParameters
        Status
    end

    methods
        function obj = DataClass(path)
            % Create an empty object of this class
            obj.RawData = [];
            obj.ProData = [];
            obj.ProDataStim = [];
            obj.ProDataRest = [];
            obj.ProDataFull = [];
            obj.Status.IsRawData = 0;
            obj.Status.IsProData = 0;
            obj.Status.IsDICOM = 0;
            obj.Status.Is2DSeq = 0;
            obj.Status.IsCSI = 0;
            obj.Status.StateRawData = 0;
            obj.Status.StateProData = 0;
            obj.AcqParameters.Sequence = '';
            obj.AcqParameters.Bandwidth = 0;
            obj.AcqParameters.Pathname = '';
            obj.AcqParameters.Filename = '';
            obj.AcqParameters.ProtocolName = '';
            obj.AcqParameters.Contrasts = 1;
            obj.AcqParameters.Averages = 1;
            obj.AcqParameters.ScanTime = 0;
            obj.AcqParameters.B0 = 0;
            obj.AcqParameters.Nulceus = '1H';
            obj.AcqParameters.Frequency = 0;
            obj.AcqParameters.ReferenceVoltage = 0;
            obj.AcqParameters.DwellTime = 0;
            obj.AcqParameters.AcqDur = 0;
            obj.AcqParameters.TR = 0;
            obj.AcqParameters.TE = 0;
            obj.AcqParameters.SliceThick = 0;
            obj.AcqParameters.FOV = [0,0,0];
            obj.AcqParameters.NDims = 1;
            obj.AcqParameters.FlipAngle = 0;
            obj.AcqParameters.Position = [0,0,0];
            obj.AcqParameters.Orientation = [0,0,0];
            obj.AcqParameters.SliceOrient = '';
            obj.AcqParameters.MatrixSize = [0,0,0,0];
            obj.AcqParameters.Channels = 0;
            obj.AcqParameters.weighting.weight = [];
            obj.AcqParameters.weighting.corr = [];
            obj.AcqParameters.NoiseCov = [];
            obj.AcqParameters.kspace = [];
            obj.AcqParameters.psf = [];
            obj.AcqParameters.Units = 'px';
            obj.RecoParameters=[];
            if nargin >= 1
                % Read data
                % First, find out what datatype
                if isfile(path)         % A file name
                    [p,n,ext] = fileparts(path);
                    if strcmp(ext,'.ima')
                        obj = obj.SiReadDICOM(path);
                    elseif strcmp(ext,'.dat')
                        obj = obj.SiReadRaw(path);
                    elseif strcmp(n,'2dseq')
                        obj = obj.BrReadPro(path);
                    elseif strcmp(n,'rawdata')
                        obj = obj.BrReadRaw(path);
                    else
                        fprintf('%s is not a valid data file!\n',path);
                        return
                    end
                elseif isfolder(path)   % A folder
                    imafiles = dir([path,filesep,'*.ima']);
                    if numel(imafiles) > 0
                        obj = obj.SiReadDICOM(path);
                    else
                        datfiles = dir([path,filesep,'*.dat']);     
                        if numel(datfiles)>0            % This is not yet finished
                        elseif exist([path,filesep,'rawdata.job0'],'file') == 2
                            if exist([path,filesep,'pdata',filesep,'1',filesep,'2dseq'],'file') == 2
                                obj = obj.BrReadPro([path,filesep,'pdata',filesep,'1',filesep,'2dseq']);
                            else
                                % Here, the raw data part is still missing
                            end
                        end
                            
                    end
                else
                    % may be a file name with *  (This works only for DICOM files)
                    a = dir(path);
                    numima = 0;
                    if numel(a) > 0
                        for cnt = 1:numel(a)
                            [p,n,ext] = fileparts(a(cnt).name);
                            if strcmpi(ext,'.ima')
                                numima = numima + 1;
                            end
                        end
                    end
                    if numima > 0
                        obj = obj.SiReadDICOM(path);
                    else
                        fprintf('File or folder %s does not exist!\n',path);                        
                        return
                    end
                end
            end

        end

        function r = plus(obj1,obj2)
            if size(obj1.RawData) ~= size(obj2.RawData)
                fprintf('Objects can not be added since they have different sizes!');
                return;
            end
            r = obj1;
            r.RawData = r.RawData + obj2.RawData;
            if size(obj1.ProData) == size(obj2.ProData)
                r.ProData = r.ProData + obj2.ProData;
            else
                r.ProData = [];
            end
        end


        function obj = SiReadDICOM(obj, filename)
            % Reads image data from a series of .ima files
            % filename can be:
            %                               - a single DICOM file to be read
            %                               - not implemented:  a cell array containing several DICOM files. These files have to be from the same scan, otherwise reading will fail.
            %                               - a path containing DICOM files. In that case, a selection dialog will appear.
            %                               - not immplemented: a path containing the unique part of the filenames of the image
            %                               - nothing: a selection dialog will appear
            images = SiReadProData(filename);
            nim = numel(images);
            if nim == 0
                return;
            end
            slice = zeros(nim,1);
            orient = zeros(nim,1);
            contrast = zeros(nim,1);
            position = zeros(nim,3);
            imsize = size(images(1).Image);
            for cnt = 1:numel(images)   % First count images
                is = size(images(cnt).Image);
                if is(1) ~= imsize(1) || is(2) ~= imsize(2)
                    fprintf('Cannot read data: Data contains images with different sizes!\n');
                end
                position(cnt,:) = images(cnt).Info.Coord;
                contrast(cnt) = images(cnt).Info.Contrast;
            end
            [position,order] = sort(position,1);
            obj.ProData = zeros(imsize(1),imsize(2),numel(unique(position(:,1))),numel(unique(contrast)), numel(unique(orient)));  % Here are some changes necessary here!!!!
            posind = 1;
            conind = 1;
            orientind = 1;
            repind = 1;
            for cnt = 1:numel(images)
                obj.ProData(:,:,posind, conind, orientind) = images(cnt).Image;
                posind = posind + 1;
            end
            numpos = posind - 1;
            obj.AcqParameters.Sequence = images(1).Info.Sequence;
            %obj.AcqParameters.Bandwidth = 0;
            obj.AcqParameters.Pathname = images(1).Info.path;
            [~,b,c] = fileparts(images(1).Info.filename);
            obj.AcqParameters.Filename = [b,c];     
            obj.AcqParameters.ProtocolName = images(1).Info.Protocol;
            obj.AcqParameters.Contrasts = images(1).Info.NContrasts;
            obj.AcqParameters.Averages = images(1).Info.Averages;
            obj.AcqParameters.ScanTime = 0;
            obj.AcqParameters.B0 = 0;
            obj.AcqParameters.Nulceus = images(1).Info.Nucleus;
            obj.AcqParameters.Frequency = 0;
            obj.AcqParameters.ReferenceVoltage = 0;
            obj.AcqParameters.DwellTime = 0;
            obj.AcqParameters.AcqDur = 0;
            obj.AcqParameters.TR = images(1).Info.TR;
            obj.AcqParameters.TE = images(1).Info.TE;
            obj.AcqParameters.SliceThick = images(1).Info.Slice;
            obj.AcqParameters.FOV = [images(1).Info.FOV,images(1).Info.Slice];
            obj.AcqParameters.NDims = numel(size(obj.ProData));
            obj.AcqParameters.Orientation = images(1).Info.Orientation;
            obj.AcqParameters.SliceOrient = '';
            obj.AcqParameters.PatientPos = images(1).Info.PatientPos;
            obj.AcqParameters.Channels = 0;     
            obj.AcqParameters.Units = 'mm';
            obj = obj.DoOrientations;
            if numpos > 1
                sliceind = find(abs(obj.RecoParameters.RearrangeDims)==3);
                indx = find(abs(obj.RecoParameters.RearrangeDims)==1);
                indy = find(abs(obj.RecoParameters.RearrangeDims)==2);
               % lowpos = min(position(:,abs(obj.RecoParameters.RearrangeDims(3))));
               % highpos = max(position(:,abs(obj.RecoParameters.RearrangeDims(3))));
                lowpos = min(position(:,sliceind));
                highpos = max(position(:,sliceind));
                slicefov = highpos-lowpos;                                                                       % This is the distance between the centers of the outer slices; just used for calculating the position (since the slice position also is the center in contrast to DICOM position defininitions)
                obj.AcqParameters.FOV(3)  = (highpos-lowpos)/(numpos-1)*numpos;   % This is the real FOV between the outer values of the outer slices
            end
            obj.AcqParameters.MatrixSize(1:2) = images(1).Info.ImageSize;
            obj.AcqParameters.MatrixSize(3) = numpos;
%             obj.AcqParameters.Position = [position(1,indx) + sign(obj.RecoParameters.RearrangeDims(indx))*obj.AcqParameters.FOV(1)/2, ...
%                                                               position(1,indy) + sign(obj.RecoParameters.RearrangeDims(indy))*obj.AcqParameters.FOV(2)/2,  ...
%                                                               (position(end,abs(obj.RecoParameters.RearrangeDims(3)))-position(1,abs(obj.RecoParameters.RearrangeDims(3))))/2+position(1,abs(obj.RecoParameters.RearrangeDims(3)))];
%             obj.AcqParameters.Position = [position(1,abs(obj.RecoParameters.RearrangeDims(1))) + sign(obj.RecoParameters.RearrangeDims(1))*obj.AcqParameters.FOV(1)/2, ...
%                                                               position(1,abs(obj.RecoParameters.RearrangeDims(2))) + sign(obj.RecoParameters.RearrangeDims(2))*obj.AcqParameters.FOV(2)/2,  ...
%                                                               (position(end,abs(obj.RecoParameters.RearrangeDims(3)))-position(1,abs(obj.RecoParameters.RearrangeDims(3))))/2+position(1,abs(obj.RecoParameters.RearrangeDims(3)))];
            obj.RecoParameters.FOV = obj.AcqParameters.FOV(abs(obj.RecoParameters.RearrangeDims(1:3)));
            obj.RecoParameters.MatrixSize = obj.AcqParameters.MatrixSize(abs(obj.RecoParameters.RearrangeDims(1:3)));
            fov = obj.AcqParameters.FOV;
            fov(3) = slicefov;
            fov = fov(abs(obj.RecoParameters.RearrangeDims(1:3)));
%            obj.RecoParameters.Position = [-(position(1,2) + sign(obj.RecoParameters.RearrangeDims(2))*fov(2)/2), ...
%                                                              -(position(1,1) + sign(obj.RecoParameters.RearrangeDims(1))*fov(1)/2),  ...
%                                                              position(1,3) + sign(obj.RecoParameters.RearrangeDims(3))*fov(3)/2];
            obj.RecoParameters.Position = [position(1,1) + sign(obj.RecoParameters.RearrangeDims(1))*obj.RecoParameters.FOV(1)/2, ...
                                                              position(1,2) + sign(obj.RecoParameters.RearrangeDims(2))*obj.RecoParameters.FOV(2)/2,  ...
                                                              position(1,3) + sign(obj.RecoParameters.RearrangeDims(3))*obj.RecoParameters.FOV(3)/2];
                                                              %(position(end,abs(obj.RecoParameters.RearrangeDims(3)))-position(1,abs(obj.RecoParameters.RearrangeDims(3))))/2+position(1,abs(obj.RecoParameters.RearrangeDims(3)))];
            %obj.RecoParameters.Position =  obj.AcqParameters.Position(abs(obj.RecoParameters.RearrangeDims(1:3)));
        end

        function obj = SiReadRaw(obj, filename,average,brokendata)
            % Reads raw data from Siemens .dat file
           
            if nargin < 4
                brokendata = 0;
            end
            if nargin < 3
                average = 1;
            end
            if exist(filename, "file")== 0
                fprintf('File %s not found!\n',filename);
                return;
            else
                fprintf("Reading file %s.\n",filename);
                obj = obj.SiReadRawParams(filename);
                if brokendata == 0
                    [obj.RawData,corr,~,mdhs] = SiReadRawData(filename,[],average,0);
                else
                    [obj.RawData,corr,~,mdhs] = SiReadRawDataBroken(filename,[],average,0);
                end
                if isfield(corr,'Noiseadj')
                    n = size(corr.Noiseadj.data);
                    corr.Noiseadj.data = reshape(corr.Noiseadj.data,[n(1)*n(2),n(3)]);
                    obj.AcqParameters.NoiseCov = cov(corr.Noiseadj.data);
                    if numel(find(SiFlags(mdhs(1).flags) == 25)) > 0                 % Usually the first scan is the noise adjustment -> remove it from the mdhs
                        mdhs = mdhs(2:end);
                    end
                end
                if numel(obj.AcqParameters.Channels) == 0 || obj.AcqParameters.Channels == 0
                    s = size(obj.RawData);
                    obj.AcqParameters.Channels = s(end);
                end
                % Some special processing for csi data
                if contains(obj.AcqParameters.Sequence,'csi','IgnoreCase',true)
                    % Read k-space weighting
                    dim2 = [mdhs.ILine];
                    dim1 = [mdhs.ISegment];
                    %dim3 = [mdhs.IPhase];
                    dim3 = [mdhs.IPartition];
                    min1 = min(dim1);
                    num1 = max(dim1)-min1+1;
                    min2 = min(dim2);
                    num2 = max(dim2)-min2+1;
                    min3 = min(dim3);
                    num3 = max(dim3)-min3+1;
                    min1 = 0;min2 = 0;min3 = 0;

                    %kspace = zeros(num1, num2, num3);
                    kspace = zeros(max(dim1)+1, max(dim2)+1, max(dim3)+1);
                    correctweighting = 0;
                    for cnt = 1:numel(dim1)
                        kspace(dim1(cnt)-min1+1, dim2(cnt)-min2+1, dim3(cnt)-min3+1) = kspace(dim1(cnt)-min1+1, dim2(cnt)-min2+1, dim3(cnt)-min3+1)+1;
                        if mdhs(cnt).Free4 >0
                            obj.AcqParameters.weighting.weight(dim1(cnt)-min1+1, dim2(cnt)-min2+1, dim3(cnt)-min3+1) = mdhs(cnt).Free3;
                            obj.AcqParameters.weighting.corr(dim1(cnt)-min1+1, dim2(cnt)-min2+1, dim3(cnt)-min3+1) = mdhs(cnt).Free4/10000;
                            correctweighting = 1;
                        end
                    end
                    obj.AcqParameters.kspace = kspace;
                    % CSI data is in the order [Line, Partition, ....., Segment], but should be in the order [Baseline, Line, Part], which means [Segment, Line, Parittion]. This does the reordering
                    %obj.RawData = flip(flip(permute(obj.RawData,[1,2,8,3,4,5,6,7,9]),3),2); 
                    obj.RawData = permute(obj.RawData,[1,8,2,3,4,5,6,7,9]);
                    if correctweighting == 1
                        sw = size( obj.AcqParameters.weighting.corr);
                        s = size(obj.RawData);
                        corr = reshape(obj.AcqParameters.weighting.corr,[1,sw]);
                        %tic
                        for cnt1 = 1:s(1)
                            %fprintf('Correction: step %d of %d: %f s\n',cnt1,s(1),toc);
                            for cnt5 = 1:s(5)
                                for cnt6 = 1:s(6)
                                    for cnt7 = 1:s(7)
                                        for cnt8 = 1:s(8)
                                            for cnt9 = 1:s(9)
                                                obj.RawData(cnt1, s(2)-sw(1)+1:s(2),s(3)-sw(2)+1:s(3),s(4)-sw(3)+1:s(4),cnt5, cnt6, cnt7, cnt8, cnt9) =  obj.RawData(cnt1, s(2)-sw(1)+1:s(2),s(3)-sw(2)+1:s(3),s(4)-sw(3)+1:s(4),cnt5, cnt6, cnt7, cnt8, cnt9).* corr;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end                        

                    % For weighted, the raw data should be filled to reach the initially given matrix size
                    s = size(obj.RawData);
                    if s(2)<obj.AcqParameters.MatrixSize(1) || s(3)<obj.AcqParameters.MatrixSize(2) || s(4)<obj.AcqParameters.MatrixSize(3)
                        s(2:4) = obj.AcqParameters.MatrixSize;
                        obj.RawData(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9)) = 0;
                    end
                    % The Matrix Size may differ from the data size, because it contains empty lines or columns
                    %s = size(obj.RawData);
                    %if max(obj.AcqParameters.MatrixSize ~= s(2:4)) == 1
                    %    obj.AcqParameters.MatrixSize = s(2:4);
                    %end
                    % Now we have to make sure the order is (x,y,z) and the RearangeDims parameter is correct
                    [~,ind] = max(obj.AcqParameters.Orientation);
                    if ind == 3          % transverse
                        obj.RecoParameters.RearrangeDims = [-1,-2,3];
                    elseif ind == 2     % coronal
                        obj.RecoParameters.RearrangeDims = [-2,3,-1];
                    elseif ind == 1     % sagittal
                        obj.RecoParameters.RearrangeDims = [3,-2,1];
                    end
                    obj.AcqParameters.CSIPosition = [0,0,0];
                end
            end
        end

        function obj = BrReadRaw(obj, filename)
            SpecCutOff = 1;
            % Reads raw data from Bruker rawdata.job0 file
            if exist(filename, "file")== 0
                fprintf('File %s not found!\n',filename);
                return;
            else
                obj = obj.BrReadRawParams(filename);
                obj.RawData = BrReadFid([filename,filesep,'rawdata.job0']);
                if obj.AcqParameters.SpectSize > 0
                    obj.RawData = reshape(obj.RawData,obj.AcqParameters.SpectSize,[]);
                else
                    obj.RawData = reshape(obj.RawData,obj.AcqParameters.MatrixSize(1)*obj.AcqParameters.Channels,[]);
                end
                s = size(obj.RawData);
                if numel(obj.AcqParameters.kspace) > 1
                    ind = 1;
                    for cnt = 1:numel(obj.AcqParameters.kspace)
                        for avecnt = 1:obj.AcqParameters.kspace(cnt)
                            if avecnt == 1
                                if ind <= s(2)
                                    obj.RawData(:,cnt) = obj.RawData(:,ind);
                                else
                                    obj.RawData(:,cnt) = 0;     % This ind <= s(2) clause was inserted to allow for reading a scan that was accidentally stopped shortly before finishing
                                end
                            else
                                if ind <= s(2)
                                    obj.RawData(:,cnt) = obj.RawData(:,cnt) +obj.RawData(:,ind);
                                end
                            end
                            ind = ind+1;
                        end
                    end
                    obj.RawData = obj.RawData(:,1:cnt);
                end
                if obj.AcqParameters.Channels >1
                    obj.RawData = reshape(obj.RawData,obj.AcqParameters.MatrixSize(1),obj.AcqParameters.Channels,[]);
                   obj.RawData = permute(obj.RawData,[1,3,2]);
                end
                % Some special processing for csi data
                if obj.AcqParameters.SpectSize > 1 && obj.AcqParameters.MatrixSize(1) > 1
                    if numel(obj.AcqParameters.kspace) > 1
                        obj.AcqParameters.kspace = reshape(obj.AcqParameters.kspace,obj.AcqParameters.MatrixSize);
                    end
                    obj.RawData = reshape(obj.RawData,[obj.AcqParameters.SpectSize,obj.AcqParameters.MatrixSize,obj.AcqParameters.Channels]);
                    % Now remove the initial filter delay in the raw data
                    % Since the delay may vary, determine it yourself
                    s = size(obj.RawData);
                    ss = floor(s/2);
                    testdat = sum(sum(sum(abs(obj.RawData(:,ss(2)-3:ss(2)+3,ss(3)-3:ss(3)+3,ss(4)-3:ss(4)+3)),4),3),2);
                    mtest = max(testdat);
                    high = find(testdat > mtest/5);
                    filtdel = high(1)-1;
                    fprintf('Removing filter delay of %d points.\n',filtdel);
                    obj.RawData = circshift(obj.RawData,[-filtdel,0,0,0,0]);
                    % What to do with the cut-off data? We can either cut it off, set it to zero or replace it with extrapolated values
                    if SpecCutOff == 0             % Cut it off
                        obj.RawData = obj.RawData(1:end-filtdel,:,:,:,:);
                    elseif SpecCutOff == 1      % Replace it with zeros (default)
                        obj.RawData(end-filtdel+1:end,:,:,:,:) = 0;
                    elseif SpecCutOff == 2    % Replace it with extrapolated values
                        obj = obj.Extrapolate();
                    end
                    obj.AcqParameters.CSIPosition = obj.AcqParameters.Position;
                    obj.AcqParameters.Position = [0,0,0];             % In Bruker data, CSI is always acquired around the center, and should be shifted during recontruction. Since this is not yet done here, the real position is 0, the shifted position in stored in CSIPosition.
                    % Now the data has to be arranged to be (x,y,z)
                    gradmat = obj.AcqParameters.GradMatrix(:,:,1);
                    [~,ind] = max(abs(gradmat'));
                    obj.RecoParameters.RearrangeDims = ind;                    
                    %obj.RecoParameters.RearrangeDims = round(gradmat')*[1;2;3];
%                     if find(ind==3) == 3
%                         obj.RecoParameters.RearrangeDims(obj.RecoParameters.RearrangeDims==1) = -obj.RecoParameters.RearrangeDims(obj.RecoParameters.RearrangeDims==1);
%                         obj.RecoParameters.RearrangeDims(obj.RecoParameters.RearrangeDims==3) = -obj.RecoParameters.RearrangeDims(obj.RecoParameters.RearrangeDims==3);
%                     elseif find(ind==3) == 2
%                         if obj.RecoParameters.RearrangeDims(1) == 1
%                             obj.RecoParameters.RearrangeDims(1) = -obj.RecoParameters.RearrangeDims(1);
%                             obj.RecoParameters.RearrangeDims(3) = -obj.RecoParameters.RearrangeDims(3);
%                         elseif obj.RecoParameters.RearrangeDims(1) == 2
%                             obj.RecoParameters.RearrangeDims(2) = -obj.RecoParameters.RearrangeDims(2);
%                             obj.RecoParameters.RearrangeDims(3) = -obj.RecoParameters.RearrangeDims(3);
%                         end
%                     end
                    for cnt = 1:3
                        if gradmat(cnt, ind(cnt))<0
                            obj.RecoParameters.RearrangeDims(cnt) = -obj.RecoParameters.RearrangeDims(cnt);
                        end
                        if cnt == 3 && ind(cnt)<3 
                            obj.RecoParameters.RearrangeDims(cnt) = -obj.RecoParameters.RearrangeDims(cnt);
                        end
                        if cnt < 3 && ind(cnt)==3 
                            obj.RecoParameters.RearrangeDims(cnt) = -obj.RecoParameters.RearrangeDims(cnt);
                        end                       
                    end
                    obj.RecoParameters.RearrangeDims(1) = -obj.RecoParameters.RearrangeDims(1);   % This may only be correct for Head_Supine
                    obj.RecoParameters.RearrangeDims(3) = -obj.RecoParameters.RearrangeDims(3);
                else
                    obj.RawData = reshape(obj.RawData,[obj.AcqParameters.MatrixSize,obj.AcqParameters.Constrasts,obj.AcqParameters.Channels]);
                end
            end
        end

        function obj = BrReadPro(obj, pathname, reconum)
            % pathname can be either the path up to the 2dseq-file, or just the scan directory. Then, the reconum is the number of the reconstruction.
            if nargin < 3
                reconum = 1;
            end
            if exist(pathname,"dir") == 7             % pathname is a directory
                if exist([pathname,filesep,'2dseq'],'file') == 2     
                    pathname = [pathname,filesep,'2dseq'];
                elseif exist([pathname,filesep,'pdata',filesep,num2str(reconum),filesep,'2dseq'],'file') == 2  
                    pathname = [pathname,filesep,'pdata',filesep,num2str(reconum),filesep,'2dseq'];
                else
                    fprintf('%s is not a valid 2dseq directory or file!\n',pathname);
                    return
                end
            elseif exist(pathname, 'file') == 2 && endsWith(pathname, '2dseq')
            else
                fprintf('%s is not a valid 2dseq directory or file!\n',pathname);
                return
            end
            scanpath = fileparts(pathname);
            obj = obj.BrReadRawParams([scanpath,filesep,'..',filesep,'..']);
            [obj.ProData, acqppars] = BrReadProImage(pathname);
            obj.ProData = squeeze(flipud(permute(obj.ProData,[1,2,3,5,4])));   % order should be [x,y,slices,contrasts]. What is dimension 5?
            obj = obj.DoOrientations;
            % Here, the real RearrangeDims also depends on the Read direction
             [~,i] = max(abs(acqppars.ACQ_GradientMatrix(:,:,1))');
             i1 = find(i==1);
             i2 = find(i==2);
             %i(i1) = 2;
             %i(i2) = 1;
             obj.RecoParameters.RearrangeDims = i;
%             obj.RecoParameters.RearrangeDims = [fliplr(i(find(i~=3))),3];
            obj.RecoParameters.MatrixSize = [acqppars.recosize(2),acqppars.recosize(1),acqppars.NSLICES];
            if numel(acqppars.fov) == 2    % slice selective 2D data
                obj.RecoParameters.FOV = [10*acqppars.fov(1),10*acqppars.fov(2),max(acqppars.ACQ_slice_offset) - min(acqppars.ACQ_slice_offset)];
            elseif numel(acqppars.fov) == 3                     % 3D data
                obj.RecoParameters.FOV = [10*acqppars.fov(1),10*acqppars.fov(2),10*acqppars.fov(3)];
            end
            obj.RecoParameters.FOV = obj.RecoParameters.FOV(abs(obj.RecoParameters.RearrangeDims(1:3)));
            s = size(obj.ProData);
            obj.RecoParameters.MatrixSize = s([2,1,3]);  %obj.RecoParameters.MatrixSize(abs(obj.RecoParameters.RearrangeDims(1:3)));
            obj.AcqParameters.Units = 'mm';
            obj.RecoParameters.Position = [-obj.AcqParameters.Position(obj.RecoParameters.RearrangeDims(1)),-obj.AcqParameters.Position(obj.RecoParameters.RearrangeDims(2)),obj.AcqParameters.Position(obj.RecoParameters.RearrangeDims(3))];
        end

        function obj = BrReadRawParams(obj, pathname)
            % Reads parameters from a raw data file
            if exist([pathname,filesep,'method'], "file")== 0 || exist([pathname,filesep,'acqp'], "file")== 0
                fprintf('%s not a valid Bruker raw data directory!\n',pathname);
                return;
            else
                paramnames = {  'Method', ...
                                    'PVM_NEchoImages', ...
                                    'PVM_NAverages',...
                                    'PVM_ScanTime',...
                                    'PVM_Nucleus1', ...
                                    'PVM_FrqWork', ...
                                    'PVM_RefPowCh1', ...
                                    'PVM_SpecSWH', ...
                                    'PVM_EffSWh', ...
                                    'PVM_RepetitionTime', ...
                                    'PVM_EchoTime', ...
                                    'PVM_SliceThick',...
                                    'PVM_Fov',...
                                    'PVM_SPackArrReadOffset',...
                                    'PVM_SPackArrPhase1Offset',...
                                    'PVM_SPackArrSliceOffset',...
                                    'PVM_SPackArrGradOrient',...
                                    'PVM_SliceGeo',...
                                    'PVM_SPackArrSliceOrient', ...
                                    'PVM_Matrix', ...
                                    'PVM_SpecMatrix',...
                                    'PVM_EncActReceivers', ...
                                    'AverageList'};
                                                    %'tProtocolName', ...???
                                    %'sProtConsistencyInfo.flNominalB0', ...???
                                    %'sSliceArray.asSlice[0].dInPlaneRot',...???

                p = BrReadParams(paramnames,[pathname,filesep,'method']);
                obj.AcqParameters.Filename = pathname;
                obj.AcqParameters.Sequence = p.Method;
                %obj.AcqParameters.ProtocolName = p.tProtocolName;
                obj.AcqParameters.Contrasts = p.PVM_NEchoImages;;
                obj.AcqParameters.Averages = p.PVM_NAverages;
                obj.AcqParameters.ScanTime = p.PVM_ScanTime/1000;
                %obj.AcqParameters.B0 = p.sProtConsistencyInfo_flNominalB0;
                obj.AcqParameters.Nucleus = p.PVM_Nucleus1;
                obj.AcqParameters.Frequency = p.PVM_FrqWork(1);
                obj.AcqParameters.ReferenceVoltage = p.PVM_RefPowCh1;   % Probably not voltage but watts
                if p.PVM_SpecSWH > 0
                    obj.AcqParameters.Bandwidth = p.PVM_SpecSWH;
                else
                    obj.AcqParameters.Bandwidth = p.PVM_EffSWh;
                end
                obj.AcqParameters.DwellTime = 1/obj.AcqParameters.Bandwidth;
                obj.AcqParameters.AcqDur = 0;
                obj.AcqParameters.TR = p.PVM_RepetitionTime;
                obj.AcqParameters.TE = p.PVM_EchoTime;
                obj.AcqParameters.SliceThick = p.PVM_SliceThick;
                obj.AcqParameters.FOV = p.PVM_Fov;
                obj.AcqParameters.NDims = numel(p.PVM_Matrix);
                obj.AcqParameters.Orientation = p.PVM_SPackArrGradOrient*[0;0;1];   % This is supposed to find the normal vector. Correct?
                obj.AcqParameters.Position = [p.PVM_SPackArrReadOffset,p.PVM_SPackArrPhase1Offset,p.PVM_SPackArrSliceOffset];
                if strcmp(p.PVM_SPackArrSliceOrient,'axial')
                    obj.AcqParameters.SliceOrient = 'xy';
                elseif strcmp(p.PVM_SPackArrSliceOrient,'coronal')
                    obj.AcqParameters.SliceOrient = 'xz';
                elseif strcmp(p.PVM_SPackArrSliceOrient,'sagittal')
                    obj.AcqParameters.SliceOrient = 'yz';
                end
                obj.AcqParameters.MatrixSize =p.PVM_Matrix';
                obj.AcqParameters.Channels = sum(contains(p.PVM_EncActReceivers,'On'));
                obj.AcqParameters.SpectSize = p.PVM_SpecMatrix;
                obj.AcqParameters.kspace = p.AverageList;
                paramnames = {  'BF5', ...
                                             'ACQ_protocol_name', ...
                                             'ACQ_GradientMatrix', ...
                                             'ACQ_grad_matrix'};
                ap = BrReadParams(paramnames,[pathname,filesep,'acqp']);
                obj.AcqParameters.ProtocolName = ap.ACQ_protocol_name;
                obj.AcqParameters.B0 = ap.BF5/42.5756;
                obj.AcqParameters.GradMatrix = p.PVM_SPackArrGradOrient;
                % The following would be valid for xyz coordinate system, but only works for PV360, not for PV6, because in PV6 ACQ_grad_matrix is equal to PVM-SPackArrGradOrient 
%                 if numel(ap.ACQ_GradientMatrix) >= 3
%                     obj.AcqParameters.GradMatrix = ap.ACQ_GradientMatrix;
%                 else
%                         obj.AcqParameters.GradMatrix = ap.ACQ_grad_matrix;
%                 end
            end
        end

        function obj = SiReadRawParams(obj, filename)
            % Reads parameters from a raw data file
            if exist(filename, "file")== 0
                fprintf('File %s not found!\n',filename);
                return;
            else
                paramnames = {  'tSequenceFileName', ...
                                    'tProtocolName', ...
                                    'lContrasts', ...
                                    'lAverages',...
                                    'lTotalScanTimeSec',...
                                    'sProtConsistencyInfo.flNominalB0', ...
                                    'sTXSPEC.asNucleusInfo[0].tNucleus', ...
                                    'sTXSPEC.asNucleusInfo[0].lFrequency', ...
                                    'sTXSPEC.asNucleusInfo[0].flReferenceAmplitude', ...
                                    'sRXSPEC.alDwellTime[0]', ...
                                    'alTR[0]', ...
                                    'alTE[0]', ...
                                    'sSliceArray.asSlice[0].dThickness',...
                                    'sSliceArray.asSlice[0].dPhaseFOV',...
                                    'sSliceArray.asSlice[0].dReadoutFOV',...
                                    'sSliceArray.asSlice[0].dInPlaneRot',...
                                    'sSliceArray.asSlice[0].sPosition.dTra',...
                                    'sSliceArray.asSlice[0].sPosition.dSag',...
                                    'sSliceArray.asSlice[0].sPosition.dCor',...
                                    'sSliceArray.asSlice[0].sNormal.dTra',...
                                    'sSliceArray.asSlice[0].sNormal.dSag',...
                                    'sSliceArray.asSlice[0].sNormal.dCor',...
                                    'sKSpace.lBaseResolution', ...
                                    'sKSpace.lPhaseEncodingLines', ...
                                    'sKSpace.lPartitions', ...
                                    'sSpecPara.lVectorSize',...
                                    'adFlipAngleDegree[0]',...
                                    'sCoilSelectMeas.aRxCoilSelectData[0].asList.__attribute__.size'};
                p = SiReadRawParams(filename,paramnames,[],1);
                obj.AcqParameters.Filename = filename;
                obj.AcqParameters.Sequence = p.tSequenceFileName;
                obj.AcqParameters.ProtocolName = p.tProtocolName;
                obj.AcqParameters.Contrasts = p.lContrasts;
                obj.AcqParameters.Averages = p.lAverages;
                obj.AcqParameters.ScanTime = p.lTotalScanTimeSec;
                obj.AcqParameters.B0 = p.sProtConsistencyInfo_flNominalB0;
                obj.AcqParameters.Nucleus = p.sTXSPEC_asNucleusInfo_0__tNucleus;
                obj.AcqParameters.Frequency = p.sTXSPEC_asNucleusInfo_0__lFrequency/1e6;
                obj.AcqParameters.ReferenceVoltage = p.sTXSPEC_asNucleusInfo_0__flReferenceAmplitude;
                obj.AcqParameters.DwellTime = p.sRXSPEC_alDwellTime_0_/1e6;
                obj.AcqParameters.Bandwidth = 1/p.sRXSPEC_alDwellTime_0_*1e9;
                obj.AcqParameters.AcqDur = 0;
                obj.AcqParameters.TR = p.alTR_0_/1000;
                obj.AcqParameters.TE = p.alTE_0_/1000;
                obj.AcqParameters.SliceThick = p.sSliceArray_asSlice_0__dThickness;
                obj.AcqParameters.FOV = [p.sSliceArray_asSlice_0__dReadoutFOV,p.sSliceArray_asSlice_0__dPhaseFOV, p.sSliceArray_asSlice_0__dThickness];
                obj.AcqParameters.NDims = 0;
                obj.AcqParameters.Orientation = [0,0,0];
                obj.AcqParameters.FlipAngle = p.adFlipAngleDegree_0_;
                if numel(p.sSliceArray_asSlice_0__sNormal_dSag) >0
                    obj.AcqParameters.Orientation(1) = p.sSliceArray_asSlice_0__sNormal_dSag;
                end
                if numel(p.sSliceArray_asSlice_0__sNormal_dCor) >0
                    obj.AcqParameters.Orientation(2) = p.sSliceArray_asSlice_0__sNormal_dCor;
                end
                if numel(p.sSliceArray_asSlice_0__sNormal_dTra) >0
                    obj.AcqParameters.Orientation(3) = p.sSliceArray_asSlice_0__sNormal_dTra;
                end
                if numel(p.sSliceArray_asSlice_0__sPosition_dSag) >0
                    obj.AcqParameters.Position(1) = p.sSliceArray_asSlice_0__sPosition_dSag;
                end
                if numel(p.sSliceArray_asSlice_0__sPosition_dCor) >0
                    obj.AcqParameters.Position(2) = p.sSliceArray_asSlice_0__sPosition_dCor;
                end
                if numel(p.sSliceArray_asSlice_0__sPosition_dTra) >0
                    obj.AcqParameters.Position(3) = p.sSliceArray_asSlice_0__sPosition_dTra;
                end
                obj.AcqParameters.SliceOrient = '';
                obj.AcqParameters.MatrixSize = [p.sKSpace_lBaseResolution,p.sKSpace_lPhaseEncodingLines,p.sKSpace_lPartitions];
                obj.AcqParameters.Channels = p.sCoilSelectMeas_aRxCoilSelectData_0__asList___attribute___size;
                obj.AcqParameters.SpectSize = p.sSpecPara_lVectorSize;
                obj.AcqParameters.Units = 'mm';
            end
        end

        function obj = DoOrientations(obj)
            if numel(obj.AcqParameters.Orientation) == 6              % DICOM orientation is: direction cosines of first row, direction cosines of first column
                if isfield(obj.AcqParameters,'PatientPos') == 0 || numel(obj.AcqParameters.PatientPos) == 0 || strcmp(obj.AcqParameters.PatientPos,'HFS')    % Head first supine is default
                    if obj.AcqParameters.Orientation(3) == 0 && obj.AcqParameters.Orientation(6) == 0      % transverse
                        obj.AcqParameters.SliceOrient = 'xy';
                        obj.RecoParameters.RearrangeDims = [2,1,3];
                        fprintf('transverse\n');
                    elseif obj.AcqParameters.Orientation(1) == 0 && obj.AcqParameters.Orientation(4) == 0      % sagittal
                        obj.AcqParameters.SliceOrient = 'yz';
                        obj.ProData = flip(permute(obj.ProData,[2,3,1]),3);
                        %obj.RecoParameters.RearrangeDims = [-3,2,1];
                        obj.RecoParameters.RearrangeDims = [3,2,-1];
                        fprintf('sagittal\n');
                    elseif obj.AcqParameters.Orientation(2) == 0 && obj.AcqParameters.Orientation(5) == 0      % coronal
                        obj.AcqParameters.SliceOrient = 'xz';
                        obj.ProData = flip(permute(obj.ProData,[3,2,1]),3);
                        obj.RecoParameters.RearrangeDims = [2,3,-1];
                        fprintf('coronal\n');
                    else
                        fprintf('Oblique slice orientations are not supported!\n');
                    end
                else
                    fprintf('Non-implemented patient position: %s\n',obj.AcqParameters.PatientPos);
                end
            elseif numel(obj.AcqParameters.Orientation) == 3        % Bruker Pro data orientation is the normal vector
                mainvec = find(max(obj.AcqParameters.Orientation) == obj.AcqParameters.Orientation);
                mainvec = mainvec(1);
                switch mainvec
                    case 1
                            obj.AcqParameters.SliceOrient = 'yz';
                            obj.ProData = flip(rot90(permute(obj.ProData,[1,3,2]),2),3);
                            obj.RecoParameters.RearrangeDims = [3,2,1];
                            fprintf('sagittal\n');                    
                    case 2
                            obj.AcqParameters.SliceOrient = 'xz';
                            obj.ProData = flip(rot90(permute(obj.ProData,[3,1,2]),2),3);
                            obj.RecoParameters.RearrangeDims = [2,3,1];
                            fprintf('coronal\n');                    
                    case 3
                            obj.AcqParameters.SliceOrient = 'xy';
                            obj.ProData = rot90(obj.ProData,-1);
                            obj.RecoParameters.RearrangeDims = [2,1,3];
                            fprintf('transverse\n');                    
                end
            end
        end

        function obj = CoilCombine(obj)
            if obj.AcqParameters.Channels <= 1
                %obj.ProData = obj.RawData;
                return
            end
            obj.RawData = squeeze(obj.RawData);
            s = size(obj.RawData);
            
            if numel(obj.AcqParameters.NoiseCov)>0
                wsvdopt.noiseCov = obj.AcqParameters.NoiseCov;
            else
                wsvdopt.noiseCov = eye(obj.AcqParameters.Channels)*0.5;
            end
           dat = fft(obj.RawData,[],2);
           dat = fft(dat,[],3);
           dat = fft(dat,[],4);

           comb = zeros(s(1:end-1));
            %comb = obj.RawData(:,:,:,:,1);
            for cnt1 = 1:s(2)
                for cnt2 = 1:s(3)
                    for cnt3 = 1:s(4)
                        if numel(s) == 5
                            [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(squeeze(dat(:,cnt1,cnt2,cnt3,:)), [], wsvdopt);
                            comb(:,cnt1,cnt2,cnt3) = wsvdCombination;
                        elseif numel(s) == 6
                            for cnt4 = 1:s(5)
                                [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(squeeze(dat(:,cnt1,cnt2,cnt3,cnt4,:)), [], wsvdopt);
                                comb(:,cnt1,cnt2,cnt3,cnt4) = wsvdCombination;
                            end
                        end
                    end
                end
            end
            obj.RawData = ifft(comb,[],4);
            obj.RawData = ifft(obj.RawData,[],3);
            obj.RawData = ifft(obj.RawData,[],2);
        end

        function obj = ReadWeighting(obj, filename)
            % Reads the weighting,txt file
            if nargin == 2
                [p,n,e] = fileparts(filename);
                if numel(p) == 0
                    [path,~,~] = fileparts(obj.AcqParameters.Filename);
                    filename = fullfile(path,n,e);
                end
            else
                [path,~,~] = fileparts(obj.AcqParameters.Filename);
                filename = fullfile(path,'weighting.txt');
            end
            if exist(filename, "file")== 0
                fprintf('File %s not found!\n',filename);
                return;
            end
            file = fopen(filename);
            vals = fscanf(file,"%d: %d, %d, %d:%d, %f\n");
            fclose(file);
            vals = reshape(vals,6,[]);
            x = vals(2,:);
            y = vals(3,:);
            z = vals(4,:);
            weight = vals(5,:);
            corr = vals(6,:);
            obj.AcqParameters.weighting.weight = zeros([max(x)-min(x)+1,max(y)-min(y)+1, max(z)-min(z)+1]);
            obj.AcqParameters.weighting.corr = zeros([max(x)-min(x)+1,max(y)-min(y)+1, max(z)-min(z)+1]);
            for cnt = 1:numel(x)
                obj.AcqParameters.weighting.weight(x(cnt)-min(x)+1,y(cnt)-min(y)+1,z(cnt)-min(z)+1) = weight(cnt);
                obj.AcqParameters.weighting.corr(x(cnt)-min(x)+1,y(cnt)-min(y)+1,z(cnt)-min(z)+1) = corr(cnt);
            end
        end

        function obj = CSIRearrange(obj)
            s = size(obj.RawData);
            % line, partition, slice, echo, phase, meas, segment, (acq, set), (channel)
            % This is for 3D CSI!!!
            obj.RawData = permute(obj.RawData, [1,2,3,8,7,4,5,6]);    
            obj.AcqParameters.NDims = numel(size(obj.RawData));
        end

        function obj = CSIApplyWeighting(obj)
            s = size(obj.ProData);
            if numel(s) >= 5
                numdat = s(5);
            else
                numdat = 1;
            end
            c = reshape(obj.AcqParameters.weighting.corr,[1,size(obj.AcqParameters.weighting.corr)]);
            for cnt = 1:numdat
                for cnt2 = 1:s(1)
                    obj.ProData(cnt2,:,:,:,cnt1) = obj.ProData(cnt2,:,:,:,cnt1).*c;
                end
            end
        end

        function obj = CSIFullFT(obj, dir)
            if nargin < 2 || dir ~= -1
                dir = 1;
            end
            % Simply does FFT in all dimensions. dir is 1 for fft and -1 for ifft
            s = size(obj.RawData,1,2,3,4,5);
            if dir == 1
                obj.RawData = circshift(fft(fft(fft(fft(obj.RawData,[],1),[],2),[],3),[],4),[-floor(s(1)/2),-floor(s(2)/2),-floor(s(3)/2),-floor(s(4)/2),0]);
            else
                obj.RawData = ifft(ifft(ifft(ifft(circshift(obj.RawData,[floor(s(1)/2),floor(s(2)/2),floor(s(3)/2),floor(s(4)/2),0]),[],1),[],2),[],3),[],4);
             end

        end
        
        function obj = CSISpatialFT(obj, zf, useProData)
            if nargin < 3
                useProData = 0;
            end            
            if useProData == 0
                s = size(obj.RawData);
            else
                s = size(obj.ProData);
            end
            if nargin < 2 || numel(zf) == 0
                zf = s(2:4);
            end

            % After Fourier transform, the data size has to be even in all directions to make the correct shift possible
            if zf(1)/2~=floor(zf(1)/2)
                zf(1) = zf(1)+1;
            end
            if zf(2)/2~=floor(zf(2)/2)
                zf(2) = zf(2)+1;
            end
            if zf(3)/2~=floor(zf(3)/2)
                zf(3) = zf(3)+1;
            end
            
            res = zeros([s(1),zf,s(5:end)]);
            if useProData == 0
                res(:,1:s(2),1:s(3),1:s(4),:) = obj.RawData;
            else
                res(:,1:s(2),1:s(3),1:s(4),:) = obj.ProData;
            end
            res = circshift(res,[0,-floor(s(2)/2),-floor(s(3)/2),-floor(s(4)/2),0]);
            res = fft(fft(fft(res,[],4),[],3),[],2);
            obj.ProData = circshift(res,[0,-floor(zf(1)/2),-floor(zf(2)/2),-floor(zf(3)/2),0]);
            % Orientation of the image:
            % The Gradient matrix (in xyz coordinates: ACQ_grad_matrix in PV6 and ACQ_GradientMatrix in PV360 - in PV 360 ACQ_grad_matrix also exists and is updated, but not used -
            % in patient coordinates PVM_SPackArrGradOrient) has columns (rps) and rows (xyz), which define the order of the dimensions.
            % minus signs change the direction of the corresponding gradient and such the image orientation.
            % A diagonal matrix just leads to a RearrangeDims of [1,2,3]. The different directions are, however, not treated equally: For a positive value in the gradient matrix, 
            % read and phase gradients go from minus to plus, but slice gradient from plus to minus. That means that the z-gradient changes orientation when going in a different direction,
            % as to x or y if they are in slice direction:
            obj.ProData = permute(obj.ProData,[1,abs(obj.RecoParameters.RearrangeDims([2,1,3]))+1,5]);
            if obj.RecoParameters.RearrangeDims(2) < 0
                obj.ProData = flip(obj.ProData,2);
            end
            if obj.RecoParameters.RearrangeDims(1) < 0
                obj.ProData = flip(obj.ProData,3);
            end
            if obj.RecoParameters.RearrangeDims(3) < 0
                obj.ProData = flip(obj.ProData,4);
            end

            % Now set the RecoParameters
            obj.RecoParameters.Zerofill = zf(abs(obj.RecoParameters.RearrangeDims(1:3)));
            obj.RecoParameters.MatrixSize = zf(abs(obj.RecoParameters.RearrangeDims(1:3)));
            obj.RecoParameters.FOV = obj.AcqParameters.FOV(abs(obj.RecoParameters.RearrangeDims(1:3)));
            if max(abs(obj.AcqParameters.CSIPosition)) > 0   % This is the case for Bruker data
                obj.RecoParameters.CSIPosition = obj.AcqParameters.CSIPosition(abs(obj.RecoParameters.RearrangeDims(1:3)));
                obj.RecoParameters.Position = [0,0,0];
            else
                %obj.RecoParameters.Position = obj.AcqParameters.Position(abs(obj.RecoParameters.RearrangeDims(1:3)));  % In Bruker data, CSI is always acquired around the center, and should be shifted during recontruction. Since this is not yet done here, the real position is 0, the shifted position in stored in CSIPosition.
                obj.RecoParameters.Position = obj.AcqParameters.Position;
                obj.RecoParameters.CSIPosition = [0,0,0];
            end
        end

        function fid = CSIGetVoxelData(obj, voxelind, zf)
            % Returns data of a selected voxel. If zf is given or the data has not been Fouriertransformed yet, it first fourier-transforms the data with the given zerofilling.
            if nargin > 2 || numel(obj.ProData) == 0
                if nargin > 2
                    obj = obj.CSISpatialFT(zf);
                else
                    obj = obj.CSISpatialFT();
                end
            end
            fid.raw = squeeze(obj.ProData(:,voxelind(1),voxelind(2),voxelind(3),:));
            fid.params.frequency = obj.AcqParameters.Frequency;
            fid.params.bandwidth = obj.AcqParameters.Bandwidth;
            fid.params.field = obj.AcqParameters.B0;
            return;
        end


        function obj = CSIPhaseAndCenter(obj,SearchRegion, Zerofreq, ppm)
            dophase = 1;
            s = ones(1,5);
            siz = size(obj.ProData);
            s(1:numel(siz)) = siz;
            N = s(1)*4;
            if nargin < 4
                ppm = 0;
            end
            if nargin < 3
                Zerofreq = 0;
            elseif ppm == 1
                Zerofreq = Zerofreq * obj.AcqParameters.Frequency;
            end            
            if nargin < 2 || numel(SearchRegion) == 1
                SearchRegion = [-obj.AcqParameters.Bandwidth/10,obj.AcqParameters.Bandwidth/10];
            end
            SearchRegionInd = floor(SearchRegion*(N-1)/obj.AcqParameters.Bandwidth + (N+1)/2);
            ZeroInd = Zerofreq*N/obj.AcqParameters.Bandwidth + (N+1)/2;
            if s(5) > 1
                fprintf('Adjusting repetition %2d of %2d',0,s(5));
            end
            for cnt5 = 1 : s(5)
            if s(5) > 1
                fprintf('\b\b\b\b\b\b\b\b%2d of %2d',cnt5,s(5));
            end
                
                for cnt4 = 1:s(4)
                    for cnt3 = 1:s(3)
                        for cnt2 = 1:s(2)
                            spec = circshift(fft(obj.ProData(:,cnt2, cnt3, cnt4, cnt5),s(1)*4),s(1)*2);
                            aspec = abs(spec(SearchRegionInd(1):SearchRegionInd(2)));
                            if max(abs(spec(SearchRegionInd(1):SearchRegionInd(2))))/std(abs(spec))> 5   %Only do this if there is a peak
                                mm = find(aspec == max(aspec))+SearchRegionInd(1)-1;
                                if dophase == 1
                                    ph = mean(angle(spec(mm-2:mm+2)));
                                else
                                    ph = 0;
                                end
    
                                diff = mm - ZeroInd;
                                obj.ProData(:,cnt2, cnt3, cnt4, cnt5) = obj.ProData(:,cnt2, cnt3, cnt4, cnt5).*exp(-complex(0,1)*(linspace(0,s(1)-1,s(1))'*diff/N*2*pi+ph));
                            end
                        end
                    end               
                end
            end
            if s(5) > 1
              fprintf('\n');
            end

        end

        function obj = CSICalcPSF(obj)
            finalresol = 1024;
            w = obj.AcqParameters.kspace;
            %w = obj.AcqParameters.weighting.weight; %.* obj.AcqParameters.weighting.corr; 
            s = size(w);
            w = padarray(w,[finalresol-s(1),0,0],0,'post');
            psf = circshift(fft(circshift(w,[-ceil(s(1)/2),0,0])),[-finalresol/2,0,0]);
            %psf = psf(floor(3/8*finalresol)+1:floor(5/8*finalresol),:,:);
            psf = padarray(psf,[0,finalresol-s(2),0],0,'post');
            psf = circshift(fft(circshift(psf,[0,-ceil(s(2)/2),0]),[],2),[0,-finalresol/2,0]);
            %psf = psf(:,floor(3/8*finalresol)+1:floor(5/8*finalresol),:);
            psf = padarray(psf,[0,0,finalresol-s(3)],0,'post');
            psf = circshift(fft(circshift(psf,[0,0,-ceil(s(3)/2)]),[],3),[0,0,-finalresol/2]);
            %psf = psf(:,:,floor(3/8*finalresol)+1:floor(5/8*finalresol));    
            psf = psf/max(max(max(psf)));
            % Detect width
            width = [numel(find(psf(:,round(finalresol/2),round(finalresol/2))>0.64)),numel(find(psf(round(finalresol/2),:,round(finalresol/2))>0.64)),numel(find(psf(round(finalresol/2),round(finalresol/2),:)>0.64))];
            width = width/finalresol;
            fprintf('Psf-width: %f, %f, %f\nResolution: %f, %f, %f\n',width(1),width(2),width(3),width(1)*obj.AcqParameters.FOV(1),width(2)*obj.AcqParameters.FOV(2),width(3)*obj.AcqParameters.FOV(3))
            rim(psf)
            obj.AcqParameters.psf = psf;
        end

        function obj = CSIRotate(obj)
            if max(obj.AcqParameters.Orientation) == obj.AcqParameters.Orientation(1)  % sagittal
                obj.RawData = flip(permute(obj.RawData,[1,2,4,3]),4);
            elseif max(obj.AcqParameters.Orientation) == obj.AcqParameters.Orientation(2)  % coronal
                obj.RawData = flip(flip(permute(obj.RawData,[1,3,4,2]),4),3);
            end
        end

        function obj = Average(obj,sets,useraw)
        %% Averages data over sets (Last dimension)
        % If called without argument: simply average over all sets
        % sets can be an array with indices of all sets to be added or a cell array with several arrays
        % normally averages over the processed data (ProData) if available, otherwise over raw data.
        % If useraw == 1, averages raw data even if processed data is available
        % Careful: all sets that do not appear in the sets-array is lost!
        if nargin < 3
            useraw = 0;
        end
            if numel(obj.ProData) == 0
                useraw = 1;
            end
            if useraw == 1 && numel(obj.RawData) == 0
                fprintf('Averaging: No raw data available\n');
                return;
            end
            if useraw == 1
                s = size(obj.RawData);
            else
                s = size(obj.ProData);
            end
            if numel(s)<5 || s(5) == 1
                fprintf('Averaging: 5th dimension nonexistant or size 1: No sets available to average\n');
            end
            if nargin == 1 || numel(sets) == 0
                if useraw == 1
                    obj.RawData = mean(obj.RawData,5);
                else
                    obj.ProData = mean(obj.ProData,5);
                end
                return
            else   % sets array or cell array given
                if ~iscell(sets)
                    sets = {sets};
                end
                averaged = zeros([s(1:4),numel(sets)]);
                for cnt = 1:numel(sets)
                    if useraw == 1
                        averaged(:,:,:,:,cnt,:,:) = mean(obj.RawData(:,:,:,:,sets{cnt},:,:),5);
                    else
                        averaged(:,:,:,:,cnt,:,:) = mean(obj.ProData(:,:,:,:,sets{cnt},:,:),5);                    
                    end
                end
                if useraw == 1
                    obj.RawData = averaged;
                else
                    obj.ProData = averaged;
                end
            end
        end

        function obj = CSIGlobalFreqCorr(obj,center,peakRegion)
        %% For CSI with several sets (for averaging), before Fourier transformation (works on raw data):
        % Tries to adjust all set to the same frequency: Takes global spectrum from k-space center, searches the largest peak 
        % and shifts the global frequency of all further sets to be equal to the first one.
        % If center is set to one, the biggest peak will be moved to the center of the spectrum.
        % If peakRegion is set to a 2-element array, the highest peak of this spectral region is used 
            if nargin == 1
                center = 0;
            end
            if numel(obj.RawData) == 0
                fprintf('Global frequency correction: requires raw data!\n');
                return
            end
            s = size(obj.RawData);
            if nargin > 2 && numel(peakRegion) == 2     % Peak Region given
                reg = sort(round(s(1)*8-peakRegion/obj.AcqParameters.Bandwidth*s(1)*8));
            else
                reg = [];
            end
            if center == 0 && (numel(s)<5 || s(5)==1)
                fprintf('Global frequency correction only possible for several sets!\n');
                return
            end
            if numel(s) < 5
                s(5) = 1;
            end
            % Find center of k-space
            kcent = ceil(s(2:4)/2);
            centdat = squeeze(mean(obj.RawData(:,kcent(1)-1:kcent(1)+1,kcent(2)-1:kcent(2)+1,kcent(3)-1:kcent(3)+1,:),[2,3,4]));
            N = s(1)*8;
            centdatft = fftshift(fft(centdat,N),1);
            for cnt = 1:s(5)
                if numel(reg) == 0
                    ind = find(abs(centdatft(:,cnt)) == max(abs(centdatft(:,cnt))));                
                else
                    ind = find(abs(centdatft(reg(1):reg(2),cnt)) == max(abs(centdatft(reg(1):reg(2),cnt)))) + reg(1)-1;  
                end
                if cnt == 1
                    if center == 1
                        targ = (N+1)/2;
                    else
                        targ = ind;
                    end
                end
                currphase = mean(angle(centdatft(ind-2:ind+2,cnt)));
                phc = exp(complex(0,1)*(linspace(0,s(1)-1,s(1))*(ind-targ)/N*2*pi+currphase));
                obj.RawData(:,:,:,:,cnt) = obj.RawData(:,:,:,:,cnt).*repmat(phc',[1,s(2:4)]);
            end
        end

        function obj = ChooseSets(obj)
            str1 = '\\mrz10\Upload9T\USERS\Rolf\fMRS\Subject13\meas_MID00522_FID48988_rpcsi_ICE_OFF_9_15_13.dat';
            str2 = '\\mrz10\Upload9T\USERS\Rolf\fMRS\Subject14\meas_MID00735_FID50925_rpcsi_ICE_OFF_9_15_13.dat';
            if isequal(obj.AcqParameters.Filename,str1) == true %requires also change in choosing of sets for ProDataRest and proDataStim (seen in spectraSum.m)
                set = 1:40;
                set(1:4:end) = [];
                set(1:3:end) = [];
                obj.ProData = obj.ProData(:,:,:,:,set);
                a = set(1:2:end);
                cnt = 1;

                for i = 1:length(a)
                    obj.ProData(:,:,:,:,i) =obj.ProData(:,:,:,:,cnt) + obj.ProData(:,:,:,:,cnt+1);%+ csi.RawData(:,:,:,:,cnt+2);  
                    cnt = cnt+2; %=3
                end
                %     if cnt == 31
                %         break;
                %     end            
                obj.ProData = obj.ProData(:,:,:,:,1:10);   
            elseif isequal(obj.AcqParameters.Filename,str2) == true %requires also change in choosing of sets for ProDataRest and proDataStim (seen in spectraSum.m)
                set = 1:32;
                set(1:4:end) = [];
                set(1:3:end) = [];
                obj.ProData = obj.ProData(:,:,:,:,set);
                a = set(1:2:end);
                cnt = 1;

                for i = 1:length(a)
                    obj.ProData(:,:,:,:,i) =obj.ProData(:,:,:,:,cnt) + obj.ProData(:,:,:,:,cnt+1);%+ csi.RawData(:,:,:,:,cnt+2);  
                    cnt = cnt+2; %=3
                end
                %     if cnt == 31
                %         break;
                %     end            
                obj.ProData = obj.ProData(:,:,:,:,1:8);   
            else           
            set = 1:40;
            set(1:4:end) = [];
            set(1:3:end) = [];
            obj.ProData = obj.ProData(:,:,:,:,set);
            a = set(1:2:end);
            cnt = 1;
            
            for i = 1:length(a)
                obj.ProData(:,:,:,:,i) =obj.ProData(:,:,:,:,cnt) + obj.ProData(:,:,:,:,cnt+1);%+ csi.RawData(:,:,:,:,cnt+2);  
                cnt = cnt+2; %=3
            end
            %     if cnt == 31
            %         break;
            %     end            
            obj.ProData = obj.ProData(:,:,:,:,1:10); 
			end
        end

        function obj = SpecShift(obj, index)
            s = size(obj.ProData);
            if nargin == 1
                for cnt2 = 1:s(2)
                    for cnt3 = 1:s(3)
                        for cnt4 = 1:s(4)
                            for cnt5 = 1:s(5)
                                temp = fftshift(fft(obj.ProData(:,cnt2,cnt3,cnt4,cnt5),256*4)); %ask Rolf if its necessary to iterate over all sets
                                xIndex = find(abs(temp) == max(abs(temp)), 1, 'first');
                                phaseramp = exp(1i*linspace(0,255,256)*(xIndex-512)/4/256*2*pi); 
                                obj.ProData(:,cnt2,cnt3,cnt4,cnt5) = obj.ProData(:,cnt2,cnt3,cnt4,cnt5).*phaseramp';
                            end
                        end
                    end
                end
            else
                if isvector(index) == true 
                    for cnt5 = 1:s(5)
                        temp = fftshift(fft(obj.ProData(:,index(1),index(2),index(3),cnt5),256*4)); %ask Rolf if its necessary to iterate over all sets
                        xIndex = find(abs(temp) == max(abs(temp)), 1, 'first');
                        phaseramp = exp(1i*linspace(0,255,256)*(xIndex-512)/4/256*2*pi); 
                        obj.ProData(:,index(1),index(2),index(3),cnt5) = obj.ProData(:,index(1),index(2),index(3),cnt5).*phaseramp';
                    end
                else
                    disp("No Voxel specified or Voxel not in form [x y z]");
                    return;
                end
            end
        end

        function obj = PhaseCor(obj, index)
            s = size(obj.ProData);
            if nargin == 1
                for cnt2 = 1:s(2)
                    for cnt3 = 1:s(3)
                        for cnt4 = 1:s(4)
                            for cnt5 = 1:s(5)
                                temp = fft(obj.ProData(:,cnt2,cnt3,cnt4,cnt5)); 
                                xIndex = find(abs(temp) == max(abs(temp)), 1, 'first');
                                theta = angle(obj.ProData(xIndex,cnt2,cnt3,cnt4,cnt5)); 
                                obj.ProData(:,cnt2,cnt3,cnt4,cnt5) = obj.ProData(:,cnt2,cnt3,cnt4,cnt5)*exp(-1i*theta);
                            end
                        end
                    end
                end
            else
                if isvector(index) == true 
                    for cnt5 = 1:s(5)
                        temp = fft(obj.ProData(:,index(1),index(2),index(3),cnt5)); 
                        xIndex = find(abs(temp) == max(abs(temp)), 1, 'first');
                        theta = angle(obj.ProData(xIndex,index(1),index(2),index(3),cnt5)); 
                        obj.ProData(:,index(1),index(2),index(3),cnt5) = obj.ProData(:,index(1),index(2),index(3),cnt5)*exp(-1i*theta);
                    end
                else
                    disp("No Voxel specified or Voxel not in form of [x y z]");
                    return;
                end
            end
        end

        function obj = CropFID(obj,cropPoints)
                obj.ProData = obj.ProData(cropPoints:end-1,:,:,:,:);
                obj.ProData = padarray(obj.ProData,[cropPoints 0],0,'post'); 
        end

        function obj = CSIDCCorrect(obj,numpoints, RawOrPro)
            % numpoints: number of FID points at the end of the FID used for the correction. If not given or <= 0: 1/32 of fid size
            % RawOrPro: 0: apply to Raw, 1: Pro data , default: Raw
            if nargin < 3 || numel(obj.ProData) == 0
                RawOrPro = 0;
            end
            if nargin < 2
                numpoints = 0;
            end
            s = size(obj.RawData,1);
            if numpoints <= 0 
                numpoints = round(s/32);
            end
            if numpoints > s/2
                fprintf('CSIDCCorrect: number of points %d too high for a number of datapoints of %d.\n',numpoints, s );
                numpoints = round(s/32);
                fprintf('Reduced to %d\n',numpoints);
            end
            if RawOrPro == 0
                corr = mean(obj.RawData(end-12:end,:,:,:,:));
                obj.RawData = obj.RawData - corr;
            else
                corr = mean(obj.ProData(end-12:end,:,:,:,:));
                obj.ProData = obj.ProData - corr;
            end
        end

        function obj = CSIExtrapolate(obj,ext_size, ProOrIndex,FrontOrBack, keepsize)
            %--------------------------------------------------------------------------
            % Linear predicion of the missing points of FID signal measured with
            % CSI-FID sequence.
            %--------------------------------------------------------------------------
            % Usage:
            %       obj = CSIExtrapolate(fid,ext_size)
            % where:
            %              fid - original FID signal;
            %         ext_size - number of points to be predicted. This depends from
            %                    acquisition duration, vector size and delay time of 
            %                    the CSI-FID sequence;
            %       fid_extrap - final fid signal with begin points extrapolated
            %--------------------------------------------------------------------------
            % Function calculates the missing points in the FID signal measured with
            % CSI-FID sequence. Function should be used with non-complex data.
            % Prediction of the missing points can be done separately for the real 
            % and imaginary part of the FID.
            %
            % Adopted from the code on: 
            %                                     https://gist.github.com/tobin/2843661
            %
            % Thread about linear prediction is available under:
            %
            %                                     http://dsp.stackexchange.com/a/110/64
            %
            % Adopted by:                                                     GCH, 2013
            % Last changes:
            % 15.04.2015 - GCH - adopted to ToolboxCSI v2.
            % 
            % modified for 31P -fMRS project
            %--------------------------------------------------------------------------
            if nargin < 5
                keepsize = 1;                 % what is added at the beginning is cut off at the end
            end
            if nargin < 4                        % if FrontOrBack isn't given, add data at the front
                FrontOrBack = 0;
            end
            if FrontOrBack == 1
                keepsize = 0;                   % if something is added at the end, it shouldn't be cut off again
            end
            % if no ext_size is given, calculate from acquisition delay and dwell time
            if (nargin == 1 || ext_size<0) 
                if isfield(obj.AcqParameters,'AcquisitionDelay') 
                    ext_size = ceil(obj.AcqParameters.AcquisitionDelay/obj.AcqParameters.DwellTime);
                else
                    fprintf('Acquisition Delay is not known. Cannot calculate number of extrapolated points.\n');
                end
            end
            s = ones(1,5);
            siz = size(obj.RawData);
            s(1:numel(siz)) = siz;
            useProData = 0;
            singlevoxel = 0;
            if nargin >= 3     % the third argument can be either 0 or 1 for use RawData or ProData, or the (3D-) index of the voxel to be treated
                if numel(ProOrIndex) == 1
                    if ProOrIndex == 1 && numel(size(obj.ProData))>3
                        useProData = 1;
                        s = ones(1,5);
                        siz = size(obj.ProData);
                        s(1:numel(siz)) = siz;
                    end
                else
                    index = ProOrIndex;
                    singlevoxel = 1;
                end
            end
            if keepsize == 0
                    newsiz = siz;
                    newsiz(1) = newsiz(1)+ext_size;
                    newdata = zeros(newsiz);               
            end
            if singlevoxel == 0
                if s(5) > 1
                        fprintf('Extrapolating repetition %2d of %2d',0,s(5));
                end
                for cnt5 = 1:s(5)
                    if s(5) > 1
                        fprintf('\b\b\b\b\b\b\b\b%2d of %2d',cnt5,s(5));
                    end
                    for cnt2 = 1:s(2)
                        for cnt3 = 1:s(3)
                            for cnt4 = 1:s(4)
                                %for cnt5 = 1:s(5)
                                    %fid = obj.ProData(:,cnt2,cnt3,cnt4,cnt5);
                                    if useProData == 1
                                        fid = obj.ProData(:,cnt2,cnt3,cnt4,cnt5);
                                    else
                                        fid = obj.RawData(:,cnt2,cnt3,cnt4,cnt5);
                                    end
                                    if fid(10) ~= 0.0
                                        LPC_start = length(fid);          % Starting point for LPC 
                                        LPC_order = length(fid)-1;        % Order of LPC autoregression
                                        LPC_datax = length(fid)+ext_size; % Size of expadned data set
                                        
                                        % FID signal is rotated so that the begin goes to the end and the LPC
                                        % algorithm can predict the end
                                        if FrontOrBack == 0
                                            fid_rot = rot90(fid,2);
                                        else
                                            fid_rot = fid;
                                        end
                                        
                                        % Autoregression with Burg method 
                                        % from: Signal Processing Toolbox / Parametric Modelling
                                        LPC_arb = arburg(fid_rot,LPC_order);
                                        %LPC_arb = lpc(fid_rot,LPC_order);
                                        
                                        % Below is the expanded dataset defined:
                                        fid_new = zeros(LPC_datax,1);
                                        fid_new(1:LPC_start) = fid_rot;
                                        
                                        % Initial FID is run through the filter to get the LPC coeficients
                                        [~, zf] = filter(-[0 LPC_arb(2:end)], 1, fid_rot(1:LPC_start));
                                        
                                        % Filter is used as an IIR (infinite impulse response) to extrapolate the
                                        % missing data points
                                        fid_new((LPC_start+1):LPC_datax) = filter([0 0], -LPC_arb, zeros(LPC_datax-LPC_start,1), zf);
                                        
                                        % Final result must be rotated back in order to match the imput time series
                                        if FrontOrBack == 0
                                            fid_0 = rot90(fid_new,2);
                                        else
                                            fid_0 = fid_new;
                                        end
                                        if keepsize == 1
                                            fid_extrap = fid_0(1:LPC_start);
                                            if useProData == 1
                                                obj.ProData(:,cnt2,cnt3,cnt4,cnt5) = fid_extrap;
                                            else
                                                obj.RawData(:,cnt2,cnt3,cnt4,cnt5) = fid_extrap;
                                            end
                                        else
                                            newdata(:,cnt2,cnt3,cnt4,cnt5) = fid_0;
                                        end
                                    end
                                %end
                            end
                        end
                    end
                end
                if s(5) > 1
                    fprintf('\n');
                end
                if keepsize == 0
                    if useProData == 1
                        obj.ProData = newdata;
                    else
                        obj.RawData = newdata;
                    end
                else
                    obj.RecoParameters.ExtrapolatedPoints = ext_size;
                end
            else
                if isvector(index) == true
                    for cnt5 = 1:s(5)
                        fid = obj.ProData(:,index(1),index(2), index(3),cnt5);
                                    
                        LPC_start = length(fid);          
                        LPC_order = length(fid)-1;       
                        LPC_datax = length(fid)+ext_size; 
                        fid_rot = rot90(fid,2);
                        LPC_arb = arburg(fid_rot,LPC_order);
                        fid_new = zeros(LPC_datax,1);
                        fid_new(1:LPC_start) = fid_rot;                 
                        [~, zf] = filter(-[0 LPC_arb(2:end)], 1, fid_rot(1:LPC_start));                    
                        fid_new((LPC_start+1):LPC_datax) = filter([0 0], -LPC_arb, zeros(LPC_datax-LPC_start,1), zf);                    
                        fid_0 = rot90(fid_new,2);
                        fid_extrap = fid_0(1:LPC_start);
                        obj.ProData(:,index(1),index(2), index(3),cnt5) = fid_extrap;
                    end
                else
                    disp("No Voxel specified or Voxel not in Form [x y z]");
                    return;
                end
            end
        end	

        function obj = CSIDenoisePCASingleVoxel(obj,voxel,strength,patchsize)
            s = size(obj.ProData,1,2,3,4,5);
            patchpos0 = [max([1,voxel(1)-floor(patchsize(1)/2)]),max([1,voxel(2)-floor(patchsize(2)/2)]),max([1,voxel(3)-floor(patchsize(3)/2)]),1];
            patchpos1 = [min([s(2),voxel(1)+floor(patchsize(1)/2)]),min([s(3),voxel(2)+floor(patchsize(2)/2)]),min([s(4),voxel(3)+floor(patchsize(3)/2)]),s(5)];
            patch = obj.ProData(:,patchpos0(1):patchpos1(1),patchpos0(2):patchpos1(2),patchpos0(3):patchpos1(3),patchpos0(4):patchpos1(4));
            ps = size(patch);
            M = [reshape(real(patch),s(1),[]),reshape(imag(patch),s(1),[])];
            [U,S,V] = svd(M);
            S(round(s(1)/strength):end,:) = 0;
            Mnew = U*S*V';
            pnew = complex(reshape(Mnew(:,1:size(Mnew,2)/2),ps),reshape(Mnew(:,size(Mnew,2)/2+1:end),ps));
            obj.ProData(:,patchpos0(1):patchpos1(1),patchpos0(2):patchpos1(2),patchpos0(3):patchpos1(3),patchpos0(4):patchpos1(4)) = reshape(pnew,ps);
        end

        function obj = CSIDenoisePCAFull(obj,strength,patchsize)
            if nargin < 4
                doFT = 0;
            end
            if nargin < 3 || numel(patchsize)==0
                patchsize = [5,5,5];
            end
            s = size(obj.RawData,1,2,3,4,5);
            datnew = complex(zeros(s),zeros(s));
            tic;
            totnum = s(2)*s(3)*s(4);
            outputstr = sprintf('Denoising: %5.1f %% done (voxel %2d, %2d, %2d), within %5.0f s',0,1,1,1,toc);
            outputlen = numel(outputstr);
            back = sprintf('\b');
            back = repmat(back,outputlen,1);
            fprintf('%s',outputstr);
            step = 0;
            for cnt2 = 1:s(2)
                for cnt3 = 1:s(3)
                    for cnt4 = 1:s(4)
                        patchpos0 = [max([1,cnt2-floor(patchsize(1)/2)]),max([1,cnt3-floor(patchsize(2)/2)]),max([1,cnt4-floor(patchsize(3)/2)]),1];
                        patchpos1 = [min([s(2),cnt2+floor(patchsize(1)/2)]),min([s(3),cnt3+floor(patchsize(2)/2)]),min([s(4),cnt4+floor(patchsize(3)/2)]),s(5)];
                        patch = obj.RawData(:,patchpos0(1):patchpos1(1),patchpos0(2):patchpos1(2),patchpos0(3):patchpos1(3),patchpos0(4):patchpos1(4));
                        ps = size(patch);
                        M = [reshape(real(patch),s(1),[]),reshape(imag(patch),s(1),[])];
                        [U,S,V] = svd(M);
                        % The following section can be used to calculate the number of singular values to be retained: In how many points is lefts > rights
%                         if cnt2 >= patchsize(1) && cnt2 < s(2)-patchsize(1) && cnt3 >= patchsize(2) && cnt3 < s(3)-patchsize(2) && cnt4 >= patchsize(3) && cnt4 < s(4)-patchsize(3)
%                             sM = size(M);
%                             lambda = diag(S).^2/sM(2);
%                             for p = 1:sM(1)-1
%                                 lefts(p,cnt2-patchsize(1)+1,cnt3-patchsize(2)+1,cnt4-patchsize(3)+1) = sum(lambda(p+1:end));
%                                 rights(p,cnt2-patchsize(1)+1,cnt3-patchsize(2)+1,cnt4-patchsize(3)+1) = (sM(1)-p)*(lambda(p+1)-lambda(end))/4/sqrt((sM(1)-p)/sM(2));
%                             end
%                         end
                        S(round(s(1)/strength+1):end,:) = 0;
                        Mnew = U*S*V';
                        pnew = complex(reshape(Mnew(:,1:size(Mnew,2)/2),ps),reshape(Mnew(:,size(Mnew,2)/2+1:end),ps));
                        datnew(:,patchpos0(1):patchpos1(1),patchpos0(2):patchpos1(2),patchpos0(3):patchpos1(3),patchpos0(4):patchpos1(4)) = datnew(:,patchpos0(1):patchpos1(1),patchpos0(2):patchpos1(2),patchpos0(3):patchpos1(3),patchpos0(4):patchpos1(4))  + reshape(pnew,ps);
                        step = step + 1;
                        fprintf('%sDenoising: %5.1f %% done (voxel %2d, %2d, %2d), within %5.0f s',back,step/totnum*100,cnt2,cnt3,cnt4,toc);
                    end
                end
            end
            fprintf('\n');

            obj.RawData = datnew;
        end

        function obj = CSIDeconvolve(obj, type, param,spec)
            % Type can be 'Wiener', 'Reg','Lucy'
            % param is an additional parameter, often around 10000 for Wiener and Reg, and the number of iterations for Lucy (default would be 10)
            % if spec is 1, the data is Fourier-transformed before deconvolution, otherwise the FIDs are used.
            s = size(obj.ProData,1,2,3,4,5);
            sk = size(obj.AcqParameters.kspace);
            if nargin < 4
                spec = 0;
            end
            if spec == 1
                obj.ProData = fft(obj.ProData);
            end
            kspace = zeros(s(2:4));
            kspace(1:sk(1), 1:sk(2), 1:sk(3)) = obj.AcqParameters.kspace;
            kspace = circshift(kspace,[-floor(sk(1)/2),-floor(sk(2)/2),-floor(sk(3)/2),0]);
            kspace = fft(fft(fft(kspace,[],1),[],2),[],3);
            kspace = abs(circshift(kspace,[-floor(s(2)/2),-floor(s(3)/2),-floor(s(4)/2)]));
            kspace = permute(kspace,abs(obj.RecoParameters.RearrangeDims([2,1,3])));
            kspace = kspace/max(kspace,[],'all');
            mag = abs(obj.ProData);
            noisevar = var(mag(end-3:end,:,:,:,:),0,"all");
            siz = size(mag);
            ss = siz(2)*siz(3)*siz(4);
            for cnt1 = 1:s(1)
                for cnt5 = 1:s(5)
                    % First option: Wiener deconvolution. Doesn't work, bacause either noise is high (for the given nsr) or the PSF is broader (for very high nsr)
                    if strcmp(type,'Wiener')
                    %imvar = var(squeeze(mag(cnt1,:,:,:,cnt5)),0,"all");
                    %nsr = noisevar/imvar;
                        mag(cnt1,:,:,:,cnt5) = deconvwnr(squeeze(mag(cnt1,:,:,:,cnt5)),kspace,param);
                    elseif strcmp(type,'Reg')
                    % Second option: regularized deconvolution
                        mag(cnt1,:,:,:,cnt5) = deconvreg(squeeze(mag(cnt1,:,:,:,cnt5)),kspace,param);
                    elseif strcmp(type,'Lucy')

                    % Third option: 
                        mag(cnt1,:,:,:,cnt5) = deconvlucy(squeeze(mag(cnt1,:,:,:,cnt5)),kspace,10,param);  %,0.01,[],noisevar,2);
                    else
                        printf('Wrong deconvolution type: %s!\n',type);
                        return
                    end
                end
            end
            obj.ProData = mag.*exp(complex(0,1)*angle(obj.ProData));
            if spec == 1
                obj.ProData = ifft(obj.ProData);
            end

        end

    end
end