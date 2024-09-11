classdef MetCon_CSI<matlab.mixin.Copyable
    properties
        time  %[s]

        FieldMap %B0 spatial [rad/s]
        mask

        twix
        DMIPara
        flags

        filename

        metabolites
        sig % raw data
        img % reconstructed image [CHAxCOLxLINxPARxSLCxREP]
        coilSens %%[CHAxCOLxLINxPARxSLC]
        coilNormMat %%[COLxLINxPARxSLC]


        SolverObj
        Metcon
        Experimental

        D % noise decoraltion matrix

    end
    methods

        %constructor: get full path of dat file or no arguments
        function obj=MetCon_CSI(varargin)
            if(nargin==0)
                %path='D:\Data\Spiral\20200918_B0test_sameres';
                path='M:\Subject_data\20201020_subject4847';
                [fn, pathname, ~] = uigetfile(strcat(path,'\*.dat'), 'Pick a DATA file');
                obj.filename=fullfile(pathname,fn);
                obj.twix = mapVBVD(obj.filename,'rmos');
            else
                if(isstruct(varargin{1})||iscell(varargin{1}))
                    obj.twix =varargin{1};
                elseif(isfile(varargin{1}))
                    obj.filename=varargin{1};
                    obj.twix = mapVBVD(obj.filename, 'rmos');
                elseif(isfolder(varargin{1}))
                    path=varargin{1};
                    [fn, pathname, ~] = uigetfile(strcat(path,'\*.dat'), 'Pick a DATA file');
                    obj.filename=fullfile(pathname,fn);
                    obj.twix = mapVBVD(obj.filename, 'rmos');
                else
                    error('wrong input or file not found')
                end


            end

            if(iscell(obj.twix)) %VE
                obj.DMIPara=getDMIPara(obj.twix);
                TW1=obj.twix{1};
                obj.twix=obj.twix{end};
                obj.getflags(varargin{2:end});
                % multi-raid twix files has noise data in the first twix
                obj.calcNoiseDecorrMatrix(TW1);
            else
                obj.DMIPara=getDMIPara(obj.twix);
                obj.getflags(varargin{2:end});
            end

            obj.performRecon();
            obj.getMask(0); % select all voxels

        end
        function getflags(obj,varargin)
            switch nargin
                case 2
                    obj.flags=varargin{1};
                otherwise %parse name valur pairs
                    p=inputParser;
                    p.KeepUnmatched=1;

                    addParameter(p,'doCoilCombine','adapt1',@(x) any(strcmp(x,{'none','sos','adapt1','wsvd'})));

                    addParameter(p,'doNoiseDecorr',true,@(x) islogical(x));
                    addParameter(p,'NormNoiseData',false,@(x) islogical(x));

                    addParameter(p,'CoilSel',1:obj.twix.image.NCha, @(x) isvector(x));
                    addParameter(p,'PCSel',1:obj.twix.image.NRep, @(x) (isvector(x)&& all(x<=1:obj.twix.image.NRep)) );
                    addParameter(p,'EchoSel',1:obj.twix.hdr.MeasYaps.sSpecPara.lVectorSize, @(x) (isvector(x)&& all(x<=obj.twix.hdr.MeasYaps.sSpecPara.lVectorSize)) );
                    addParameter(p,'is3D',(obj.twix.image.NPar>1),@(x)islogical(x));
                    addParameter(p,'doPhaseCorr','none',@(x) any(strcmp(x,{'none','Manual','Burg'})));
                    %                   addParameter(p,'precision','single',@(x) any(strcmp(x,{'single','double'})));
                    %-1 for debug
                    addParameter(p,'doDenosing',0,@(x) isscalar(x));
                    addParameter(p,'Solver','IDEAL',@(x) any(strcmp(x,{'phaseonly','IDEAL','AMARES','LorentzFit'})));
                    addParameter(p,'IdealSetttings',{},@(x)iscell(x) && mod(length(x),2)==0 );
                    %                     addParameter(p,'maxit',10,@(x)isscalar(x));
                    %                     addParameter(p,'tol',1e-6,@(x)isscalar(x));
                    %                     addParameter(p,'reg','none',@(x) any(strcmp(x,{'none','Tikhonov'})));
                    %                     addParameter(p,'reg_lambda',0,@(x)isscalar(x));
                    addParameter(p,'parfor',true,@(x) islogical(x));
                    addParameter(p,'doZeroPad',[1 1 1 0], @(x) isvector(x)); % [3 physical x 1 time]
                    parse(p,varargin{:});

                    obj.flags=p.Results;
                    if(isfield(p.Unmatched,'csm')) %[CHAxCOLxLINxPARxSLC]
                        obj.coilSens=p.Unmatched.csm;
                        obj.flags.doCoilCombine='none';
                        %                         obj.SpiralPara.R_PE=2;
                        %                         obj.flags.doPAT='CGSENSE';
                    end
                    if(isfield(p.Unmatched,'fm'))
                        obj.FieldMap=p.Unmatched.fm;
                    else
                        obj.FieldMap=[];
                    end
                    if(isfield(p.Unmatched,'mask'))
                        obj.mask=p.Unmatched.mask;
                    end
                    if(isfield(p.Unmatched,'metabolites'))
                        obj.metabolites=p.Unmatched.metabolites;
                    end
                    if(isfield(p.Unmatched,'phaseoffset'))
                        obj.flags.phaseoffset=p.Unmatched.phaseoffset;
                    else
                        Acqdelay=obj.twix.hdr.Phoenix.sSpecPara.lAcquisitionDelay*1e-6; %s
                        obj.flags.phaseoffset=[0 2*pi*Acqdelay];
                    end


            end

            if( strcmp(obj.flags.Solver,'AMARES'))
%                 assert(strcmp(obj.flags.doPhaseCorr,'Burg'),'use ''doPhaseCorr'' flag with ''Burg'' method for ''AMARES''');
                assert(~isempty(which('AMARES.amaresFit')),'add AMARES to path: <a href="https://github.com/OXSAtoolbox/OXSA.git">OXSA git link</a>');
            end

        end


        function performNoiseDecorr(obj)

            if(obj.flags.doNoiseDecorr)
                if(isempty(obj.D))
                    obj.calcNoiseDecorrMatrix(obj.twix);
                end
                if(ismatrix(obj.D) )
                    sz     = size(obj.sig);
                    obj.sig   = obj.D*obj.sig(:,:);
                    obj.sig    = reshape(obj.sig,sz);
                end
            end
        end

        function calcNoiseDecorrMatrix(obj,twix)
            if (~isempty(obj.D))
                %skip if decorrelation matrix is set by getflags()
                return;
            elseif(isfield(twix,'noise'))
                noise                = permute(twix.noise(:,obj.flags.CoilSel,:),[2,1,3]);
                noise                = noise(:,:)';
                R                    = cov(noise);
                R(eye(size(R,1))==1) = abs(diag(R));
                if(obj.flags.NormNoiseData)
                    R= R./mean(abs(diag(R)));
                    obj.D               = sqrtm(inv(R));
                else
                    scale_factor=1; %dwell time are the same
                    Rinv = inv(chol(R,'lower'));
                    obj.D = Rinv*sqrt(2)*sqrt(scale_factor);

                end
            else
                warning('no Noise data');
                obj.D = 1;
            end
        end

        function performRecon(obj)
            tic;
            print_str='';
            fprintf('starting reco\n');
            obj.getSig();
            obj.performNoiseDecorr();

            obj.performZeroPaddding();
            %scale
            obj.img=myfft(obj.sig,[2 3 4]).*sqrt(numel(obj.sig));%./(sqrt(sum(obj.sig(:)>0)./numel(obj.sig))*sqrt((1685+667+252+49)./(1685)));
            obj.performCoilCombination();
            obj.performPhaseCorr();
            obj.performSVDdenosing();
            print_str = sprintf( 'reco  time = %6.1f s\n', toc);
            fprintf(print_str);
            obj.getMask(0);
            if(~isempty(obj.metabolites))
                obj.performMetCon();
            end
            fprintf(sprintf( 'Metabolite mapping time = %6.1f s\n', toc));


        end

        function getSig(obj)
            obj.sig = obj.twix.image(obj.flags.EchoSel,obj.flags.CoilSel,:,:,:,:,:,:,:,:,:,:,:,:,:);
            % scale factor for SNR units reco : step 1
            obj.sig=obj.sig./(sqrt(sum(abs(obj.sig(:))>0)));

            % Sum averages and repetitions(if not already done)
            obj.sig  = sum( obj.sig ,[6 9]);

            % Permute the data
            % COL x CHAx LIN x PAR x SLC x AVE x PxExMxSetsx SEG             ->  LIN  xSEG x PAR CHA xCOL
            obj.sig = permute( obj.sig ,[2 3 11 4  1 [5:6 7:10]]);

            ksp=obj.twix.hdr.Phoenix.sKSpace;
            MissingVOXLIN  = double(ksp.lBaseResolution -obj.twix.image.NLin);
            MissingVOXPAR = double(ksp.lPhaseEncodingLines-  obj.twix.image.NPar);
            MissingVOXSEG = double(ksp.lPartitions- obj.twix.image.NSeg);
            obj.sig  = padarray( obj.sig ,[0,MissingVOXLIN,MissingVOXSEG,MissingVOXPAR,0,0],'post');
        end
        function performZeroPaddding(obj)

            disp(['initial CSI data size:           ', num2str(size(obj.sig))])
            zp_PRS=obj.flags.doZeroPad(1:3);
            pad_size=[0 round(size(obj.sig,2)*zp_PRS(1)) round(size(obj.sig,3)*zp_PRS(2))  round(size(obj.sig,4)*zp_PRS(3)) 0];

            obj.sig=padarray(obj.sig,pad_size,0,'both'); % spatial zeropad
            obj.sig=padarray(obj.sig,[zeros(1,4) round(size(obj.sig,5)*obj.flags.doZeroPad(4))],0,'post'); % spectral zerpoad

            disp(['final CSI data size:           ', num2str(size(obj.sig))])
        end
        function performPhaseCorr(obj)
            if(strcmpi(obj.flags.doPhaseCorr,'none'))
                return;
            elseif (strcmpi(obj.flags.doPhaseCorr,'Burg'))
                fprintf('Starting Phase Correction using Burg method \n')
                %use extrapolation
                dw=obj.DMIPara.dwell; % s
                %  Acqdelay=mcobj.twix.hdr.MeasYaps.alTE{1}*1e-6; %s
                Acqdelay=obj.twix.hdr.Phoenix.sSpecPara.lAcquisitionDelay*1e-6; %s
                ext_size=round(Acqdelay/dw);

                [obj.img] = fidExtrp(permute(obj.img,[5 1 2 3 4]),ext_size);
                obj.img=ipermute(obj.img,[5 1 2 3 4]);
                phi0=obj.flags.phaseoffset(1);
                obj.img=obj.img.*exp(-1i*phi0);
                fprintf('Done .... \n')
            elseif (strcmpi(obj.flags.doPhaseCorr,'manual'))
                dw=obj.DMIPara.dwell; % s
                n=size(obj.img,5);
                faxis=(-n/2:(n/2-1) )*(1/(dw*n));
                obj.img=fftshift(fft(obj.img,[],5),5);
                phi0=obj.flags.phaseoffset(1);
                phi1=obj.flags.phaseoffset(2);
                obj.img=obj.img.*reshape(exp(-1i*(phi0+faxis(:).*phi1)),1,1,1,1,[]);
                obj.img=ifftshift(ifft(obj.img,[],5),5);
            end
        end

        function performCoilCombination(obj)

            if(~(size(obj.img,1)>1))
                obj.flags.doCoilCombine='none';
            end
            %coil combination
            switch(obj.flags.doCoilCombine) %{'none','sos','adapt','adapt2'}
                case 'sos'
                    obj.img(1,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = sqrt(sum(abs(obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep)).^2,1));
                case 'adapt1' %preserves phase better?

                    [~,obj.coilSens,obj.coilNormMat]=adaptiveCombine(sum(obj.img(:,:,:,:,:,:),[5 6]));

                    obj.img=(sum(obj.coilSens.*obj.img,1));
                case 'adapt2'
                    error('not implemented')
                case 'wsvd'


                    %% Combine coil data
                    CSI_DataSize=size(obj.img);
                    %we already have noise decorrelated data!
                    CSI_wsvdOption.noiseCov         =0.5*eye(CSI_DataSize(1));
                    CSI_Data=reshape(obj.img,size(obj.img,1),[],size(obj.img,5));
                    CSI_Data=permute(CSI_Data,[2 3 1]);
                    CSI_Combined=zeros([prod(CSI_DataSize(2:4)),CSI_DataSize(5)]);
                    CoilWeights = zeros([prod(CSI_DataSize(1:3)),CSI_DataSize(1)]);
                    if(obj.flags.parfor)
                        parfor vx=1:size(CSI_Data,1)
                            rawSpectra=squeeze(CSI_Data(vx,:,:,1));  % Spectrum x Coil
                            [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(rawSpectra, [], CSI_wsvdOption);
                            CoilWeights(vx,:) = wsvdWeights;
                            CSI_Combined(vx,:)=wsvdCombination;
                        end
                    else
                        for vx=1:size(CSI_Data,1)
                            rawSpectra=squeeze(CSI_Data(vx,:,:,1));  % Spectrum x Coil
                            [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(rawSpectra, [], CSI_wsvdOption);
                            CoilWeights(vx,:) = wsvdWeights;
                            CSI_Combined(vx,:)=wsvdCombination;
                        end
                    end

                    % Combine later?
                    %  CSI = sum(bsxfun (@times,CSI_Data_Filtered,permute(CoilWeights,[1 3 2])),3);
                    obj.coilSens=reshape(CoilWeights.',[CSI_DataSize(1:4)]);
                    obj.img=reshape(CSI_Combined,[1 CSI_DataSize(2:5)]);


            end

        end

        function performSVDdenosing(obj)
            if(all(obj.flags.doDenosing==0))
                return;
            end
            Ncomp=obj.flags.doDenosing;
            imsz=size(obj.img);
            mat=reshape(squeeze(obj.img),[],prod(imsz(5:end)));
            mat_mean=mean(mat,'all','omitnan');
            [U,S,V]=svd(mat-mat_mean,'econ');
            S=diag(S);
            if(Ncomp<0)

                figure,plot(S);
                Ncomp=input('give Number of component to keep');
            end
            S((Ncomp+1):end)=0;

            S=diag(S);
            % figure,plot(S);
            obj.img= reshape(U*S*V',imsz)+mat_mean;
        end



        function cMask=getMask(obj,thres_prctl)
            if(nargin==1)
                thres_prctl=95;
            end
            %for initialization
            if(thres_prctl==0 && isempty(obj.mask) )
                cMask=ones(size(obj.img,2),size(obj.img,3),size(obj.img,4),'logical');
            elseif(thres_prctl>0)

                im_abs=abs(squeeze(sum(obj.img,[5 6 7])));
                cMask=im_abs>prctile(im_abs(:),thres_prctl);
                voxel_size=1e3*mean(obj.DMIPara.resolution/(obj.flags.doZeroPad(1:3)+1)); %mm
                cMask=imerode(cMask,strel('sphere',ceil(15/voxel_size))); % 15 mm
                cMask=imdilate(cMask,strel('sphere',ceil(30/voxel_size))); %30 mm
            elseif(isscalar(obj.mask))
                thres_prctl=obj.mask;
                im_abs=abs(squeeze(sum(obj.img,[5 6 7])));
                cMask=im_abs>prctile(im_abs(:),thres_prctl);
                voxel_size=1e3*mean(obj.DMIPara.resolution/(obj.flags.doZeroPad(1:3)+1)); %mm
                cMask=imerode(cMask,strel('sphere',ceil(15/voxel_size))); % 15 mm
                cMask=imdilate(cMask,strel('sphere',ceil(30/voxel_size))); %30 mm

            else 
                return;
            end

            if(nargout==0)
                obj.mask=cMask;
            end

        end
        function performMetCon(obj)

            freq=[obj.metabolites.freq_shift_Hz];
            fids=MetCon_CSI.mat2col(permute(obj.img,[2 3 4 5 1]),obj.mask);
            nMet=length(freq);
            if(strcmpi(obj.flags.Solver,'AMARES'))

                fprintf('Performing AMARES fitting on %d voxels\n', size(fids,1));

                  if(sum(obj.FieldMap,'all')==0 || isempty(obj.FieldMap))
                    fm_est= zeros(size(fids,1),1);
                else
                    fm_est=-1*obj.FieldMap;
                    fm_est=obj.mat2col(fm_est,obj.mask);
                    fm_est=fm_est./(2*pi)*(6.536 /42.567); % [Hz]   
                    fm_est=fm_est./obj.DMIPara.ImagingFrequency_MHz; %ppm
                  end

                [expParams,pk]=getAMARES_structs(obj);
                metabol_con=zeros(size(fids,1),nMet*3+3);
                %             output format (4th dim) : use xFit order [chemicalshift x N] [linewidth xN] [amplitude xN] [phase x1] fitStatus.relativeNorm fitStatus.resNormSq]
                 if(obj.flags.parfor)
                    parfevalOnAll(@warning,0,'off','all');
                    parfor i=1:size(fids,1)
                        expParams_2=expParams;
                        expParams_2.offset=fm_est(i);
                        [fitResults, fitStatus,~,CRBResults] = AMARES.amaresFit(double(fids(i,:).'), expParams_2, pk, 0,'quiet',true);
                         metabol_con(i,:)=[fitResults.chemShift,fitResults.linewidth,fitResults.amplitude,...
                            fitResults.phase(1),fitStatus.relativeNorm(1),fitStatus.resNormSq(1)];
                    end
                    parfevalOnAll(@warning,0,'on','all');
                 else
                     warning('off','all');
                   for i=1:size(fids,1)
                        expParams.offset=fm_est(i);
                        [fitResults, fitStatus,~,CRBResults] = AMARES.amaresFit(double(fids(i,:).'), expParams, pk, 0,'quiet',true);
                        metabol_con(i,:)=[fitResults.chemShift,fitResults.linewidth,fitResults.amplitude,...
                            fitResults.phase(1),fitStatus.relativeNorm(1),fitStatus.resNormSq(1)];
                   end
                   warning('on','all');
                 end

                % convert to 3D matrix and give resonable names
                metabol_con=MetCon_CSI.col2mat(metabol_con,obj.mask);

                obj.Experimental.chemicalshift=metabol_con(:,:,:,((1:nMet)+nMet*(1-1)))*obj.DMIPara.ImagingFrequency_MHz; %chemical shift in Hz
                obj.Metcon=metabol_con(:,:,:,((1:nMet)+nMet*(3-1)));  %amplitudes
                obj.Experimental.linewidth=metabol_con(:,:,:,((1:nMet)+nMet*(2-1)));
                obj.Experimental.phase=metabol_con(:,:,:,((1)+nMet*(4-1))); %phase
                obj.Experimental.relativeNorm=metabol_con(:,:,:,((2)+nMet*(4-1)));
                obj.Experimental.residue=metabol_con(:,:,:,((3)+nMet*(4-1)));
                obj.Experimental.fm_est=obj.Experimental.chemicalshift(:,:,:,1)*obj.DMIPara.ImagingFrequency_MHz; %water peak in Hz

            elseif(strcmpi(obj.flags.Solver,'LorentzFit'))
                dw=obj.DMIPara.dwell;
                faxis=linspace(-0.5/dw,0.5/dw,size(obj.img,5));
                % nlorentz fitting

                % make the fitting faster by limiting the range +-200 Hz
                [~,minIdx]=min(abs(faxis-(min(freq)-200)));
                [~,maxIdx]=min(abs(faxis-(min(freq)+200)));
                freqSel=minIdx:maxIdx; %1:length(faxis)

                xdata=faxis(freqSel);
                spectrum_All=fftshift(fft(fids,[],2),2);
                spectrum_All=1*double(abs(spectrum_All(:,freqSel)));

                wbhandle = waitbar(0,'Performing nLorentzian fitting');
                metabol_con=zeros(size(spectrum_All,1),nMet);
                chemicalshift=zeros(size(spectrum_All,1),nMet);
                gamma=zeros(size(spectrum_All,1),nMet);
                rSQ=zeros(size(spectrum_All,1),1);
                sRMSE=zeros(size(spectrum_All,1),1);
                %             output format (4th dim) : use xFit order [chemicalshift x N] [linewidth xN] [amplitude xN] [phase x1] fitStatus.relativeNorm fitStatus.resNormSq]

                for i=1:size(spectrum_All,1)
                    %                 if(mod(i,1000)==0), fprintf('%.0f %d done\n',i/size(fids,1)*100); drawnow(); end
                    waitbar(i/size(fids,1),wbhandle,'Performing nLorentzian fitting');

                    [fitf,gof,fitoptions]=NLorentzFit(xdata(:),spectrum_All(i,:).',freq);


                    rSQ(i)=gof.adjrsquare;
                    sRMSE(i)=gof.rmse;
                    chemicalshift(i,:)=[fitf.center1; fitf.center2; fitf.center3;];

                    coeff=coeffvalues(fitf{1});
                    coeff=reshape(coeff(1:end-1),nMet,[]);
                    metabol_con(i,:)= coeff(:,1);
                    chemicalshift(i,:)= coeff(:,2);
                    gamma(i,:)= coeff(:,2);



                end
                close(wbhandle);

                % convert to 3D matrix and give resonable names
                obj.Metcon=MetCon_CSI.col2mat(metabol_con,obj.mask);
                obj.Experimental.rSQ=MetCon_CSI.col2mat(rSQ,obj.mask);
                obj.Experimental.residue=MetCon_CSI.col2mat(sRMSE,obj.mask);
                obj.Experimental.chemicalshift=MetCon_CSI.col2mat(chemicalshift,obj.mask);
                obj.Experimental.gamma=MetCon_CSI.col2mat(gamma,obj.mask);

            elseif(strcmpi(obj.flags.Solver,'IDEAL'))
                obj.SolverObj=IDEAL(obj.metabolites,obj.DMIPara.TE,'fm',obj.FieldMap,'solver',obj.flags.Solver, ...
                    'mask',obj.mask,'parfor',obj.flags.parfor,obj.flags.IdealSetttings{:});
                obj.Metcon=obj.SolverObj'*squeeze(obj.img);
                obj.Experimental.fm_est=obj.SolverObj.experimental.fm_est;
                obj.Experimental.residue=sum(obj.SolverObj.experimental.residue.^2,[4 5]); %squared sum of residual
                obj.SolverObj.experimental=[];
            elseif(strcmpi(obj.flags.Solver,'phaseonly'))
                obj.SolverObj=IDEAL(obj.metabolites,obj.DMIPara.TE,'fm',obj.FieldMap,'solver',obj.flags.Solver,'maxit',10,'mask',obj.mask,'parfor',obj.flags.parfor);
                obj.Metcon=obj.SolverObj'*squeeze(obj.img);
                obj.Experimental.residue=sum((obj.SolverObj.experimental.residue).^2,[4 5]); %Squared norm of residual
                obj.SolverObj.experimental=[];
            else
                error('Unknown Solver : use {''phaseonly'',''IDEAL'',''AMARES'',''LorentzFit''} \n')
            end
        end

        function demoFit(obj,pxl_idx)
            %             demoFit(obj,pxl_idx)


            freq=[obj.metabolites.freq_shift_Hz];
            fid=squeeze(obj.img(1,pxl_idx(1),pxl_idx(2),pxl_idx(3),:));
            fid=padarray(fid,[size(fid,1) 0],0,'post');
            fprintf('Performing Nlorentz fit\n')

            dw=obj.DMIPara.dwell;
            faxis=linspace(-0.5/dw,0.5/dw,length(fid));
            [~,minIdx]=min(abs(faxis-(min(freq)-200)));
            [~,maxIdx]=min(abs(faxis-(max(freq)+200)));
            freqSel=minIdx:maxIdx; %1:length(faxis)

            xdata=faxis(freqSel);
            spect=(fftshift(fft(fid*exp(-1i*rad2deg(42.94)))));
            [fitf,gof,fitoptions]=NLorentzFit(xdata(:),abs(spect(freqSel)),freq);
            figure(34),clf,plot(xdata,abs(spect(freqSel)));
            hold on
            plot(xdata,real(spect(freqSel)));
            plot(xdata,fitf(xdata));


            %test AMARES fit
            fid=squeeze(obj.img(1,pxl_idx(1),pxl_idx(2),pxl_idx(3),:));
            [expParams,pk]=getAMARES_structs(mcobj);
            if(ismatrix(obj.FieldMap))
                fm_est=-1*obj.FieldMap(pxl_idx(1),pxl_idx(2),pxl_idx(3));
                fm_est=obj.mat2col(fm_est,obj.mask);
                fm_est=fm_est./(2*pi)*(6.536 /42.567); % [Hz]
                fm_est=fm_est./obj.DMIPara.ImagingFrequency_MHz;
            else
                fm_est=0;
            end
            
            expParams.offset=fm_est;
            [fitResults, fitStatus,~,CRBResults] = AMARES.amaresFit(double(fid(:)), expParams, pk,true,'quiet',false);
            disp(fitResults)

            % try to reconstruct fitted spectrum
            taxis=expParams.timeAxis;
            sig1=zeros(size(taxis));
            a=fitResults.amplitude;
            f=fitResults.chemShift*amares_struct.imagingFrequency;
            phi=0.*deg2rad(fitResults.phase);
            d=fitResults.linewidth; % [Hz]
            for i=1:length(a)
                sig1=sig1+a(i)*exp(-1i*phi(i)).*exp(-(pi*d(i))*taxis).*exp(1i*2*pi*f(i)*taxis);
            end
            spec_amares=fftshift(fft(sig1(:)));
            figure(34),
            plot(faxis,spec_amares,'LineWidth',2);
            xlim([min(freq)-200 max(freq)+200])
            legend({'Data-abs','Data-real','Nlorentzian','AMRARES'})
            xlabel('Frequency [Hz]')
            title(sprintf('Fit results of the voxel at (%d,%d,%d)',pxl_idx))
        end

        function WriteImages(obj)

            if(nargin==1)
                fPath=pwd;
                fn=sprintf('M%d_%s.nii',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription);
                niiFileName=fullfile(fPath,fn);
            elseif(isfolder(niiFileName))
                fPath=niiFileName;
                fn=sprintf('M%d_%s.nii',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription);
                niiFileName=fullfile(fPath,fn);
            else
                [fPath,~]=fileparts(niiFileName);
            end

            fPath=pwd;
            vol_PRS=squeeze(sos(flip(obj.img,3),[5 6])); % 9.4T specific
            description='averaged image across echo and phase cycle';
            MyNIFTIWrite_CSI(squeeze(single(abs(vol_PRS))),obj.twix,niiFileName,description,-1*[0 0 0]*1.5);
            if(~isempty(obj.Metcon))
                fn=sprintf('Metcon_m%d_%s_%s_%s.nii',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription,obj.flags.Solver,obj.flags.doPhaseCorr);
                vol_PRS=flip(obj.getNormalized,2); %9.4T specific
                description=sprintf('dim4_%s_%s_%s_%s_',obj.metabolites.name);
                MyNIFTIWrite_CSI(squeeze(single(abs(vol_PRS))),obj.twix,fullfile(fPath,fn),description,1*[1 2 4]*8.3);
            end

        end
        function SaveResults(obj)
            % save results to a MAT file

            OutFile=sprintf('m%d_%s_%s_%s.mat',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription,obj.flags.Solver,obj.flags.doPhaseCorr);
            mcobj_copy=copy(obj);
            mcobj_copy.twix=[];
            mcobj_copy.sig=[];
            save(OutFile,"mcobj_copy");



        end

        function PlotResults(obj,fh)
            if(~exist('fh','var'))
                fh=figure;
            end


            tt=tiledlayout(fh,2,ceil(length(obj.metabolites)/2+1),'TileSpacing','compact','Padding','compact');
            slcFac=0.3;slcdim=2;
            slcSel=round(slcFac*size(obj.Metcon,slcdim));
            slcSel=slcSel:(size(obj.Metcon,slcdim)-slcSel);

            imtransFunc=@(x) flip(flip(permute(x(:,slcSel,:),[setdiff(1:3,slcdim) slcdim 4]),10),30);
            %get metcon in normalized SNR units
            [im_snr,sf]=obj.getNormalized();

%             [im_snr,sf]=obj.getmM();
%             im_snr(:,:,:,1)=im_snr1(:,:,:,1);% replace water with SNR 
            for i=1:length(obj.metabolites)
                nexttile(tt)
                im_curr=createImMontage(imtransFunc(abs(im_snr(:,:,:,i))));
                imagesc(im_curr);
                colorbar,
                axis image

                cax_im=[0 prctile(im_curr(:),99)];
                clim(cax_im);
                xticks([]),yticks([]),title(obj.metabolites(i).name)

            end

            % display residue
            nexttile(tt)
            residue_norm=sos(imtransFunc(obj.Experimental.residue)./reshape(sf,1,1,1,[]),4);
            im_curr=createImMontage(residue_norm);
            imagesc(im_curr);
            colorbar,

            axis image

            cax_im=[0 prctile(im_curr(:),95)];
            clim(cax_im);
            xticks([]),yticks([]),title('residue')

            %display fieldmaps

            if(isfield(obj.Experimental,'fm_est'))
                fm_Hz=createImMontage(imtransFunc(obj.Experimental.fm_est));
            elseif(~isempty(obj.FieldMap))
                fm_Hz=createImMontage(imtransFunc(obj.FieldMap))./(2*pi)*(6.536 /42.567);
            end
            if(exist('fm_Hz','var'))
                nexttile()
                imagesc(fm_Hz);
                title('input 2H fieldmap [Hz]')
                cax_im=[-1*prctile(fm_Hz(:),95) prctile(fm_Hz(:),95)+1];
                clim(cax_im);
                colorbar,
                axis image,colormap(gca,'jet');
                xticks([]),yticks([]),
            end

            temp_str=sprintf('%s|%s',obj.DMIPara.ShortDescription,obj.flags.Solver);

            if(isfield(obj.DMIPara,'IntakeTime'))
                temp_str=sprintf('%s|%d min',temp_str,obj.getMinutesAfterIntake());
            end

            annotation('textbox',[0.5 0.9 0.1 0.1],'String',temp_str,'HorizontalAlignment','center','FontSize',12)

            OutfigFile=sprintf('%s.fig',strrep(temp_str,'|','_'));
            OutfigFile=strrep(OutfigFile,' ','');
            savefig(fh,OutfigFile)

        end
        function [Metcon_norm,scale_fac]= getNormalized(obj)
            % normalize metabolite maps to get SNR units hoprfully!
            if(any(strcmp(obj.flags.Solver,{'IDEAL','IDEAL-modes','phaseonly'})))
                %more analytical
                Ai=pinv(  obj.SolverObj.getA() );
                scale_fac=(sum(abs(Ai).^2,2).^(1/2))/sqrt(2);
                 scale_fac=diag(abs(Ai*Ai'))*sqrt(size(Ai,2))*sqrt(2);
                Metcon_norm=((obj.Metcon)./reshape(scale_fac,1,1,1,[]));

            else % scale by std of noise region
                Metcon_norm=zeros(size(obj.Metcon));
                noiseMask=~obj.getMask(30);
                %for noise measurements
                if(sum(noiseMask,'all')==0),noiseMask=ones(size(noiseMask),'logical');end
                calcStd=@(x) std(obj.mat2col(x,noiseMask),[],'all','omitnan');
                calcMean=@(x) mean(obj.mat2col(x,noiseMask),'all','omitnan');
                for i=1:length(obj.metabolites)
                    scale_fac(i)=calcStd(abs(obj.Metcon(:,:,:,i)));
                    Metcon_norm(:,:,:,i)=abs(obj.Metcon(:,:,:,i))./calcStd(abs(obj.Metcon(:,:,:,i))); %normalize to SNR units


                    noise_SNR= abs(calcMean((obj.Metcon(:,:,:,i))))./calcStd(abs(obj.Metcon(:,:,:,i)));
                    if(noise_SNR>0.1), warning('Noise mask probably has signal and your SNR metric is bullshit!'); end
                end
            end
        end
        function [Metcon_mM,scale_fac]= getmM(obj)
            % convert metabolite amplitudes into mM
            % assuming water concentration is 10 mM 
            %reference: peters et al DOI: 10.1002/mrm.28906 , equation 6
            
            TR=obj.DMIPara.TR;
            TE=obj.DMIPara.TE(1);
            DC=obj.DMIPara.DutyCycle;
            %Actual reference voltage is between 500 and 550
            FA_rad=obj.DMIPara.FlipAngle*(obj.DMIPara.RefVoltage/500);

            if(contains(obj.twix.hdr.Config.SequenceFileName,'fid'))
                [Msig_all]=MetSignalModel(obj.metabolites,TE,0,TR,0,FA_rad,'GRE-peters',DC);
            else
                [Msig_all]=MetSignalModel(obj.metabolites,TE,0,TR,0,FA_rad,'bSSFP',DC);
            end

            Sig_theory=mean(abs(squeeze(Msig_all)),[2 3 4]);
            Sig_theory=Sig_theory(1)./Sig_theory;
            scale_fac=10*Sig_theory./[1;2;2;2]; % mM
            metcon_w=obj.Metcon./obj.Metcon(:,:,:,1);
            Metcon_mM=metcon_w.*permute(scale_fac(:),[2 3 4 1]);
%             as(Metcon_mM,'select',':,:,25,3','windowing',[1.5 3])
        end



        function MinuteElapsed=getMinutesAfterIntake(obj,IntakeTime)
            %mcobj.getMinutesAfterIntake('hh:mm')
            if(exist("IntakeTime",'var'))
                IntakeTime=duration(IntakeTime,'InputFormat','hh:mm');
                obj.DMIPara.IntakeTime=IntakeTime;
            elseif(isfield(obj.DMIPara,'IntakeTime'))
                fprintf('Using glucose intake time : %s\n',obj.DMIPara.IntakeTime);
            else
                error('Glucose Intake in hh:mm')
            end

            MeasStartTime=duration(string(obj.DMIPara.SeriesDateTime_start,'hh:mm:ss'),'InputFormat','hh:mm:ss');
            timeElapsed=(MeasStartTime-obj.DMIPara.IntakeTime)+seconds(obj.DMIPara.acq_duration_s*0.5);
            MinuteElapsed=round(minutes(timeElapsed));
        end
    function ppm=Hz2ppm(freq)
        imagingFrequency=obj.twix.hdr.Dicom.lFrequency*1e-6; %MHz
        ppm=freq./imagingFrequency;
    end
    end

    methods(Static)
        function col_vec=mat2col(mat,mask)
            sz=[size(mat),1,1,1];
            mat=reshape(mat,[prod(sz(1:3)) sz(4:end)]);
            col_vec=mat(mask(:),:,:,:);

        end
        function mat=col2mat(col_vec,mask)
            sz=size(col_vec);
            col_vec=reshape(col_vec,sz(1),[]);
            mat=zeros([numel(mask) prod(sz(2:end))]);
            for i=1:size(mat,2)
                mat(mask,i)=col_vec(:,i);
            end
            mat=reshape(mat,[size(mask) sz(2:end)]);

        end

    end
end


