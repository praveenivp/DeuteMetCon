classdef MetCon_CSI<matlab.mixin.Copyable
    properties
        time  %[s]

        B0Map %spatial [rad/s]
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
                    addParameter(p,'Solver','AMARES',@(x) any(strcmp(x,{'pinv','IDEAL','AMARES','LorentzFit'})));
                    %                     addParameter(p,'maxit',10,@(x)isscalar(x));
                    %                     addParameter(p,'tol',1e-6,@(x)isscalar(x));
                    %                     addParameter(p,'reg','none',@(x) any(strcmp(x,{'none','Tikhonov'})));
                    %                     addParameter(p,'reg_lambda',0,@(x)isscalar(x));
                    
                    addParameter(p,'ZeroPadSize',[], @(x) isvector(x)); % [3 physical x 1 time]
                    parse(p,varargin{:});

                    obj.flags=p.Results;
                    if(isfield(p.Unmatched,'csm')) %[CHAxCOLxLINxPARxSLC]
                        obj.coilSens=p.Unmatched.csm;
                        obj.flags.doCoilCombine='none';
                        %                         obj.SpiralPara.R_PE=2;
                        %                         obj.flags.doPAT='CGSENSE';
                    end
                    if(isfield(p.Unmatched,'fm'))
                        obj.B0Map=p.Unmatched.fm;
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
                        Acqdelay=mcobj.twix.hdr.Phoenix.sSpecPara.lAcquisitionDelay*1e-6; %s
                        obj.flags.phaseoffet=[0 2*pi*Acqdelay];
                    end


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
            if isfield(twix,'noise')
                noise                = permute(twix.noise(:,obj.flags.CoilSel,:),[2,1,3]);
                noise                = noise(:,:)';
                R                    = cov(noise);
                R(eye(size(R,1))==1) = abs(diag(R));
                if(obj.flags.NormNoiseData)
                    R= R./mean(abs(diag(R)));
                    obj.D               = sqrtm(inv(R)).';
                else
                    scale_factor=1; %dwell time are the same
                    Rinv = inv(chol(R,'lower')).';
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
            obj.img=myfft(obj.sig,[2 3 4]);
            obj.performCoilCombination();
            obj.performPhaseCorr();
            obj.performSVDdenosing();
            print_str = sprintf( 'reco  time = %6.1f s\n', toc);
            fprintf(print_str);

            
        end

        function getSig(obj)
            obj.sig = obj.twix.image(obj.flags.EchoSel,obj.flags.CoilSel,:,:,:,:,:,:,:,:,:,:,:,:,:);
            % Sum averages (if not already done)
            obj.sig  = sum( obj.sig ,6);

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
            % zeropad data
            if(isempty(obj.flags.ZeroPadSize))
                pad_size=[0 0 0 0];
                sSpecPara=obj.twix.hdr.MeasYaps.sSpecPara;
                % Add missing voxels
                MissingVOXRead  = double(sSpecPara.lFinalMatrixSizeRead-size(obj.sig,2));
                MissingVOXPhase = double(sSpecPara.lFinalMatrixSizePhase-size(obj.sig,3));
                MissingVOXSlice = double(sSpecPara.lFinalMatrixSizeSlice-size(obj.sig,4));
                disp(['initial CSI data size:           ', num2str(size(obj.sig))])
                obj.sig= padarray(obj.sig,floor([0,MissingVOXRead,MissingVOXPhase,MissingVOXSlice,0,0]./2),'pre');
                obj.sig = padarray(obj.sig,ceil([0,MissingVOXRead,MissingVOXPhase,MissingVOXSlice,0,0]./2),'post');
                disp(['final CSI data size:           ', num2str(size(obj.sig))])
            else
                zp_PRS=obj.flags.ZeroPadSize(1:3);
                pad_size=[0 round(size(obj.sig,2)*zp_PRS(1)) round(size(obj.sig,3)*zp_PRS(2))  round(size(obj.sig,4)*zp_PRS(3))];
            end
            obj.sig=padarray(obj.sig,pad_size,0,'both'); % spatial zeropad
            obj.sig=padarray(obj.sig,[zeros(1,4) round(size(obj.sig,5)*obj.flags.ZeroPadSize(4))],0,'post'); % spectral zerpoad
            
        
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
                n=size(size(obj.img,5));
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
                    for vx=1:size(CSI_Data,1)
                        rawSpectra=squeeze(CSI_Data(vx,:,:,1));  % Spectrum x Coil
                        [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(rawSpectra, [], CSI_wsvdOption);
                        CoilWeights(vx,:) = wsvdWeights;
                        CSI_Combined(vx,:)=wsvdCombination;
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


        function performMetaboliteMapping(obj)
            obj.SolverObj=LeastSquares(obj.metabolites,obj.DMIPara.TE,'fm',obj.B0Map);
            obj.Metcon=obj.SolverObj'*squeeze(obj.img);


        end

        function peformB0mapReslice(obj,regfiles)
            % we need spm12 in path
            if(exist('spm','file'))

                % write average volume
                bSSFP_file='S:\Deuterium\20230920_xpulseq\processeddata\meas_MID00079_FID62457_pvrh_trufi_mc_2H_15mm_150PC.nii';

                % B0 map
                fm_file='S:\Deuterium\20230920_xpulseq\M80_B0map_rad_s.nii';
                im_file='S:\Deuterium\20230920_xpulseq\gre_WIP1441B_shim_2mm_R2_e1_mag.nii';

                % write average volume
                bSSFP_file=regfiles{1};

                % B0 map
                fm_file=regfiles{2};
                im_file=regfiles{3};


                %spm reslice
                % Define the input NIfTI files and output directory
                [input_dir,~,~] = fileparts(fm_file);
                [output_dir,~,~] = fileparts(bSSFP_file);

                % Define the reslicing parameters (e.g., target space)
                matlabbatch{1}.spm.spatial.coreg.write.ref = {bSSFP_file};
                matlabbatch{1}.spm.spatial.coreg.write.source = {fm_file;im_file};
                matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
                matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r_';


                % Run the reslicing job
                spm_jobman('run', matlabbatch);

%                 % Move the resliced files to the output directory
%             
%                     [~, niifile, ext] = fileparts(matlabbatch{1}.spm.spatial.coreg.write.source{i});
%                     movefile(fullfile(input_dir,sprintf('r_*')), output_dir);


            else
                error('Needs spm for reslicing/registration')
            end
        end



        function getMask(obj,thres_prctl)
            if(nargin==1)
                thres_prctl=95;
            end
            im_abs=abs(squeeze(sum(obj.img,5)));
            obj.mask=im_abs>0.1*prctile(im_abs(:),thres_prctl);  
 
            obj.mask=imerode(obj.mask,strel('sphere',2));
            obj.mask=imdilate(obj.mask,strel('sphere',3));

        end
        function performMetcon(obj,freq)

            fids=MetCon_CSI.mat2col(permute(obj.img,[2 3 4 5 1]),obj.mask);
            nMet=3;
            if(strcmpi(obj.flags.Solver,'AMARES'))







                wbhandle = waitbar(0,'Performing AMARES fitting');


                            DW= obj.DMIPara.dwell;
            BW=1/DW; %Hz
            samples=size(obj.img,5);
            timeAxis=0:DW:(samples-1)*DW;
            imagingFrequency=obj.twix.hdr.Dicom.lFrequency*1e-6;
            offset=0;
            ppmAxis=linspace(-BW/2,BW/2-BW/samples,samples)/imagingFrequency;
            nMet=3;
            chemShift=freq./imagingFrequency;
            phase=ones(nMet,1).*0;
            amplitude=[1 1 1];
            linewidth=[12 12  12];

            amares_struct=struct('chemShift',chemShift, ...
                'phase', phase,...
                'amplitude', amplitude,...
                'linewidth',linewidth, ...
                'imagingFrequency', imagingFrequency,...
                'BW', BW,...
                'timeAxis', timeAxis(:), ...
                'dwellTime', DW,...
                'ppmAxis',ppmAxis(:), ...
                'beginTime',0, ...
                'offset',offset,...
                'samples',samples, ...
                'peakName',["water","Glc","Glx","Lac"]);

            pk=PriorKnowledge_DMI(amares_struct);


                metabol_con=zeros(size(fids,1),nMet*3+3);
                %             output format (4th dim) : use xFit order [chemicalshift x N] [linewidth xN] [amplitude xN] [phase x1] fitStatus.relativeNorm fitStatus.resNormSq]

                parfor i=1:size(fids,1)
                    %                 if(mod(i,1000)==0), fprintf('%.0f %d done\n',i/size(fids,1)*100); drawnow(); end
%                     waitbar(i/size(fids,1),wbhandle,'Performing AMARES fitting');
                    [fitResults, fitStatus, figureHandle, CRBResults] = AMARES.amaresFit(double(fids(i,:).'), amares_struct, pk, 0,'quiet',true);
                    metabol_con(i,:)= [fitStatus.xFit fitStatus.relativeNorm fitStatus.resNormSq]';

                end
                close(wbhandle);

                % convert to 3D matrix and give resonable names
                metabol_con=MetCon_CSI.col2mat(metabol_con,obj.mask);

                obj.Experimental.chemicalshift=metabol_con(:,:,:,((1:nMet)+nMet*(1-1))); %chemical shift ppm
                obj.Metcon=metabol_con(:,:,:,((1:nMet)+nMet*(3-1)));  %amplitudes
                obj.Experimental.linewidth=metabol_con(:,:,:,((1:nMet)+nMet*(2-1)));
                obj.Experimental.phase=metabol_con(:,:,:,((1)+nMet*(4-1))); %phase
                obj.Experimental.relativeNorm=metabol_con(:,:,:,((2)+nMet*(4-1)));
                obj.Experimental.resNormSq=metabol_con(:,:,:,((3)+nMet*(4-1)));

            elseif(strcmpi(obj.flags.Solver,'LorentzFit'))
                dw=obj.DMIPara.dwell;
                faxis=linspace(-0.5/dw,0.5/dw,size(obj.img,5));
%                 freq=[-148.7 -54.79 7.82]'; %Hz
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
                rSQ=zeros(size(spectrum_All,1),1);
                sRMSE=zeros(size(spectrum_All,1),1);
                %             output format (4th dim) : use xFit order [chemicalshift x N] [linewidth xN] [amplitude xN] [phase x1] fitStatus.relativeNorm fitStatus.resNormSq]
                
                for i=1:size(spectrum_All,1)
                    %                 if(mod(i,1000)==0), fprintf('%.0f %d done\n',i/size(fids,1)*100); drawnow(); end
                    waitbar(i/size(fids,1),wbhandle,'Performing nLorentzian fitting');

                    [fitf,gof,fitoptions]=NLorentzFit(xdata(:),spectrum_All(i,:).',freq);
                    metabol_con(i,:)= [fitf.A1; fitf.A2; fitf.A3;];

                    rSQ(i)=gof.adjrsquare;
                    sRMSE(i)=gof.rmse;

                end
                close(wbhandle);

                % convert to 3D matrix and give resonable names
                obj.Metcon=MetCon_CSI.col2mat(metabol_con,obj.mask);
                obj.Experimental.rSQ=MetCon_CSI.col2mat(rSQ,obj.mask);
                obj.Experimental.sRMSE=MetCon_CSI.col2mat(sRMSE,obj.mask);


            end
        end

        function demoFit(obj,pxl_idx,freq)
            

%              freq=[-142.7 -51.79 7.82]'; %Hz

            fid=squeeze(obj.img(1,pxl_idx(1),pxl_idx(2),pxl_idx(3),:));
            fid=padarray(fid,[size(fid,1) 0],0,'post');
            fprintf('Performing Nlorentz fit\n')
           
            dw=obj.DMIPara.dwell;
            faxis=linspace(-0.5/dw,0.5/dw,length(fid));
            [~,minIdx]=min(abs(faxis-(min(freq)-200)));
            [~,maxIdx]=min(abs(faxis-(max(freq)+200)));
            freqSel=minIdx:maxIdx; %1:length(faxis)

            xdata=faxis(freqSel);
            spect=1*abs(fftshift(fft(fid*exp(-1i*rad2deg(42.94)))));
            [fitf,gof,fitoptions]=NLorentzFit(xdata(:),spect(freqSel),freq);
            figure(34),clf,plot(xdata,spect(freqSel));
            hold on
            plot(xdata,fitf(xdata));

            fprintf('Performing AMARES fit\n')


            DW= obj.DMIPara.dwell;
            BW=1/DW; %Hz
            timeAxis=0:DW:(length(fid)-1)*DW;
            imagingFrequency=obj.twix.hdr.Dicom.lFrequency*1e-6;
            offset=0;
            ppmAxis=linspace(-BW/2,BW/2-BW/length(fid),length(fid))/imagingFrequency;
            nMet=3;
            chemShift=freq./imagingFrequency;
            phase=ones(nMet,1).*0;
            amplitude=[1 1 1];
            linewidth=[12 12  12];

            samples=length(fid);

            amares_struct=struct('chemShift',chemShift, ...
                'phase', phase,...
                'amplitude', amplitude,...
                'linewidth',linewidth, ...
                'imagingFrequency', imagingFrequency,...
                'BW', BW,...
                'timeAxis', timeAxis(:), ...
                'dwellTime', DW,...
                'ppmAxis',ppmAxis(:), ...
                'beginTime',0, ...
                'offset',offset,...
                'samples',samples);

            pk=PriorKnowledge_DMI(amares_struct);




            [fitResults, fitStatus, figureHandle, CRBResults] = AMARES.amaresFit(double(fid), amares_struct, pk, 0,'quiet',true);

            taxis=0:dw:dw*(size(fid,1)-1);
            sig1=zeros(size(taxis));
            a=fitResults.amplitude;
            f=fitResults.chemShift*amares_struct.imagingFrequency;
            phi=0.*deg2rad(fitResults.phase);
            d=fitResults.linewidth; % [Hz]

            for i=1:length(a)
                sig1=sig1+a(i)*exp(-1i*phi(i)).*exp(-(pi*d(i))*taxis).*exp(1i*2*pi*f(i)*taxis);
            end
            spec_amares=fftshift(fft(sig1(:)));
            plot(faxis,spec_amares,'LineWidth',2);
            xlim([min(freq)-200 max(freq)+200])
            legend({'Data','Nlorentzian','AMRARES'})
        end

        function WriteImages(obj)
            [fPath,fn,~]=fileparts(obj.filename);
%             try
%                 if(~isfolder(fullfile(fPath,'processeddata')))
%                     mkdir(fullfile(fPath,'processeddata'))
%                 end
%             catch
                fPath=pwd;
                if(~isfolder(fullfile(fPath,'processeddata')))
                    mkdir(fullfile(fPath,'processeddata'))
                end
%             end
            vol_PRS=squeeze(sos(flip(obj.img,3),[5 6])); % 9.4T specific
            description='averaged image across echo and phase cycle';
            MyNIFTIWrite_CSI(squeeze(single(abs(vol_PRS))),obj.twix,fullfile(fPath,'processeddata',fn),description,-1*[9.8+4;0;9.8]*1.5);


        end
        function SaveResults(obj)
            % save results toa MAT file
            im=squeeze(obj.img);

            flags=obj.flags;

            OutFile=sprintf('m%d_B0%s_PhaseC%s.mat',obj.twix.hdr.Config.MeasUID,obj.flags.doMet,obj.flags.doPhaseCorr);
            sp=obj.DMIPara;
            fn=obj.filename;
%             ro=(2*sp.ADCLength*sp.DwellTime)/1e6; % ms
%             vTR=(sp.TR*sp.Ninterleaves*sp.NPartitions)/(sp.R_PE*sp.R_3D*1e6); %s
% %             descrip=(sprintf('R%dx%dC%d TR=%.1fms RO=%.2fms vTR=%.1fs',sp.R_PE,sp.R_3D,sp.CAIPIShift,sp.TR/1e3,ro,vTR));
%             descrip_reco=sprintf('%s PAT=%s coilcomb=%s B0=%s DCF=%s CompMode=%s',flags.CompMode,flags.doPAT, flags.doCoilCombine, flags.doB0Corr,flags.doDCF,flags.CompMode);
            save(OutFile,'im','sp','flags','descrip','descrip_reco','fn','-v7.3')


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


