classdef MetCon_CSI<matlab.mixin.Copyable
    properties
        time  %[s]

        B0Map %spatial [rad/s]

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
                obj.twix=obj.twix{2};
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

                    addParameter(p,'doCoilCombine','adapt1',@(x) any(strcmp(x,{'none','sos','adapt1','adapt2'})));

                    addParameter(p,'doNoiseDecorr',true,@(x) islogical(x));
                    addParameter(p,'NormNoiseData',false,@(x) islogical(x));

                    addParameter(p,'CoilSel',1:obj.twix.image.NCha, @(x) isvector(x));
                    addParameter(p,'PCSel',1:obj.twix.image.NRep, @(x) (isvector(x)&& all(x<=1:obj.twix.image.NRep)) );
                    addParameter(p,'EchoSel',1:obj.twix.hdr.Phoenix.sSpecPara.lVectorSize, @(x) (isvector(x)&& all(x<=obj.twix.hdr.Phoenix.lContrasts)) );
                    addParameter(p,'is3D',(obj.twix.image.NPar>1),@(x)islogical(x));
                    %                     addParameter(p,'doB0Corr','none',@(x) any(strcmp(x,{'none','MTI','MFI'})));
                    %                   addParameter(p,'precision','single',@(x) any(strcmp(x,{'single','double'})));
                    addParameter(p,'Solver','pinv',@(x) any(strcmp(x,{'pinv','IDEAL'})));
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
                        %                         obj.flags.doB0Corr='MTI';
                    end

                    if(isfield(p.Unmatched,'metabolites'))
                        obj.metabolites=p.Unmatched.metabolites;
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
                noise                = noise(:,:).';
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
                zp_PRS=[1 1 1];
                pad_size=[0 round(size(obj.sig,2)*zp_PRS(1)) round(size(obj.sig,3)*zp_PRS(2))  round(size(obj.sig,4)*zp_PRS(3))];
            end
            obj.sig=padarray(obj.sig,pad_size,0,'both');
        
        end
        function performPhaseCorr(obj)
            % remove phase of first echo!
            obj.img=bsxfun(@times,obj.img,exp(-1i*angle(obj.img(:,:,:,:,1,1))));
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

                    [~,obj.coilSens,obj.coilNormMat]=adaptiveCombine(sum(obj.img(:,:,:,:,1,:),6));
                    obj.img=(sum(obj.coilSens.*obj.img,1));
                case 'adapt2'
                    if(obj.LoopCounter.cRep==1 && ~any(obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc),'all') )
                        [obj.img(1,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep), ...
                            obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc),obj.coilNormMat]...
                            = adaptiveCombine2(obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep),[],false);
                        obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc)=conj(obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc));
                    else
                        obj.img(1,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep)=sum(obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep).*...
                            conj(obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc)),1);
                    end
            end

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

            OutFile=sprintf('m%d_B0%s_DCF%s.mat',obj.twix.hdr.Config.MeasUID,obj.flags.doB0Corr,obj.flags.doDCF);
            sp=obj.SpiralPara;
            fn=obj.filename;
            ro=(2*sp.ADCLength*sp.DwellTime)/1e6; % ms
            vTR=(sp.TR*sp.Ninterleaves*sp.NPartitions)/(sp.R_PE*sp.R_3D*1e6); %s
            descrip=(sprintf('R%dx%dC%d TR=%.1fms RO=%.2fms vTR=%.1fs',sp.R_PE,sp.R_3D,sp.CAIPIShift,sp.TR/1e3,ro,vTR));
            descrip_reco=sprintf('%s PAT=%s coilcomb=%s B0=%s DCF=%s CompMode=%s',flags.CompMode,flags.doPAT, flags.doCoilCombine, flags.doB0Corr,flags.doDCF,flags.CompMode);
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


