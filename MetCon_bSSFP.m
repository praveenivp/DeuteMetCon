classdef MetCon_bSSFP<matlab.mixin.Copyable
    properties
        time  %[s]

        FieldMap %spatial [rad/s]
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
        function obj=MetCon_bSSFP(varargin)
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

                    addParameter(p,'doCoilCombine','adapt1',@(x) any(strcmp(x,{'none','sos','adapt1','adapt2'})));
                    addParameter(p,'doKspaceFilter','none',@(x) any(strcmp(x,{'none','hann','hamm'})));

                    addParameter(p,'doNoiseDecorr',true,@(x) islogical(x));
                    addParameter(p,'doPhaseCorr',true,@(x) islogical(x));
                    addParameter(p,'NormNoiseData',false,@(x) islogical(x));
                    addParameter(p,'doAverage',true,@(x) islogical(x));


                    addParameter(p,'doZeroPad',[1,1,1],@(x) isvector(x));
                    %-1 for debug
                    addParameter(p,'doDenosing',0,@(x) isscalar(x));

                    addParameter(p,'CoilSel',1:obj.twix.image.NCha, @(x) isvector(x));
                    addParameter(p,'PCSel',1:obj.twix.image.NRep, @(x) (isvector(x)&& all(x<=1:obj.twix.image.NRep)) );
                    addParameter(p,'EchoSel',1:obj.twix.hdr.Phoenix.lContrasts, @(x) (isvector(x)&& all(x<=obj.twix.hdr.Phoenix.lContrasts)) );
                    addParameter(p,'is3D',(obj.twix.image.NPar>1),@(x)islogical(x));
                    %                     addParameter(p,'doB0Corr','none',@(x) any(strcmp(x,{'none','MTI','MFI'})));
                    %                   addParameter(p,'precision','single',@(x) any(strcmp(x,{'single','double'})));
                    addParameter(p,'Solver','pinv',@(x) any(strcmp(x,{'pinv','IDEAL'})));
                    %                     addParameter(p,'maxit',10,@(x)isscalar(x));
                    %                     addParameter(p,'tol',1e-6,@(x)isscalar(x));
                    %                     addParameter(p,'reg','none',@(x) any(strcmp(x,{'none','Tikhonov'})));
                    %                     addParameter(p,'reg_lambda',0,@(x)isscalar(x));


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
                        %                         obj.flags.doB0Corr='MTI';
                    end
                    if(isfield(p.Unmatched,'mask'))
                        obj.mask=p.Unmatched.mask;
                        %                         obj.flags.doB0Corr='MTI';
                    end
                    if(isfield(p.Unmatched,'metabolites'))
                        obj.metabolites=p.Unmatched.metabolites;
                        %                         obj.flags.doB0Corr='MTI';
                    end
                    if(isfield(p.Unmatched,'NoiseDecorr'))
                        obj.D=p.Unmatched.NoiseDecorr;
                        %                         obj.flags.doB0Corr='MTI';
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
            if (isfield(twix,'noise') &&  isempty(obj.D))
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
            obj.performKspaceFiltering();
            obj.img=myfft(obj.sig,[2 3 4]);
            obj.performCoilCombination();
            obj.performSVDdenosing();
            obj.performPhaseCorr();
            
            print_str = sprintf( 'reco  time = %6.1f s\n', toc);
            fprintf(print_str);

            
        end

        function getSig(obj)
            
        if(~obj.DMIPara.isCSI)
                                    %             {'Col','Cha','Lin','Par','Sli','Ave','Phs','Eco','Rep',
            %     'Set','Seg','Ida','Idb','Idc','Idd','Ide'}
            obj.sig = obj.twix.image(:, obj.flags.CoilSel,:,:, 1,:,1,obj.flags.EchoSel,obj.flags.PCSel,1,1);
            % Sum averages (if not already done)
            if(obj.flags.doAverage)
            obj.sig  = mean( obj.sig ,6);
            end
             %Coil x Phase x Read xSlice x echo x PC
            obj.sig=permute(obj.sig,[2 3 1 4 8 9 6 5 7]);
        else

            error('not implemented');

        end

        
        end
        function performKspaceFiltering(obj)
       

            if strcmp(obj.flags.doKspaceFilter,'hann')
            alpha = 0.5;
            elseif strcmp(obj.flags.doKspaceFilter,'hamm')
            alpha = 0.54;
            else 
                return;
            end


            sz=size(obj.sig);%iso tropic voxels

            radius = 0.;
            for ax=2:4
                r = (-floor(sz(ax)/2):ceil(sz(ax)/2)-1)/floor(sz(ax)/2);
                radius = bsxfun(@plus,radius,shiftdim(r(:).^2,1-ax));
            end
            radius = sqrt(radius); % *0.75; %reduce filtering

 
            func = @(x) alpha - (1-alpha) * cos(pi*(x+1.));
            filter=func(radius);
            obj.sig=obj.sig.*(filter);
        end
        function performZeroPaddding(obj)
            
            if(all(obj.flags.doZeroPad==0))
                return;
            end
            % zeropad data
            zp_PRS=obj.flags.doZeroPad;
            pad_size=[0 round(size(obj.sig,2)*zp_PRS(1)) round(size(obj.sig,3)*zp_PRS(2))  round(size(obj.sig,4)*zp_PRS(3))];
            obj.sig=padarray(obj.sig,pad_size,0,'both');
        
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
        function performPhaseCorr(obj)
            if(obj.flags.doPhaseCorr)
            % remove phase of first echo!
            obj.img=bsxfun(@times,obj.img,exp(-1i*angle(obj.img(:,:,:,:,1,1))));
            end
        end

        

        function performCoilCombination(obj)
         

            if(~(size(obj.img,1)>1))
                obj.flags.doCoilCombine='none';
            end
            %coil combination
            switch(obj.flags.doCoilCombine) %{'none','sos','adapt','adapt2'}
                case 'sos'
                    obj.img= sqrt(sum(abs(obj.img).^2,1));
                case 'adapt1' %preserves phase better?

                    [~,obj.coilSens,obj.coilNormMat]=adaptiveCombine(sum(obj.img(:,:,:,:,1,:),6));
                    obj.img=(sum(obj.coilSens.*obj.img,1));
                     obj.img=obj.img.*permute(obj.coilNormMat,[5 1 2 3 4]); %norm
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

        function performMetCon(obj)
            PhSel=obj.flags.PCSel;
            EchoSel=obj.flags.EchoSel;
            FA=obj.DMIPara.FlipAngle; %rad
            TE=obj.DMIPara.TE; %s
            TR=obj.DMIPara.TR; %s
            PC=obj.DMIPara.PhaseCycles; %rad
           
            if(strcmpi(obj.flags.Solver,'pinv'))

            tic;
            B0=obj.FieldMap(obj.mask)./(2*pi)*(6.536 /42.567);%+1e6*(Spectroscopy_para.PVM_FrqWork(1)- ExpPara.PVM_FrqWork(1)); %Hz
            im1=reshape(obj.img(:,:,:,:,obj.flags.EchoSel,obj.flags.PCSel), ...
                [],length(obj.flags.EchoSel)*length(obj.flags.PCSel));
            im1=im1(obj.mask(:),:);
            %         [Msig_all]=bSSFP_sim_analytical(metabolites,TE(EchoSel),PC(PhSel),TR,zeros(size(B0)),FA);
            [Msig_all]=bSSFP_sim_analytical(obj.metabolites,TE(EchoSel),PC(PhSel),TR,B0,FA);
            Msig_all=bsxfun(@times,Msig_all,exp(-1i*angle(Msig_all(:,:,1,:,:,:,:,:,:,:))));
            metabol_con=zeros(size(im1,1),length(obj.metabolites));
            resi=zeros(size(im1));
            condnm=zeros(size(im1,1),1);
            
            for i=1:size(im1,1)
                if(mod(i,1000)==0), fprintf('%.0f %d done\n',i/size(im1,1)*100); drawnow(); end
                
                A= padarray(reshape(squeeze(Msig_all(1,:,:,:,1,i)),length(obj.metabolites),length(PhSel)*length(EchoSel)),[0 0],1,'pre').';
                b=[ im1(i,:)];
                %     b=b./max(abs(b(:)));
                metabol_con(i,:)=A\b(:);
                resi(i,:)=b(:) -A*metabol_con(i,:).';
                condnm(i)=cond(A);
            end
            
            obj.Metcon=obj.col2mat(metabol_con,obj.mask);
            obj.Experimental.residue=obj.col2mat(resi,obj.mask);
            obj.Experimental.condition=obj.col2mat(condnm(:),obj.mask);

            fprintf('\n Metabolite fitting done in %0.1f s ! \n',toc)
            elseif(strcmpi(obj.flags.Solver,'IDEAL'))

                fm_meas_Hz=obj.FieldMap./(2*pi)*(6.536 /42.567); % 2H field map in Hz
                im_me=squeeze(sum(obj.img,6)).*obj.mask; % sum phase cycles to get FISP contrast
                obj.SolverObj=IDEAL(obj.metabolites,TE,'fm',0.*fm_meas_Hz,'solver','IDEAL','maxit',5,'mask',obj.mask);
                obj.Metcon=obj.SolverObj'*im_me;

            end

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

        function WriteImages(obj,niiFileName)

            if(nargin==1)

            [fPath,fn,~]=fileparts(obj.filename);
            try
                if(~isfolder(fullfile(fPath,'processeddata')))
                    mkdir(fullfile(fPath,'processeddata'))
                end
            catch
                fPath=pwd;
                if(~isfolder(fullfile(fPath,'processeddata')))
                    mkdir(fullfile(fPath,'processeddata'))
                end
            end
            niiFileName=fullfile(fPath,'processeddata',fn);
            end

            vol_PRS=squeeze(sos(flip(obj.img,3),[5 6])); % 9.4T specific
            description='averaged image across echo and phase cycle';
            MyNIFTIWrite_bSSFP2(squeeze(single(abs(vol_PRS))),obj.twix,niiFileName,description);


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
        function OutMat=col2mat(col_vec,inMask)
            sz=size(col_vec);
            col_vec=reshape(col_vec,sz(1),[]);
            OutMat=zeros([numel(inMask) prod(sz(2:end))]);
            for i=1:size(OutMat,2)
            OutMat(inMask,i)=col_vec(:,i);
            end
            OutMat=reshape(OutMat,[size(inMask) sz(2:end)]);

        end
    end
end


