classdef MetCon_ME<matlab.mixin.Copyable

% MetCon_ME
% A class to process bSSFP multi-echo spectral-spatial data. The constructor accepts
% raw/twix filename and other optional arguments as name value pairs as
% listed below: 
%
% | name                          | description                                                | default   | possible options                                                 |
% |-------------------------------|------------------------------------------------------------|-----------|------------------------------------------------------------------|
% | 'metabolites'                 | struct array with definition of metabolites                | []        | see getMetaboliteStruct.m function                               |
% | 'NoiseDecorr'                 | Noise Decorrelation matrix [Coilxcoil]                                      | []        | 2D numeric matrix                                     |
% | 'fm'                          | 1H fieldmap in rad/s                                       | []        | 3D numeric matrix or 'IDEAL'                                     |
% | 'csm'                         | coil maps                                                  | []        | 3D numeric matrix                                                |
% | 'mask'                        | mask for spectral separation                               | []        | 3D logical matrix                                                |
% | 'doDenosing'                  | SVD denoising                                              | 0         | scalar No of components, -1 for debug                            |
% | 'Solver'                      | spectral separation method                                 | 'IDEAL'   |{'pinv','IDEAL','IDEAL-modes','phaseonly'}          |
% |                               |                                                            |           | 'phaseonly'- linear method with only phase evolution             |
% |                               |                                                            |           | 'pinv'- linear method with full signal model                     |
% |                               |                                                            |           | 'IDEAL'-  IDEAL algorithm with data averaged along phase cycles  |
% |                               |                                                            |           | 'IDEAL-modes'-IDEAL algorithm for phase cycled data              |
% | 'parfor'                      | flag to use parfor                                         | true      | boolean                                                          |
% | 'doZeroPad'                   | zero pad factor                                            | [1 1 1]   | positive scalar array [3 physical axis]                          |
% | 'doSmoothFM','maxit'          | IDEAL flags: fieldmap smooth factor and maximum iterations | 1,10      | scalar(+ve: gaussian, -ve: median),postive scalar                |
% | 'doPhaseCorr'                 | phase correction(subtract phase of all image with 1st echo)| false     | boolean                                                          |
% | 'CoilSel','PCSel','EchoSel'   | arrays to picks some of coils, echoes and phasecyles. | 1:max()   | positive integer array                                                |
% | 'doNoiseDecorr'               | flag to perform noise decorrelation                        | true      | boolean                                                          |
% | 'doKspaceFilter'              | apply restrospective skpace filter                         | 'none'    | {'none','hann','hamm'})                                          |
% | 'doCoilCombine'               | coil combine mode                                          | 'adapt1'  | {'none','sos','adapt1','adapt2'}                                 |
%
% Example:
% mcobj_me=MetCon_ME(ME_filename,'metabolites',getMetaboliteStruct('invivo'),'doZeropad',[0.5 0.5 0.5 0],'Solver','IDEAL');
%
% author: praveen

    properties

        FieldMap %spatial 1H B0 map [rad/s]
        mask % binary mask for metabolite estimation

        twix %mapVBVD object
        DMIPara % parsed sequence parameters
        flags

        filename % of raw data file

        metabolites % metabolite struct see getMetaboliteStruct.m
        sig % raw data
        img % reconstructed image [CHAxCOLxLINxPARxSLCxREP]
        coilSens %%[CHAxCOLxLINxPARxSLC]
        coilNormMat %%[COLxLINxPARxSLC]

        SolverObj % IDEAL solver object
        Metcon % metbolite amplitudes

        Experimental % experimental outputs like residue, other fit parameters , fieldmap ,etc
        D % noise decorrelation matrix

    end
    methods

        %constructor: get full path of dat file or no arguments
        function obj=MetCon_ME(varargin)
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
                    addParameter(p,'parfor',true,@(x) islogical(x));


                    addParameter(p,'doZeroPad',[1,1,1],@(x) isvector(x));
                    %-1 for debug
                    addParameter(p,'doDenoising',0,@(x) isscalar(x));
                    % posistive value for imgaulsfilt3 or negative value
                    % for medfilt3
                    addParameter(p,'doSmoothFM',1,@(x) isscalar(x));
                    addParameter(p,'maxit',10,@(x) isscalar(x));
                    %                     addParameter(p,'doParFor',inf,@(x) isscalar(x));

                    addParameter(p,'CoilSel',1:obj.twix.image.NCha, @(x) isvector(x));
                    addParameter(p,'PCSel',1:obj.twix.image.NRep, @(x) (isvector(x)&& all(x<=1:obj.twix.image.NRep)) );
                    addParameter(p,'EchoSel',1:obj.twix.hdr.Phoenix.lContrasts, @(x) (isvector(x)&& all(x<=obj.twix.hdr.Phoenix.lContrasts)) );
                    addParameter(p,'is3D',(obj.twix.image.NPar>1),@(x)islogical(x));
                    %                   addParameter(p,'precision','single',@(x) any(strcmp(x,{'single','double'})));
                    addParameter(p,'Solver','pinv',@(x) any(strcmp(x,{'pinv','IDEAL','IDEAL-modes','IDEAL-modes2','phaseonly'})));
                   addParameter(p, 'PxlShiftPerformed',true,@(x)islogical(x));
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
                        obj.mask=p.Unmatched.mask>0;
                        %                         obj.flags.doB0Corr='MTI';
                    end
                    if(isfield(p.Unmatched,'metabolites'))
                        obj.metabolites=p.Unmatched.metabolites;
                        %                         obj.flags.doB0Corr='MTI';
                    else
                        obj.metabolites=[];
                    end
                    if(isfield(p.Unmatched,'NoiseDecorr'))
                        obj.D=p.Unmatched.NoiseDecorr;
                        %                         obj.flags.doB0Corr='MTI';
                    end

                    if (obj.twix.image.NRep<4 && strcmp(obj.flags.Solver,'IDEAL-modes'))
                        obj.flags.Solver='IDEAL';
                        warning('Not enough phase cycles for IDEAL-modes: Using IDEAL mode')
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
            elseif(isempty(obj.D))
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
            obj.img=myfft(obj.sig,[2 3 4]).*sqrt(numel(obj.sig));
            obj.performCoilCombination();
            obj.performSVDdenoising();
            obj.performPhaseCorr();

            obj.getMask(0);


            obj.getFieldmap();
            fprintf(sprintf( 'reco  time = %6.1f s\n', toc));
            if(~isempty(obj.metabolites))
                obj.performMetCon();
                obj.performPxlShiftCorrection();
            end
            fprintf(sprintf( 'Metabolite mapping time = %6.1f s\n', toc));

        end

        function getSig(obj)
            %             {'Col','Cha','Lin','Par','Sli','Ave','Phs','Eco','Rep',
            %     'Set','Seg','Ida','Idb','Idc','Idd','Ide'}
            obj.sig = obj.twix.image(:, obj.flags.CoilSel,:,:, 1,:,1,obj.flags.EchoSel,obj.flags.PCSel,1,1);
            % SNR scaling for unit SD : part 1

            obj.sig=obj.sig./(sqrt(sum(abs(obj.sig(:))>0)));
            % Sum averages (if not already done)
            if(obj.flags.doAverage)
                obj.sig  = sum(obj.sig,6);
            end


            %Coil x Phase x Read xSlice x echo x PC
            obj.sig=permute(obj.sig,[2 3 1 4 8 9 6 5 7]);
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
        function performSVDdenoising(obj)
            if(all(obj.flags.doDenoising==0))
                return;
            end
            Ncomp=obj.flags.doDenoising;
            imsz=size(obj.img);
            mat=reshape(squeeze(obj.img),[],prod(imsz(5:end)));
            mat_mean=mean(mat,'all','omitnan');
            [U,S,V]=svd(mat-mat_mean,'econ');
            S=diag(S);
            if(Ncomp<0)
                figure,plot(S);
                Ncomp=input('give Number of component to keep');
                obj.flags.doDenoising=Ncomp;
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
                obj.DMIPara.TE=obj.DMIPara.TE-obj.DMIPara.TE(1);
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

                    [~,obj.coilSens,obj.coilNormMat]=adaptiveCombine(sum(obj.img(:,:,:,:,1,:),[5,6]),[1 1 1]*7);
                    obj.img=(sum(obj.coilSens.*obj.img,1));
                    %                     obj.img=obj.img.*permute(obj.coilNormMat,[5 1 2 3 4]); %norm
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

        function getFieldmap(obj)
            if(strcmpi(obj.FieldMap,'IDEAL') && ~any(strcmpi(obj.flags.Solver,{'IDEAL','IDEAL-modes'}) ))
                 %                 fm_meas_Hz=obj.FieldMap./(2*pi)*(6.536 /42.567); % 2H field map in Hz
                im_me=squeeze(sum(obj.img,6))./sqrt(size(obj.img,6)); % sum phase cycles to get FISP contrast
                IdealObj=IDEAL(obj.metabolites,obj.DMIPara.TE,'fm',[],'solver','IDEAL', ...
                    'maxit',obj.flags.maxit,'mask',obj.mask,'SmoothFM',obj.flags.doSmoothFM,'parfor',obj.flags.parfor);
                Metcon_temp=IdealObj'*im_me;
                obj.FieldMap=IdealObj.experimental.fm_est*(-2*pi)/(6.536 /42.567);
            elseif(ischar(obj.FieldMap)&&isfile(obj.FieldMap))
                %field map should be 1H fielmap in rad/s
                assert(exist('spm','file'),'Need spm for registering the field map\n');
                if(isunix)
                    im_space=obj.WriteImages('/tmp/im.nii',{'image'});
                else
                     im_space=obj.WriteImages(fullfile(getenv('temp'),'im.nii'),{'image'});
                end
                obj.Experimental.fm_file=obj.FieldMap;
                obj.FieldMap=myspm_reslice(dir(im_space),dir(obj.Experimental.fm_file), 'linear','r');    
                obj.FieldMap=obj.FieldMap{1};
            elseif(strcmpi(obj.FieldMap,'IDEAL'))
                obj.FieldMap=[];
            end
        end

        function performMetCon(obj)
            PCSel=obj.flags.PCSel;
            EchoSel=obj.flags.EchoSel;
            FA=obj.DMIPara.FlipAngle*obj.DMIPara.pulseCorrectionFactor; %rad
            TE=obj.DMIPara.TE(EchoSel); %s
            TR=obj.DMIPara.TR; %s
            PC=obj.DMIPara.PhaseCycles(PCSel); %rad
            met_struct=obj.metabolites;

            if(strcmpi(obj.flags.Solver,'pinv'))

                tic;
                if(isempty(obj.FieldMap)),obj.FieldMap=zeros(size(obj.mask));end
                B0=obj.FieldMap(obj.mask)./((2*pi)/(6.536 /42.567));%+1e6*(Spectroscopy_para.PVM_FrqWork(1)- ExpPara.PVM_FrqWork(1)); %Hz
                im1=reshape(obj.img(:,:,:,:,EchoSel,PCSel),[],length(EchoSel)*length(PCSel));
                im1=im1(obj.mask(:),:);
                %         [Msig_all]=bSSFP_sim_analytical(metabolites,TE(EchoSel),PC(PhSel),TR,zeros(size(B0)),FA);

                metabol_con=zeros(size(im1,1),length(obj.metabolites));
                resi=zeros(size(im1));
                condnm=zeros(size(im1,1),1);
                sclfac=zeros(size(im1,1),length(obj.metabolites));
                fprintf('Calculating bSSFP profile basis\n');
                Msig_all=MetSignalModel(met_struct,TE,PC,TR,B0,FA,'bSSFP');
                if(TE(1)==0)
                    Msig_all=bsxfun(@times,Msig_all,exp(-1i*angle(Msig_all(:,1,:,:,:,:,:,:,:))));
                end
                fprintf('done.....\n');

                if(obj.flags.parfor)
                    parfor i=1:size(im1,1)
                        A= reshape(Msig_all(:,:,:,1,i,1),length(met_struct),length(PCSel)*length(EchoSel)).';
                        b=[ im1(i,:)];
                        metabol_con(i,:)=A\b(:);
                        resi(i,:)=b(:) -A*metabol_con(i,:).';
                        Ai=pinv(A);
                        sclfac(i,:)=(sum(abs(Ai).^2,2).^(1/2))/sqrt(2);
                        condnm(i)=cond(A);
                    end
                else
                    for i=1:size(im1,1)
                        %low mem mode!
                        %                   [Msig_all]=bSSFP_sim_analytical(met_struct,TE,PC,TR,B0(i),FA);
                        %                   A= reshape(Msig_all,length(met_struct),length(PCSel)*length(EchoSel)).';

                        A= reshape(Msig_all(:,:,:,1,i,1),length(met_struct),length(PCSel)*length(EchoSel)).';
                        b=[ im1(i,:)];
                        metabol_con(i,:)=A\b(:);
                        resi(i,:)=b(:) -A*metabol_con(i,:).';
                        Ai=pinv(A);
                        sclfac(i,:)=(sum(abs(Ai).^2,2).^(1/2))/sqrt(2);
                        condnm(i)=cond(A);
                    end
                end


                obj.Metcon=obj.col2mat(metabol_con,obj.mask);
                obj.Experimental.residue=sum(obj.col2mat(single(resi),obj.mask).^2,[4 5]);
                obj.Experimental.condition=obj.col2mat(single(condnm(:)),obj.mask);
                obj.Experimental.sclfac=obj.col2mat(sclfac,obj.mask);

                fprintf('\n Metabolite fitting done in %0.1f s ! \n',toc)
            elseif(strcmpi(obj.flags.Solver,'IDEAL'))

                %                 fm_meas_Hz=obj.FieldMap./(2*pi)*(6.536 /42.567); % 2H field map in Hz
                im_me=squeeze(sum(obj.img,6))./sqrt(size(obj.img,6)); % sum phase cycles to get FISP contrast
                obj.SolverObj=IDEAL(obj.metabolites,TE,'fm',obj.FieldMap,'solver','IDEAL', ...
                    'maxit',obj.flags.maxit,'mask',obj.mask,'SmoothFM',obj.flags.doSmoothFM,'parfor',obj.flags.parfor);
                obj.Metcon=obj.SolverObj'*im_me;
                obj.Experimental.fm_est=obj.SolverObj.experimental.fm_est;
                obj.Experimental.residue=sum(obj.SolverObj.experimental.residue.^2,[4 5]);
                obj.Experimental.residue=sum(abs(obj.SolverObj.experimental.residue),[4 5 6])./sqrt(2*prod(size(obj.SolverObj.experimental.residue,4:6)));
                obj.SolverObj.experimental=[];


            elseif(strcmpi(obj.flags.Solver,'IDEAL-modes'))


                obj.SolverObj=IDEAL(obj.metabolites,TE,'fm',obj.FieldMap,'solver','IDEAL', ...
                    'maxit',obj.flags.maxit,'mask',obj.mask,'SmoothFM',obj.flags.doSmoothFM,'parfor',obj.flags.parfor);

                Np=floor((length(obj.DMIPara.PhaseCycles)-1)/2);
                if(Np>10), Np=10; end % higher order doesn't hold that much signal
                %calculate SSFP configuration modes
                Fn=calc_Fn2(squeeze(obj.img),obj.DMIPara.PhaseCycles,Np);
                %sqrt(length(obj.DMIPara.PhaseCycles)) scaling for SNR units
                Fn=Fn.*sqrt(length(obj.DMIPara.PhaseCycles));

                %estimate fieldmap from F0
                im_me=Fn(:,:,:,:,Np+1); % F0
                obj.Metcon=obj.SolverObj'*im_me;
                fm_ideal=smooth3(obj.SolverObj.experimental.fm_est)*(-2*pi)/(6.536 /42.567); % scaled to 1H field map in rad/s


                %estimate metabolites concentration from all modes
                Metcon_modes=[];res_all=[];
                for i=1:size(Fn,5)
                    TEs=obj.DMIPara.TE +TR*(i-(Np+1));
                    SolverObj_pinv=IDEAL(obj.metabolites,TEs,'fm',fm_ideal,'solver','phaseonly','mask',obj.mask,'parfor',obj.flags.parfor);
                    Metcon_temp= SolverObj_pinv'*Fn(:,:,:,:,i);

                    Metcon_modes=cat(5,Metcon_modes,Metcon_temp);
                    res_all=cat(6,res_all,SolverObj_pinv.experimental.residue);
                end

                %SVD combine  metcon_all
                sz=size(Metcon_modes);
                Metcon_comb=zeros(size(Metcon_modes));
                obj.Experimental.sing_val=cell(length(obj.metabolites),1);
                obj.Experimental.sing_vec=cell(length(obj.metabolites),1);
                for cMet=1:length(obj.metabolites)
                    [~,S,V]=svd(reshape(Metcon_modes(:,:,:,cMet,:),[],size(Fn,5)),'econ');
                    Metcon_comb(:,:,:,cMet,:)=reshape(reshape(Metcon_modes(:,:,:,cMet,:),[],size(Fn,5))*V,[sz(1:3),1,sz(5)]);
                    obj.Experimental.sing_val{cMet}=diag(S);
                    obj.Experimental.sing_vec{cMet}=V(:,1);
                end
                %pick the first component
                obj.Metcon=Metcon_comb(:,:,:,:,1);
                
                obj.Experimental.Metcon_modes=Metcon_modes;             
                obj.Experimental.Metcon_modes_svd=Metcon_comb;
                obj.Experimental.fm_est=obj.SolverObj.experimental.fm_est;
                obj.Experimental.residue=sum(res_all.^2,[4 5 6]);
                obj.Experimental.Rsq=1-sum(abs(res_all).^2,[4 5 6])./sum(abs(Fn-mean(Fn,[4 5 6])).^2,[4 5 6]);
                obj.SolverObj.experimental=[];
            elseif(strcmpi(obj.flags.Solver,'IDEAL-modes2'))
                %we try to fit all modes together

                obj.SolverObj=IDEAL(obj.metabolites,TE,'fm',obj.FieldMap,'solver','IDEAL', ...
                    'maxit',obj.flags.maxit,'mask',obj.mask,'SmoothFM',1,'parfor',obj.flags.parfor);

                Np=round(length(obj.DMIPara.PhaseCycles)/2);
                if(Np>10), Np=10; end % higher order doesn't hold that much signal
                %calculate SSFP configuration modes
                Fn=calc_Fn2(squeeze(obj.img),obj.DMIPara.PhaseCycles,Np);
                Fn=Fn.*sqrt(length(obj.DMIPara.PhaseCycles)); %SNR scaling

                %estimate fieldmap from F0
                im_me=Fn(:,:,:,:,Np+1); % F0
                obj.Metcon=obj.SolverObj'*im_me;
                fm_ideal=obj.SolverObj.experimental.fm_est*(-2*pi)/(6.536 /42.567); % scaled to 1H field map in rad/s



                TEs=[];
                for i=1:size(Fn,5)
                    TEs=[TEs;(obj.DMIPara.TE(:) +TR*(i-(Np+1)))];
                end
                    SolverObj_pinv=IDEAL(obj.metabolites,TEs','fm',fm_ideal,'solver','phaseonly','mask',obj.mask,'parfor',obj.flags.parfor);
                    obj.Metcon= SolverObj_pinv'*Fn(:,:,:,:);
                obj.Experimental.fm_est=obj.SolverObj.experimental.fm_est;
                obj.Experimental.residue=sos(SolverObj_pinv.experimental.residue,[4 5]);
                obj.SolverObj.experimental=[];


            elseif(strcmpi(obj.flags.Solver,'phaseonly'))

                %                 fm_meas_Hz=obj.FieldMap./(2*pi)*(6.536 /42.567); % 2H field map in Hz
                im_me=squeeze(sum(obj.img,6)); % sum phase cycles to get FISP contrast
                obj.SolverObj=IDEAL(obj.metabolites,TE,'fm',obj.FieldMap,'solver','phaseonly','mask',obj.mask);
                obj.Metcon= obj.SolverObj'*im_me;
                obj.Experimental.fm_est=[];
                obj.Experimental.residue=obj.SolverObj.experimental.residue;
                obj.SolverObj.experimental=[];

            end
            obj.flags.PxlShiftPerformed=false;

        end

        function pxlShift=performPxlShiftCorrection(obj)
            % correct the pixel shift along read direction due to
            % off-resonance
            if(~isempty(obj.Metcon) && ~obj.flags.PxlShiftPerformed )
            
            dx=obj.DMIPara.resolution(1); %m
            Nx=obj.DMIPara.MatrixSize(1);
            dt=obj.DMIPara.dwell; %s
            gammaH2=obj.DMIPara.gammaH2*1e6; %Hz/T

            Gread=1/(dt*gammaH2*Nx*dx); %T/m
            freq=[obj.metabolites.freq_shift_Hz]; % Hz
            pxlShift=freq./(Gread*gammaH2); %m
            Metcon_temp=permute(obj.Metcon,[2 1 3 4 5]); % first dim read
            for cMet=1:length(obj.metabolites)
               readAxis=linspace(0,dx*(Nx-1),size(Metcon_temp,1));
               Metcon_temp(:,:,:,cMet)=interp1(readAxis,Metcon_temp(:,:,:,cMet),readAxis(:)+pxlShift(cMet),'linear',0);
            
            end
            obj.Metcon=ipermute(Metcon_temp,[2 1 3 4 5]);
            fprintf('Performed pixel shift along read: (%.1f,%.1f,%.1f ,%.1f) mm\n',pxlShift*1e3);
            obj.flags.PxlShiftPerformed=true;
            end
            
        end
        function ShiftMetcon(obj,Shift_PRS_mm)
        % crappy motion correction: only translation
        Metcon_temp=abs(obj.Metcon);
        FOV_RPS=obj.DMIPara.MatrixSize.*obj.DMIPara.resolution(1:3)*1e3; %mm
        imSize_PRS=size(obj.Metcon);
        res_PRS=FOV_RPS([2 1 3])./imSize_PRS(1:3);
        Phase_mm=linspace(0,FOV_RPS(2)-res_PRS(1),imSize_PRS(1));
        Metcon_temp=interp1(Phase_mm,Metcon_temp,Phase_mm+Shift_PRS_mm(1),'linear','extrap');

        Metcon_temp=permute(Metcon_temp,[2 1 3 4 5]); % first dim read
        Read_mm=linspace(0,FOV_RPS(1)-res_PRS(2),imSize_PRS(2));
        Metcon_temp=interp1(Read_mm,Metcon_temp,Read_mm+Shift_PRS_mm(2),'linear','extrap');
        Metcon_temp=ipermute(Metcon_temp,[2 1 3 4 5]); % first dim read

        Metcon_temp=permute(Metcon_temp,[3 1 2 4 5]); % first dim read
        Slice_mm=linspace(0,FOV_RPS(3)-res_PRS(3),imSize_PRS(3));
        Metcon_temp=interp1(Slice_mm,Metcon_temp,Slice_mm+Shift_PRS_mm(3),'linear','extrap');
        Metcon_temp=ipermute(Metcon_temp,[3 1 2 4 5]); % first dim read
        obj.Metcon=Metcon_temp;
        
        end


        function demoFit(obj,voxel_idx)
            % plot the fittting results for the voxel with index for debuggig
            % (voxel_idx: 3x1 array)
            assert(strcmpi(obj.flags.Solver,'pinv'),'only implemented for pinv mode');
            vMask=zeros(size(obj.mask),'logical');
            vMask(voxel_idx(1),voxel_idx(2),voxel_idx(3))=true;
            Mask_baskup=obj.mask;
            obj.mask=vMask;
            obj.Metcon();

            PCSel=obj.flags.PCSel;
            EchoSel=obj.flags.EchoSel;
            FA=obj.DMIPara.FlipAngle*obj.DMIPara.pulseCorrectionFactor; %rad
            TE=obj.DMIPara.TE(EchoSel); %s
            TR=obj.DMIPara.TR; %s
            PC=obj.DMIPara.PhaseCycles(PCSel); %rad
            met_struct=obj.metabolites;
            B0=obj.FieldMap(obj.mask)./(2*pi)*(6.536 /42.567);
            %             [Msig_all]=bSSFP_sim_analytical(met_struct,TE,PC,TR,B0,FA);

            Msig_all=MetSignalModel(met_struct,TE,PC,TR,B0,FA,'bSSFP');
            Msig_all = reshape(Msig_all,length(met_struct),length(EchoSel),length(PCSel));
            imdata=squeeze(obj.img(1,voxel_idx(1),voxel_idx(2),voxel_idx(3),:,:));
            figure(17),clf
            tt=tiledlayout(2,3,'TileSpacing','tight');
            Met_amp=squeeze(obj.Metcon(voxel_idx(1),voxel_idx(2),voxel_idx(3),:));
            disp(abs(Met_amp))
            Msig_scaled=Msig_all(:,:).'*Met_amp(:);


            nexttile();
            plot(rad2deg(PC),abs(reshape(Msig_scaled,[],length(PC)).'))
            hold on
            plot(rad2deg(PC),abs(imdata).','.')
            ylabel('amplitude [a.u]')
            xlabel('RF phase increment [deg]')
            title('magnitude: data vs fit')
            set(gca,'colororder',jet(size(imdata,1)))

            %plot phase wrt PC
            nexttile();
            plot(rad2deg(PC),rad2deg(angle(reshape(Msig_scaled,[],length(PC)).')))
            hold on
            plot(rad2deg(PC),rad2deg(angle(imdata).'),'.')
            set(gca,'colororder',jet(size(imdata,1)))
            ylabel('phase [deg]')
            xlabel('RF phase increment [deg]')
            title('Phase: data vs fit')


            for i=1:length(met_struct)
                nexttile();
                plot(rad2deg(PC),abs(reshape(Msig_all(i,:,:),[],length(PC)).'))
                yyaxis right
                plot(rad2deg(PC),rad2deg(angle(reshape(Msig_all(i,:,:),[],length(PC)).')),'--')
                set(gca,'colororder',jet(size(imdata,1)))
                title([obj.metabolites(i).name,' basis'])
            end

            obj.mask=Mask_baskup;


        end
     
        function niiFileName=WriteImages(obj,niiFileName,Select,norm_mat)

            % write nifti files of avaearged ME-images, Metabolite
            % amplitudes in SNR units and metabolite concentrations in mM
            % eg: mcobj.WriteImages('',{'image','snr','mm','fm'})
            if(nargin==1 || isempty(niiFileName))
                fPath=pwd;
                fn=sprintf('M%05d_%s.nii',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription);
                niiFileName=fullfile(fPath,fn);
            elseif(isfolder(niiFileName))
                fPath=niiFileName;
                fn=sprintf('M%05d_%s.nii',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription);
                niiFileName=fullfile(fPath,fn);
            else
                [fPath,~]=fileparts(niiFileName);
            end
            if(~exist("Select",'var')),Select={'image','SNR','mM'};end            
            if(~exist("norm_mat",'var')),norm_mat=[];end

            if(any(strcmpi(Select,'image')))
            vol_PRS=single(squeeze(sos(flip(obj.img,3),[5 6]))); % 9.4T specific read flip
            description='averaged image across echo and phase cycle';
            MyNIFTIWrite_ME(vol_PRS,obj.twix,niiFileName,description);
            end
            if(~isempty(obj.Metcon)&&any(strcmpi(Select,'SNR')))
                fn2=sprintf('Metcon_SNR_m%05d_%s_%s.nii',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription,obj.flags.Solver);
                vol_PRS=single(abs(flip(obj.getNormalized,2))); %9.4T specific read flip
                description=sprintf('dim4_%s_%s_%s_%s_',obj.metabolites.name);
                MyNIFTIWrite_ME(squeeze(single(abs(vol_PRS))),obj.twix,fullfile(fPath,fn2),description);
            end
            if(~isempty(obj.Metcon)&&any(strcmpi(Select,'mM')))
                fn3=sprintf('Metcon_mM_m%05d_%s_%s.nii',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription,obj.flags.Solver);
                vol_PRS=flip(obj.getmM(norm_mat),2); %9.4T specific
                description=sprintf('dim4_%s_%s_%s_%s',obj.metabolites.name);
                MyNIFTIWrite_ME(squeeze(single(abs(vol_PRS))),obj.twix,fullfile(fPath,fn3),description);
            end
            %export fir
            if(isfield(obj.Experimental,'fm_est') && ~isempty(obj.Experimental.fm_est)&&any(strcmpi(Select,'fm')))
                fn3=sprintf('fm_Hz_m%05d_%s_%s.nii',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription,obj.flags.Solver);
                vol_PRS=flip(obj.Experimental.fm_est,2); %9.4T specific
                description='fieldmap estimated with IDEAL';
                MyNIFTIWrite_ME(squeeze(single((vol_PRS))),obj.twix,fullfile(fPath,fn3),description);
            end

        end

        function SaveResults(obj)
            % save results toa MAT file


            OutFile=sprintf('m%d_%s_%s.mat',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription,obj.flags.Solver);
            mcobj_copy=copy(obj);
            mcobj_copy.twix=[];
            mcobj_copy.sig=[];
%             varName=sprintf('mcobj_m%d_%s',obj.twix.hdr.Config.MeasUID,obj.flags.Solver);
%             evalc(sprintf('%s=mcobj_copy;',varName));
            save(OutFile,"mcobj_copy",'-v7.3');


        end
        function PlotResults(obj,fh)
            if(~exist('fh','var'))
                fh=figure;
            end


            tt=tiledlayout(fh,2,ceil(length(obj.metabolites)/2+1),'TileSpacing','compact','Padding','compact');
            slcFac=0.3;slcdim=2;
            slcSel=round(slcFac*size(obj.Metcon,slcdim));
            slcSel=slcSel:(size(obj.Metcon,slcdim)-slcSel);

            imtransFunc=@(x) flip(flip(permute(x(:,slcSel,:),[1 3 slcdim 4]),10),30);
            %get metcon in normalized SNR units
            [im_snr,sf]=obj.getNormalized();
            for i=1:length(obj.metabolites)
                nexttile(tt)
                im_curr=createImMontage(imtransFunc(abs(im_snr(:,:,:,i))));
                im_curr(isnan(im_curr))=0;
                imagesc(im_curr);
                colorbar,axis image
                cax_im=[0 prctile(im_curr(:),99)];
                clim(cax_im);
                xticks([]),yticks([]),title(obj.metabolites(i).name)
            end

            % display residue
            nexttile(tt)
            residue_norm=sos(obj.Experimental.residue,[4 5])./mean(sf);
            im_curr=createImMontage(imtransFunc(abs(residue_norm)));
            imagesc(im_curr);
            colorbar,axis image

            cax_im=[0 prctile(im_curr(:),95)];
            clim(cax_im);
            xticks([]),yticks([]),title('residue')

            %display fieldmaps
            nexttile(tt)
            if(isfield(obj.Experimental,'fm_est') && ~isempty(obj.Experimental.fm_est))
                fm_Hz=createImMontage(imtransFunc(obj.Experimental.fm_est(:,:,:,:)));
                imagesc(fm_Hz);
                title('estimated 2H fieldmap [Hz]')
            elseif(~isempty(obj.FieldMap))
                fm_Hz=createImMontage(imtransFunc(obj.FieldMap(:,:,:,:)))./(2*pi)*(6.536 /42.567);
                imagesc(fm_Hz);
                title('input 2H fieldmap [Hz]')
            else
                title('no estimated/input fieldmap');
            end
        

            colorbar,
            axis image,colormap(gca,'jet');
            cax_im=[-1*prctile(fm_Hz(:),95) prctile(fm_Hz(:),95)+1];
            clim(cax_im);
            xticks([]),yticks([]),

            str=sprintf('%s|%s',obj.DMIPara.ShortDescription,obj.flags.Solver);
            if(isfield(obj.DMIPara,'IntakeTime'))
                str=sprintf('%s|%d min',str,obj.getMinutesAfterIntake());
            end

            annotation('textbox',[0.5 0.9 0.1 0.1],'String',str,'HorizontalAlignment','center','FontSize',12)

            OutFile=sprintf('%s.fig',strrep(str,'|','_'));
            OutFile=strrep(OutFile,' ','');
            savefig(fh,OutFile)

        end
        function [Metcon_norm,scale_fac]= getNormalized(obj)
            % normalize metabolite maps same as plotresults()
            if(any(strcmp(obj.flags.Solver,{'IDEAL','IDEAL-modes','phaseonly'}))) 
                %analytical SNR
                Ai=pinv(obj.SolverObj.getA() );
                scale_fac=(sum(abs(Ai).^2,2).^(1/2))/sqrt(2);
                Metcon_norm=abs(obj.Metcon)./reshape(scale_fac,1,1,1,[]);
            elseif(any(strcmp(obj.flags.Solver,{'pinv'}))) 
                Metcon_norm=abs(obj.Metcon)./obj.Experimental.sclfac;
                scale_fac=squeeze(mean(obj.Experimental.sclfac,[1 2 3]));
            else % scale by std of noise region
                Metcon_norm=zeros(size(obj.Metcon));
                noiseMask=~obj.getMask(90);
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

        function [Metcon_mM,scale_fac]= getmM(obj,norm_mat)
%           %[Metcon_mM,scale_fac]= getmM(obj,norm_mat)
            %norm_mat is the smoothed water image for scaling other metabolites
            % converts metabolite amplitudes into mM
            % assuming water concentration is 10 mM 
            %reference: peters et al DOI: 10.1002/mrm.28906 , equation 6
            
            TR=obj.DMIPara.TR;
            TE=obj.DMIPara.TE;
            DC=obj.DMIPara.DutyCycle;
            %Actual reference voltage is higher than set refvoltage (500-550 V) 
            FA_rad=obj.DMIPara.FlipAngle*(obj.DMIPara.pulseCorrectionFactor);

            if(exist('norm_mat','var')&&~isempty(norm_mat))
            assert(isequal(size(norm_mat),[size(obj.Metcon,1),size(obj.Metcon,2),size(obj.Metcon,3)]),...
                'Size of input norm_mat doesn''t match the obj.MetCon size');
            assert(isreal(norm_mat),'input norm_mat should not be complex');
            else
                mc=obj.getNormalized;
                norm_mat=1./imgaussfilt3(abs(mc(:,:,:,1)),2);
            end
            

            if(~contains(obj.twix.hdr.Config.SequenceFileName,'trufi'))
                [Msig_all,dc_fac]=MetSignalModel(obj.metabolites,TE,0,TR,0,FA_rad,'FISP',DC);
                Msig_all=Msig_all.*dc_fac;
            else
                [Msig_all,dc_fac]=MetSignalModel(obj.metabolites,TE,obj.DMIPara.PhaseCycles,TR,0,FA_rad,'bSSFP',DC);
                Msig_all=Msig_all.*dc_fac;
            end

            Sig_theory=mean(abs(squeeze(Msig_all)),[2 3 4]);
            Sig_theory=Sig_theory(1)./Sig_theory;
            % 111 M H20 *0.0115% 2H *80% water in brain= 10.12 mM
            % 1.33 Glx label loss in TCA cycle De Feyter et al 2018
            Average_deuterons=[1;2;1.33;2];
            scale_fac=10.12*Sig_theory./Average_deuterons(1:length(obj.metabolites)); % mM
            metcon_w=abs(obj.getNormalized).*norm_mat;
            Metcon_mM=metcon_w.*permute(scale_fac(:),[2 3 4 1]);
%             as(Metcon_mM,'select',':,:,25,3','windowing',[1.5 3])
            % all higher values are probably noise
            Metcon_mM(Metcon_mM>11)=0;
            Mask50=obj.getMask(50);
            Metcon_mM(:,:,:,1)=Metcon_mM(:,:,:,1).*Mask50;
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
        function setRefereneceVoltage(obj,RefVoltage)
        % as we cannot go beyond 447 V as reference voltage, we correct the
        % flip angles retrospectively with the actual reference voltage 
        % (550 V is default for pulseCorrectionFactor calculation)
        obj.DMIPara.pulseCorrectionFactor=obj.DMIPara.RefVoltage/RefVoltage;
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


