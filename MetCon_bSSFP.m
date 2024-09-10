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
                    addParameter(p,'parfor',true,@(x) islogical(x));


                    addParameter(p,'doZeroPad',[1,1,1],@(x) isvector(x));
                    %-1 for debug
                    addParameter(p,'doDenoising',0,@(x) isscalar(x));
                    %                     addParameter(p,'doParFor',inf,@(x) isscalar(x));

                    addParameter(p,'CoilSel',1:obj.twix.image.NCha, @(x) isvector(x));
                    addParameter(p,'PCSel',1:obj.twix.image.NRep, @(x) (isvector(x)&& all(x<=1:obj.twix.image.NRep)) );
                    addParameter(p,'EchoSel',1:obj.twix.hdr.Phoenix.lContrasts, @(x) (isvector(x)&& all(x<=obj.twix.hdr.Phoenix.lContrasts)) );
                    addParameter(p,'is3D',(obj.twix.image.NPar>1),@(x)islogical(x));
                    %                     addParameter(p,'doB0Corr','none',@(x) any(strcmp(x,{'none','MTI','MFI'})));
                    %                   addParameter(p,'precision','single',@(x) any(strcmp(x,{'single','double'})));
                    addParameter(p,'Solver','pinv',@(x) any(strcmp(x,{'pinv','IDEAL','IDEAL-modes','IDEAL-modes2','phaseonly'})));
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
                        warning('Not enough phase cycles for IDEAL-modes: Using IDEAL')
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

            fprintf(sprintf( 'reco  time = %6.1f s\n', toc));
            if(~isempty(obj.metabolites))
                obj.performMetCon();
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

                    [~,obj.coilSens,obj.coilNormMat]=adaptiveCombine(sum(obj.img(:,:,:,:,1,:),6),[1 1 1]*7);
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

        function performMetCon(obj)
            PCSel=obj.flags.PCSel;
            EchoSel=obj.flags.EchoSel;
            FA=obj.DMIPara.FlipAngle; %rad
            TE=obj.DMIPara.TE(EchoSel); %s
            TR=obj.DMIPara.TR; %s
            PC=obj.DMIPara.PhaseCycles(PCSel); %rad
            met_struct=obj.metabolites;

            if(strcmpi(obj.flags.Solver,'pinv'))

                tic;
                if(isempty(obj.FieldMap)),obj.FieldMap=zeros(size(obj.mask));end
                B0=obj.FieldMap(obj.mask)./(2*pi)*(6.536 /42.567);%+1e6*(Spectroscopy_para.PVM_FrqWork(1)- ExpPara.PVM_FrqWork(1)); %Hz
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

                %                 pb = parwaitbar(size(im1,1),'WaitMessage','Least squares with bSSFP model:estimating concentrations');
                for i=1:size(im1,1)
                    %low mem mode!
                    %                   [Msig_all]=bSSFP_sim_analytical(met_struct,TE,PC,TR,B0(i),FA);
                    %                   A= reshape(Msig_all,length(met_struct),length(PCSel)*length(EchoSel)).';

                    A= reshape(Msig_all(:,:,:,1,i,1),length(met_struct),length(PCSel)*length(EchoSel)).';
                    b=[ im1(i,:)];
                    metabol_con(i,:)=A\b(:);
                    resi(i,:)=b(:) -A*metabol_con(i,:).';
                    condnm(i)=cond(A);
                    %                     pb.progress();
                    Ai=pinv(A);
                   sclfac(i,:)=(sum(abs(Ai).^2,2).^(1/2))/sqrt(2);
                end


                obj.Metcon=obj.col2mat(metabol_con,obj.mask);
                obj.Experimental.residue=obj.col2mat(single(resi),obj.mask);
                obj.Experimental.condition=obj.col2mat(single(condnm(:)),obj.mask);
                obj.Experimental.sclfac=obj.col2mat(sclfac,obj.mask);

                fprintf('\n Metabolite fitting done in %0.1f s ! \n',toc)
            elseif(strcmpi(obj.flags.Solver,'IDEAL'))

                %                 fm_meas_Hz=obj.FieldMap./(2*pi)*(6.536 /42.567); % 2H field map in Hz
                im_me=squeeze(sum(obj.img,6)); % sum phase cycles to get FISP contrast
                obj.SolverObj=IDEAL(obj.metabolites,TE,'fm',obj.FieldMap,'solver','IDEAL', ...
                    'maxit',10,'mask',obj.mask,'SmoothFM',1,'parfor',obj.flags.parfor);
                obj.Metcon=obj.SolverObj'*im_me;
                obj.Experimental.fm_est=obj.SolverObj.experimental.fm_est;
                obj.Experimental.residue=sos(obj.SolverObj.experimental.residue,[4 5]);
                obj.SolverObj.experimental=[];


            elseif(strcmpi(obj.flags.Solver,'IDEAL-modes'))


                obj.SolverObj=IDEAL(obj.metabolites,TE,'fm',obj.FieldMap,'solver','IDEAL', ...
                    'maxit',10,'mask',obj.mask,'SmoothFM',1,'parfor',obj.flags.parfor);

                Np=round(length(obj.DMIPara.PhaseCycles)/2);
                if(Np>10), Np=10; end % higher order doesn't hold that much signal
                %calculate SSFP configuration modes
                Fn=calc_Fn2(squeeze(obj.img),obj.DMIPara.PhaseCycles,Np);
                %sqrt(Np*2+1) scaling for SNR units
                Fn=Fn.*sqrt(Np*2+1);

                %estimate fieldmap from F0
                im_me=Fn(:,:,:,:,Np+1); % F0
                obj.Metcon=obj.SolverObj'*im_me;
                fm_ideal=obj.SolverObj.experimental.fm_est*(-2*pi)/(6.536 /42.567); % scaled to 1H field map in rad/s


                %estimate metabolites concentration from all modes
                Metcon_modes=[];res_all=[];
                for i=1:size(Fn,5)
                    TEs=obj.DMIPara.TE +TR*(i-(Np+1));
                    SolverObj_pinv=IDEAL(obj.metabolites,TEs,'fm',fm_ideal,'solver','phaseonly','mask',obj.mask,'parfor',obj.flags.parfor);
                    Metcon_temp= SolverObj_pinv'*Fn(:,:,:,:,i);

                    Metcon_modes=cat(5,Metcon_modes,Metcon_temp);
                    res_all=cat(5,res_all,sos(SolverObj_pinv.experimental.residue,[4 5]));
                end

                %SVD combine  metcon_all

                [~,S,V]=svd(reshape(Metcon_modes,[],size(Fn,5)),'econ');
                sz=size(Metcon_modes);
                Metcon_comb=reshape(reshape(Metcon_modes,[],size(Fn,5))*V,sz(1:5));
                %pick the first
                obj.Metcon=Metcon_comb(:,:,:,:,1);

                obj.Experimental.Metcon_modes=Metcon_modes;
                obj.Experimental.sing_val=S;
                obj.Experimental.Metcon_modes_svd=Metcon_comb;
                obj.Experimental.fm_est=obj.SolverObj.experimental.fm_est;
                obj.Experimental.residue=sos(res_all,[4 5]);
                obj.SolverObj.experimental=[];
            elseif(strcmpi(obj.flags.Solver,'IDEAL-modes2'))
                %we try to fit all modes togther

                obj.SolverObj=IDEAL(obj.metabolites,TE,'fm',obj.FieldMap,'solver','IDEAL', ...
                    'maxit',10,'mask',obj.mask,'SmoothFM',1,'parfor',obj.flags.parfor);

                Np=round(length(obj.DMIPara.PhaseCycles)/2);
                if(Np>10), Np=10; end % higher order doesn't hold that much signal
                %calculate SSFP configuration modes
                Fn=calc_Fn2(squeeze(obj.img),obj.DMIPara.PhaseCycles,Np);

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

        end

        function plotFit(obj,voxel_idx)
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
            FA=obj.DMIPara.FlipAngle; %rad
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
     
        function niiFileName=WriteImages(obj,niiFileName)
            %average echo and phase cycling dimension and write nifti
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

            vol_PRS=single(squeeze(sos(flip(obj.img,3),[5 6]))); % 9.4T specific read flip
            description='averaged image across echo and phase cycle';
            MyNIFTIWrite_bSSFP3(vol_PRS,obj.twix,niiFileName,description);

            if(~isempty(obj.Metcon))
                fn2=sprintf('Metcon_m%d_%s_%s.nii',obj.twix.hdr.Config.MeasUID,obj.twix.hdr.Config.SequenceDescription,obj.flags.Solver);
                vol_PRS=single(abs(flip(obj.getNormalized,2))); %9.4T specific read flip
                description=sprintf('dim4_%s_%s_%s_%s_',obj.metabolites.name);
                MyNIFTIWrite_bSSFP2(squeeze(single(abs(vol_PRS))),obj.twix,fullfile(fPath,fn2),description);
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

        function [Metcon_mM,scale_fac]= getmM(obj)
            % convert metabolite amplitudes into mM
            % assuming water concentration is 10 mM 
            %reference: peters et al DOI: 10.1002/mrm.28906 , equation 6
            
            TR=obj.DMIPara.TR;
            TE=obj.DMIPara.TE(1);
            DC=obj.DMIPara.DutyCycle;
            %Actual reference voltage is between 500 and 550
            FA_rad=obj.DMIPara.FlipAngle*(obj.DMIPara.RefVoltage/500);

            if(~contains(obj.twix.hdr.Config.SequenceFileName,'trufi'))
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


