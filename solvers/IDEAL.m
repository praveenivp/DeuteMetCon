classdef IDEAL < matlab.mixin.Copyable
    %     IDEALsolver(metabolties,TE_s,FieldMap_Hz,Mode{'pinv','IDEAL'},nIteration)
    properties
        metabolites % struvt array struct('freq_shift_Hz',79,'con',randi([2,10],1)/10,'name','1.3 ppm','T1_s',297e-3,'T2_s',61e-3)
        FieldMap_Hz % Initial field map in Hz
        TE_s % Echo time in seconds
        transp=false; % transpose flag
        experimental% experimental outputs
        flags
        mask
    end
    methods
        function obj=IDEAL(varargin)

            if(nargin==0)
                obj.TE_s=[2.2:2.9:15]*1e-3; %s
                obj.metabolites=[struct('freq_shift_Hz',79,'con',randi([2,10],1)/10,'name','1.3 ppm','T1_s',297e-3,'T2_s',61e-3),...
                    struct('freq_shift_Hz',226,'con',randi([5,10],1)/10,'name','3.7 ppm ','T1_s',64e-3,'T2_s',32e-3),...
                    struct('freq_shift_Hz',293,'con',randi([8,10],1)/10,'name','4.8 ppm','T1_s',320e-3,'T2_s',12e-3)];
                obj.FieldMap_Hz=0;
            elseif(nargin==2)
                obj.metabolites=varargin{1}; %struct
                obj.TE_s=varargin{2}; %s
                obj.getFlag();
            elseif (nargin>2)
                obj.metabolites=varargin{1}; %struct
                obj.TE_s=varargin{2}; %s
                obj.getFlags(varargin{3:end})
            else
                error('Input parameter: IDEALsolver(metabolties,TE_s)');
            end


        end
        function [data]=performPhaseCorr(obj,data)
            % Stupid fucntion to subtract the phase of first echo from all
            % other echo.
            if(obj.flags.PhaseCorr &&  obj.TE_s(1)>0)
                obj.TE_s= obj.TE_s-obj.TE_s(1);
                data=data.*exp(-1i*angle(data(:,:,:,1,1)));
            end
        end

        function res = mtimes(obj,inp)
            %             bb=[nFE,nIntlv,nPar,nCha]% adj case
            %             bb=[ImX,Imy,Imz,nCha] % forward case
            res=0;
            if (obj.transp)
                A=getA(obj);
                Ainv=pinv(A);
                if(isempty(obj.mask))
                    obj.mask=ones(size(inp,1),size(inp,2),size(inp,3))>0;
                end

                inp= obj.performPhaseCorr(inp);
                for cPC=1:size(inp,5)
                    inp_col=reshape(obj.mat2col(inp(:,:,:,:,cPC),obj.mask),[],length(obj.TE_s));

                    res=zeros(size(inp_col,1),size(A,2));
                    if(all(obj.FieldMap_Hz(:)==0) || isempty(obj.FieldMap_Hz) )
                        fm_est= zeros(size(inp_col,1),1);
                    else
                        fm_est=obj.mat2col(obj.FieldMap_Hz,obj.mask);
                    end
                    fm_est(isnan(fm_est))=0;

                    if(strcmpi(obj.flags.solver,'pinv'))
                        for i=1:size(inp_col,1)
                            inp_corr=inp_col(i,:).*exp(-1i*2*pi*fm_est(i)*obj.TE_s);
                            res(i,:)=Ainv*[real(inp_corr) imag(inp_corr)]';
                        end

                    else % IDEAL algorithm

                        pb = waitbar(0, 'estimating field map(1/2)');
                        %                     fm_est_delta=zeros(numel(fm_Hz),1);
                        for jj=1 % number of fieldmap refinements
                            for i=1:size(inp_col,1)
                                all_B0=[];
                                for it=1:obj.flags.maxit

                                    Sn_cap=inp_col(i,:).*exp(-1i*2*pi*fm_est(i)*obj.TE_s); % remove Known b0 offresonance
                                    metabol_con_cv=Ainv*[real(Sn_cap) imag(Sn_cap)]';
                                    %                             metabol_con(i,:)=metabol_con_cv

                                    metabol_est=complex(metabol_con_cv(1:length(obj.metabolites)),metabol_con_cv(length(obj.metabolites)+1:end));
                                    Sn_doublecap=getS_doublehat(obj,metabol_est,Sn_cap);
                                    B=getB(obj,metabol_est);
                                    delta_= pinv(B)*Sn_doublecap;
                                    %                             fprintf('%3.4f, ',delta_(1))
                                    %                             all_B0=[all_B0 delta_(1)];
                                    if(abs(delta_(1))<obj.flags.tol) %figure,plot(all_B0);
                                        break; end
                                    fm_est(i)=fm_est(i)+delta_(1);

                                end
                                waitbar(i/size(inp_col,1), pb);
                            end
                            fm_est=obj.col2mat(fm_est,obj.mask);
                            %                     fm_est=imfilter(fm_est,ones(2,2)/4,'replicate');
                            fm_est=medfilt3(fm_est,[1 1 1]*(2*(3-jj)+1));
                            fm_est=obj.mat2col(fm_est,obj.mask);
                            %     fprintf('iter :%d fm_delta_norm: %.4f \n',it,norm(fm_est_delta));
                        end
                        fm_est=obj.col2mat(fm_est,obj.mask);
                        fm_est=medfilt3(fm_est,[1 1 1]*3);
                        %     as(fm_est.*(sum(metbol_mask,4)>0),'title','estimated field map(Hz)','colormap','jet')
                        fm_est=obj.mat2col(fm_est,obj.mask);
                        close(pb)

                        pb = waitbar(0, 'estimating concentrations(2/2)');
                        for i=1:size(inp_col,1)
                            Sn_cap=inp_col(i,:).*exp(-1i*2*pi*fm_est(i)*obj.TE_s); % remove Known b0 offresonance
                            metabol_con_cv=Ainv*[real(Sn_cap) imag(Sn_cap)]';
                            res(i,:)=metabol_con_cv;
                            waitbar(i/size(inp_col,1), pb);
                        end
                        close(pb)

                        obj.experimental.fm_est{cPC}=obj.col2mat(fm_est,obj.mask);
                    end
                    res=reshape(res,size(res,1),length(obj.metabolites),2);
                    res=complex(res(:,:,1),res(:,:,2));
                    res=obj.col2mat(res,obj.mask);
                    res_all{cPC}=res;
                end
                res=cat(5,res_all{:});
            else
                error('forward operation not implemented')
            end
            obj.transp = false;
        end



        function res = ctranspose(obj)
            obj.transp = true;%xor(obj.transp,1);
            res = obj;
        end
        function res = transpose(obj)
            obj.transp = true;%xor(obj.transp,1);
            res = obj;
        end



        %% supporting functions
        function A=getA(obj)
            CDij=@(freq,tn) exp(1i*2*pi*freq*tn); % Cij i->chemical shift(Hz) tn-> time(s)
            A=zeros(2*length(obj.TE_s), 2*length(obj.metabolites));
            for cMet=1:length(obj.metabolites)
                CD_cMet=CDij(obj.metabolites(cMet).freq_shift_Hz,obj.TE_s(:));
                A(:,((1:2)+2*(cMet-1)))=  [[real(CD_cMet(:)); imag(CD_cMet(:))]  ,[ imag(conj(CD_cMet(:))); real(conj(CD_cMet(:)));] ];
            end
        end
        function B=getB(obj,metabol_est)
            % metabol_est should be in complex
            CDij=@(freq,tn) exp(1i*2*pi*freq*tn); % Cij i->chemical shift(Hz) tn-> time(s)
            A=getA(obj);
            Gjn_real=zeros(length(obj.TE_s),1);
            Gjn_imag=zeros(length(obj.TE_s),1);
            for cMet=1:length(obj.metabolites)
                CD_cMet=CDij(obj.metabolites(cMet).freq_shift_Hz,obj.TE_s(:));
                Cjn=real(CD_cMet(:));
                Djn=imag(CD_cMet(:));
                Gjn_real= Gjn_real+(-1*real(metabol_est(cMet)).*Djn-1*imag(metabol_est(cMet)).*Cjn);
                Gjn_imag= Gjn_imag+(real(metabol_est(cMet)).*Cjn-1*imag(metabol_est(cMet)).*Djn);
            end
            B=[ [2*pi*obj.TE_s(:).*Gjn_real; 2*pi*obj.TE_s(:).*Gjn_imag], A];
        end

        function S_hathat=getS_doublehat(obj,metabol_est,sig)
            % metabol_est should be in complex
            %sig is echo images
            sig=sig(:);
            CDij=@(freq,tn) exp(1i*2*pi*freq*tn); % Cij i->chemical shift(Hz) tn-> time(s)

            metabol_est_real=real(metabol_est);
            metabol_est_imag=imag(metabol_est);


            sum_real=zeros(size(sig));
            sum_imag=zeros(size(sig));
            for cMet=1:length(obj.metabolites)
                CD_cMet=CDij(obj.metabolites(cMet).freq_shift_Hz,obj.TE_s(:));
                Cjn=real(CD_cMet);
                Djn=imag(CD_cMet);
                sum_real= sum_real+metabol_est_real(cMet).*Cjn-metabol_est_imag(cMet).*Djn;
                sum_imag= sum_imag+metabol_est_real(cMet).*Djn+metabol_est_imag(cMet).*Cjn;
            end
            S_hathat=[real(sig)-sum_real(:); imag(sig)-sum_imag(:)];
        end





        function sig=simSig(obj,MatSz,fm_Scale,noise_scale)
            % generate some circles and call it pahnatom
            obj.FieldMap_Hz= fm_Scale*getRandomFmap(obj,[MatSz(:); 1]);
            Phantom=zeros([MatSz(:)' 1 length(obj.metabolites)]);
            pos=[0.35 0.35; 0.65 0.65; 0.35 0.65; 0.65 0.35;  ].*MatSz(:)';
            r2=(MatSz(1)*0.3).^2;
            for k=1:length(obj.metabolites)
                for i=1:MatSz(1)
                    for j=1:MatSz(2)
                        Phantom(i,j,1,k)= (( (i-pos(k,2)).^2+ (j-pos(k,1)).^2 ) <r2);
                    end
                end
            end
            obj.experimental.Phantom=Phantom;
            initial_phase=0;
            T2star=30e-3; %ms
            %relaxation is ignored
            sig=zeros(size(Phantom,1),size(Phantom,2),1,length(obj.TE_s));
            for k=1:length(obj.metabolites)
                for t=1:length(obj.TE_s)
                    sig(:,:,:,t)=sig(:,:,:,t)+ obj.metabolites(k).con.*Phantom(:,:,:,k).* ...
                        exp(1i*(initial_phase+2*pi*obj.FieldMap_Hz*obj.TE_s(t)+2*pi*obj.metabolites(k).freq_shift_Hz *obj.TE_s(t)));
                end

            end
            sig=sig.*exp(-permute(obj.TE_s(:),[2 3 4 1])/T2star);
            sig=sig+noise_scale *complex(randn(size(sig)),randn(size(sig)));

        end

        function B0map= getRandomFmap(obj,imSize)
            % generate random low frequency maps in range -0.5 to 0.5
            %eg: fmap=getRandomFmap([100 130 1])
            [x,y,z]=meshgrid(linspace(-1,1,imSize(1)),linspace(-1,1,imSize(2)),linspace(-1,1,imSize(3)));
            Ncomp=100;
            B0map=zeros(size(x));
            meanG=(rand([3, Ncomp])-0.5)*1.3;
            varG=randi([40,70],[3, Ncomp])./1000;

            for comp=1:Ncomp
                fac=double(rand(1)>0.5)*2 -1;
                B0map=B0map+fac*exp(-1*( ((x-meanG(1,comp)).^2)./(2*varG(1,comp)) +...
                    double(imSize(2)>1)* ((y-meanG(2,comp)).^2)./(2*varG(2,comp)) +...
                    double(imSize(3)>1)*((z-meanG(3,comp)).^2)./(2*varG(3,comp))  ));
            end
            B0map=imfilter(B0map,ones(ceil(size(B0map)./10)),'replicate');
            %make it in range -0.5 to
            B0map=(B0map-min(B0map(:)))./(max(B0map(:))-min(B0map(:)))-0.5;

        end


        function [h]=poltFrame(obj,metabol_con,metabol_con_pinv,metbol_mask,fm_Hz,fm_est)



            %calculate mean and standard deviation
            orig_mean=col([[obj.metabolites.con]',[zeros(size([obj.metabolites.con]))]']');
            temp=metbol_mask;
            temp(temp==0)=nan;
            temp=cat(3,temp,temp);
            temp=reshape(temp,size(metabol_con));
            con_stat_mean=squeeze(mean(metabol_con.*temp,[1 2],'omitnan'));
            con_stat_std=squeeze(std(metabol_con.*temp,[],[1 2],'omitnan'));

            con_stat_meanp=squeeze(mean(metabol_con_pinv.*temp,[1 2],'omitnan'));
            con_stat_stdp=squeeze(std(metabol_con_pinv.*temp,[],[1 2],'omitnan'));

            figure,clf

            im_orig=metbol_mask.*reshape([obj.metabolites.con],1,1,1,[]);
            im=reshape(cat(5,im_orig,metabol_con_pinv(:,:,:,1:2:end),metabol_con(:,:,:,1:2:end)),size(metabol_con,1),size(metabol_con,2),[]);
            subplot(2,2,1), imagesc(createImMontage(im,0.5*size(metabol_con,4))), axis image, colormap jet,title('real'),colorbar

            im=reshape(cat(5,0*im_orig,metabol_con_pinv(:,:,:,2:2:end),metabol_con(:,:,:,2:2:end)),size(metabol_con,1),size(metabol_con,2),[]);
            subplot(2,2,2), imagesc(createImMontage(im,0.5*size(metabol_con,4))), axis image, colormap jet,title('imag'),colorbar

            subplot(2,2,3),imagesc([fm_Hz fm_est]),caxis([min(fm_Hz(:)),max(fm_Hz(:))]),axis image, colormap jet,title('field map(orig/est) Hz'),colorbar

            subplot(2,2,4)
            bar([ con_stat_meanp(:) con_stat_mean(:) orig_mean(:) ])
            hold on
            er = errorbar(0.78:2*length(obj.metabolites),con_stat_meanp(:),-1*con_stat_stdp(:),con_stat_stdp(:));
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            er.LineWidth=2;
            er = errorbar(1:2*length(obj.metabolites),con_stat_mean(:),-1*con_stat_std(:),con_stat_std(:));
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            er.LineWidth=2;
            % xlim([0.5,4.5])
            legend('pinv','IDEAL','Original')
            xticklabels({[obj.metabolites(1).name,'(real)'],'imag',[obj.metabolites(2).name,'(real)'],'imag',...
                [obj.metabolites(3).name,'(real)'],'imag'})
            ylabel('Concentrations (AU)'),title('Metabolite quantification @ 9.4T')
        end

        function flags =getFlags(obj,varargin)

            p=inputParser;
            p.KeepUnmatched=true;
            addParameter(p,'solver','pinv',@(x) any(validatestring(x,{'pinv','IDEAL'})));
            % only for IDEAL solver
            addParameter(p,'maxit',20,@(x) isscalar(x));
            addParameter(p,'tol',0.1,@(x) isscalar(x)); % convergence criteria in Hz
            addParameter(p,'PhaseCorr',true,@(x) islogical(x));
            %                     addParameter(p,'doMasking',true,@(x) islogical(x)); %just during resampling
            %                     addParameter(p,'Interpmode','linear',@(x) any(validatestring(x,{'linear','pchip','spline'})));
            %                     addParameter(p,'doRegistration',false,@(x)islogical(x));

            %                     addParameter(p,'is3D',(obj.reco_obj.twix.image.NPar>1),@(x)islogical(x));

            parse(p,varargin{:});
            obj.flags=p.Results;

            if(isfield(p.Unmatched,'fm'))
                obj.FieldMap_Hz=p.Unmatched.fm;
            end
            if(isfield(p.Unmatched,'mask'))
                obj.mask=p.Unmatched.mask;
            end
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