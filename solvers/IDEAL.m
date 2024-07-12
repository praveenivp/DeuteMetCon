classdef IDEAL < matlab.mixin.Copyable
    %     IDEALsolver(metabolties,TE_s,FieldMap_Hz,Mode{'phaseonly','IDEAL'},nIteration)
    properties
        metabolites % struvt array struct('freq_shift_Hz',79,'con',randi([2,10],1)/10,'name','1.3 ppm','T1_s',297e-3,'T2_s',61e-3)
        TE_s % Echo time in seconds
        transp=false; % transpose flag
        experimental% experimental outputs
        flags
        mask
        FieldMap_Hz
    end
    methods
        function obj=IDEAL(varargin)
            if (nargin>=2)
                obj.metabolites=varargin{1}; %struct
                obj.TE_s=varargin{2}; %s
                obj.getFlags(varargin{3:end})
            else
                error('Input parameter: IDEALsolver(metabolties,TE_s)');
            end
            obj.TE_s=obj.TE_s(:)';% make it row vector
        end
        function [data]=performPhaseCorr(obj,data)
            % Stupid fucntion to subtract the phase of first echo from all
            % other echo.
            if(obj.flags.PhaseCorr &&  obj.TE_s(1)>0)
                obj.TE_s= obj.TE_s-obj.TE_s(1);
                data=data.*exp(-1i*angle(data(:,:,:,1,1)));
            end
        end

        function res= mtimes(obj,inp)
            %             bb=[nFE,nIntlv,nPar,nCha]% adj case
            %             bb=[ImX,Imy,Imz,nCha] % forward case

            if(obj.flags.parfor) 
                res= mtimes_parallel(obj,inp);
                return;
            end

            res=0;
            if (obj.transp)
                A=getA(obj);
%                 Ainv=pinv(A);
                if(isempty(obj.mask))
                    obj.mask=ones(size(inp,1),size(inp,2),size(inp,3))>0;
                end

                inp_col=reshape(obj.mat2col(inp(:,:,:,:,1),obj.mask),[],length(obj.TE_s));

                res=zeros(size(inp_col,1),size(A,2));
%                 res2=zeros(size(inp_col,1),size(A,2));
                residue=zeros(size(inp_col));
                if(sum(obj.FieldMap_Hz,'all')==0)
                    fm_est= zeros(size(inp_col,1),1);
                else
                    fm_est=-1*obj.FieldMap_Hz;
                    fm_est=obj.mat2col(fm_est,obj.mask);
                end


                if(strcmpi(obj.flags.solver,'phaseonly'))
                    for i=1:size(inp_col,1)
                        inp_corr=inp_col(i,:).*exp(1i*2*pi*fm_est(i)*obj.TE_s);
                        %                             res(i,:)=Ainv*[real(inp_corr) imag(inp_corr)]';
                        res(i,:)=A\inp_corr(:);
                        residue(i,:)=single(inp_corr(:)-A*res(i,:).');
                    end

                else % IDEAL algorithm

                    pb = waitbar(0, 'estimating field map(1/2)');
                    for i=1:size(inp_col,1)
                        for it=1:obj.flags.maxit
                            A3=getA(obj);
                            Sn_cap=inp_col(i,:).*exp(1i*2*pi*fm_est(i)*obj.TE_s); % remove Known b0 offresonance
                            metabol_est=A3\Sn_cap(:); % linear prediction
                            % calc residue
                            residue_iter=Sn_cap(:)-A*metabol_est;
                            B=getB(obj,metabol_est);
                            delta_2= B\residue_iter;

                            fm_est(i)=fm_est(i)-real(delta_2(1));
                           if(abs(real(delta_2(1)))<obj.flags.tol) %figure,plot(all_B0);
                                break; end

                        end
                        waitbar(i/size(inp_col,1), pb);
                     end

                    fm_est=obj.col2mat(fm_est,obj.mask);

                    % smooth estimated fieldmap(deafult=0)
                    if(obj.flags.SmoothFM>0)
                        fm_est=imgaussfilt3(fm_est,obj.flags.SmoothFM);
                    elseif((obj.flags.SmoothFM<0))
                        fm_est=medfilt3(fm_est,[1 1 1]*obj.flags.SmoothFM*-1);
                    end


                    fm_est=obj.mat2col(fm_est,obj.mask);
                    close(pb)

                    pb = waitbar(0, 'estimating concentrations(2/2)');
                    for i=1:size(inp_col,1)
                        Sn_cap=inp_col(i,:).*exp(1i*2*pi*fm_est(i)*obj.TE_s); % remove Known b0 offresonance
                        res(i,:)=A\Sn_cap(:);%[real(Sn_cap) imag(Sn_cap)]';              
                        residue(i,:)=single(Sn_cap(:)-A*res(i,:).');
                        waitbar(i/size(inp_col,1), pb);
                    end
                    close(pb)

                    obj.experimental.fm_est=obj.col2mat(fm_est,obj.mask);
                end
                %                     res=reshape(res,size(res,1),length(obj.metabolites),2);
                %                     res=complex(res(:,:,1),res(:,:,2));
                res=obj.col2mat(res,obj.mask);
                obj.experimental.residue=obj.col2mat(single(residue),obj.mask);
            else
                error('forward operation not implemented')
            end
            obj.transp = false;
        end
                function res= mtimes_parallel(obj,inp)
            %             bb=[nFE,nIntlv,nPar,nCha]% adj case
            %             bb=[ImX,Imy,Imz,nCha] % forward case

      
            res=0;
            if (obj.transp)
                A=getA(obj);
                TE= obj.TE_s;
%                 Ainv=pinv(A);
                if(isempty(obj.mask))
                    obj.mask=ones(size(inp,1),size(inp,2),size(inp,3))>0;
                end

                inp_col=reshape(obj.mat2col(inp(:,:,:,:,1),obj.mask),[],length(obj.TE_s));

                res=zeros(size(inp_col,1),size(A,2));
%                 res2=zeros(size(inp_col,1),size(A,2));
                residue=zeros(size(inp_col));
                if(sum(obj.FieldMap_Hz,'all')==0)
                    fm_est= zeros(size(inp_col,1),1);
                else
                    fm_est=-1*obj.FieldMap_Hz;
                    fm_est=obj.mat2col(fm_est,obj.mask);
                end


                if(strcmpi(obj.flags.solver,'phaseonly'))
                    parfor i=1:size(inp_col,1)
                        inp_corr=inp_col(i,:).*exp(1i*2*pi*fm_est(i)*TE);
                        %res(i,:)=Ainv*[real(inp_corr) imag(inp_corr)]';
                        res(i,:)=A\inp_corr(:);
                        residue(i,:)=single(inp_corr(:)-A*res(i,:).');
                    end

                else % IDEAL algorithm

                   fprintf('estimating field map(1/2) \n');
               
                    maxit=obj.flags.maxit;
                    parfor i=1:size(inp_col,1)
                        for it=1:maxit
                            Sn_cap=inp_col(i,:).*exp(1i*2*pi*fm_est(i)*TE); % remove Known b0 offresonance
                            metabol_est=A\Sn_cap(:); % linear prediction
                            % calc residue
                            residue_iter=Sn_cap(:)-A*metabol_est;
                            B=getB(obj,metabol_est);
                            delta_2= B\residue_iter;

                            fm_est(i)=fm_est(i)-real(delta_2(1));
                           if(abs(real(delta_2(1)))<obj.flags.tol) %figure,plot(all_B0);
                                break; end

                        end
                     end

                    fm_est=obj.col2mat(fm_est,obj.mask);

                    % smooth estimated fieldmap(deafult=0)
                    if(obj.flags.SmoothFM>0)
                        fm_est=imgaussfilt3(fm_est,obj.flags.SmoothFM);
                    elseif((obj.flags.SmoothFM<0))
                        fm_est=medfilt3(fm_est,[1 1 1]*obj.flags.SmoothFM*-1);
                    end


                    fm_est=obj.mat2col(fm_est,obj.mask);
                    fprintf('estimating metabolities(2/2) \n');
                    parfor i=1:size(inp_col,1)
                        Sn_cap=inp_col(i,:).*exp(1i*2*pi*fm_est(i)*TE); % remove Known b0 offresonance
                        res(i,:)=A\Sn_cap(:);%[real(Sn_cap) imag(Sn_cap)]';              
                        residue(i,:)=single(Sn_cap(:)-A*res(i,:).');
                    end

                    obj.experimental.fm_est=obj.col2mat(fm_est,obj.mask);
                end
                %res=reshape(res,size(res,1),length(obj.metabolites),2);
                %res=complex(res(:,:,1),res(:,:,2));
                res=obj.col2mat(res,obj.mask);
                obj.experimental.residue=obj.col2mat(single(residue),obj.mask);
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
            A=zeros(length(obj.TE_s), length(obj.metabolites));
            for cMet=1:length(obj.metabolites)
                CD_cMet=CDij(obj.metabolites(cMet).freq_shift_Hz,obj.TE_s(:));
                A(:,cMet)=CD_cMet(:);
%                 A(:,((1:2)+2*(cMet-1)))=  [[real(CD_cMet(:)); imag(CD_cMet(:))]  ,[ imag(conj(CD_cMet(:))); real(conj(CD_cMet(:)));] ];
            end
        end
        function B=getB(obj,metabol_est)
            A=obj.getA();
            Gjn=1i*(A*metabol_est);
            B=[ 2*pi*obj.TE_s(:).*Gjn, A];
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
            legend('phaseonly','IDEAL','Original')
            xticklabels({[obj.metabolites(1).name,'(real)'],'imag',[obj.metabolites(2).name,'(real)'],'imag',...
                [obj.metabolites(3).name,'(real)'],'imag'})
            ylabel('Concentrations (AU)'),title('Metabolite quantification @ 9.4T')
        end

        function flags =getFlags(obj,varargin)

            p=inputParser;
            p.KeepUnmatched=true;
            addParameter(p,'solver','IDEAL',@(x) any(validatestring(x,{'phaseonly','IDEAL'})));
            % only for IDEAL solver
            addParameter(p,'maxit',20,@(x) isscalar(x));
            addParameter(p,'tol',0.1,@(x) isscalar(x)); % convergence criteria in Hz
            addParameter(p,'PhaseCorr',false,@(x) islogical(x));
            addParameter(p,'SmoothFM',0,@(x)isscalar(x));
            addParameter(p,'parfor',false,@(x)islogical(x));
%             addParameter(p,'doMasking',true,@(x) islogical(x)); %just during resampling
%             addParameter(p,'Interpmode','linear',@(x) any(validatestring(x,{'linear','pchip','spline'})));
%             addParameter(p,'doRegistration',false,@(x)islogical(x));
%             addParameter(p,'is3D',(obj.reco_obj.twix.image.NPar>1),@(x)islogical(x));

            parse(p,varargin{:});
            obj.flags=p.Results;

            if(isfield(p.Unmatched,'fm'))
                obj.FieldMap_Hz=p.Unmatched.fm./(2*pi)*(6.536 /42.567);
            else
                obj.FieldMap_Hz=0;
            end
            if(isfield(p.Unmatched,'mask'))
                obj.mask=p.Unmatched.mask>0;
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
            for i=1:size(col_vec,2)
                mat(mask,i)=col_vec(:,i);
            end
            mat=reshape(mat,[size(mask) sz(2:end)]);

        end
    end
end