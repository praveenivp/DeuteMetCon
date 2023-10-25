classdef LeastSquares < matlab.mixin.Copyable
    %     LeastSquares(metabolties,TE_s,FieldMap_Hz,Mode{'pinv','IDEAL'},nIteration)
    properties
        metabolites % struvt array struct('freq_shift_Hz',79,'con',randi([2,10],1)/10,'name','1.3 ppm','T1_s',297e-3,'T2_s',61e-3)
        FieldMap % Initial field map in rad/s
        DMIPara % struct('TE',[],'TR',[],'PhaseCycles',[],FlipAngle,offset,'')
        transp=false; % transpose flag
        experimental% experimental outputs
        flags
        mask
    end
    methods
        function obj=LeastSquares(varargin)

            if(nargin==0)
                obj.DMIPara=struct('TE',[2.2:2.9:15]*1e-3); %s
                obj.metabolites=[struct('freq_shift_Hz',79,'con',randi([2,10],1)/10,'name','1.3 ppm','T1_s',297e-3,'T2_s',61e-3),...
                    struct('freq_shift_Hz',226,'con',randi([5,10],1)/10,'name','3.7 ppm ','T1_s',64e-3,'T2_s',32e-3),...
                    struct('freq_shift_Hz',293,'con',randi([8,10],1)/10,'name','4.8 ppm','T1_s',320e-3,'T2_s',12e-3)];
                obj.FieldMap=0;
            elseif(nargin==2)
                obj.metabolites=varargin{1}; %struct
                obj.DMIPara=varargin{2}; %s
                obj.getFlags();
            elseif (nargin>2)
                obj.metabolites=varargin{1}; %struct
                obj.DMIPara=varargin{2}; %s
                obj.getFlags(varargin{3:end})
            else
                error('Input parameter: LeastSquares(metabolties,TE_s)');
            end


        end
        function [data]=performPhaseCorr(obj,data)
            % Stupid fucntion to subtract the phase of first echo from all
            % other echo.
            if(obj.flags.PhaseCorr &&  obj.DMIPara.TE(1)>0)
                obj.DMIPara.TE= obj.DMIPara.TE-obj.DMIPara.TE(1);
                data=data.*exp(-1i*angle(data(:,:,:,1,1)));
            end
        end

        function outp = mtimes(obj,inp)
            %             bb=[ImX x Imy xImz x Echo x PC x av]% adj case
            %             bb=[tbd] % forward case
            outp=0;
            datasz=size(inp);
            if (obj.transp)

                %get signal model
                A=getA(obj);
                % average get the same columns
                A=repmat(A,[size(inp,6), 1]);


                if(isempty(obj.mask))
                    obj.mask=ones(size(inp,1),size(inp,2),size(inp,3),'logical');
                end

                inp= obj.performPhaseCorr(inp);
                inp_col=reshape(obj.mat2col(inp,obj.mask),[],size(A,1));

                outp=zeros(size(inp_col,1),size(A,2));
                Residue= zeros([size(inp_col),size(inp,5)]);
                obj.experimental.CondNumber= cond(A);

                if(all(obj.FieldMap(:)==0) || isempty(obj.FieldMap) )
                    fm_est= zeros(size(inp_col,1),1);
                else
                    fm_est=obj.mat2col(obj.FieldMap,obj.mask);
                end
                fm_est(isnan(fm_est))=0;

                for i=1:size(inp_col,1)
                    inp_corr=double(inp_col(i,:));
                    A_B0=A.*repmat(exp(-1i*fm_est(i)*obj.DMIPara.TE(:)),[size(inp,6), 1]);
                    outp(i,:)=A_B0\inp_corr(:);
                    Residue(i,:)=inp_corr(:) -A_B0*outp(i,:).';
                end

                outp=obj.col2mat(outp,obj.mask);
                obj.experimental.residue=obj.col2mat(Residue,obj.mask);

                if(obj.flags.plotresults)
                    obj.plotresults(outp)
                end
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
        function SigModel=getA(obj)
            if(strcmp(obj.flags.SignalModel,'phase'))
                CDij=@(freq,tn) exp(-1i*2*pi*freq*tn); % Cij i->chemical shift(Hz) tn-> time(s)
                SigModel=zeros(length(obj.DMIPara.TE), length(obj.metabolites));
                for cMet=1:length(obj.metabolites)
                    SigModel(:,cMet)=CDij(obj.metabolites(cMet).freq_shift_Hz,obj.DMIPara.TE(:));
                end
            elseif (strcmp(obj.flags.SignalModel,'bSSFP-Analytical'))
                [SigModel]=bSSFP_sim_analytical(obj.metabolites,obj.DMIPara.TE,obj.DMIPara.PhaseCycles,obj.DMIPara.TR,obj.FieldMap(obj.mask),obj.DMIPara.FlipAngle);
                SigModel=bsxfun(@times,SigModel,exp(-1i*angle(SigModel(:,:,1,:,:,:,:,:,:,:))));
            else
                error('Signal model : %s not implemented yet',obj.flags.SignalModel)
            end
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

        function plotresults(obj,res)



            middle_slice=round(size(res)/2);
            figure,
            tt=tiledlayout(3,size(res,4)+1+(numel(obj.FieldMap)>1),'TileSpacing','compact','Padding','compact');
            titlestr={obj.metabolites.name,'residue','2H B0 map [Hz]'};
            for cView=1:3
                switch(cView)
                    case 1
                        temp=cat(4,abs(res(:,:,middle_slice(3),:)),...
                            sos(obj.experimental.residue(:,:,middle_slice(3),:,:),[4 5]));
                        if(numel(obj.FieldMap)>1),  temp=cat(4,temp,obj.FieldMap(:,:,middle_slice(3))/(2*pi)); end



                    case 2
                        temp=cat(4,abs(res(:,middle_slice(2),:,:)), ...
                            sos(obj.experimental.residue(:,middle_slice(2),:,:,:),[4 5]));
                        if(numel(obj.FieldMap)>1),  temp=cat(4,temp,obj.FieldMap(:,middle_slice(2),:)/(2*pi)); end
                    case 3
                        temp=cat(4,abs(res(middle_slice(1),:,:,:)), ...
                            sos(obj.experimental.residue(middle_slice(1),:,:,:,:),[4 5]));
                        if(numel(obj.FieldMap)>1),  temp=cat(4,temp,obj.FieldMap(middle_slice(1),:,:)/(2*pi)); end
                end



                temp=squeeze(temp);


                for i=1:size(temp,3)
                    nexttile()
                    imagesc(temp(:,:,i)),colorbar
                    title(titlestr{i})
                end
            end



        end
      

        function getFlags(obj,varargin)

            p=inputParser;
            p.KeepUnmatched=true;
            %                     addParameter(p,'solver','pinv',@(x) any(validatestring(x,{'pinv','IDEAL'})));
            addParameter(p,'SignalModel','phase',@(x) any(validatestring(x,{'phase','bSSFP-Analytical','bSSFP-Numerical'})))
            % only for IDEAL solver
            addParameter(p,'maxit',20,@(x) isscalar(x));
            addParameter(p,'tol',0.1,@(x) isscalar(x)); % convergence criteria in Hz
            addParameter(p,'PhaseCorr',true,@(x) islogical(x));

            addParameter(p,'CalcResidual',true,@(x) islogical(x));
            addParameter(p,'CalcCondNumber',true,@(x) islogical(x));
            addParameter(p,'plotresults',true,@(x) islogical(x));
            %                     addParameter(p,'doMasking',true,@(x) islogical(x)); %just during resampling
            %                     addParameter(p,'Interpmode','linear',@(x) any(validatestring(x,{'linear','pchip','spline'})));
            %                     addParameter(p,'doRegistration',false,@(x)islogical(x));

            %                     addParameter(p,'is3D',(obj.reco_obj.twix.image.NPar>1),@(x)islogical(x));

            parse(p,varargin{:});
            obj.flags=p.Results;

            if(isfield(p.Unmatched,'fm'))
                obj.FieldMap=p.Unmatched.fm;
            end
            if(isfield(p.Unmatched,'mask'))
                obj.mask=p.Unmatched.mask;
            else
                obj.mask=ones(size(obj.FieldMap),'logical');
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


% function sig=simSig(obj,MatSz,fm_Scale,noise_scale)
%             % generate some circles and call it pahnatom
%             obj.FieldMap= fm_Scale*getRandomFmap(obj,[MatSz(:); 1]);
%             Phantom=zeros([MatSz(:)' 1 length(obj.metabolites)]);
%             pos=[0.35 0.35; 0.65 0.65; 0.35 0.65; 0.65 0.35;  ].*MatSz(:)';
%             r2=(MatSz(1)*0.3).^2;
%             for k=1:length(obj.metabolites)
%                 for i=1:MatSz(1)
%                     for j=1:MatSz(2)
%                         Phantom(i,j,1,k)= (( (i-pos(k,2)).^2+ (j-pos(k,1)).^2 ) <r2);
%                     end
%                 end
%             end
%             obj.experimental.Phantom=Phantom;
%             initial_phase=0;
%             T2star=30e-3; %ms
%             %relaxation is ignored
%             sig=zeros(size(Phantom,1),size(Phantom,2),1,length(obj.DMIPara.TE));
%             for k=1:length(obj.metabolites)
%                 for t=1:length(obj.DMIPara.TE)
%                     sig(:,:,:,t)=sig(:,:,:,t)+ obj.metabolites(k).con.*Phantom(:,:,:,k).* ...
%                         exp(1i*(initial_phase+obj.FieldMap*obj.DMIPara.TE(t)+2*pi*obj.metabolites(k).freq_shift_Hz *obj.DMIPara.TE(t)));
%                 end
% 
%             end
%             sig=sig.*exp(-permute(obj.DMIPara.TE(:),[2 3 4 1])/T2star);
%             sig=sig+noise_scale *complex(randn(size(sig)),randn(size(sig)));
% 
%         end


% 
%   function [h]=poltFrame(obj,metabol_con,metabol_con_pinv,metbol_mask,fm_Hz,fm_est)
% 
% 
% 
%             %calculate mean and standard deviation
%             orig_mean=col([[obj.metabolites.con]',[zeros(size([obj.metabolites.con]))]']');
%             temp=metbol_mask;
%             temp(temp==0)=nan;
%             temp=cat(3,temp,temp);
%             temp=reshape(temp,size(metabol_con));
%             con_stat_mean=squeeze(mean(metabol_con.*temp,[1 2],'omitnan'));
%             con_stat_std=squeeze(std(metabol_con.*temp,[],[1 2],'omitnan'));
% 
%             con_stat_meanp=squeeze(mean(metabol_con_pinv.*temp,[1 2],'omitnan'));
%             con_stat_stdp=squeeze(std(metabol_con_pinv.*temp,[],[1 2],'omitnan'));
% 
%             figure,clf
% 
%             im_orig=metbol_mask.*reshape([obj.metabolites.con],1,1,1,[]);
%             im=reshape(cat(5,im_orig,metabol_con_pinv(:,:,:,1:2:end),metabol_con(:,:,:,1:2:end)),size(metabol_con,1),size(metabol_con,2),[]);
%             subplot(2,2,1), imagesc(createImMontage(im,0.5*size(metabol_con,4))), axis image, colormap jet,title('real'),colorbar
% 
%             im=reshape(cat(5,0*im_orig,metabol_con_pinv(:,:,:,2:2:end),metabol_con(:,:,:,2:2:end)),size(metabol_con,1),size(metabol_con,2),[]);
%             subplot(2,2,2), imagesc(createImMontage(im,0.5*size(metabol_con,4))), axis image, colormap jet,title('imag'),colorbar
% 
%             subplot(2,2,3),imagesc([fm_Hz fm_est]),caxis([min(fm_Hz(:)),max(fm_Hz(:))]),axis image, colormap jet,title('field map(orig/est) Hz'),colorbar
% 
%             subplot(2,2,4)
%             bar([ con_stat_meanp(:) con_stat_mean(:) orig_mean(:) ])
%             hold on
%             er = errorbar(0.78:2*length(obj.metabolites),con_stat_meanp(:),-1*con_stat_stdp(:),con_stat_stdp(:));
%             er.Color = [0 0 0];
%             er.LineStyle = 'none';
%             er.LineWidth=2;
%             er = errorbar(1:2*length(obj.metabolites),con_stat_mean(:),-1*con_stat_std(:),con_stat_std(:));
%             er.Color = [0 0 0];
%             er.LineStyle = 'none';
%             er.LineWidth=2;
%             % xlim([0.5,4.5])
%             legend('pinv','IDEAL','Original')
%             xticklabels({[obj.metabolites(1).name,'(real)'],'imag',[obj.metabolites(2).name,'(real)'],'imag',...
%                 [obj.metabolites(3).name,'(real)'],'imag'})
%             ylabel('Concentrations (AU)'),title('Metabolite quantification @ 9.4T')
%         end