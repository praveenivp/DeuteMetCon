classdef MetSim < matlab.mixin.Copyable
    %     IDEALsolver(metabolties,TE_s,FieldMap_Hz,Mode{'pinv','IDEAL'},nIteration)
    properties
       
        experimental% experimental outputs
        flags
        
        %experimental parameters
        TE_s % Echo time in seconds
        seqType % {'GRE','bSSFP'}
        Noisefac
        TR_s % repetition
        PhaseCylces % radians
        FA %radians

        %phantom things
        mask
        FieldMap_Hz
        metabolites % struvt array struct('freq_shift_Hz',79,'con',randi([2,10],1)/10,'name','1.3 ppm','T1_s',297e-3,'T2_s',61e-3)
        sig% simulated sig


        
    end
    methods
        function obj=MetSim(Metabolites,TE,PhaseCyles,TR,FA,seqType,varargin)
            % obj=MetSim(Metabolites,TE,PhaseCyles,TR,FA,seqType,varargin)  
            obj.metabolites=Metabolites;
               obj.TE_s=TE;
               obj.PhaseCylces=PhaseCyles;
               obj.TR_s=TR;
               obj.FA=FA;
               obj.seqType=seqType;
               obj.getFlags(varargin{:})
        end
        function [data]=performPhaseCorr(obj,data)
            % Stupid fucntion to subtract the phase of first echo from all
            % other echo.
            if(obj.flags.PhaseCorr &&  obj.TE_s(1)>0)
                obj.TE_s= obj.TE_s-obj.TE_s(1);
                data=data.*exp(-1i*angle(data(:,:,:,1,1)));
            end
        end

       

  
        
        function getSig(obj)

            % simualte the signal
            noise_scale=obj.flags.NoiseFac;
%             amp=randi
            [obj.experimental.Phantom,obj.experimental.mask_labels,obj.experimental.id_labels]=...
                getMetabolMask(obj.flags.MatSize,length(obj.metabolites)) ;

            % generate some circles and call it pahnatom
            obj.FieldMap_Hz= zeros(size(getRandomFmap(obj,[obj.flags.MatSize; obj.flags.MatSize;  1])));
%              -0.5 to 0.5 Hz
             obj.FieldMap_Hz= 10*((getRandomFmap(obj,[obj.flags.MatSize; obj.flags.MatSize;  1])));
           
%             obj.FieldMap_Hz=ones(size(obj.FieldMap_Hz))*10.567;
            obj.sig=MetSignalModel(obj.metabolites,obj.TE_s,obj.PhaseCylces,obj.TR_s,obj.FieldMap_Hz(:),obj.FA,obj.seqType);
            
            obj.sig=permute(obj.sig,[5 2 3 4 6 1]);
            obj.sig=reshape(obj.sig,obj.flags.MatSize,obj.flags.MatSize,size(obj.sig,2),size(obj.sig,3),size(obj.sig,4),[]);
            obj.sig=obj.sig.*permute(obj.experimental.Phantom,[1 2 4 5 6 3]);
             obj.sig=sum(obj.sig,6);
            obj.sig=obj.sig+noise_scale *complex(randn(size(obj.sig)),randn(size(obj.sig)));
            
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
%             temp=cat(3,temp,temp);
            temp=reshape(temp,size(metabol_con));
            con_stat_mean=squeeze(mean(metabol_con.*temp,[1 2],'omitnan'));
            con_stat_std=squeeze(std(metabol_con.*temp,[],[1 2],'omitnan'));
            
            con_stat_meanp=squeeze(mean(metabol_con_pinv.*temp,[1 2],'omitnan'));
            con_stat_stdp=squeeze(std(metabol_con_pinv.*temp,[],[1 2],'omitnan'));
            
            figure,clf
            
            im_orig=metbol_mask.*reshape([obj.metabolites.con],1,1,1,[]);
            im=reshape(cat(3,squeeze(im_orig),metabol_con_pinv(:,:,:,1:end),metabol_con(:,:,:,1:end)),size(metabol_con,1),size(metabol_con,2),[]);
            subplot(2,2,1), imagesc(createImMontage(real(im),3)), axis image, colormap jet,title('real'),colorbar
            
           
            subplot(2,2,2), imagesc(createImMontage(imag(im),3)), axis image, colormap jet,title('imag'),colorbar
            
            subplot(2,2,3),imagesc([fm_Hz fm_est]),caxis([min(fm_Hz(:)),max(fm_Hz(:))]),axis image, colormap jet,title('field map(orig/est) Hz'),colorbar
            
            subplot(2,2,4)
            bar([ [real(con_stat_meanp(:));imag(con_stat_meanp(:))],  [real(con_stat_mean(:));imag(con_stat_mean(:))] orig_mean(:) ])
            hold on
%             er = errorbar(0.78:2:1*length(obj.metabolites),real(con_stat_meanp(:)),-1*con_stat_stdp(:),con_stat_stdp(:));
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            er.LineWidth=2;
%             er = errorbar(1:2:1*length(obj.metabolites),real(con_stat_mean(:)),-1*con_stat_std(:),con_stat_std(:));
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            er.LineWidth=2;
            % xlim([0.5,4.5])
            legend('pinv','IDEAL','Original')
            xticklabels({[obj.metabolites(1).name,'(real)'],'imag',[obj.metabolites(2).name,'(real)'],'imag',...
                [obj.metabolites(3).name,'(real)'],'imag'})
            ylabel('Concentrations (AU)'),title('Metabolite quantification @ 9.4T')
        end

        function getFlags(obj,varargin)

            p=inputParser;
            p.KeepUnmatched=true;
            addParameter(p,'model','bSSFP',@(x) any(validatestring(x,{'bSSFP','GRE'})));
            addParameter(p,'PhaseCorr',false,@(x) islogical(x));

            addParameter(p,'MatSize',16,@(x) isscalar(x));
            addParameter(p,'NoiseFac',0.05,@(x) isscalar(x));

            
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