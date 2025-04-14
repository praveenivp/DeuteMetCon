function [cb,cax]=overlayplot(anat,metabol,varargin)

% metabol : cell array of filenames
%S{tra,timepoints}


%tra
if(isstruct(metabol))
 anat_vol=double(niftiread(fullfile(anat(1).folder,anat(1).name)));
else
    anat_vol=anat;
end
if(isstruct(metabol))
met_vol=[];
for i=1:length(metabol)
met_vol=cat(5,met_vol, ...
    niftiread(fullfile(metabol(i).folder,metabol(i).name)));
end
else
    met_vol=metabol;
    
end
nTimePoints=size(met_vol,5);

st=Parseinput(anat_vol,varargin{:});

% transform
anat_vol=st.transform(anat_vol);
anat_vol=double(anat_vol)./double(max(anat_vol));
met_vol=st.transform(met_vol);
imSize=size(anat_vol);

% 5 slice average to improve noise 
% met_vol=met_vol+circshift(met_vol,1,3)++circshift(met_vol,-1,3)...
%     +circshift(met_vol,2,3)+circshift(met_vol,-2,3);
%   met_vol=met_vol.*5;
% met_vol=met_vol+circshift(met_vol,1,3)++circshift(met_vol,-1,3);%...
%     +circshift(met_vol,2,3)+circshift(met_vol,-2,3);
if(st.norm)
%norm
normat=met_vol(:,:,:,1,1); %water 1st time point
normat(normat<prctile(normat(:),10))=prctile(normat(:),10);
 normat=1./(normat);
%  normat(normat>prctile(normat(:),80))=prctile(normat(:),80);
%  normat=ndCircShift(normat,[0 -5 0],[1 2 3]);
        met_vol=met_vol.*normat;
end




doMedianFilt=false;
if(doMedianFilt)
    for cMet=1:2
        for tp=1:size(met_vol,5)
            inp_vol=met_vol(:,:,st.SlcSel+(1:3),cMet,tp);

                 inp_mask=imerode(inp_vol>0,strel('sphere',5));

            fac=10/tp;
            medFilt_mask=...
                            inp_vol>prctile(inp_vol(inp_mask),100-fac);
            medFilt_mask=medFilt_mask&inp_mask>0;
             inp_vol_medfilt=imgaussfilt3(inp_vol,1);
              inp_vol_medfilt=medfilt3(inp_vol_medfilt,[51 51 3]);
            
            inp_vol(medFilt_mask)=inp_vol_medfilt(medFilt_mask);
            met_vol(:,:,st.SlcSel+(1:3),cMet,tp)=inp_vol;
       end
    end
end

%     
%remove first time points?
%   met_vol(:,:,:,:,1)=[];
% nTimePoints=nTimePoints-1;

if(isempty(st.Mask))
    st.Mask=anat_vol>0.01;
end


slcSel=st.SlcSel;
anat_vol=repmat(anat_vol(:,:,slcSel),[1 1 1 nTimePoints]);
mask2=repmat(st.Mask(:,:,slcSel),[1 1 1 nTimePoints]);
% mask1=imdilate(anat_vol(:,:,1)>0,strel('sphere',10));
% mask1=imdilate(mask1,strel('sphere',100)); %extra
% mask1=ones(size(mask1),'logical');
met_vol=met_vol(:,:,slcSel,st.MetIdx,:);%.*mask1;

if(nTimePoints>1)
anat_vol=createImMontage(squeeze(anat_vol),1);
met_vol=createImMontage(squeeze(met_vol),1);
mask2=createImMontage(squeeze(mask2),1);
end
anat_vol=double(anat_vol)./max(double(anat_vol),[],'all');

% figure,
h1=imagesc(anat_vol);colormap(gca,'gray'),caxis(st.cax_im)
axis image
hold on
%  mask2=repmat(mask(:,:,slcSel),[2 1]);
%  mask3=imerode(mask2,strel('disk',6));
% dat2=createImMontage(squeeze(CSIdata_norm(:,:,slcSel,i,:)),1);
% dat2(~mask3)=0;
% dat2=imgaussfilt(dat2,5);
% axis('image')

xticks([]),

yticks(round(linspace(0.5/nTimePoints,1-0.5/nTimePoints,nTimePoints)*size(anat_vol,1)))
yticklabels(1:nTimePoints)

set(gca,'FontSize',12)
set(gca,'FontWeight','bold')

if(isempty(st.cax))
 cax=[0 prctile(met_vol(:),st.prctile)];
else
    cax=st.cax;
end


if(isempty(st.Mask))
    mask2=anat_vol>0.01;
end


 h2=image(gca,ind2rgb(uint16((2^8)*(met_vol./max(cax(2)))),st.cmap),'AlphaData',double(mask2).*st.alpha_overlay);


 if(st.colorbar)
cb=colorbar;
 else
     cb=[];
 end




end


function st=Parseinput(im,varargin)
    
    Nslc=size(im,3);
    
    

    p=inputParser;
    p.KeepUnmatched=1;
    addParameter(p,'cax',[],@(x) isvector(x)||isempty(x));
    addParameter(p,'cax_im',[0 0.5],@(x) isvector(x)||isempty(x));
    addParameter(p,'SlcSel',floor(0.1*Nslc):ceil(0.9*Nslc),@(x) isvector(x));
    addParameter(p,'caxis_im',[0,2.5e-4],@(x) isvector(x));
    addParameter(p,'alpha_overlay',0.6,@(x) isscalar(x));
    addParameter(p,'cmap',turbo(256),@(x) ismatrix(x));
    addParameter(p,'im_horz',1,@(x) isscalar(x));
    addParameter(p,'transform',@(x)x);
    addParameter(p,'title_im','',@(x) ischar(x))
    addParameter(p,'title_im_format',{'Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1]},@(x) iscell(x))
    addParameter(p,'norm',false,@(x)islogical(x))
    addParameter(p,'colorbar',false,@(x)islogical(x))
    addParameter(p,'prctile',90,@(x)isscalar(x))
    addParameter(p,'Mask',[],@(x) ismatrix(x))
    addParameter(p,'MetIdx',1,@(x)isscalar(x))
    addParameter(p,'AllAx',{})
    
   
    
    
%     addParameter(p,'doDCF','Jackson',@(x) any(strcmp(x,{'none','Jackson','voronoi'})));
    
    parse(p,varargin{:});
    
    st=p.Results;   
end