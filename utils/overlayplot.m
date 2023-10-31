function overlayplot(anat,metabol,varargin)

% metabol : cell array of filenames
%S{tra,timepoints}
nTimePoints=size(metabol,2);

%tra
anat_vol=double(niftiread(anat{1}));
met_vol=[];
for i=1:size(metabol,2)
met_vol=cat(5,met_vol,niftiread(metabol{1,i}));
end

st=Parseinput(anat_vol,varargin{:});

% transform
anat_vol=st.transform(anat_vol);
met_vol=st.transform(met_vol);
imSize=size(anat_vol);
%norm
normat=met_vol(:,:,:,3,1); %water 1st time point
% normat=1./imgaussian(normat,5);
% normat(normat>25)=25;
%   met_vol=met_vol.*normat;


slcSel=st.SlcSel;
anat_vol=repmat(anat_vol(:,:,slcSel),[1 1 1 nTimePoints]);
mask1=imdilate(anat_vol(:,:,1)>0.1,strel('sphere',10));
met_vol=met_vol(:,:,slcSel,st.MetIdx,:).*mask1;

anat_vol=createImMontage(squeeze(anat_vol),1);
anat_vol=anat_vol./max(anat_vol,[],'all');
met_vol=createImMontage(squeeze(met_vol),1);

% figure,
h1=imagesc(anat_vol);colormap('gray'),caxis([0 1])
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
cax=[0 prctile(met_vol(:),st.prctile)];
mask2=anat_vol>0.1;
 h2=image(gca,ind2rgb(uint16((2^8)*(met_vol./max(cax(2)))),jet),'AlphaData',double(mask2).*0.4);






end


function st=Parseinput(im,varargin)
    
    Nslc=size(im,3);
    
    

    p=inputParser;
    p.KeepUnmatched=1;
    addParameter(p,'Thres',[3.1 15],@(x) isvector(x));
    addParameter(p,'SlcSel',floor(0.1*Nslc):ceil(0.9*Nslc),@(x) isvector(x));
    addParameter(p,'caxis_im',[0,2.5e-4],@(x) isvector(x));
    addParameter(p,'alpha_blobs',0.95,@(x) isscalar(x));
    addParameter(p,'cmap',hot(4096),@(x) ismatrix(x));
    addParameter(p,'im_horz',1,@(x) isscalar(x));
    addParameter(p,'transform',@(x)x);
    addParameter(p,'title_im','',@(x) ischar(x))
    addParameter(p,'title_im_format',{'Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1]},@(x) iscell(x))
    addParameter(p,'negBold',false,@(x)islogical(x))
    addParameter(p,'colorbar',false,@(x)islogical(x))
    addParameter(p,'prctile',90,@(x)isscalar(x))
    addParameter(p,'MetIdx',1,@(x)isscalar(x))
    addParameter(p,'AllAx',{})
    
   
    
    
%     addParameter(p,'doDCF','Jackson',@(x) any(strcmp(x,{'none','Jackson','voronoi'})));
    
    parse(p,varargin{:});
    
    st=p.Results;   
end