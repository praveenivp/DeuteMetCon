function imMontage = createImMontage(imarray, varargin)
% This function creates a grid (montage) of images from imarray to display or save
% ----------------------------------------------------
% Input arguments:
% 1st arg: imarray is input array of dimensions H x W x NumImages or H x W x C x NumImages where C 
%          is the number of channels in each image
% 2nd arg: gridcols = Number of images to appear horizontally
%                     Default = 3
% 3rd arg: dsf = Down Sampling Factor (to downsample or reduce the size of the input images)
%                range [0 1]. Default = 1 (no downsampling) 
% 4th arg: border = Thickness of border between two images (in pixels)
%                   Default = 0;
% ----------------------------------------------------
% Output argument:
% imMontage = output image montage or grid  
% ----------------------------------------------------
% Examples:
% imMontage = createImMontage(imarray); % automatically select grid size and populate without any 
% downsampling with 2 pixels border between consecutive images (default)
%
% imMontage = createImMontage(imarray,5,0.25); % 5 images horizontally and every
% image downsampled to quarter with 2 pixels border (default) between
% consecutive images
% 
% imMontage = createImMontage(imarray,7,0.5,5); % 7 images horizontally and every 
% image downsampled to half with 5 pixels border between consecutive images  
%
% NOTE: Input image array can be created by: 
% imarray = cat(3,imgray1,imgray2,imgray3,imgray4); % for gray scale images
% imarray = cat(4,imcolor1, imcolor2, imcolor4, imcolor5); for color images
%
% ----------------------------------------------------
% Author: Wajahat Kazmi
% Aston University, Birmingham, UK
% email: s.kazmi3@aston.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array_dims = numel(size(imarray));
if array_dims == 3 
    imn = size(imarray,3);
elseif array_dims == 4
    imn = size(imarray,4);
else
    error('Unknown image array shape. It should be  H x W x NumImages or H x W x C x NumImages where C is the number of channels in each image')
end
if ~isempty(varargin)
    if numel(varargin)>2 
        gridcols = varargin{1};
        dsf = varargin{2};
        border = varargin{3};   
    elseif numel(varargin) == 2
        gridcols = varargin{1};
        dsf = varargin{2};
        border = 0;
    else
        gridcols = varargin{1};
        dsf = 1;
        border = 0;
    end
else
    if imn<10
        gridcols = 3;
    else
        gridcols = floor(imn/3);
    end
        dsf = 1;
        border = 0;
end
blankImgs = 0;
imd = [];
%% for imarray dims: W x H x numImages
switch array_dims
    case 3
        if dsf ~= 1
            imd = imresize(imarray, dsf);
        else
            imd = imarray;
        end
        if rem(imn,gridcols)
            blankImgs = gridcols-rem(imn,gridcols);
        end
        % add blank(white images) to the imagearray
        for dex = 1:blankImgs
            imd(:,:,end+1) = NaN(size(imd,1),size(imd,2));
        end
        gridrows = size(imd,3)/gridcols;
        % populate the image grid to create montage
        imMontage=[];
        for rdex = 1:gridrows
            temp_mon=[];
            for cdex = 1:gridcols
                temp_mon = [temp_mon 255*ones(size(imd,1),border) imd(:,:,(rdex-1)*gridcols+cdex)];
            end
            temp_mon = [temp_mon 255*ones(size(imd,1),border)];
            imMontage = [imMontage; 255*ones(border,size(temp_mon,2)); temp_mon];
        end
        imMontage = [imMontage; 255*ones(border,size(temp_mon,2))];
    
%% for imarray dims: W x H x C x numImages
    case 4
        if dsf ~= 1
            for ii = 1:imn
                imd(:,:,:,ii) = imresize(imarray(:,:,:,ii), dsf);
            end
        else
            imd = imarray;
        end
        
        % number of blank images to complete the grid
        if rem(imn,gridcols)
            blankImgs = gridcols-rem(imn,gridcols);
        end
        % add blank(white images) to the imagearray
        for dex = 1:blankImgs
            imd(:,:,:,end+1) = repmat(255*ones(size(imd,1),size(imd,2)),1,1,size(imd,3));
        end
        gridrows = size(imd,4)/gridcols;
        % populate the image grid to create montage
        imMontage=[];
        for rdex = 1:gridrows
            temp_mon=[];
            for cdex = 1:gridcols
                temp_mon = [temp_mon repmat(255*ones(size(imd,1),border),1,1,size(imd,3)) imd(:,:,:,(rdex-1)*gridcols+cdex)];
            end
            temp_mon = [temp_mon repmat(255*ones(size(imd,1),border),1,1,size(imd,3))];
            imMontage = [imMontage; repmat(255*ones(border,size(temp_mon,2)),1,1,size(imd,3)); temp_mon];
        end
            imMontage = [imMontage; repmat(255*ones(border,size(temp_mon,2)),1,1,size(imd,3))];    
end
