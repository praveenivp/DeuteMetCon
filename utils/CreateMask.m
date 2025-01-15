function [mask,roi] = CreateMask(MaskImage,ROIshape,vertices,varargin)
% [mask] = CreateMask(Masksize,ROIshape,vertices)
% Creates a binary mask of the given size with polygon points
% MaskImage: image to draw mask on
% ROIshape: 'circle','polygon','freehand'
%size - [m n] a vector with two element size vector
%vertices- vertices of the plygonal mask
%praveen.ivp@gmail.com
switch nargin
    case 0
        error('need atleast one parameter')
    case {1,2,4}
        if(~exist('ROIshape','var'))
        ROIshape= questdlg('Select ROI Shape', ...
            'ROI shape?', ...
            'circle','polygon','freehand','polygon');
        end
        imagesc(MaskImage)
        MaskImage=size(MaskImage);
        [x,y]=meshgrid(1:MaskImage(2),1:MaskImage(1));
        
        switch(ROIshape)
            case 'circle'
                circle=drawcircle;%('Radius',6.1,'InteractionsAllowed','translate','Center',varargin{1});
%                 pos=customWait(circle);
                mask=(((x-circle.Center(1)).^2+(y-circle.Center(2)).^2)<circle.Radius.^2);
                roi=circle.Center;
                
            case 'polygon'
                poly=drawpolygon;
                [in,on] =inpolygon(x(:),y(:),poly.Position(:,1),poly.Position(:,2));
                mask=reshape(in|on,MaskImage);
                
            case 'freehand'
                fh=drawfreehand;
                if(size(fh.Position,1)==1)
                    disp(ceil(fh.Position))
                    mask=zeros(MaskImage);
                    mask(ceil(fh.Position(2)),ceil(fh.Position(1)))=1;
                else
                    [in,on] =inpolygon(x(:),y(:),fh.Position(:,1),fh.Position(:,2));
                    mask=reshape(in|on,MaskImage);
                end
        end
        
        
    case 3
        [x,y]=meshgrid(1:MaskImage(2),1:MaskImage(1));
        % boundary = x< min(vertices(:,1))|| x>max(vertices(:,1))||y< min(vertices(:,2))|| y>max(vertices(:,2));
        [in,on] =inpolygon(x(:),y(:),vertices(:,1),vertices(:,2));
        mask=reshape(in|on,MaskImage);
end

clf,imagesc(mask);
end


function pos = customWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);

% Return the current position
pos = hROI.Position;

end


function clickCallback(~,evt)

if strcmp(evt.SelectionType,'double')
    uiresume;
end

end
