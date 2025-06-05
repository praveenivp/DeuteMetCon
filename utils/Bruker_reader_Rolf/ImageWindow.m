classdef ImageWindow < handle
    % Window for watching 2D or 3D images
    % Possible values for varargin:
    %   'OverlayPosPath'
    %   'Mask'
    %   'MaskLevel'
    %   'MaskColor'
    %   'Parent'
    %   A program that implements this object as child of a parent figure, has to implement the following function:
    %       PlotROI
    properties
        Parent
        Image
        WindowSize
        WindowPos
        widgets
        mode
        mask
        rois
        Overlay
        paths
        
    end
    
    methods
        function obj = ImageWindow(ImageObject, WindowSize, baseWidget, varargin)
            obj.mode.Click = 0;
            obj.mode.Orientation = 0;   %0: xy, 1: xz, 2: yz
            obj.mode.Marker = [];
            obj.mode.ROIactive = 0;
            obj.mode.SelectColLim = 0;
            obj.mode.ZoomRect = 0;
            obj.mode.Cross = [];
            obj.mode.roiselect = 1;
            obj.mode.roitype = 1;  % 1: circular, 2: arbitrary, 3: markers
            obj.mode.plotmask = 0;
            obj.mode.UseRealCoords = 1;
            obj.mode.colors = [[1,0,0];[0,1,0];[0,0,1];[0.7,0.7,0];[1,0,1];[0,1,1];[0.1,0.1,0.1];[0.0,0.7,0.5]];
            obj.rois.nrois = 0;
            obj.rois.rois = []; 
            obj.mask.mask = [];
            obj.mask.maskLevel = 0.5;
            obj.mode.maskcolor = 'red';
            obj.mode.Colorbar = 0;
            obj.mode.name = inputname(1);
            obj.mode.input = 0;
            obj.mode.CurrentPoint = [1,1];
            obj.mode.ROIDataOpen = [];
            Windowpos = [0,0];
            if numel(obj.mode.name) == 0
                obj.mode.name = 'ImageWindow';
            end
            if nargin < 2 || numel(WindowSize) == 0
                WindowSize = [580,500];
            end
            if numel(WindowSize)==1
                WindowSize = [WindowSize,WindowSize];
            elseif numel(WindowSize) == 4
                Windowpos = WindowSize(1:2);
                WindowSize = WindowSize(3:4);
            end
            obj.WindowSize = [WindowSize(1),WindowSize(2)];
            obj.WindowPos = Windowpos;
            % Set the image object
            if isobject(ImageObject) 
                if isa(ImageObject,'ImageClass')
                    obj.Image = ImageObject;
                elseif isa(ImageObject, 'DataClass')
                    obj.Image = ImageClass(ImageObject);
                    obj.Image.Params(ImageObject.AcqParameters, ImageObject.RecoParameters);
                end
            else                                       % ImageObject is not an ImageObject, but a simple 3D array
                if ndims(ImageObject) > 3
                    ImageObject = squeeze(ImageObject);
                    if ndims(ImageObject) > 3
                        disp([num2str(ndims(ImageObject)),'-dimensional array - removing excessive dimensions!']);
                        ImageObject = ImageObject(:,:,:,1,1,1,1,1,1,1,1,1);
                    end
                end
                if ~isreal(ImageObject) 
                    ImageObject = abs(ImageObject);
                end
                obj.Image = ImageClass(ImageObject);
                if isobject(obj.Image) == 0 || numel(obj.Image.NDims) == 0
                    return
                end
            end
            
            % Now set the ImageParams           
            if obj.mode.UseRealCoords == 1 && ~isa(obj.Image,'ImageClass')
                obj.Image.Params();
            end
            % Interpret varargin
            if nargin > 3 
                obj.Setvararginparams(varargin);
            end
            
            % Generate the widget
            if nargin < 3 || numel(baseWidget)~=1
                obj.widgets.base = figure('Name',obj.mode.name,'NumberTitle','off','DockControls','off','Pointer','crosshair',...
                    'Toolbar','none','Menubar','none','Unit','Pixels','Position',[100,100,obj.WindowSize(1), obj.WindowSize(2)],'Color',[0.8,0.8,0.8],'Tag','rimWindow');%,...
            else
                obj.widgets.base = baseWidget;
            end
            obj.widgets.main = uipanel(obj.widgets.base,'Unit','Pixels','Position',[obj.WindowPos(1),obj.WindowPos(2),obj.WindowSize(1), obj.WindowSize(2)],'BackgroundColor',[0.8,0.8,0.8]);
            if obj.mode.Colorbar == 0
                axespos = [50,30,obj.WindowSize(1)-80, obj.WindowSize(2)-50];
            else
                axespos = [50,30,obj.WindowSize(1)-140, obj.WindowSize(2)-50];
            end
            obj.widgets.ImageAxes = axes('Parent',obj.widgets.main,'Unit','Pixels','Position',axespos,'Box','On','XTickLabelMode','manual','YTickLabelMode','manual', 'Colormap',obj.Image.ColorMap.cm);
            if obj.mode.Colorbar ~= 0
                obj.mode.Colorbar = colorbar(obj.widgets.ImageAxes,'Units','Pixels','Position',[obj.WindowSize(1)-82,30,15,obj.WindowSize(2)-60]);
            end
            if obj.Image.NDims == 3
                pos = [30,axespos(2),15,axespos(4)];
                obj.Image.CurrentSlice = round(obj.Image.MatrixSize(3)/2);
                obj.widgets.slider = uicontrol('Parent',obj.widgets.main,'Style','slider','Min',1,'Max',obj.Image.MatrixSize(3),'SliderStep',[1.0/double(obj.Image.MatrixSize(3)-1.0),1.0/double(obj.Image.MatrixSize(3)-1.)],'Value',obj.Image.CurrentSlice,'BackgroundColor',[0.8,0.8,0.8],'Unit','Pixels','Position',pos);
                obj.widgets.sliderlabel = uicontrol('Parent',obj.widgets.main, 'Style','text','Position',[1,pos(4)/2-5,25,15],'String',num2str(obj.Image.CurrentSlice),'BackgroundColor',[0.8,0.8,0.8]);
                obj.widgets.orientbuttons = uibuttongroup(obj.widgets.main,'Unit','Pixels','Position',[obj.WindowSize(1)/2-150, obj.WindowSize(2)-20,300,20],'BackgroundColor',[0.8,0.8,0.8],'BorderType','none','SelectionChangedFcn',{@obj.ImageWindowCallback,'Orientation'});
                r1 = uicontrol(obj.widgets.orientbuttons,'Style','radiobutton','Unit','Pixels','Position',[0, 0,100,20],'String','xy','BackgroundColor',[0.8,0.8,0.8]);
                r2 = uicontrol(obj.widgets.orientbuttons,'Style','radiobutton','Unit','Pixels','Position',[100, 0,100,20],'String','xz','BackgroundColor',[0.8,0.8,0.8]);
                r3 = uicontrol(obj.widgets.orientbuttons,'Style','radiobutton','Unit','Pixels','Position',[200, 0,100,20],'String','yz','BackgroundColor',[0.8,0.8,0.8]);
            else
                obj.Image.CurrentSlice = 1;
            end
            obj.widgets.rightslider = uicontrol('Parent',obj.widgets.main,'Style','slider','Min',0,'Max',1,'SliderStep',[0.1,0.2],'Value',0.5,'BackgroundColor',[0.8,0.8,0.8],'Unit','Pixels','Position',[obj.WindowSize(1)-20,axespos(2),15,axespos(4)],'Visible',0);
            %obj.widgets.rightsliderlabel = uicontrol('Parent',obj.widgets.main, 'Style','text','Position',[1,pos(4)/2-5,25,15],'String',num2str(obj.Image.CurrentSlice),'BackgroundColor',[0.8,0.8,0.8]);
            obj.widgets.Image = image('Parent',obj.widgets.ImageAxes,'CDataMapping','scaled');%,...
            obj.widgets.Overlay = image('Parent',obj.widgets.ImageAxes, 'Visible','off'); 
            obj.widgets.mask = image('Parent',obj.widgets.ImageAxes, 'Visible','off'); 
            obj.widgets.coords = uicontrol('Parent',obj.widgets.main,'Style','text','Position',[axespos(1)+axespos(3)/10-50,15,100,15],'String','  ','BackgroundColor',[0.8,0.8,0.8]);
            obj.widgets.value = uicontrol('Parent',obj.widgets.main,'Style','text','Position',[axespos(1)+axespos(3)/2-60,15,150,15],'String','  ','BackgroundColor',[0.8,0.8,0.8]);
            obj.widgets.ovalue = uicontrol('Parent',obj.widgets.main,'Style','text','Position',[axespos(1)+axespos(3)/2+100,15,150,15],'String','  ','BackgroundColor',[0.8,0.8,0.8]);
            obj.widgets.realcoords = uicontrol('Parent',obj.widgets.main,'Style','text','Position',[axespos(1)+axespos(3)/10-60,3,110,15],'String','     ','BackgroundColor',[0.8,0.8,0.8]);
            obj.widgets.inputlabel = uicontrol('Parent',obj.widgets.main,'Style','text','Position',[axespos(1)+axespos(3)-120,8,50,15],'String','Shift x','BackgroundColor',[0.8,0.8,0.8],'Visible',0);
            obj.widgets.inputfield = uicontrol('Parent',obj.widgets.main,'Style','edit','Position',[axespos(1)+axespos(3)-60,8,50,15],'String','Shift x','BackgroundColor',[0.8,0.8,0.8],'Visible',0);
            obj.widgets.inputok = uicontrol('Parent',obj.widgets.main,'Style','pushbutton','Position',[axespos(1)+axespos(3),8,20,15],'String','OK','BackgroundColor',[0.8,0.8,0.8],'Visible',0);
            set(obj.widgets.ImageAxes, 'YLim',[1,obj.Image.MatrixSize(1)],'XLim',[1,obj.Image.MatrixSize(2)]);
            cm = uicontextmenu;
            cmZoom = uimenu(cm,'Label','Zoom');
            uimenu(cmZoom,'Label','Zoom Out','Callback',{@obj.ImageWindowCallback,'ZoomOut'});
            uimenu(cmZoom,'Label','Zoom','Callback',{@obj.ImageWindowCallback,'Zoom'});
            cmModify = uimenu(cm,'Label','Modify');
            uimenu(cmModify,'Label','Shift x','Callback',{@obj.ImageWindowCallback,'ShiftX'});
            uimenu(cmModify,'Label','Shift y','Callback',{@obj.ImageWindowCallback,'ShiftY'});
            cmColors = uimenu(cm,'Label','Colors');
            uimenu(cmColors,'Label','Colormap','Callback',{@obj.ImageWindowCallback,'Colormap'});
            menuScaleGlobal = uimenu(cmColors,'Label','Scale globally','Callback',{@obj.ImageWindowCallback,'ScaleGlobal'});
            uimenu(cmColors,'Label','Adjust max','Callback',{@obj.ImageWindowCallback,'AdjustMax'})
            cmColorsRange = uimenu(cmColors,'Label','Color Range');            
            uimenu(cmColorsRange,'Label','Select min/max','Callback',{@obj.ImageWindowCallback,'ColorMinMax',menuScaleGlobal});
            uimenu(cmColorsRange,'Label','Point to min','Callback',{@obj.ImageWindowCallback,'ColorPointMin'});
            uimenu(cmColorsRange,'Label','Point to max','Callback',{@obj.ImageWindowCallback,'ColorPointMax'});
            uimenu(cmColors,'Label','Toggle colorbar','Callback',{@obj.ImageWindowCallback,'ColorBar'});
            uimenu(cmColors,'Label','Set Lowest to White','Callback',{@obj.ImageWindowCallback,'LowestWhite'});
            cmOverlay = uimenu(cm,'Label','Overlay');
            %uimenu(cm3,'Label','Import Overlay','Callback',{@obj.ImageWindowCallback,'ImportOverlay'});
            uimenu(cmOverlay,'Label','Store Overlay Position','Callback',{@obj.ImageWindowCallback,'OverlayStorePos'});
            uimenu(cmOverlay,'Label','Load Overlay Position','Callback',{@obj.ImageWindowCallback,'OverlayLoadPos'});
            obj.widgets.widgetOverlayHide = uimenu(cmOverlay,'Label','  Hide Overlay','Callback',{@obj.ImageWindowCallback,'OverlayHide'});
            cmMask = uimenu(cm,'Label','Mask');
            uimenu(cmMask,'Label','Import mask','Callback',{@obj.ImageWindowCallback,'ImportMask'});
            uimenu(cmMask,'Label','Toggle mask','Callback',{@obj.ImageWindowCallback,'ToggleMask'});
            uimenu(cmMask,'Label','mask outline','Callback',{@obj.ImageWindowCallback,'MaskOutline'});
            uimenu(cmMask,'Label','plot mask','Callback',{@obj.ImageWindowCallback,'MaskPlot'});
            cmRoi = uimenu(cm,'Label','ROI');
            obj.widgets.roitype(1) = uimenu(cmRoi,'Label','   Circular Rois','Callback',{@obj.ImageWindowCallback,'RoiCircular'});
            obj.widgets.roitype(2) = uimenu(cmRoi,'Label','   Arbitrary Rois','Callback',{@obj.ImageWindowCallback,'RoiArbitrary'});
            obj.widgets.roitype(3) = uimenu(cmRoi,'Label','   Marker','Callback',{@obj.ImageWindowCallback,'RoiMarker'});
            uimenu(cmRoi,'Label','Inverted ROI','Callback',{@obj.ImageWindowCallback,'RoiInvert'});
            uimenu(cmRoi,'Label','Store ROIs','Callback',{@obj.ImageWindowCallback,'RoiStore'});
            uimenu(cmRoi,'Label','Load ROIs','Callback',{@obj.ImageWindowCallback,'RoiLoad'});
            uimenu(cmRoi,'Label','Show ROI data','Callback',{@obj.ImageWindowCallback,'RoiData'});
            uimenu(cmRoi,'Label','Calculate SNR','Callback',{@obj.ImageWindowCallback,'CalcSNR'});
            set(obj.widgets.roitype(obj.mode.roitype), 'Label', ['+ ',strtrim(get(obj.widgets.roitype(obj.mode.roitype), 'Label'))]);
            cmFile = uimenu(cm,'Label','File');
            uimenu(cmFile,'Label','Export','Callback',{@obj.ImageWindowCallback,'FileExport'});
            uimenu(cmFile,'Label','Clipboard','Callback',{@obj.ImageWindowCallback,'FileClipboard'});
            set(obj.widgets.main,'UIContextMenu',cm);
            set(obj.widgets.Image, 'UIContextMenu',cm);
            set(obj.widgets.Overlay, 'UIContextMenu',cm);
            set(obj.widgets.mask, 'UIContextMenu',cm);
            obj.DrawImage;
            if nargin < 3  || numel(baseWidget)~=1  %These functions are only set if the figure is generated locally. Otherwise the calling figure has to handle this
                set(obj.widgets.base, 'CurrentAxes',obj.widgets.ImageAxes, 'WindowButtonMotionFcn',{@obj.ImageWindowCallback,'Motion'},'KeyPressFcn',{@obj.ImageWindowCallback,'Key'}, 'WindowButtonUpFcn',{@obj.ImageWindowCallback,'ClickUp'});
                if obj.Image.NDims == 3
                    set(obj.widgets.base, 'WindowScrollWheelFcn',{@obj.ImageWindowCallback,'ScrollWheel'});
                end
            end
            set(obj.widgets.base, 'ResizeFcn',{@obj.ImageWindowCallback,'Resize'},'DeleteFcn',{@obj.ImageWindowCallback,'Exit'});
            set(obj.widgets.Image, 'ButtonDownFcn',{@obj.ImageWindowCallback,'Click'});
            set(obj.widgets.Overlay, 'ButtonDownFcn',{@obj.ImageWindowCallback,'Click'});
            set(obj.widgets.mask, 'ButtonDownFcn',{@obj.ImageWindowCallback,'Click'});
            set(obj.widgets.main,'DeleteFcn',{@obj.Destroy});
            if obj.Image.NDims == 3
                set(obj.widgets.slider,'Callback',{@obj.ImageWindowCallback,'Slider'});
            end
            set(obj.widgets.rightslider,'Callback',{@obj.ImageWindowCallback,'RightSlider'});
            set(obj.widgets.inputfield,'Callback',{@obj.ImageWindowCallback,'InputField'});
            set(obj.widgets.inputok,'Callback',{@obj.ImageWindowCallback,'InputOK'});
            obj.mode.ZoomRect = line('visible','off');
            set(gcf,'NextPlot','new');
            % Set variable paths
            if isfield(obj.paths,'OverlayPos') == 0 || numel(obj.paths.OverlayPos) == 0
                obj.paths.OverlayPos = cd;
            end
            if numel(obj.mask.mask)>1
                obj.updateMask
%                 m = obj.mask.mask;
%                 m(m<=obj.mask.maskLevel) = 0;
%                 m = m/max(max(m));
%                 coldat = repmat(m,1,1,3);
%                 alphdat = m*1;
%                 if strcmpi(obj.mode.maskcolor,'blue') == 1
%                     coldat(:,:,1:2) = 0;
%                 elseif strcmpi(obj.mode.maskcolor,'red') == 1
%                     coldat(:,:,2:3) = 0;
%                 elseif strcmpi(obj.mode.maskcolor,'green') == 1
%                     coldat(:,:,[1,3]) = 0;
%                 end
%                 set(obj.widgets.mask, 'Visible','on','CData',coldat,'AlphaData',alphdat);
            end
        end
        
        function [obj, ret] = ImageWindowCallback(obj, src, evt, action, param)
            ret = 0;
%              if strcmp(action, 'Motion') == 0
%                  disp(['IW:', action]);
%              end
            switch action
                case 'Motion'
                    p = obj.GetCurrentPoint(0);
                    ret = round(p);
                    obj.mode.CurrentPoint = p;
                    if obj.mode.Click == 0
                        return;
                    end
                    
                    if obj.mode.Click == 1    %select roi
                        if obj.mode.roitype == 1   %circular roi
                            diameter = sqrt((obj.mode.Marker(1,1)-p(1,1))^2+(obj.mode.Marker(2,1)-p(2,1))^2);
                            center = [obj.mode.Marker(1,1)+(p(1,1)-obj.mode.Marker(1,1))/2,obj.mode.Marker(2,1)+(p(2,1)-obj.mode.Marker(2,1))/2];
                            if diameter > 70
                                nroipoints = 50;
                            elseif diameter > 50 
                                nroipoints = 40;
                            elseif diameter > 30
                                nroipoints = 30;
                            else
                                nroipoints = 20;
                            end
                            x = zeros(nroipoints+1,1);
                            y = zeros(nroipoints+1,1);
                            for cnt = 0:nroipoints
                                x(cnt+1)=center(1)-diameter/2*cos(cnt*2*pi()/nroipoints);
                                y(cnt+1)=center(2)-diameter/2*sin(cnt*2*pi()/nroipoints);
                            end
                        elseif obj.mode.roitype >= 2    %arbitrary roi or marker
                            x = get(obj.rois.rois(obj.rois.nrois).line,'XData');
                            y = get(obj.rois.rois(obj.rois.nrois).line,'YData');
                            x(end) = p(1,1);
                            y(end) = p(2,1);
                            set(obj.rois.rois(obj.rois.nrois).line,'XData',x,'YData',y);
                        end
                        set(obj.rois.rois(obj.rois.nrois).line,'XData',x,'YData',y);
                    elseif obj.mode.Click == 2     %Zooming
                        x = get(obj.mode.ZoomRect,'XData');
                        y = get(obj.mode.ZoomRect,'YData');
                        x(3:4) = p(1,1);
                        y(2:3) = p(2,1);
                        set(obj.mode.ZoomRect,'XData',x,'YData',y); 
                    end
                case 'Click'
                    [p,val] = obj.GetCurrentPoint(0);
                    ret = p;
                    if evt.Button == 1
                        if obj.mode.SelectColLim >0                            
                            cl = obj.widgets.ImageAxes.CLim;
                            if obj.mode.SelectColLim == 1
                                cl(1) = val;                                
                            else
                                cl(2) = val;
                            end
                            obj.Image.ColorMap.clim = cl;
                            obj.widgets.ImageAxes.CLim = cl;
                            obj.mode.SelectColLim = 0;
                        else
                            if strcmp(get(obj.widgets.base,'SelectionType'), 'extend')         % Shift+Click: Zoom in
                                obj.mode.Click = 2;
                                set(obj.mode.ZoomRect,'visible','on','Color','red','Parent',obj.widgets.ImageAxes,'XData',[p(1,1),p(1,1),p(1,1),p(1,1),p(1,1)],'YData',[p(2,1),p(2,1),p(2,1),p(2,1),p(2,1)]);
                            else
                                if obj.mode.roiselect == 1
                                    if obj.mode.Click == 1 && obj.mode.roitype >= 2   % Continue arbitrary ROI or marker
                                        p = obj.GetCurrentPoint;
                                        x = get(obj.rois.rois(obj.rois.nrois).line,'XData');
                                        y = get(obj.rois.rois(obj.rois.nrois).line,'YData');
                                        if strcmp(get(obj.widgets.base,'SelectionType'), 'open')  % Double click: finish arbitrary ROI or marker definition                                            
                                            if obj.mode.roitype == 2
                                                x = [x(1:end-2),x(1)];  
                                                y = [y(1:end-2),y(1)];
                                            end
                                            obj.mode.Click = 0;
                                            s = double(obj.Image.MatrixSize);
                                            set(obj.rois.rois(obj.rois.nrois).line,'XData',x,'YData',y);
                                            if numel(x) == 2    % Already the second point is double clicked: Select only the first point
                                                set(obj.rois.rois(obj.rois.nrois).line,'Marker','.');
                                                obj.mode.Click = 10;
                                            elseif numel(x) == 3  % This is a line connecting two points
                                                obj.mode.Click = 10;
                                            else
                                                if obj.mode.roitype == 2
                                                    roi = poly2mask(x,y,s(1),s(2));
                                                    obj.mode.Click = 10;               % This is to mark that a new ROI was defined. The ROI is passed on at ClickUp
                                                end
                                            end

                                        else    
                                            x = [x(1:end-1),p(1,1),p(1,1)];  % Single click: continue ROI definition
                                            y = [y(1:end-1),p(2,1),p(2,1)];
                                            set(obj.rois.rois(obj.rois.nrois).line,'XData',x,'YData',y);
                                        end
                                        return;
                                    end
                                    if strcmp(get(obj.widgets.base,'SelectionType'), 'normal') && obj.mode.plotmask == 0 % No Ctrl-button pressed - remove old ROIS
                                        for cnt = 1:obj.rois.nrois
                                            set(obj.rois.rois(cnt).line,'visible','off');
                                        end
                                        obj.mode.ROIactive = 0;
                                        obj.rois.nrois = 0;
                                        obj.rois.rois = [];
                                    end
                                    obj.mode.Marker = obj.GetCurrentPoint;
                                    %set(obj.widgets.base, 'WindowButtonUpFcn',{@obj.ImageWindowCallback,'ClickUp'});
                                    obj.rois.nrois = obj.rois.nrois+1;
                                    %cols = {'red','green','blue','magenta','cyan','yellow','black','white'};
                                    obj.rois.rois(obj.rois.nrois).line = line('Color',obj.mode.colors(mod(obj.rois.nrois-1,numel(obj.mode.colors(:,1)))+1,:),'visible','on','Parent',obj.widgets.ImageAxes, 'XData',[p(1,1),p(1,1)],'YData',[p(2,1),p(2,1)],'ButtonDownFcn',{@obj.ImageWindowCallback,'Click'}); 
                                    obj.rois.rois(obj.rois.nrois).name = [];
                                    obj.mode.Click = 1;
                                end
                            end
                        end
                    elseif evt.Button == 2 %Middle mouse button: Destroy all defined ROIS
                        for cnt = 1:obj.rois.nrois
                            set(obj.rois.rois(cnt).line,'visible','off');
                        end
                        obj.mode.ROIactive = 0;
                        obj.rois.nrois = 0;
                        obj.rois.rois = [];
                        obj.mode.Click = 0;
                    end
                case 'ClickUp'                  
                    ret = obj.GetCurrentPoint(0);
                    if obj.mode.Click == 1 || obj.mode.Click == 10
                        if obj.mode.roitype >= 2 && obj.mode.Click ~= 10  % For arbitrary ROIs, the Click up should not have any effect (except at the end)
                            return;
                        end
                        XData = get(obj.widgets.Image,'XData');
                        YData = get(obj.widgets.Image,'YData');
                        x = get(obj.rois.rois(obj.rois.nrois).line,'XData');
                        y = get(obj.rois.rois(obj.rois.nrois).line,'YData');
                        if numel(x) > 2       % A ROI is selected  
%                             if obj.mode.Click ~= 10
%                                 set(obj.rois.rois(obj.rois.nrois).line,'XData',[x,x(1)]);
%                                 set(obj.rois.rois(obj.rois.nrois).line,'YData',[y,y(1)]);
%                             end
                            x = x-XData(1)+1;
                            y = y-YData(1)+1;
                            s = double(obj.Image.MatrixSize);
                            roi = poly2mask(x,y,s(1),s(2));
                            if obj.mode.plotmask == 1
                                obj.mask.mask = roi;
                                obj.mask.mask = 1.2* obj.mask.mask*obj.mask.maskLevel;
                                obj.mode.plotmask = 0;
                                obj.mode.Click = 0;
                                delete(obj.rois.rois(obj.rois.nrois).line);
                                obj.rois.nrois = obj.rois.nrois-1;
                                if obj.rois.nrois == 0
                                    obj.rois.rois = [];
                                else
                                    obj.rois.rois = obj.rois.rois(1:obj.rois.nrois);
                                end
                                obj.updateMask(0.5);
                            else
                                %obj.rois.nrois = obj.rois.nrois+1;
                                obj.mode.ROIactive = obj.rois.nrois;
                                obj.rois.rois(obj.mode.ROIactive).data = find(roi == 1);
                                obj.rois.rois(obj.mode.ROIactive).name = ['ROI',num2str(obj.mode.ROIactive)];
                            end
                            if strcmp(get(obj.widgets.mask,'Visible'),'on')
                                mm = find(obj.mask.mask > obj.mask.maskLevel);
                                maskandroi = intersect(obj.rois.rois(obj.mode.ROIactive).data,mm);
                                obj.rois.rois(obj.mode.ROIactive).data = maskandroi;
                            end
                            sl = obj.Image.ImageData(:,:,obj.Image.CurrentSlice);
                            obj.rois.rois(obj.mode.ROIactive).size = numel(sl(obj.rois.rois(obj.mode.ROIactive).data));
                            if obj.rois.rois(obj.mode.ROIactive).size > 0
                                obj.rois.rois(obj.mode.ROIactive).mean = mean(sl(obj.rois.rois(obj.mode.ROIactive).data), 'omitnan');
                                obj.rois.rois(obj.mode.ROIactive).std = std(sl(obj.rois.rois(obj.mode.ROIactive).data), 'omitnan');
                            else
                                obj.rois.rois(obj.mode.ROIactive).mean = 0;
                                obj.rois.rois(obj.mode.ROIactive).std = 0;
                            end
                            if abs(obj.rois.rois(obj.mode.ROIactive).std) > 0.1
                                set(obj.widgets.value,'String',[num2str(obj.rois.rois(obj.mode.ROIactive).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(obj.mode.ROIactive).std,'%8.2f'),' (',num2str(obj.rois.rois(obj.mode.ROIactive).size,'%d'),')']);
                            else
                                set(obj.widgets.value,'String',[num2str(obj.rois.rois(obj.mode.ROIactive).mean,'%8.2g'),' ± ',num2str(obj.rois.rois(obj.mode.ROIactive).std,'%8.2g'),' (',num2str(obj.rois.rois(obj.mode.ROIactive).size,'%d'),')']);
                            end    
                        else    %No ROI, just a click on a point
                            x = round(ret(1,1));
                            y = round(ret(2,1));
                            set(obj.rois.rois(obj.rois.nrois).line,'XData',[x,x],'YData',[y,y]);
                            set(obj.rois.rois(obj.rois.nrois).line,'Marker','.');
                            %obj.rois.nrois = obj.rois.nrois+1;
                            obj.rois.rois(obj.rois.nrois).data = sub2ind(obj.Image.MatrixSize,[x,y]);
                            obj.rois.rois(obj.rois.nrois).mean = obj.Image.ImageData(y,x,obj.Image.CurrentSlice);
                            set(obj.widgets.value,'String',num2str(obj.rois.rois(obj.rois.nrois).mean,'%8.2f'));
                            obj.mode.ROIactive = obj.rois.nrois;
                        end
                        ret = obj.rois.rois(obj.rois.nrois).data;
                    elseif obj.mode.Click == 2
                        set(obj.mode.ZoomRect,'visible','off');
                        x = get(obj.mode.ZoomRect,'XData');
                        y = get(obj.mode.ZoomRect,'YData');
                        set(obj.widgets.ImageAxes,'XLim',[min(x),max(x)],'YLim',[min(y),max(y)]);
                    else
                        ret = 0;
                    end
                    obj.mode.Click = 0;
                case 'Key'
                    ret = [];
                    if numel(obj.Overlay)>0
                        if numel(evt.Character) == 0
                            if strcmp(evt.Key,'hyphen')  % '-' with Ctrl
                                y = get(obj.widgets.Overlay,'YData');
                                set(obj.widgets.Overlay,'YData',[y(1),y(2)-1]);
                            end
                        else
                            switch evt.Character
                                case '-'
                                    x = get(obj.widgets.Overlay,'XData');
                                    y = get(obj.widgets.Overlay,'YData');
                                    set(obj.widgets.Overlay,'XData',[x(1),x(2)-1],'YData',[y(1),y(2)-1]);  
                                    ret = [x(1),x(2)-1,y(1),y(2)-1,0];
                                    %disp([x(2)-x(1)+1, y(2)-y(1)-1, (x(2)-x(1)+1)/(y(2)-y(1)-1)]);
                                case '_'    %Shift - '-'
                                    x = get(obj.widgets.Overlay,'XData');
                                    y = get(obj.widgets.Overlay,'YData');
                                    set(obj.widgets.Overlay,'XData',[x(1),x(2)-1]);
                                    ret = [x(1),x(2)-1,y(1),y(2),0];
                                    %disp([x(2)-x(1)+1, y(2)-y(1)-1, (x(2)-x(1)+1)/(y(2)-y(1)-1)]);
                                case '+'
                                    x = get(obj.widgets.Overlay,'XData');
                                    y = get(obj.widgets.Overlay,'YData');
                                    set(obj.widgets.Overlay,'XData',[x(1),x(2)+1],'YData',[y(1),y(2)+1]);
                                    ret = [x(1),x(2)+1,y(1),y(2)+1,0];
                                    %disp([x(2)-x(1)+1, y(2)-y(1)-1, (x(2)-x(1)+1)/(y(2)-y(1)-1)]);
                                case '*'    %Shift-+
                                    x = get(obj.widgets.Overlay,'XData');
                                    y = get(obj.widgets.Overlay,'YData');
                                    set(obj.widgets.Overlay,'XData',[x(1),x(2)+1]);
                                    ret = [x(1),x(2)+1,y(1),y(2),0];
                                    %disp([x(2)-x(1)+1, y(2)-y(1)-1, (x(2)-x(1)+1)/(y(2)-y(1)-1)]);
                               case char(30)  %arrow up
                                    x = get(obj.widgets.Overlay,'XData');
                                    y = get(obj.widgets.Overlay,'YData');
                                    set(obj.widgets.Overlay,'YData',[y(1)+1,y(2)+1]);
                                    ret = [x(1),x(2),y(1)+1,y(2)+1,0];
                                case char(31)  %arrow down
                                    x = get(obj.widgets.Overlay,'XData');
                                    y = get(obj.widgets.Overlay,'YData');
                                    set(obj.widgets.Overlay,'YData',[y(1)-1,y(2)-1]);
                                    ret = [x(1),x(2),y(1)-1,y(2)-1,0];
                                case char(28)   % arrow left
                                    x = get(obj.widgets.Overlay,'XData');
                                    y = get(obj.widgets.Overlay,'YData');
                                    set(obj.widgets.Overlay,'XData',[x(1)-1,x(2)-1]); 
                                    ret = [x(1)-1,x(2)-1,y(1),y(2),0];
                                case char(29)   % arrow right or Ctrl-+
                                    x = get(obj.widgets.Overlay,'XData');
                                    y = get(obj.widgets.Overlay,'YData');
                                    if strcmp(evt.Key,'rightarrow')
                                        set(obj.widgets.Overlay,'XData',[x(1)+1,x(2)+1]);  
                                        ret = [x(1)+1,x(2)+1,y(1),y(2),0];
                                    elseif strcmp(evt.Key,'0')
                                        set(obj.widgets.Overlay,'YData',[y(1),y(2)+1]); 
                                        ret = [x(1),x(2),y(1),y(2)+1,0];
                                        %disp([x(2)-x(1)+1, y(2)-y(1)-1, (x(2)-x(1)+1)/(y(2)-y(1)-1)]);
                                    end
                                case 'r'    %rotate
                                    dat = obj.Overlay.ImageData;
                                    s = obj.Overlay.MatrixSize;
                                    dat = imrotate(dat,90);
                                    obj.Overlay.ImageData = dat;
                                    dat = obj.Overlay.TrueColorData;
                                    dat = imrotate(dat,90);
                                    obj.Overlay.TrueColorData = dat;
                                    obj.Overlay.MatrixSize(1:2) = [s(2),s(1)];
                                    obj.DrawImage;  
                                    ret = [0,0,0,0,90];
                                case 'R'    %rotate
                                    s = obj.Overlay.MatrixSize;
                                    obj.Overlay.ImageData = imrotate(obj.Overlay.ImageData,-90);
                                    obj.Overlay.TrueColorData = imrotate(obj.Overlay.TrueColorData,-90);
                                    obj.Overlay.MatrixSize(1:2) = [s(2),s(1)];
                                    obj.DrawImage; 
                                    ret = [0,0,0,0,-90];
                                case '0'    % Show dimensions of overlay
                                    x = get(obj.widgets.Overlay,'XData');
                                    y = get(obj.widgets.Overlay,'YData');
                                    fprintf('Overlay dimensions: %f, %f \nRelation: %f\n',x(2)-x(1)+1, y(2)-y(1)-1, (x(2)-x(1)+1)/(y(2)-y(1)-1));
                            end
                        end
                    end
                case 'Resize'
                    p = get(obj.widgets.base,'Position');
                    s = p(3:4);
                    obj.WindowSize(1) = s(1);
                    obj.WindowSize(2) = s(2);
                    set(obj.widgets.main, 'Position',[0,0, obj.WindowSize(1),  obj.WindowSize(2)]);
%                     set(obj.widgets.ImageAxes,'Position',[50,30,obj.WindowSize(1)-60, obj.WindowSize(2)-50]);
%                     if obj.Image.NDims == 3
%                         set(obj.widgets.slider,'Position', [30,30,15,obj.WindowSize(2)-50]);
%                     end
%                     set(obj.widgets.value,'Position', [50+(obj.WindowSize(1)-60)/2-25,5,150,15]);
%                     set(obj.widgets.coords,'Position', [30+(obj.WindowSize(1)-60)/10-50,5,100,15]);
                    obj.ImageWindowCallback(src,evt,'PanelResize');
                case 'PanelResize'
                    p = get(obj.widgets.main,'Position');
                    s = p(3:4);
                    ap = get(obj.widgets.ImageAxes,'Position');
                    if obj.mode.Colorbar == 0
                        set(obj.widgets.ImageAxes,'Position',[50,30,s(1)-80, s(2)-50]);
                    else
                        set(obj.widgets.ImageAxes,'Position',[50,30,s(1)-130, s(2)-50]);
                        set(obj.mode.Colorbar,'Units','Pixels','Position',[s(1)-72,30,15,s(2)-60]);
                    end
                    if isfield(obj.widgets,'slider')
                        set(obj.widgets.slider,'Position', [30,30,15,s(2)-50]);
                    end
                    if isfield(obj.widgets,'orientbuttons')
                        set(obj.widgets.orientbuttons,'Position', [s(1)/2-150, s(2)-20,300,20]);
                    end
                    set(obj.widgets.rightslider,'Position',[s(1)-20,30,15,s(2)-50]);
                    set(obj.widgets.inputlabel,'Position',[ap(1)+ap(3)-120,8,50,15]);
                    set(obj.widgets.inputfield,'Position',[ap(1)+ap(3)-60,8,50,15]);
                    set(obj.widgets.inputok,'Position',[ap(1)+ap(3),8,20,15]);
                    set(obj.widgets.value,'Position', [50+(s(1)-60)/2-60,15,150,15]);
                    set(obj.widgets.ovalue,'Position', [50+(s(1)-60)/2+100,15,150,15]);
                    set(obj.widgets.coords,'Position', [30+(s(1)-60)/10-50,15,100,15]);
                    set(obj.widgets.realcoords,'Position', [30+(s(1)-60)/10-50,3,150,15]);
                    obj.WindowSize = s;
                case 'ScrollWheel'
                    if obj.Image.NDims == 3
                    newslice = obj.Image.CurrentSlice - evt.VerticalScrollCount;
                    if newslice > 0 & newslice <= obj.Image.MatrixSize(3)
                        obj.Image.CurrentSlice = newslice;
                        set(obj.widgets.slider,'Value',newslice);
                        obj.DrawImage;      
                        p = obj.GetCurrentPoint(0);
                        if obj.rois.nrois > 0 && obj.mode.ROIactive > 0
                            obj.CalcRoi(obj.mode.ROIactive);
                            if obj.rois.nrois > 1
                                set(obj.widgets.value,'String',[num2str(obj.rois.rois(obj.rois.nrois).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(obj.rois.nrois).std,'%8.2f'),' (',num2str(obj.rois.rois(obj.rois.nrois).size,'%d'),')']);
                            else   % obj.rois.nrois == 1: Just a single point
                                set(obj.widgets.value,'String',[num2str(obj.rois.rois(obj.rois.nrois).mean,'%8.2f')]);
                            end
                        end
                    end
                    end
                case 'Slider' 
                    newslice = round(get(src, 'Value'));
                    obj.Image.CurrentSlice = newslice;
                     if obj.mode.input == 1 || obj.mode.input == 2 || obj.mode.input == 3      
                         val = round(obj.widgets.rightslider.Value);
                         if obj.mode.input == 1
                             obj.DrawImage([val,0,0]);
                         elseif obj.mode.input == 2
                             obj.DrawImage([0,val,0]);
                         elseif obj.mode.input == 3
                             obj.DrawImage([0,0,val]);
                         end
                     else
                        obj.DrawImage; 
                     end
                     p = obj.GetCurrentPoint(0,obj.mode.CurrentPoint(1:2));
%                      if p(1,1) < 1
%                          p(1,1) = 1;
%                          p(1,2) = 1;
%                          obj.GetCurrentPoint(0,p(:,1)');
%                      end
                    if obj.rois.nrois > 0 && obj.mode.ROIactive > 0
                        obj.CalcRoi(obj.mode.ROIactive);
                        set(obj.widgets.value,'String',[num2str(obj.rois.rois(obj.mode.ROIactive).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(obj.mode.ROIactive).std,'%8.2f'),' (',num2str(obj.rois.rois(obj.mode.ROIactive).size,'%d'),')']);
                        if obj.rois.rois(obj.mode.ROIactive).size > 1
                            set(obj.widgets.value,'String',[num2str(obj.rois.rois(obj.mode.ROIactive).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(obj.mode.ROIactive).std,'%8.2f'),' (',num2str(obj.rois.rois(obj.mode.ROIactive).size,'%d'),')']);
                        else   % obj.rois.nrois == 1: Just a single point
                            set(obj.widgets.value,'String',[num2str(obj.rois.rois(obj.mode.ROIactive).mean,'%8.2f')]);
                        end
                        if numel(obj.mode.ROIDataOpen) == 1
                            a.Source.Tag = 'ROIUpdate';
                            [obj, ret] = ImageWindowCallback(obj, src, a, 'ROIDataWidget');
                        end
                    end
                case 'RightSlider' 
                    if obj.mode.input == 1      % Shift x
                        val = round(obj.widgets.rightslider.Value);
                        obj.widgets.inputfield.String = num2str(val);
                        obj.DrawImage([val,0,0]);
                    elseif obj.mode.input == 2      % Shift y
                        val = round(obj.widgets.rightslider.Value);
                        obj.widgets.inputfield.String = num2str(val);
                        obj.DrawImage([0,val,0]);
                    elseif obj.mode.input == 3      % Shift z
                        val = round(obj.widgets.rightslider.Value);
                        obj.widgets.inputfield.String = num2str(val);
                        obj.DrawImage([0,0,val]);
                    elseif obj.mode.input == 4      % Zoom
                        val = obj.widgets.rightslider.Value;
                        xzoom = obj.Image.MatrixSize(2)/(obj.widgets.ImageAxes.XLim(2)-obj.widgets.ImageAxes.XLim(1)+1)
                        yzoom = obj.Image.MatrixSize(1)/(obj.widgets.ImageAxes.YLim(2)-obj.widgets.ImageAxes.YLim(1)+1)
                        obj.widgets.ImageAxes.XLim
                        center = [obj.widgets.ImageAxes.XLim(1) + (obj.widgets.ImageAxes.XLim(2)-obj.widgets.ImageAxes.XLim(1))/2,obj.widgets.ImageAxes.YLim(1) + (obj.widgets.ImageAxes.YLim(2)-obj.widgets.ImageAxes.YLim(1))/2];
                        if xzoom >= yzoom
                            xnum = obj.Image.MatrixSize(2)/val;
                            xlim = [center(1) - xnum/2,center(1)+xnum/2];
                            if xlim(1)<1
                                xlim = xlim + (xlim(1)+0.5);
                            elseif xlim(2) > obj.Image.MatrixSize(2)
                                xlim = xlim - (xlim(2) - obj.Image.MatrixSize(2));
                            end
                            obj.widgets.ImageAxes.XLim = xlim;
                        end
                        obj.widgets.inputfield.String = num2str(obj.widgets.rightslider.Value);
                    elseif obj.mode.input == 10      % Adjust max color
                        val = obj.widgets.rightslider.Value;
                        obj.widgets.ImageAxes.CLim(2) = val;
                        obj.widgets.inputfield.String = num2str(obj.widgets.rightslider.Value);
                    end
                    
                case 'InputField' 
                    if obj.mode.input == 1      % Shift x
                        value = round(str2double(obj.widgets.inputfield.String));
                        if numel(value) ~= 1
                            obj.widgets.inputfield.String = num2str(round(obj.widgets.rightslider.Value));
                            return;
                        else
                            if value > obj.widgets.rightslider.Max
                                value = obj.widgets.rightslider.Max;
                                obj.widgets.inputfield.String = num2str(value);
                            elseif value < obj.widgets.rightslider.Min
                                value = obj.widgets.rightslider.Min;
                                obj.widgets.inputfield.String = num2str(value);
                            end
                            obj.widgets.rightslider.Value = value;
                            obj.DrawImage([value,0,0]);
                        end
                    elseif obj.mode.input == 2      % Shift y
                        value = round(str2double(obj.widgets.inputfield.String));
                        if numel(value) ~= 1
                            obj.widgets.inputfield.String = num2str(round(obj.widgets.rightslider.Value));
                            return;
                        else
                            if value > obj.widgets.rightslider.Max
                                value = obj.widgets.rightslider.Max;
                                obj.widgets.inputfield.String = num2str(value);
                            elseif value < obj.widgets.rightslider.Min
                                value = obj.widgets.rightslider.Min;
                                obj.widgets.inputfield.String = num2str(value);
                            end
                            obj.widgets.rightslider.Value = value;
                            obj.DrawImage([0,value,0]);
                        end
                    elseif obj.mode.input == 3      % Shift z
                        value = round(str2double(obj.widgets.inputfield.String));
                        if numel(value) ~= 1
                            obj.widgets.inputfield.String = num2str(round(obj.widgets.rightslider.Value));
                            return;
                        else
                            if value > obj.widgets.rightslider.Max
                                value = obj.widgets.rightslider.Max;
                                obj.widgets.inputfield.String = num2str(value);
                            elseif value < obj.widgets.rightslider.Min
                                value = obj.widgets.rightslider.Min;
                                obj.widgets.inputfield.String = num2str(value);
                            end
                            obj.widgets.rightslider.Value = value;
                            obj.DrawImage([0,0,value]);
                        end
                    elseif obj.mode.input == 10       % Adjust max
                        value = str2double(obj.widgets.inputfield.String);
                        if value <= obj.widgets.ImageAxes.CLim(1)
                            obj.widgets.inputfield.String = num2str(obj.widgets.ImageAxes.CLim(1));
                        else
                            obj.widgets.ImageAxes.CLim(2) = value;
                            if value < obj.widgets.rightslider.Max
                                obj.widgets.rightslider.Value = value;
                            else
                                obj.widgets.rightslider.Max = value*1.1;
                                obj.widgets.rightslider.Value = value;
                            end
                        end
                    end
                case 'InputOK'
                    if obj.mode.input == 1      % Shift x
                        value = round(str2double(obj.widgets.inputfield.String));
                        obj.Image.ImageData = circshift(obj.Image.ImageData,[0,value,0]);
                        obj.DrawImage;
                    elseif obj.mode.input == 2   % Shift y
                        value = round(str2double(obj.widgets.inputfield.String));
                        obj.Image.ImageData = circshift(obj.Image.ImageData,[value,0,0]);
                        obj.DrawImage;
                    elseif obj.mode.input == 3   % Shift z
                        value = round(str2double(obj.widgets.inputfield.String));
                        obj.Image.ImageData = circshift(obj.Image.ImageData,[0,0,value]);
                        obj.DrawImage;
                    end
                    obj.widgets.rightslider.Visible=0;
                    obj.widgets.inputlabel.Visible = 0;
                    obj.widgets.inputfield.Visible = 0;
                    obj.widgets.inputok.Visible = 0;
                    obj.mode.input = 0;
                    for cnt = 1:3
                        obj.widgets.orientbuttons.Children(cnt).Enable = 'on';
                    end
                case 'ZoomOut'
                    set(obj.widgets.ImageAxes,'XLim',[1,obj.Image.MatrixSize(2)],'YLim',[1,obj.Image.MatrixSize(1)]);
                case 'Zoom'
                    obj.widgets.rightslider.Visible=1;
                    obj.widgets.inputlabel.String = 'Zoom';
                    obj.widgets.inputlabel.Visible = 1;
                    obj.widgets.inputfield.Visible = 1;
                    obj.widgets.inputok.Visible = 1;
                    obj.mode.input = 4;
                    xzoom = obj.Image.MatrixSize(2)/(obj.widgets.ImageAxes.XLim(2)-obj.widgets.ImageAxes.XLim(1)+1);
                    yzoom = obj.Image.MatrixSize(1)/(obj.widgets.ImageAxes.YLim(2)-obj.widgets.ImageAxes.YLim(1)+1);
                    zf = max([xzoom,yzoom]);
                    obj.widgets.rightslider.Value = zf;
                    obj.widgets.rightslider.Min = 1;
                    obj.widgets.rightslider.Max = max([10,zf]);
                    obj.widgets.rightslider.SliderStep = [0.1,1];
                    obj.widgets.inputfield.String = num2str(zf);
                    for cnt = 1:3
                        obj.widgets.orientbuttons.Children(cnt).Enable = 'off';
                    end                    
                case 'Orientation'   
                    if strcmp(evt.OldValue.String, 'xy')         % was xy, has now changed to some other value
                        if strcmp(evt.NewValue.String, 'xz') == 1     % changed to xz
                            obj.Image.FlipOrient(3);     
                            obj.mode.Orientation = 1;
                        elseif strcmp(evt.NewValue.String, 'yz')                     % 2: changed to yz              !!!!!!!! The yz (sagittal) orientation is lr-flipped compared to medical images, but otherwise the coordinate system wouldn't fit
                            obj.Image.FlipOrient(2);
                            obj.mode.Orientation = 2;
                        end
                    elseif strcmp(evt.OldValue.String, 'xz')     % was xz
                        if strcmp(evt.NewValue.String, 'xy')    % changed to xy
                            obj.Image.FlipOrient(-3);
                            obj.mode.Orientation = 0;
                        elseif strcmp(evt.NewValue.String, 'yz')                     % 2: changed to yz
                            obj.Image.FlipOrient(5);
                            obj.mode.Orientation = 2;
                        end
                    elseif strcmp(evt.OldValue.String, 'yz')      % was yz
                        if strcmp(evt.NewValue.String, 'xy')     % changed to xy
                            obj.Image.FlipOrient(-2);
                            obj.mode.Orientation = 0;
                        elseif strcmp(evt.NewValue.String, 'xz')                     % 1: changed to xz
                            obj.Image.FlipOrient(-5);
                            obj.mode.Orientation = 1;
                        end
                    end
                    set(obj.widgets.ImageAxes, 'YLim',[1,obj.Image.MatrixSize(1)],'XLim',[1,obj.Image.MatrixSize(2)]);
                    set(obj.widgets.slider, 'Max', obj.Image.MatrixSize(3),'Value',obj.Image.CurrentSlice,'SliderStep',[1.0/double(obj.Image.MatrixSize(3)-1.0),1.0/double(obj.Image.MatrixSize(3)-1.)]);
                    obj.DrawImage;                            
                case 'ShiftX'
                    obj.widgets.rightslider.Visible=1;
                    obj.widgets.inputlabel.String = 'Shift x';
                    obj.widgets.inputlabel.Visible = 1;
                    obj.widgets.inputfield.Visible = 1;
                    obj.widgets.inputok.Visible = 1;
                    obj.mode.input = 1;
                    obj.widgets.rightslider.Value = 0;
                    obj.widgets.rightslider.Min = -obj.Image.MatrixSize(2);
                    obj.widgets.rightslider.Max = obj.Image.MatrixSize(2);
                    obj.widgets.rightslider.SliderStep = [1/(2*obj.Image.MatrixSize(2)),10/(2*obj.Image.MatrixSize(2))];
                    obj.widgets.inputfield.String = '0';
                    for cnt = 1:3
                        obj.widgets.orientbuttons.Children(cnt).Enable = 'off';
                    end
                case 'ShiftY'
                    obj.widgets.rightslider.Visible=1;
                    obj.widgets.inputlabel.String = 'Shift y';
                    obj.widgets.inputlabel.Visible = 1;
                    obj.widgets.inputfield.Visible = 1;
                    obj.widgets.inputok.Visible = 1;
                    obj.mode.input = 2;
                    obj.widgets.rightslider.Value = 0;
                    obj.widgets.rightslider.Min = -obj.Image.MatrixSize(1);
                    obj.widgets.rightslider.Max = obj.Image.MatrixSize(1);
                    obj.widgets.rightslider.SliderStep = [1/(2*obj.Image.MatrixSize(1)),10/(2*obj.Image.MatrixSize(1))];
                    obj.widgets.inputfield.String = '0';
                    for cnt = 1:3
                        obj.widgets.orientbuttons.Children(cnt).Enable = 'off';
                    end
                case 'Colormap'
                    cmlist = {'parula','jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','lines','colorcube','prism','flag','white'};
                    [indx, tf] = listdlg('ListString',cmlist,'PromptString','Colormap','SelectionMode','single');
                    if tf == 1
                        try
                        cm = eval([cmlist{indx},'(256)']);
                        catch f
                           errordlg(['Colormap ', cmlist{indx},' does not exist!']);
                           return;
                        end                        
                        obj.Image.ColorMap.cm = cm;
                        obj.widgets.ImageAxes.Colormap = cm;
                    end
                case 'ScaleGlobal'
                    if strcmp(obj.widgets.ImageAxes.CLimMode, 'auto')
                        cl = [min(min(min(obj.Image.ImageData()))),max(max(max(obj.Image.ImageData())))];
                        if cl(1)~=0 || cl(2)~=0
                            if cl(1) < cl(2)
                                obj.Image.ColorMap.clim = cl;
                                obj.widgets.ImageAxes.CLim = cl;
                                src.Text = '+ Scale globally';
                            end
                            return
                        end
                    else
                        obj.widgets.ImageAxes.CLimMode = 'auto';
                        obj.Image.ColorMap.clim = [0,0];
                        src.Text = 'Scale globally';
                    end
                case 'AdjustMax'
                    obj.widgets.rightslider.Visible=1;
                    obj.widgets.inputlabel.String = 'max';
                    obj.widgets.inputlabel.Visible = 1;
                    obj.widgets.inputfield.Visible = 1;
                    obj.widgets.inputok.Visible = 1;
                    obj.mode.input = 10;
                    currclim = obj.widgets.ImageAxes.CLim;
                    obj.widgets.rightslider.Value = currclim(2);
                    obj.widgets.rightslider.Min = currclim(1);
                    if strcmp(obj.widgets.ImageAxes.CLimMode, 'auto')
                        obj.widgets.rightslider.Max = 1.2*currclim(2);
                    else
                        obj.widgets.rightslider.Max = 2*currclim(2);
                    end
                    obj.widgets.rightslider.SliderStep = [0.01,0.1];
                    obj.widgets.inputfield.String = num2str(currclim(2));
                case 'ColorMinMax'
                    j = inputdlg({'min','max'},'Set color limits',1,strtrim(cellstr(num2str(obj.widgets.ImageAxes.CLim'))));
                    if numel(j) == 0
                        return
                    end
                    cl = cellfun(@str2num,j);
                    if cl(1)~=0 || cl(2)~=0
                        if cl(1) < cl(2)
                            obj.Image.ColorMap.clim = cl;
                            obj.widgets.ImageAxes.CLim = cl;
                            param.Text = '+ Scale globally';
                        end
                        return
                    else
                        obj.widgets.ImageAxes.CLimMode = 'auto';
                        obj.Image.ColorMap.cm = obj.widgets.ImageAxes.CLim;
                        param.Text = 'Scale globally';
                    end
                case 'ColorPointMin'
                    obj.mode.SelectColLim = 1;
                case 'ColorPointMax'
                    obj.mode.SelectColLim = 2;
                case 'LowestWhite'
                        cm = obj.Image.ColorMap.cm;
                        cm(1,:) = [1,1,1];
                        obj.Image.ColorMap.cm = cm;
                        obj.widgets.ImageAxes.Colormap = cm;
                case 'ColorBar'
                    if isa(obj.mode.Colorbar,'numeric')==1 || isvalid(obj.mode.Colorbar) == 0
                        axespos = [50,30,obj.WindowSize(1)-130, obj.WindowSize(2)-50];                        
                        obj.mode.Colorbar = colorbar(obj.widgets.ImageAxes,'Units','Pixels','Position',[obj.WindowSize(1)-72,30,15,obj.WindowSize(2)-60]);
                    else                        
                        axespos = [50,30,obj.WindowSize(1)-80, obj.WindowSize(2)-50];
                        colorbar(obj.mode.Colorbar,'off');
                    end
                    obj.widgets.ImageAxes.Position = axespos;
                case 'OverlayStorePos'
                    Position.x = get(obj.widgets.Overlay,'XData');
                    Position.y = get(obj.widgets.Overlay,'YData');
                    (Position.y(2)-Position.y(1)+1)/(Position.x(2)-Position.x(1)+1)
                    [file,path] = uiputfile('*.rop','Set file name',[obj.paths.OverlayPos,filesep,'OverlayPos.rop']);
                    if numel(file) <= 1
                        return
                    end
                    save([path, filesep,file],'Position');
                    obj.paths.OverlayPosPath = path;
                case 'OverlayLoadPos'
                    if nargin < 5 || exist(param,'file') ~= 2
                        [file,path] = uigetfile([obj.paths.OverlayPos,filesep,'*.rop'],'Select file',[obj.paths.OverlayPos,filesep,'OverlayPos.rop']);
                        if numel(file) <= 1
                            return
                        end
                        load([path, filesep,file],'-mat','Position');
                    else
                        load(param,'-mat','Position');
                    end
                    set(obj.widgets.Overlay,'XData',Position.x);
                    set(obj.widgets.Overlay,'YData',Position.y);
                    obj.paths.OverlayPosPath = path;
                case 'OverlayHide'
                    if src.Text(1) == '+'
                        src.Text(1) = ' ';
                        obj.widgets.Overlay.Visible = 'on';
                    else
                        src.Text(1) = '+';
                        obj.widgets.Overlay.Visible = 'off';
                    end
                case 'ImportMask'
                    m = ImportData(1,'mask');
                    if numel(m) == 0
                        return
                    end
                    m = imresize(m,obj.Image.MatrixSize(1:2));
                    obj.mask.mask = m;                    
                    m(m<=obj.mask.maskLevel) = 0;
                    m = m/max(max(m));
                    coldat = repmat(m,1,1,3);
                    alphdat = m*0.8;
                    if strcmpi(obj.mode.maskcolor,'blue') == 1
                        coldat(:,:,1:2) = 0;
                    elseif strcmpi(obj.mode.maskcolor,'red') == 1
                        coldat(:,:,2:3) = 0;
                    elseif strcmpi(obj.mode.maskcolor,'green') == 1
                        coldat(:,:,[1,3]) = 0;
                    end
                    set(obj.widgets.mask, 'Visible','on','CData',coldat,'AlphaData',alphdat);
                case 'ToggleMask'
                    if numel(obj.mask.mask) > 0
                        if strcmp(get(obj.widgets.mask,'Visible'), 'on')
                            set(obj.widgets.mask,'Visible','off');
                        else
                            set(obj.widgets.mask,'Visible','on');
                        end
                    end
                case 'MaskOutline'
                    if numel(obj.mask.mask) > 0
                        m = abs(obj.widgets.mask.CData(:,:,1) - circshift(obj.widgets.mask.CData(:,:,1),[1,0]));
                        m2 = abs(obj.widgets.mask.CData(:,:,1) - circshift(obj.widgets.mask.CData(:,:,1),[0,1]));
                        m = m|m2;
                        set(obj.widgets.mask, 'Visible','on','CData',repmat(m,1,1,3),'AlphaData',m);
                    end
                case 'MaskPlot' 
                    obj.mode.plotmask = 1;
                case 'RoiCircular'
                    if obj.mode.roitype == 1;
                        return
                    end
                    set(obj.widgets.roitype(obj.mode.roitype),'Label', ['   ',extractAfter(get(obj.widgets.roitype(obj.mode.roitype),'Label'),'+ ')]); 
                    obj.mode.roitype = 1;
                    set(obj.widgets.roitype(obj.mode.roitype),'Label', ['+ ',strtrim(get(obj.widgets.roitype(obj.mode.roitype),'Label'))]);                    
                case 'RoiArbitrary'
                    if obj.mode.roitype == 2;
                        return
                    end
                    set(obj.widgets.roitype(obj.mode.roitype),'Label', ['   ',extractAfter(get(obj.widgets.roitype(obj.mode.roitype),'Label'),'+ ')]); 
                    obj.mode.roitype = 2;
                    set(obj.widgets.roitype(obj.mode.roitype),'Label', ['+ ',strtrim(get(obj.widgets.roitype(obj.mode.roitype),'Label'))]);
                case 'RoiMarker'
                    if obj.mode.roitype == 3;
                        return
                    end
                    set(obj.widgets.roitype(obj.mode.roitype),'Label', ['   ',extractAfter(get(obj.widgets.roitype(obj.mode.roitype),'Label'),'+ ')]); 
                    obj.mode.roitype = 3;
                    set(obj.widgets.roitype(obj.mode.roitype),'Label', ['+ ',strtrim(get(obj.widgets.roitype(obj.mode.roitype),'Label'))]);
                case 'RoiInvert'
                    if obj.rois.nrois == 0
                        return
                    end
                    if obj.mode.roitype == 3
                        disp('Inverting ROIS doesn''t work for Markers!');
                        return
                    end
                    obj.rois.nrois = obj.rois.nrois + 1;
                    obj.rois.rois(obj.rois.nrois).line = [];
                        
                    obj.mode.ROIactive = obj.rois.nrois;
                    obj.CalcRoi(obj.mode.ROIactive);
                    set(obj.widgets.value,'String',[num2str(obj.rois.rois(obj.mode.ROIactive).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(obj.mode.ROIactive).std,'%8.2f'),' (',num2str(obj.rois.rois(obj.mode.ROIactive).size,'%d'),')']);                                  
                    if isa(obj.Parent,'matlab.ui.Figure')
                        cb = obj.Parent.WindowButtonUpFcn;
                        cb{1}(obj.Parent,[],'PlotROI',cb{3},cb{4},obj.rois.rois(obj.mode.ROIactive).data);
                    end
                case 'RoiStore'
                    if obj.rois.nrois == 0
                        return
                    end
                    [file,path] = uiputfile('*.roi','Set file name',[obj.paths.OverlayPos,filesep,'rrois.roi']);
                    if numel(file) <= 1
                        return
                    end
                    storedrois.nrois = obj.rois.nrois;
                    storedrois.dimensions = obj.Image.MatrixSize(1:2);
                    for cnt = 1:obj.rois.nrois
                        if numel(obj.rois.rois(cnt).line) >0
                            storedrois.roi(cnt).x = obj.rois.rois(cnt).line.XData;
                            storedrois.roi(cnt).y = obj.rois.rois(cnt).line.YData;
                        else
                            storedrois.roi(cnt).data = obj.rois.rois(cnt).data;
                        end
                        if isfield(obj.rois.rois(cnt),'name')
                            storedrois.roi(cnt).name = obj.rois.rois(cnt).name;
                        end
                    end
                    save([path, filesep,file],'storedrois');
                    display([num2str(obj.rois.nrois), ' ROIs stored in ', path, filesep,file,'.']);
                case 'RoiLoad'
                    % Get filename
                    if nargin < 5 || exist(param,'file') ~= 2
                        [file,path] = uigetfile([obj.paths.OverlayPos,filesep,'*.roi'],'Select file',[obj.paths.OverlayPos,filesep,'rrois.roi']);
                        if numel(file) <= 1
                            return
                        end
                        load([path, filesep,file],'-mat','storedrois');
                    else
                        load(param,'-mat','storedrois');
                    end
                    % Check: only if the current image has the same matrixsize as the image the rois were saved from can they be imported
                    if storedrois.dimensions ~= obj.Image.MatrixSize(1:2) 
                        display 'Current image has different dimensions as the one the ROIs were created in!'
                        return
                    end
                    %Now generate the ROIS
                    noldrois = obj.rois.nrois;
                    %cols = {'red','green','blue','magenta','cyan','yellow','black','white'};
                    for cnt = 1:storedrois.nrois
                        obj.rois.rois(noldrois + cnt).line = line('Color',obj.mode.colors(mod(noldrois + cnt-1,numel(obj.mode.colors(:,1)))+1,:),'visible','on','Parent',obj.widgets.ImageAxes, 'XData',[storedrois.roi(cnt).x],'YData',[storedrois.roi(cnt).y],'ButtonDownFcn',{@obj.ImageWindowCallback,'Click'}); 
                        obj.rois.rois(noldrois + cnt).data = [];
                        obj.rois.rois(noldrois + cnt).size = [];
                        obj.rois.rois(noldrois + cnt).mean = [];
                        obj.rois.rois(noldrois + cnt).std = [];
                        if (isfield(storedrois.roi(cnt),'name'))
                            obj.rois.rois(noldrois + cnt).name = storedrois.roi(cnt).name;
                        else
                            obj.rois.rois(noldrois + cnt).name = ['ROI',num2str(noldrois + cnt)];
                        end
                        obj.rois.nrois = noldrois + cnt;
                    end
                    obj.mode.ROIactive = obj.rois.nrois;
                    obj.CalcRoi(obj.mode.ROIactive);
                    set(obj.widgets.value,'String',[num2str(obj.rois.rois(obj.mode.ROIactive).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(obj.mode.ROIactive).std,'%8.2f'),' (',num2str(obj.rois.rois(obj.mode.ROIactive).size,'%d'),')']);
                case 'RoiData'
                    if obj.rois.nrois == 0
                        return
                    end
                    if numel(obj.mode.ROIDataOpen) > 0
                        return
                    end
                    roiwidgets.base = figure('Name','ROIs','NumberTitle','off','DockControls','off',...
                        'Toolbar','none','Menubar','none','Unit','Pixels','Position',[100,100,300, 30*obj.rois.nrois + 55],'Color',[0.8,0.8,0.8],'DeleteFcn',{@obj.ImageWindowCallback,'ROIDataWidget'},'Tag','Del');
                    bg = uibuttongroup('Parent',roiwidgets.base,'Units','Pixel','Position',[5,50,290,30*obj.rois.nrois],'BackgroundColor',[0.9,0.9,0.9],'Tag','ROIS','SelectionChangedFcn',{@obj.ImageWindowCallback,'ROIDataWidget'});
                    roiwidgets.DeleteButton = uicontrol(roiwidgets.base,'Style','pushbutton', 'Unit','Pixels', 'Position',[20,10,120,30], 'String','Delete ROI','Tag','RoiDelete');
                    roiwidgets.OKButton = uicontrol(roiwidgets.base,'Style','pushbutton', 'Unit','Pixels', 'Position',[150,10,120,30], 'String','OK','Tag','OK');
                    obj.mode.ROIDataOpen = roiwidgets.base;
                    colnum = 1;
                    for cnt = 1:obj.rois.nrois
                        obj.CalcRoi(cnt);
                        text = ['                     : ',num2str(obj.rois.rois(cnt).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(cnt).std,'%8.2f'),' (',num2str(obj.rois.rois(cnt).size,'%d'),')']; 
                        if colnum > numel(obj.mode.colors(:,1))
                            colnum = 1;
                        end
                        % Output also in Matlab window
                        fprintf('ROI %d: %f ±%f (%f)\n',cnt,obj.rois.rois(cnt).mean,obj.rois.rois(cnt).std, obj.rois.rois(cnt).mean/obj.rois.rois(cnt).std);
                        roiwidgets.entry(cnt) = uicontrol(bg,'Style','radiobutton','Unit','Pixel','Position',[0,30*(obj.rois.nrois-cnt),290,30],'ForegroundColor',[0,0,0],'BackgroundColor',[0.9,0.9,0.9],'FontSize',10,'String',text,'Tag',num2str(cnt));
                        roiwidgets.names(cnt) = uicontrol(roiwidgets.base,'Style','edit','Unit','Pixel','Position',[24,50+30*(obj.rois.nrois-cnt),80,28],'ForegroundColor',obj.mode.colors(colnum,:),'BackgroundColor',roiwidgets.entry(cnt).BackgroundColor,'FontSize',10,'String',obj.rois.rois(cnt).name,'HorizontalAlignment','left','Tag',num2str(cnt));
                        colnum = colnum + 1;
                    end
                    set(roiwidgets.DeleteButton,'Callback',{@obj.ImageWindowCallback,'ROIDataWidget',roiwidgets.names},'UserData',roiwidgets);
                    set(roiwidgets.OKButton,'Callback',{@obj.ImageWindowCallback,'ROIDataWidget',roiwidgets.names});
                    roiwidgets.base.UserData = roiwidgets;
                case 'CalcSNR'
                    if obj.rois.nrois < 2
                        return
                    end
                    if numel(obj.mode.ROIDataOpen) > 0
                        %return
                    end
                    roiwidgets.base = figure('Name','ROIs','NumberTitle','off','DockControls','off',...
                        'Toolbar','none','Menubar','none','Unit','Pixels','Position',[100,100,300, 30*obj.rois.nrois + 55 + 55],'Color',[0.9,0.9,0.9],'DeleteFcn',{@obj.ImageWindowCallback,'ROIDataWidget'},'Tag','Del');
                    t = uicontrol(roiwidgets.base,'Style','text','Unit','Pixels','Position',[10,80+ 30*obj.rois.nrois,30,25],'BackgroundColor',[0.9,0.9,0.9],'String','S    N');
                    roiwidgets.sigbg = uibuttongroup('Parent',roiwidgets.base,'Units','Pixel','Position',[5,90,20,30*obj.rois.nrois],'BackgroundColor',[0.9,0.9,0.9],'Tag','SigRoi','SelectionChangedFcn',{@obj.ImageWindowCallback,'ROIDataWidget'});
                    roiwidgets.noisebg = uibuttongroup('Parent',roiwidgets.base,'Units','Pixel','Position',[25,90,290,30*obj.rois.nrois],'BackgroundColor',[0.9,0.9,0.9],'Tag','NoiseRoi','SelectionChangedFcn',{@obj.ImageWindowCallback,'ROIDataWidget'});
                    roiwidgets.OKButton = uicontrol(roiwidgets.base,'Style','pushbutton', 'Unit','Pixels', 'Position',[50,10,220,30], 'String','OK','Tag','OK');
                    obj.mode.ROIDataOpen = roiwidgets.base;
                    colnum = 1;
                    for cnt = 1:obj.rois.nrois
                        obj.CalcRoi(cnt);
                        text = ['                     : ',num2str(obj.rois.rois(cnt).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(cnt).std,'%8.2f'),' (',num2str(obj.rois.rois(cnt).size,'%d'),')']; 
                        if colnum > numel(obj.mode.colors(:,1))
                            colnum = 1;
                        end
                        roiwidgets.sigbutton(cnt) = uicontrol(roiwidgets.sigbg,'Style','radiobutton','Unit','Pixel','Position',[0,30*(obj.rois.nrois-cnt),20,30],'ForegroundColor',[0,0,0],'BackgroundColor',[0.9,0.9,0.9],'FontSize',10,'String','','Tag',['s',num2str(cnt)]);
                        roiwidgets.entry(cnt) = uicontrol(roiwidgets.noisebg,'Style','radiobutton','Unit','Pixel','Position',[0,30*(obj.rois.nrois-cnt),290,30],'ForegroundColor',[0,0,0],'BackgroundColor',[0.9,0.9,0.9],'FontSize',10,'String',text,'Tag',['n',num2str(cnt)]);
                        roiwidgets.names(cnt) = uicontrol(roiwidgets.base,'Style','edit','Unit','Pixel','Position',[44,90+30*(obj.rois.nrois-cnt),80,28],'ForegroundColor',obj.mode.colors(colnum,:),'BackgroundColor',roiwidgets.entry(cnt).BackgroundColor,'FontSize',10,'String',obj.rois.rois(cnt).name,'HorizontalAlignment','left','Tag',num2str(cnt));
                        colnum = colnum + 1;
                    end
                    roiwidgets.sigbutton(1).Value = 1;
                    roiwidgets.entry(2).Value = 1;
                    text = ['SNR (mean/mean) =    ',num2str(obj.rois.rois(1).mean/obj.rois.rois(2).mean)];
                    roiwidgets.SNR1 = uicontrol(roiwidgets.base,'Style','text','Unit','Pixels','Position',[10,65,290,25],'BackgroundColor',[0.9,0.9,0.9],'String',text,'FontSize',10,'HorizontalAlignment','left');
                    text = ['SNR (mean/std) =       ',num2str(obj.rois.rois(1).mean/obj.rois.rois(2).std)];
                    roiwidgets.SNR2 = uicontrol(roiwidgets.base,'Style','text','Unit','Pixels','Position',[10,40,290,25],'BackgroundColor',[0.9,0.9,0.9],'String',text,'FontSize',10,'HorizontalAlignment','left');                    
                    set(roiwidgets.OKButton,'Callback',{@obj.ImageWindowCallback,'ROIDataWidget',roiwidgets.names});
                    roiwidgets.base.UserData = roiwidgets;
                case 'ROIDataWidget' 
                    if strcmp(evt.Source.Tag, 'OK')
                        for cnt = 1:obj.rois.nrois
                            if numel(strtrim(param(cnt).String)) > 0
                                obj.rois.rois(cnt).name = strtrim(param(cnt).String);
                            end
                        end
                        delete(evt.Source.Parent);
                        obj.mode.ROIDataOpen = [];
                    elseif strcmp(evt.Source.Tag, 'Del')
                        
                    elseif strcmp(evt.Source.Tag, 'ROIS')
                        ind = str2num(evt.NewValue.Tag);
                        obj.mode.ROIactive = ind;
                        if isa(obj.Parent,'matlab.ui.Figure')
                            cb = obj.Parent.WindowButtonUpFcn;
                            cb{1}(obj.Parent,[],'PlotROI',cb{3},cb{4},obj.rois.rois(obj.mode.ROIactive).data);
                        end
                    elseif strcmp(evt.Source.Tag, 'RoiDelete')
                        roiwidget = evt.Source.UserData;
                        
                        delete(obj.rois.rois(obj.mode.ROIactive).line);
                        if obj.mode.ROIactive == 1
                            obj.rois.rois = obj.rois.rois(2:end);
                            
                        elseif obj.mode.ROIactive == obj.rois.nrois
                            obj.rois.rois = obj.rois.rois(1:end-1);
                            obj.mode.ROIactive = obj.rois.nrois-1;
                        else
                            obj.rois.rois = obj.rois.rois([1:obj.mode.ROIactive-1,obj.mode.ROIactive+1:end]);
                        end
                        obj.rois.nrois = obj.rois.nrois-1;   
                        if obj.rois.nrois >=1
                        for cnt = 1:obj.rois.nrois
                            obj.rois.rois(cnt).line.Color = obj.mode.colors(mod(cnt-1,numel(obj.mode.colors(:,1)))+1,:);
                            set(roiwidget.entry(cnt),'String',['                     : ',num2str(obj.rois.rois(cnt).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(cnt).std,'%8.2f'),' (',num2str(obj.rois.rois(cnt).size,'%d'),')'],'ForegroundColor',obj.mode.colors(mod(cnt-1,numel(obj.mode.colors(:,1)))+1,:));
                            set(roiwidget.names(cnt),'String', obj.rois.rois(cnt).name, 'ForegroundColor',obj.mode.colors(mod(cnt-1,numel(obj.mode.colors(:,1)))+1,:));
                        end
                        delete(roiwidget.entry(obj.rois.nrois+1));
                        delete(roiwidget.names(obj.rois.nrois+1));                        
                        set(roiwidget.entry(obj.mode.ROIactive),'Value',1);
                        if isa(obj.Parent,'matlab.ui.Figure')
                            cb = obj.Parent.WindowButtonUpFcn;
                            cb{1}(obj.Parent,[],'PlotROI',cb{3},cb{4},obj.rois.rois(obj.mode.ROIactive).data);
                        end
                        else
                           delete(evt.Source.Parent);
                           obj.mode.ROIDataOpen = [];
                        end
                    elseif strcmp(evt.Source.Tag, 'SigRoi') || strcmp(evt.Source.Tag, 'NoiseRoi')
                        roiwidgets = evt.Source.Parent.UserData;
                        sigroi = sscanf(roiwidgets.sigbg.SelectedObject.Tag,'s%d');
                        noiseroi = sscanf(roiwidgets.noisebg.SelectedObject.Tag,'n%d');
                        text = ['SNR (mean/mean) =    ',num2str(obj.rois.rois(sigroi).mean/obj.rois.rois(noiseroi).mean)];
                        roiwidgets.SNR1.String = text;
                        text = ['SNR (mean/std) =        ',num2str(obj.rois.rois(sigroi).mean/obj.rois.rois(noiseroi).std)];
                        roiwidgets.SNR2.String = text;
                    elseif strcmp(evt.Source.Tag, 'ROIUpdate')
                        roiwidgets = obj.mode.ROIDataOpen.UserData;
                        for cnt = 1:obj.rois.nrois
                            obj.CalcRoi(cnt);
                            text = ['                     : ',num2str(obj.rois.rois(cnt).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(cnt).std,'%8.2f'),' (',num2str(obj.rois.rois(cnt).size,'%d'),')']; 
                            roiwidgets.entry(cnt).String = text;
                            if isfield(roiwidgets,'SNR1')
                                sigroi = sscanf(roiwidgets.sigbg.SelectedObject.Tag,'s%d');
                                noiseroi = sscanf(roiwidgets.noisebg.SelectedObject.Tag,'n%d');
                                text = ['SNR (mean/mean) =    ',num2str(obj.rois.rois(sigroi).mean/obj.rois.rois(noiseroi).mean)];
                                roiwidgets.SNR1.String = text;
                                text = ['SNR (mean/std) =        ',num2str(obj.rois.rois(sigroi).mean/obj.rois.rois(noiseroi).std)];
                                roiwidgets.SNR2.String = text;
                            end
                        end
                    end
                case 'FileClipboard' 
                    copygraphics(obj.widgets.ImageAxes);
                case 'FileExport'
                    ExportImage(obj.widgets.ImageAxes)
                case 'Exit'
                    if numel(obj.mode.ROIDataOpen) > 0
                        delete(obj.mode.ROIDataOpen);
                    end
            end
            set(gcf,'NextPlot','new');
        end
        
        function updateMask(obj,level)
            if nargin == 1
                level = 1;
            end
            m = obj.mask.mask;
            m(m<=obj.mask.maskLevel) = 0;
            m = m/max(max(m))*level;
            coldat = repmat(m,1,1,3);
            alphdat = m*1;
            if strcmpi(obj.mode.maskcolor,'blue') == 1
                coldat(:,:,1:2) = 0;
            elseif strcmpi(obj.mode.maskcolor,'red') == 1
                coldat(:,:,2:3) = 0;
            elseif strcmpi(obj.mode.maskcolor,'green') == 1
                coldat(:,:,[1,3]) = 0;
            end
            set(obj.widgets.mask, 'Visible','on','CData',coldat,'AlphaData',alphdat);
        end
        
        function obj = setOverlay(obj, ImageData, TransMax, limit)
            if nargin < 4
                limit = -1;
            end
            if nargin < 3
                TransMax = 0;
            end
            if isobject(ImageData) & strcmp(class(ImageData),'ImageClass')
                obj.Overlay = ImageData;
            else
                obj.Overlay = ImageClass(ImageData);
                if isobject(obj.Overlay) == 0
                    return
                end
                if obj.Overlay.NDims == 3
                    if obj.Overlay.MatrixSize(3) > 1 && obj.Overlay.MatrixSize(3) ~= obj.Image.MatrixSize(3)
                        display('Overlay has to have either one slice or as many as the main image!');
                    end
                end
            end
            if numel(obj.Overlay.TrueColorData) == 0
                obj.Overlay.TrueColorData = zeros([obj.Overlay.MatrixSize,4]);
                ind = find(obj.Overlay.ImageData<0);
                if limit == -1
                    pllimit = max(max(max(obj.Overlay.ImageData)));
                    milimit = min(min(min(obj.Overlay.ImageData)));
                    limit = max(abs([pllimit,milimit]));
                end
                obj.Overlay.TrueColorData(:,:,:,1) = obj.Overlay.ImageData/limit;   % all colors are red, scaled from 0 to 1
                if numel(ind)>0
                    obj.Overlay.TrueColorData = reshape(obj.Overlay.TrueColorData,[],4);
                    obj.Overlay.TrueColorData(ind,1) = 0;
                    obj.Overlay.TrueColorData(ind,3) = obj.Overlay.ImageData(ind)/(-limit);
                    obj.Overlay.TrueColorData = reshape(obj.Overlay.TrueColorData,[obj.Overlay.MatrixSize,4]);
                end
                if TransMax <= 0
                    if limit ~= -1
                        obj.Overlay.TrueColorData(:,:,:,4) = abs(obj.Overlay.ImageData/limit);   % Transparency (alpha) channel
                    else
                        obj.Overlay.TrueColorData(:,:,:,4) = abs(obj.Overlay.ImageData)/max(max(max(obj.Overlay.ImageData)));   % Transparency (alpha) channel
                    end
                else
                    obj.Overlay.TrueColorData(:,:,:,4) = abs(obj.Overlay.ImageData/TransMax);
                end
            end
            if (strcmp(obj.widgets.Overlay.Visible,'off'))
                obj.widgets.obj.widgets.widgetOverlayHide.Text(1) = ' ';
            end
            set(obj.widgets.Overlay,'Visible','on');
            set(obj.widgets.Overlay,'XData',[1,obj.Image.MatrixSize(2)],'YData',[1,obj.Image.MatrixSize(1)])
            obj.DrawImage;
        end
        
        function [p, val] = GetCurrentPoint(obj, inside, CurrPoint)  %Inside is 0 or one. If one, then -1 will be returned if the cursor is outside the image.
                                                                     % If CurrPoint is given as 2x1 array: uses not the current point, but the point given as in get(obj.widgets.ImageAxes,'CurrentPoint');
            if nargin < 2                                                                   
                inside = 1;
            end
            val = [];
            XData = get(obj.widgets.Image,'XData');
            YData = get(obj.widgets.Image,'YData');
            if nargin < 3 || numel(CurrPoint) ~= 2
                p = get(obj.widgets.ImageAxes,'CurrentPoint');
                p = squeeze(p(1,1:2));     
            else
                p = CurrPoint;                                       % Valid CurrPoint given
            end
            pp(1) = p(1)-XData(1)+1;
            pp(2) = p(2)-YData(1)+1;
            if pp(2)<YData(end) && pp(2) > 0.5 && pp(1) < XData(end) && pp(1) > 0.5
                set(obj.widgets.coords,'String',[num2str(round(p(1)),4),',' num2str(round(p(2)),4),',',num2str(obj.Image.CurrentSlice)]);
                if obj.mode.UseRealCoords == 1
                    rp = obj.Image.getrealcoords([round(p(1)),round(p(2)),obj.Image.CurrentSlice]);
                    %set(obj.widgets.realcoords,'String',[num2str(round(rp(1)),4),',' num2str(round(rp(2)),4),',',num2str(round(rp(3)),4),' ',obj.Image.Coords.unit]);
                    if obj.mode.Orientation == 0
                        set(obj.widgets.realcoords,'String',[num2str(rp(1),'%.1f'),',' num2str(rp(2),'%.1f'),',',num2str(rp(3),'%.1f'),' ',obj.Image.Coords.unit]);
                    elseif obj.mode.Orientation == 1
                        set(obj.widgets.realcoords,'String',[num2str(rp(1),'%.1f'),',' num2str(rp(3),'%.1f'),',',num2str(rp(2),'%.1f'),' ',obj.Image.Coords.unit]);
                    elseif obj.mode.Orientation == 2
                        set(obj.widgets.realcoords,'String',[num2str(rp(3),'%.1f'),',' num2str(rp(1),'%.1f'),',',num2str(rp(2),'%.1f'),' ',obj.Image.Coords.unit]);
                    end
                end
                val = obj.Image.ImageData(round(p(2)),round(p(1)),obj.Image.CurrentSlice);
                if obj.mode.ROIactive == 0
                    set(obj.widgets.value,'String',num2str(val));
                end
                p = [p',pp'];
                set(obj.widgets.base,'Pointer','crosshair');
            else
                if inside == 1
                    p = -1;
                else
                    p = [p',pp'];
                end
                set(obj.widgets.base,'Pointer','arrow');
            end
            if numel(obj.Overlay) > 0
                x = get(obj.widgets.Overlay,'XData');
                y = get(obj.widgets.Overlay,'YData');
                if pp(1)>x(1) && pp(1)<x(2) && pp(2)>y(1) && pp(2)<y(2)
                    oo(1) = p(1) - x(1);
                    oo(2) = p(2) - y(1);
                    oo = round(double(oo)./(double([x(2)-x(1),y(2)-y(1)])./double(obj.Overlay.MatrixSize(2:-1:1)))); %why is x and y exchanged?
                    if min(oo)>0
                        p = [p,oo'];
                        set(obj.widgets.ovalue,'String',num2str(obj.Overlay.ImageData(oo(2),oo(1))));
                    else
                        set(obj.widgets.ovalue,'String',' ');
                    end
                else
                    set(obj.widgets.ovalue,'String',' ');
                end
            end
            return
        end
        
        function Setvararginparams(obj, argcell)
            for cnt = 1:numel(argcell)
                if ischar(argcell{cnt})
                    if strcmpi(argcell{cnt},'OverlayPosPath')
                        if numel(argcell) > cnt
                            obj.paths.OverlayPos = argcell{cnt+1};
                        end
                    elseif strcmpi(argcell{cnt},'Mask')
                        if numel(argcell) > cnt
                            m = argcell{cnt+1};
                            ms = size(m);
                            if min(ms == obj.Image.MatrixSize(1:2)) >0
                                obj.mask.mask = m;
                            end
                        end
                    elseif strcmpi(argcell{cnt},'MaskLevel')
                        if numel(argcell) > cnt
                            obj.mask.maskLevel = argcell{cnt+1};
                        end
                    elseif strcmpi(argcell{cnt},'MaskColor')
                        if numel(argcell) > cnt
                            obj.mode.maskcolor = argcell{cnt+1};
                        end
                    elseif strcmpi(argcell{cnt},'Parent')
                        if numel(argcell) > cnt
                            obj.Parent = argcell{cnt+1};
                        end
                    elseif strcmpi(argcell{cnt},'Title')
                        if numel(argcell) > cnt
                            obj.mode.name = argcell{cnt+1};
                        end
                    end
                end
            end
        end
        
        function CalcRoi(obj, n, slice)
            % Calculates mean value and standard deviation over a ROI and outputs it.
            % Takes into account the mask (if a mask is shown, only the voxels inside the mask are taken into account
            if n == 0
                return;
            end
            if nargin <3
                slice = obj.Image.CurrentSlice;
            elseif slice < 0
                slice = obj.Image.CurrentSlice;
            end
            if numel(obj.rois.rois(n).line) > 0
                XData = get(obj.widgets.Image,'XData');
                YData = get(obj.widgets.Image,'YData');
                x = get(obj.rois.rois(n).line,'XData');
                y = get(obj.rois.rois(n).line,'YData');
                x = x-XData(1)+1;
                y = y-YData(1)+1;
                if numel(x) > 2       % A ROI is selected  
                    s = double(obj.Image.MatrixSize);
                    roi = poly2mask(x,y,s(1),s(2));
                    obj.rois.rois(n).data = find(roi == 1);
                    if strcmp(get(obj.widgets.mask,'Visible'),'on')
                        mm = find(obj.mask.mask > obj.mask.maskLevel);
                        maskandroi = intersect(obj.rois.rois(n).data,mm);
                        obj.rois.rois(n).data = maskandroi;
                    end
                    sl = obj.Image.ImageData(:,:,slice);
                    obj.rois.rois(n).size = numel(sl(obj.rois.rois(n).data));
                    if obj.rois.rois(n).size > 0
                        obj.rois.rois(n).mean = mean(sl(obj.rois.rois(n).data), 'omitnan');
                        obj.rois.rois(n).std = std(sl(obj.rois.rois(n).data), 'omitnan');
                    else
                        obj.rois.rois(n).mean = 0;
                        obj.rois.rois(n).std = 0;
                    end
                else                   %A single point
                    obj.rois.rois(n).mean = obj.Image.ImageData(y(1),x(1),slice);
                    obj.rois.rois(n).std = 0;
                    obj.rois.rois(n).size = 1;
                end
            else                      % obj.rois.rois(n).line ist leer (inverted ROI) -> neuberechnung des ROIs
                % Create a new ROI, which is one everywhere where there is no other ROI
                temproi = ones(obj.Image.MatrixSize(1:2));
                for cnt = 1:obj.rois.nrois
                    if numel(obj.rois.rois(cnt).line) > 0    % Only, if the ROI is plotted. Especially not, if the ROI is already an inverted one
                        temproi(obj.rois.rois(cnt).data) = 0;
                    end
                end
                temproidat = find(temproi == 1);
                if numel(temproidat) > 0
                    obj.rois.rois(n).data = temproidat;
                    if strcmp(get(obj.widgets.mask,'Visible'),'on')
                        mm = find(obj.mask.mask > obj.mask.maskLevel);
                        maskandroi = intersect(obj.rois.rois(n).data,mm);
                        obj.rois.rois(n).data = maskandroi;
                    end
                    sl = obj.Image.ImageData(:,:,slice);
                    obj.rois.rois(n).size = numel(sl(obj.rois.rois(n).data));
                    if obj.rois.rois(n).size > 0
                        obj.rois.rois(n).mean = mean(sl(obj.rois.rois(n).data), 'omitnan');
                        obj.rois.rois(n).std = std(sl(obj.rois.rois(n).data), 'omitnan');
                    else
                        obj.rois.rois(n).mean = 0;
                        obj.rois.rois(n).std = 0;
                    end
                else
                    obj.rois.rois(n).data = [];
                    obj.rois.rois(n).mean = 0;
                    obj.rois.rois(n).std = 0;
                    obj.rois.rois(n).size = 0;
                end
            end
            %set(obj.widgets.value,'String',[num2str(obj.rois.rois(obj.rois.nrois).mean,'%8.2f'),' ± ',num2str(obj.rois.rois(obj.rois.nrois).std,'%8.2f'),' (',num2str(obj.rois.rois(obj.rois.nrois).size,'%d'),')']);

        end

        function DrawCross(obj, coords)
            if numel(coords) == 1
                if coords == 1
                    obj.mode.Cross = line('Parent',obj.widgets.ImageAxes,'XData',[-1,1,0,0,0]*(obj.widgets.ImageAxes.XLim(2)-obj.widgets.ImageAxes.XLim(1))/20,'YData',[0,0,0,-1,1]*(obj.widgets.ImageAxes.YLim(2)-obj.widgets.ImageAxes.YLim(1))/20,'Color','red','Visible',1);
                elseif coords == 0
                    if isobject(obj.mode.Cross)
                        obj.mode.Cross.Visible = "off";
                    end
                end
            else
                if obj.mode.UseRealCoords == 1
                    coords = obj.Image.getindcoords([coords,0]);
                end
                if isobject(obj.mode.Cross)
                    set(obj.mode.Cross,'XData',[-1,1,0,0,0]*(obj.widgets.ImageAxes.XLim(2)-obj.widgets.ImageAxes.XLim(1))/20+coords(1),'YData',[0,0,0,-1,1]*(obj.widgets.ImageAxes.YLim(2)-obj.widgets.ImageAxes.YLim(1))/20+coords(2),'Color','red','Visible',1);
                else
                    obj.mode.Cross = line('Parent',obj.widgets.ImageAxes,'XData',[-1,1,0,0,0]*(obj.widgets.ImageAxes.XLim(2)-obj.widgets.ImageAxes.XLim(1))/20+coords(1),'YData',[0,0,0,-1,1]*(obj.widgets.ImageAxes.YLim(2)-obj.widgets.ImageAxes.YLim(1))/20+coords(2),'Color','red','Visible',1);
                end
            end
        end



        function DrawImage(obj, shift)
            if nargin < 2
                shift = [0,0,0];
            elseif numel(shift)<3
                shift = [shift(1),0,0];                
            end
            if obj.Image.NDims == 3
                set(obj.widgets.sliderlabel, 'String',num2str(obj.Image.CurrentSlice));
            end
            if numel(obj.Overlay) == 0   %There's no overlay
                if max(abs(shift)) == 0
                    set(obj.widgets.Image,'CData',obj.Image.ImageData(:,:,obj.Image.CurrentSlice));
                else
                    set(obj.widgets.Image,'CData',circshift(obj.Image.ImageData(:,:,obj.Image.CurrentSlice+shift(3)),[shift(2), shift(1)]));
                end
            else   %There exists an overlay
                if max(abs(shift)) == 0
                    set(obj.widgets.Image,'CData',obj.Image.ImageData(:,:,obj.Image.CurrentSlice));
                else
                    set(obj.widgets.Image,'CData',circshift(obj.Image.ImageData(:,:,obj.Image.CurrentSlice+shift(3)),[shift(2), shift(1)]));
                end
                % Define the combined colormap
                %cm = [obj.Image.ColorMap; obj.Overlay.ColorMap];
                
                if obj.Overlay.NDims == 3
                    set(obj.widgets.Overlay,'CData',squeeze(obj.Overlay.TrueColorData(:,:,obj.Image.CurrentSlice,1:3)),'AlphaData',obj.Overlay.TrueColorData(:,:,obj.Image.CurrentSlice,4));
                else
                    set(obj.widgets.Overlay,'CData',squeeze(obj.Overlay.TrueColorData(:,:,1,1:3)),'AlphaData',obj.Overlay.TrueColorData(:,:,1,4).^1.0);
                end
            end
        end
        
        function Destroy(obj, src, evt)
            delete(obj);
        end
        

    end
    
    
end

