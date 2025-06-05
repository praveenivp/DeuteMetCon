function out = rCSIView(data, image,repetition)
% rCSIView: Function to view CSI data 
% data is a DataClass object containing CSI data
% image is an optional DataClass object containing the reference images

out=0;
if nargin == 0
    return;
end
if nargin < 3
    repetition = 1;
end
if nargin == 1 || numel(image) ==0
    mode.RefImage = 0;
else
    if isa(image, 'DataClass')                  %Currently, only 'DataClass' reference images are allowed
        mode.RefImage = 1;
        RefImage = ImageClass(flipud(image.ProData));   % The flipud is to make sure the orientation is the same as in a simple imagesc  
        RefImage.Params(image.AcqParameters, image.RecoParameters);
    else
        mode.refImage = 0;
    end
end
if isa(data,'DataClass')
    mode.size = size(data.ProData,1,2,3,4,5);
    acqparams = data.AcqParameters;
    recoparams = data.RecoParameters;
    if mode.size(5) > 1 
        if mode.size(5) < repetition
            repetition = 1;
        end
        data = data.ProData(:,:,:,:,repetition);
    else
        data = data.ProData;
    end
else
    if isstruct(data)                   % Here, there should be some code to transform data in a different format (struct, just array) into a DataClass object
        if isfield(data,'raw')
            if isfield(data,'pars')
                acqparams = data.pars;
            end
            data = data.raw;
        else
            disp('Structure doesnt have a raw field!');
            return
        end
    end
    mode.size = size(data);
end
if numel(mode.size) == 4 || numel(mode.size) == 5
    mode.dims = 3;
    mode.slice = round(mode.size(4)/2);
else
    if numel(mode.size) == 3
        mode.dims = 2;
        mode.size(4) = 1;
        mode.slice = 1;
    else
        disp('This is not CSI data!');
        return;
    end
end
mode.select = 0;
 
if mode.dims == 2
    CSIim = ImageClass(squeeze(sum(abs(data(:,:,:,mode.slice)))));
elseif mode.dims == 3
    CSIim = ImageClass(squeeze(sum(abs(data))));
end
if nargin == 2
    wsize = [0,0,1800,500];
    refwpos = [0,0,wsize(3)/4,wsize(4)];
    imwpos = [wsize(3)/4,0,wsize(3)/4,wsize(4)];
    specwpos = [wsize(3)/2,0,wsize(3)/2,wsize(4)];
else
    wsize = [100,100,1500,500];
    imwpos = [0,0,500,500];
    specwpos = [500,0,1000,500];
end
widgets.base = figure('Name',['CSIView: ',inputname(1)],'NumberTitle','off','DockControls','off','Pointer','crosshair',...
                    'Toolbar','none','Menubar','none','Unit','Pixels','Position',wsize,'Color',[0.8,0.8,0.8],'UserData',mode);
if nargin == 2
    widgets.ref = ImageWindow(RefImage,refwpos,widgets.base);
end
widgets.image = ImageWindow(CSIim, imwpos, widgets.base);
CSIim.Params(acqparams,recoparams)
%widgets.image.SetRealCoords()
spec.raw = squeeze(data(:,round(mode.size(2)/2),round(mode.size(3)/2),round(mode.size(4)/2)));
spec.params.field = acqparams.B0;
spec.params.frequency = acqparams.Frequency;
spec.params.bandwidth = acqparams.Bandwidth;
widgets.plot = SpecWin(spec,specwpos, widgets.base);
set(widgets.base, 'WindowButtonMotionFcn',{@CSIViewCallback,'Motion',widgets,data},'KeyPressFcn',{@CSIViewCallback,'Key',widgets},'WindowScrollWheelFcn',{@CSIViewCallback,'ScrollWheel',widgets}, 'ResizeFcn',{@CSIViewCallback,'Resize',widgets}, 'WindowButtonUpFcn',{@CSIViewCallback,'ClickUp',widgets, data});
set(widgets.plot.widgets.PlotAxes, 'ButtonDownFcn',{@CSIViewCallback,'Click', widgets});
set(widgets.image.widgets.orientbuttons, 'SelectionChangedFcn',{@CSIViewCallback,'ChangeImOrient', widgets});
set(widgets.image.widgets.slider, 'Callback',{@CSIViewCallback,'ImSlider', widgets});
if mode.RefImage == 1
    set(widgets.ref.widgets.orientbuttons, 'SelectionChangedFcn',{@CSIViewCallback,'ChangeRefOrient', widgets});
    set(widgets.ref.widgets.slider, 'Callback',{@CSIViewCallback,'RefSlider', widgets});
end
end

function res=CSIViewCallback(src, evt, action, widgets, data, slicepos)
%% Callback Routines for the CSIView widget
res = 0;
switch action
    case 'Motion'        
        mode = get(widgets.base,'UserData');
        p = widgets.base.CurrentPoint;
        if mode.RefImage == 1
            rwpos = widgets.ref.widgets.main.Position;
            rpos = widgets.ref.widgets.ImageAxes.Position;
        end
        wpos = widgets.image.widgets.main.Position;
        pos = widgets.image.widgets.ImageAxes.Position;
        plotpos = widgets.plot.widgets.main.Position;
        if p(1)>wpos(1)+pos(1) && p(1) < wpos(1)+pos(1)+pos(3) && p(2)>wpos(2)+pos(2) && p(2) < wpos(2)+pos(2)+pos(4)   % We're in the csi image 
            [~,pp] = widgets.image.ImageWindowCallback(src,evt,action);
            coords = widgets.image.Image.getrealcoords([pp(1:2),1]);
            widgets.image.DrawCross(0);
            if mode.RefImage == 1
                widgets.ref.DrawCross(coords(1:2));
            end
            if widgets.image.mode.ROIactive ~= 1 && numel(pp)>1   %Not if a ROI is selected
                if widgets.image.mode.Orientation == 0
                   % [pp(2),pp(1),widgets.image.Image.CurrentSlice]
                    widgets.plot.SetSpec(data(:,pp(2),pp(1),widgets.image.Image.CurrentSlice));
                elseif widgets.image.mode.Orientation == 1
                    widgets.plot.SetSpec(data(:,widgets.image.Image.CurrentSlice,pp(1),pp(2)));
                elseif widgets.image.mode.Orientation == 2
                    widgets.plot.SetSpec(data(:,pp(1),widgets.image.Image.CurrentSlice,pp(2)));
                end
            end
        elseif mode.RefImage == 1 && p(1)>rpos(1) && p(1) < rpos(1)+rpos(3) && p(2)>rpos(2) && p(2) < rpos(2)+rpos(4)   % We're in the reference image
            [~,pp] = widgets.ref.ImageWindowCallback(src,evt,action);
            widgets.ref.DrawCross(0);
            coords = widgets.ref.Image.getrealcoords([pp(1:2),1]);
            widgets.image.DrawCross(coords(1:2));
            if widgets.image.mode.ROIactive ~= 1 && numel(pp)>1   %Not if a ROI is selected
               % widgets.plot.SetSpec(data(:,pp(2),pp(1),widgets.image.Image.CurrentSlice));
            end
        elseif p(1)>plotpos(1) && p(1) < plotpos(1)+plotpos(3) && p(2)>plotpos(2) && p(2) < plotpos(2)+plotpos(4)
            widgets.plot.PlotWindowCallback(src,evt,action);
            if mode.RefImage == 1
                widgets.ref.DrawCross(0);
            end
        else
            if mode.RefImage == 1
                widgets.ref.DrawCross(0);
            end
        end
    case 'Click'
        if evt.Button == 1
            mode = get(widgets.base,'UserData');
            widgets.plot.DefineXPoint(2);
            mode.select = 2;
            set(widgets.base,'UserData',mode);
        end
        widgets.plot.PlotWindowCallback(src, evt, action);
    case 'Key'
    case 'ScrollWheel'
    case 'Resize'
            p = get(widgets.base,'Position');
            mode = get(widgets.base,'UserData');
            if mode.RefImage == 1
                set(widgets.image.widgets.main, 'Position',[floor(p(3)/4),0, floor(p(3)/4),  p(4)]);
                widgets.image.ImageWindowCallback(src,evt,'PanelResize');
                set(widgets.ref.widgets.main, 'Position',[0,0, floor(p(3)/4),  p(4)]);
                widgets.ref.ImageWindowCallback(src,evt,'PanelResize');
                set(widgets.plot.widgets.main,'Position',[ floor(p(3)/2),0, floor(p(3)/2),  p(4)]);
            else
                set(widgets.image.widgets.main, 'Position',[0,0, floor(p(3)/3),  p(4)]);
                widgets.image.ImageWindowCallback(src,evt,'PanelResize');
                set(widgets.plot.widgets.main,'Position',[ floor(p(3)/3),0, 2*floor(p(3)/3),  p(4)]);
            end
            widgets.plot.SpecWinCallback(widgets.plot.widgets.main,evt,'Resize',widgets.plot.specwidgets)
%             w.phase0slider = findobj('Tag','Phase0Slider');
%             w.phase1slider = findobj('Tag','Phase1Slider');
%             w.phase1value = findobj('Tag','Phase1Value');
%             w.phase0value = findobj('Tag','Phase0Value');
%             w.phase1label = findobj('Tag','Phase1Label');
%             w.phase0label = findobj('Tag','Phase0Label');
%             w.buttongroup1 = findobj('Tag','bgFIDSpec');
%             w.buttongroup2 = findobj('Tag','bgCompl');
%             w.buttongroup3 = findobj('Tag','bgHzppm');            
%             pos = [widgets.image.WindowSize(1),0,s(1)-widgets.image.WindowSize(1),s(2)];
%             widgets.plot.SpecWinCallback(src,evt,'Resize',w,pos);
    case 'ClickUp'
        mode = get(widgets.base,'UserData');
        p = widgets.base.CurrentPoint;
        impos = widgets.image.widgets.ImageAxes.Position;
        imwpos = widgets.image.widgets.main.Position;
        impos(1:2) = impos(1:2) + imwpos(1:2);
        plotpos = widgets.plot.widgets.main.Position;
        if mode.RefImage == 1
            rwpos = widgets.ref.widgets.main.Position;
            rpos = widgets.ref.widgets.ImageAxes.Position;
        end
        if p(1)>impos(1) && p(1) < impos(1)+impos(3) && p(2)>impos(2) && p(2) < impos(2)+impos(4)
            [~,pp] = widgets.image.ImageWindowCallback(src,evt,action);
            sp = size(pp);
            if sp(1) == 1 && sp(2) == 2   % Just one point selected
            elseif numel(pp)>1
                if mode.dims == 2   % only works for 2D data
                    widgets.plot.Plot.yDat(fftshift(fft(sum(data(:,pp),2))));
                    widgets.plot.Plot.plot;
                end
            end
        elseif mode.RefImage == 1 && p(1)>rpos(1) && p(1) < rpos(1)+rpos(3) && p(2)>rpos(2) && p(2) < rpos(2)+rpos(4)
            [~,pp] = widgets.ref.ImageWindowCallback(src,evt,action);
            sp = size(pp);
            if sp(1) == 1 && sp(2) == 2   % Just one point selected
            elseif numel(pp)>1
                if mode.dims == 2   % only works for 2D data
                    widgets.plot.Plot.yDat(fftshift(fft(sum(data(:,pp),2))));
                    widgets.plot.Plot.plot;
                end
            end
        elseif p(1)>plotpos(1) && p(1) < plotpos(1)+plotpos(3) && p(2)>plotpos(2) && p(2) < plotpos(2)+plotpos(4)
            widgets.plot.SpecWinCallback(src,evt,action);
            if mode.select > 0            % Definition of a spectral region
                mode.select = mode.select-1;
                set(widgets.base,'UserData',mode);
                if mode.select == 0
                    p = widgets.plot.mode.clicked([1, 3])';
                    s = find(widgets.plot.Plot.Dat.x >= min(p));
                    s1 = s(1);
                    s = find(widgets.plot.Plot.Dat.x <= max(p));
                    s2 = s(end);
                    if mode.dims == 2
                        widgets.image.Image.Data(squeeze(sum(abs(data(s1:s2,:,:,mode.slice)))));
                    elseif mode.dims == 3
                        fprintf('Spectral region: %f to %f (%d to %d)\n',p(1),p(2),s1,s2);
                        siz = size(data);
                        im = abs(circshift(fft(data,widgets.plot.SpecMode.Zerofill),[widgets.plot.SpecMode.Zerofill/2,0,0,0]));
                        widgets.image.Image.ImageData = squeeze(sum(im(s1:s2,:,:,:)));
                        widgets.image.DrawImage;
                    end
                end
            end
        end
    case 'ChangeImOrient'
        widgets.image.ImageWindowCallback(src, evt, 'Orientation');
        mode = get(widgets.base,'UserData');
        if mode.RefImage == 1
            if widgets.image.widgets.orientbuttons.SelectedObject ==  widgets.image.widgets.orientbuttons.Children(3)
                widgets.ref.widgets.orientbuttons.SelectedObject =  widgets.ref.widgets.orientbuttons.Children(3);
            elseif widgets.image.widgets.orientbuttons.SelectedObject ==  widgets.image.widgets.orientbuttons.Children(2)
                widgets.ref.widgets.orientbuttons.SelectedObject =  widgets.ref.widgets.orientbuttons.Children(2);
            elseif widgets.image.widgets.orientbuttons.SelectedObject ==  widgets.image.widgets.orientbuttons.Children(1)
                widgets.ref.widgets.orientbuttons.SelectedObject =  widgets.ref.widgets.orientbuttons.Children(1);
            end
            widgets.ref.ImageWindowCallback(src, evt, 'Orientation');
            CSIViewCallback( widgets.image.widgets.slider,evt,'ImSlider',widgets);
        end
    case 'ChangeRefOrient'
        widgets.ref.ImageWindowCallback(src, evt, 'Orientation');
        mode = get(widgets.base,'UserData');
        if widgets.ref.widgets.orientbuttons.SelectedObject ==  widgets.ref.widgets.orientbuttons.Children(3)
            widgets.image.widgets.orientbuttons.SelectedObject =  widgets.image.widgets.orientbuttons.Children(3);
        elseif widgets.ref.widgets.orientbuttons.SelectedObject ==  widgets.ref.widgets.orientbuttons.Children(2)
            widgets.image.widgets.orientbuttons.SelectedObject =  widgets.image.widgets.orientbuttons.Children(2);
        elseif widgets.ref.widgets.orientbuttons.SelectedObject ==  widgets.ref.widgets.orientbuttons.Children(1)
            widgets.image.widgets.orientbuttons.SelectedObject =  widgets.image.widgets.orientbuttons.Children(1);
        end
        widgets.image.ImageWindowCallback(src, evt, 'Orientation');
        CSIViewCallback( widgets.ref.widgets.slider,evt,'RefSlider',widgets);
    case 'ImSlider'
        widgets.image.ImageWindowCallback(src, evt, 'Slider');
        mode = get(widgets.base,'UserData');
        if mode.RefImage == 1
            slicecoord = widgets.image.Image.Coords.z(widgets.image.Image.CurrentSlice);
            refcoord = abs(widgets.ref.Image.Coords.z - slicecoord);
            refslice = find(refcoord == min(refcoord));
            refslice = refslice(1);
            widgets.ref.widgets.slider.Value = refslice;
            widgets.ref.ImageWindowCallback(widgets.ref.widgets.slider, evt, 'Slider');
        end
    case 'RefSlider'
        widgets.ref.ImageWindowCallback(src, evt, 'Slider');
        mode = get(widgets.base,'UserData');
            slicecoord = widgets.ref.Image.Coords.z(widgets.ref.Image.CurrentSlice);
            imcoord = abs(widgets.image.Image.Coords.z - slicecoord);
            imslice = find(imcoord == min(imcoord));
            imslice = imslice(1);
            widgets.image.widgets.slider.Value = imslice;
            widgets.image.ImageWindowCallback(widgets.image.widgets.slider, evt, 'Slider');
end
end
