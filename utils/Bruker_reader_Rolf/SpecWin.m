classdef SpecWin < PlotWindow 
% SpecWin: class to view spectroscopy data   
properties
    rawFID
    proFID
    spectrum
    nspectra
    rawFIDSize
    proFIDSize
    SpecWinSize
    SpecMode
    datapars
    specwidgets
end

methods
    
function [obj, ret] = SpecWin(data, WindowSize, baseWidget, varargin)
    datapars.frequency = -1;
    datapars.bandwidth = -1;
    datapars.field = -1;
    ret = [];
    SpecWinWait = 0;
    % Default parameters
    defval.Zerofill = 0;
    defval.Spectrum = 0;
    defval.ppm = 0;
    defval.SpecRange = [];

    if (isstruct(data))    % it is assumed that data has a certain structure, with data.raw being the 
                            % FID data and data.params being the parameters
        datapars = data.params;
        data = data.raw;
    end                  % data is just the fid -> parameters have to come from varargin or have to be asked
    plotnames = [];
    for cnt = 1:numel(varargin)
        if ischar(varargin{cnt})
            if contains(varargin{cnt},'params','IgnoreCase',true)
                if numel(varargin) > cnt
                    datapars = varargin{cnt+1};
                end
            elseif contains(varargin{cnt},'Frequency','IgnoreCase',true)
                if numel(varargin) > cnt
                    datapars.frequency = varargin{cnt+1};
                end
            elseif contains(varargin{cnt},'Bandwidth','IgnoreCase',true)
                if numel(varargin) > cnt
                    datapars.bandwidth = varargin{cnt+1};
                end
            elseif contains(varargin{cnt},'Field','IgnoreCase',true)
                if numel(varargin) > cnt
                    datapars.field = varargin{cnt+1};
                end
            elseif contains(varargin(cnt),'Zerofill','IgnoreCase',true)
                if numel(varargin) > cnt
                    defval.Zerofill = varargin{cnt+1};
                end
            elseif contains(varargin(cnt),'Title','IgnoreCase',true)
                if numel(varargin) > cnt
                    wintitle = varargin{cnt+1};
                end
            elseif contains(varargin(cnt),'SpecNames','IgnoreCase',true)
                if numel(varargin) > cnt
                    plotnames = varargin{cnt+1};
                end
            elseif contains(varargin(cnt),'SpecZoom','IgnoreCase',true)
                if numel(varargin) > cnt
                    defval.SpecRange = varargin{cnt+1};
                end
            elseif contains(varargin(cnt),'Spectrum','IgnoreCase',true)
                defval.Spectrum = 1;
            elseif contains(varargin(cnt),'ppm','IgnoreCase',true)
                defval.ppm = 1;
            elseif contains(varargin{cnt},'wait','IgnoreCase',true)
                    SpecWinWait = 1;
            end
        end
    end
    if exist('wintitle','var') == 0
        wintitle = inputname(1);
    end
    wintitle = ['SpecWin: ',wintitle];
    if datapars.frequency == -1 || datapars.bandwidth == -1 || datapars.field == -1
        prompt = {'Bandwidth:','Fieldstrength:','Frequency:'};
        dlg_title = 'Measurement Parameters';
        def = {'10000','9.4','161.807825'};
        if datapars.frequency > -1
            def{3} = num2str(datapars.frequency);
        end
        if datapars.bandwidth > -1
            def{1} = num2str(datapars.bandwidth);
        end
        if datapars.field > -1
            def{2} = num2str(datapars.field);
        end
        answer = inputdlg(prompt,dlg_title,1,def);
        datapars.frequency = str2double(answer{3});
        datapars.bandwidth = str2double(answer{1});
        datapars.field = str2double(answer{2});
    end   
    if nargin < 2 || numel(WindowSize) == 0
        WindowSize = [1000,500];
    end
    if numel(WindowSize)==1
        WindowSize = [WindowSize,WindowSize];
    end
    if numel(WindowSize)<4
        Pos = [100,100,WindowSize(1), WindowSize(2)];
    else
        Pos = WindowSize;
    end
    WindowSize = [Pos(3),Pos(4)];
    data = squeeze(data);
    if nargin < 3 || numel(baseWidget) == 0
        widgets.base = figure('Name',wintitle,'NumberTitle','off','DockControls','off','Pointer','crosshair',...
                              'Toolbar','none','Menubar','none','Unit','Pixels','Position',Pos,'Color',[0.8,0.8,0.8]);
        wpos = [0,0,1000,500];
    else
        widgets.base = baseWidget;
        wpos = Pos;
    end
    PlotWinSize = [wpos(1)+82,wpos(2)+1,WindowSize(1)-82-100,WindowSize(2)-1];
    obj@PlotWindow(data,PlotWinSize,widgets.base,'PlotNames',plotnames);
    obj.SpecMode.FidSpec = defval.Spectrum;   % 0: FID, 1: Spectrum
    obj.SpecMode.Compl = 0;     % 0: Real, 1: Imaginary, 2: Phase, 3: Abs, 4: Both
    obj.SpecMode.Phase0 = 0;
    obj.SpecMode.Phase1 = 0;
    obj.SpecMode.Phase1Ref = 0;
    obj.SpecMode.SigRegion = [];
    obj.SpecMode.NoiseRegion = [];
    obj.SpecMode.Filter = 0;
    obj.SpecMode.Scaling = 1;
    obj.SpecMode.HzPpm = defval.ppm;     % 0: Hz, 1: ppm
    obj.SpecMode.SliderPhase = 0;  % 0: Phase, 1: Filter
    obj.SpecMode.Zerofill = defval.Zerofill;
    obj.SpecMode.Zero = 0;
    obj.SpecMode.Def = 0;         % 1: define zero, 2: define phase1ref, 3: define signal region, 4: define noise region
    if obj.SpecMode.HzPpm == 0
        obj.SpecMode.SpecRange = defval.SpecRange;
    else
        obj.SpecMode.SpecRange = defval.SpecRange*datapars.frequency;
    end
    obj.SpecMode.FidRange = [];
    
    % Additional widget elements
    widgets.buttongroup1 = uibuttongroup('Parent',widgets.base,'Units','Pixels','Position',[wpos(1)+1,wpos(2)+WindowSize(2)-41,80,40],'BackgroundColor',[0.8,0.8,0.8],'Tag','bgFIDSpec');
    widgets.buttons(1) = uicontrol('Parent',widgets.buttongroup1,'Style','radiobutton','Position',[5,22,70,13],'String','FID','BackgroundColor',[0.8,0.8,0.8],'Value',1);
    widgets.buttons(2) = uicontrol('Parent',widgets.buttongroup1,'Style','radiobutton','Position',[5,5,70,13],'String','Spectrum','BackgroundColor',[0.8,0.8,0.8],'Value',obj.SpecMode.FidSpec);
    widgets.buttongroup2 = uibuttongroup('Parent',widgets.base,'Units','Pixels','Position',[wpos(1)+1,wpos(2)+WindowSize(2)-41-98,80,97],'BackgroundColor',[0.8,0.8,0.8],'Tag','bgCompl');
    widgets.buttons(3) = uicontrol('Parent',widgets.buttongroup2,'Style','radiobutton','Position',[5,73,70,13],'String','Real','BackgroundColor',[0.8,0.8,0.8]);
    widgets.buttons(4) = uicontrol('Parent',widgets.buttongroup2,'Style','radiobutton','Position',[5,56,70,13],'String','Imaginary','BackgroundColor',[0.8,0.8,0.8]);
    widgets.buttons(5) = uicontrol('Parent',widgets.buttongroup2,'Style','radiobutton','Position',[5,39,70,13],'String','Abs','BackgroundColor',[0.8,0.8,0.8]);
    widgets.buttons(6) = uicontrol('Parent',widgets.buttongroup2,'Style','radiobutton','Position',[5,22,70,13],'String','Phase','BackgroundColor',[0.8,0.8,0.8]);
    widgets.buttons(7) = uicontrol('Parent',widgets.buttongroup2,'Style','radiobutton','Position',[5,5,70,13],'String','Both','BackgroundColor',[0.8,0.8,0.8]);
    widgets.buttongroup3 = uibuttongroup('Parent',widgets.base,'Units','Pixels','Position',[wpos(1)+1,wpos(2)+WindowSize(2)-41-98-41,80,40],'BackgroundColor',[0.8,0.8,0.8],'Tag','bgHzppm');
    widgets.buttons(8) = uicontrol('Parent',widgets.buttongroup3,'Style','radiobutton','Position',[5,22,70,13],'String','Hz','BackgroundColor',[0.8,0.8,0.8],'Value',1);
    widgets.buttons(9) = uicontrol('Parent',widgets.buttongroup3,'Style','radiobutton','Position',[5,5,70,13],'String','ppm','BackgroundColor',[0.8,0.8,0.8],'Value',obj.SpecMode.HzPpm);

    widgets.phase0slider =  uicontrol('Parent',widgets.base,'Units','Pixels','Style','slider','Position',[PlotWinSize(1)+PlotWinSize(3)+15,PlotWinSize(2)+50,20,PlotWinSize(4)-60],'Value',obj.SpecMode.Phase0, 'Min',-180,'Max',180, 'SliderStep',[1/360,5/360],'Tag','Phase0Slider');
    widgets.phase0value = uicontrol('Parent',widgets.base,'Units','Pixels','Style','text','Position',[PlotWinSize(1)+PlotWinSize(3),PlotWinSize(2)+30,50,15],'String', [num2str(obj.SpecMode.Phase0),'°'],'BackgroundColor',[0.8,0.8,0.8],'Tag','Phase0Value');
    widgets.phase0label = uicontrol('Parent',widgets.base,'Units','Pixels','Style','text','Position',[PlotWinSize(1)+PlotWinSize(3),PlotWinSize(2)+10,50,15],'String', 'Phase 0','BackgroundColor',[0.8,0.8,0.8],'Tag','Phase0Label');
    widgets.phase1slider =  uicontrol('Parent',widgets.base,'Units','Pixels','Style','slider','Position',[PlotWinSize(1)+PlotWinSize(3)+65,PlotWinSize(2)+50,20,PlotWinSize(4)-60],'Value',obj.SpecMode.Phase1, 'Min',-180,'Max',180, 'SliderStep',[1/360,5/360],'Tag','Phase1Slider');
    widgets.phase1value = uicontrol('Parent',widgets.base,'Units','Pixels','Style','text','Position',[PlotWinSize(1)+PlotWinSize(3)+50,PlotWinSize(2)+30,50,15],'String', [num2str(obj.SpecMode.Phase1),'°/ppm'],'BackgroundColor',[0.8,0.8,0.8],'Tag','Phase1Value');
    widgets.phase1label = uicontrol('Parent',widgets.base,'Units','Pixels','Style','text','Position',[PlotWinSize(1)+PlotWinSize(3)+50,PlotWinSize(2)+10,50,15],'String', 'Phase 1','BackgroundColor',[0.8,0.8,0.8],'Tag','Phase1Label');

    % Expand the Plotwidget menus
    widgets.specmenu1 = uimenu(obj.widgets.cm,'Label','Define','Visible',0);
    widgets.SetZero = uimenu(widgets.specmenu1,'Label','Zero','Callback',{@obj.SpecWinCallback,'SetZero',widgets}, 'Visible',0);
    widgets.Phase1Ref = uimenu(widgets.specmenu1,'Label','Phase Reference','Callback',{@obj.SpecWinCallback,'Phase1Ref',widgets}, 'Visible',0);
    widgets.SetSig = uimenu(widgets.specmenu1,'Label','Signal Region','Callback',{@obj.SpecWinCallback,'SigReg',widgets}, 'Visible',0);
    widgets.SetNoise = uimenu(widgets.specmenu1,'Label','Noise Region','Callback',{@obj.SpecWinCallback,'NoiseReg',widgets}, 'Visible',0);

    widgets.specmenu2 = uimenu(obj.widgets.cm,'Label','Process');
    uimenu(widgets.specmenu2,'Label','Phase','Callback',{@obj.SpecWinCallback,'Phase',widgets});
    uimenu(widgets.specmenu2,'Label','Filter','Callback',{@obj.SpecWinCallback,'Filter',widgets});
    uimenu(widgets.specmenu2,'Label','Zerofill','Callback',{@obj.SpecWinCallback,'Zerofill',widgets});
    uimenu(widgets.specmenu2,'Label','DC offset','Callback',{@obj.SpecWinCallback,'DCOffset',widgets});
    uimenu(widgets.specmenu2,'Label','Normalize','Callback',{@obj.SpecWinCallback,'Normalize',widgets});
    uimenu(obj.widgets.cm2,'Label','Export for jMRUI','Callback',{@obj.SpecWinCallback,'ExportjMRUI'});
    uimenu(obj.widgets.cm2,'Label','Fast Export','Callback',{@obj.SpecWinCallback,'FastExport'});
    widgets.specmenu2 = uimenu(obj.widgets.cm,'Label','Analyse');
    uimenu(widgets.specmenu2,'Label','Linewidth','Callback',{@obj.SpecWinCallback,'Linewidth',widgets});

    % Callback functions
    set(widgets.buttongroup1, 'SelectionChangeFcn',{@obj.SpecWinCallback,'ChooseFidSpec',widgets});
    set(widgets.buttongroup2, 'SelectionChangeFcn',{@obj.SpecWinCallback,'ChooseComplex',widgets});
    set(widgets.buttongroup3, 'SelectionChangeFcn',{@obj.SpecWinCallback,'ChooseHzPpm',widgets});
    set(widgets.phase0slider,'Callback',{@obj.SpecWinCallback,'Phase0',widgets});
    set(widgets.phase1slider,'Callback',{@obj.SpecWinCallback,'Phase1',widgets});
    % Capture callbacks from the PlotWindow
    set(obj.widgets.cm0,'Callback',{@obj.SpecWinCallback,'zoomOut'});

    set(widgets.base, 'WindowButtonUpFcn',{@obj.SpecWinCallback,'ClickUp', widgets},'ResizeFcn',{@obj.SpecWinCallback,'Resize',widgets});
    obj.specwidgets = widgets;
    obj.datapars = datapars;
    obj.rawFID = data;
    s = size(data);
    obj.rawFIDSize = s(1);
    obj.nspectra = s(2);
    obj.proFID = obj.rawFID;
    obj.proFIDSize = obj.rawFIDSize;
    obj.SpecMode.FidRange = [0,(obj.rawFIDSize-1)/datapars.bandwidth*1000];
    if numel(obj.SpecMode.SpecRange) == 0
        obj.SpecMode.SpecRange = [-obj.datapars.bandwidth/2,obj.datapars.bandwidth/2] - obj.SpecMode.Zero;
    end
    obj.SpecMode.Zerofill = max([obj.proFIDSize, obj.SpecMode.Zerofill]);
    obj.CalcCurrentSpec;
    obj.DrawSpec;
    if SpecWinWait == 1
       waitfor(obj.specwidgets.base); 
       ret.proFID = obj.proFID;
       ret.spectrum = obj.spectrum;
       ret.SpecMode = obj.SpecMode;
    end
end

function obj=SpecWinCallback(obj, src, evt, action, widgets, data, slicepos)
%% Callback Routines for the SpecWin widget
switch action
    case 'Key'
    case 'ScrollWheel'
    case 'Resize'
        if nargin < 6
            WindowSize = get(src,'Position');
            %WindowSize(1:2) = 0;
        else                    % There exists a 'data'-argument, that signifies the position of the plot window
            WindowSize = data;
        end
        PlotWinSize = [82,1,WindowSize(3)-82-100,WindowSize(4)-1];
        obj.WindowSize = PlotWinSize(3:4);
        set(widgets.phase0slider,'Position',[PlotWinSize(1)+PlotWinSize(3)+15,50,20,PlotWinSize(4)-60]);
        set(widgets.phase0value, 'Position',[PlotWinSize(1)+PlotWinSize(3),30,50,15]);
        set(widgets.phase0label, 'Position',[PlotWinSize(1)+PlotWinSize(3),10,50,15]);
        set(widgets.phase1slider,'Position',[PlotWinSize(1)+PlotWinSize(3)+65,50,20,PlotWinSize(4)-60]);
        set(widgets.phase1value ,'Position',[PlotWinSize(1)+PlotWinSize(3)+50,30,50,15]);
        set(widgets.phase1label, 'Position',[PlotWinSize(1)+PlotWinSize(3)+50,10,50,15]);
        set(widgets.buttongroup1, 'Position',[1,PlotWinSize(4)-41,80,40]);
        set(widgets.buttongroup2 ,'Position',[1,PlotWinSize(4)-41-98,80,97]);
        set(widgets.buttongroup3 ,'Position',[1,PlotWinSize(4)-41-98-41,80,40]);
        set(obj.widgets.main, 'Position',PlotWinSize);
        obj.PlotWindowCallback(src,evt,'PanelResize');
    case 'ClickUp'
          % Run standard callback
          obj.PlotWindowCallback(src,evt,action);

          %Handle Hz / ppm
          if obj.SpecMode.FidSpec == 1
              if obj.SpecMode.HzPpm == 0
                  obj.SpecMode.SpecRange = obj.Plot.xRange;
              else
                  obj.SpecMode.SpecRange = obj.Plot.xRange*obj.datapars.frequency;
              end
          elseif obj.SpecMode.FidSpec == 0
              obj.SpecMode.FidRange = obj.Plot.xRange;
          end

          if obj.SpecMode.Def == 1     % define zero
              coord = obj.mode.clicked;
              obj.SpecMode.Def = 0;
              f = inputdlg('Frequency of selected point','Set frequency',1,{'0'});
              if numel(f) == 0
                  return;
              end
              frequency = str2num(f{1});
              if obj.SpecMode.HzPpm == 0
                  obj.SpecMode.Zero = obj.SpecMode.Zero + coord(1) - frequency;
                  obj.SpecMode.SpecRange = obj.SpecMode.SpecRange - coord(1) + frequency;
              elseif obj.SpecMode.HzPpm == 1
                  obj.SpecMode.Zero = obj.SpecMode.Zero + (coord(1) - frequency)*obj.datapars.frequency;
                  obj.SpecMode.SpecRange = obj.SpecMode.SpecRange - coord(1)*obj.datapars.frequency + frequency*obj.datapars.frequency;
              end
              obj.DrawSpec; 

          elseif obj.SpecMode.Def == 2 % define phase1ref
              coord = obj.mode.clicked;
              oldref = obj.SpecMode.Phase1Ref;
              if obj.SpecMode.HzPpm == 0
                  obj.SpecMode.Phase1Ref = coord(1);
              else
                  obj.SpecMode.Phase1Ref = coord(1)*obj.datapars.frequency;
              end
              obj.SpecMode.Phase0 = obj.SpecMode.Phase0 + obj.SpecMode.Phase1*(obj.SpecMode.Phase1Ref-oldref)/obj.datapars.frequency;
              if (obj.SpecMode.Phase0 > pi)
                  obj.SpecMode.Phase0 = obj.SpecMode.Phase0 - 2*pi;
              elseif (obj.SpecMode.Phase0 < -pi)
                  obj.SpecMode.Phase0 = obj.SpecMode.Phase0 + 2*pi;
              end
              set(obj.specwidgets.phase0slider,'Value',obj.SpecMode.Phase0*180/pi);
              set(obj.specwidgets.phase0value,'String',[num2str(obj.SpecMode.Phase0*180/pi),'°']);
              if obj.SpecMode.FidSpec == 0
                obj.SpecMode.FidRange = obj.Plot.xRange;  
              else
                if obj.SpecMode.HzPpm == 0
                    obj.SpecMode.SpecRange = obj.Plot.xRange;
                else
                    obj.SpecMode.SpecRange = obj.Plot.xRange*obj.datapars.frequency;
                end
              end
              obj.SpecMode.Def = 0;
              obj.CalcCurrentSpec;
              obj.DrawSpec;
           elseif obj.SpecMode.Def == 3 % define signal region 
              coord = obj.mode.clicked;
              obj.SpecMode.Def = 13;
              if obj.SpecMode.HzPpm == 0
                  obj.SpecMode.SigRegion(1) = coord(1);
              elseif obj.SpecMode.HzPpm == 1
                  obj.SpecMode.SigRegion(1) = coord(1)*obj.datapars.frequency;
              end
              obj.DefineXPoint(1);
           elseif obj.SpecMode.Def == 13    % signal region, second point
              coord = obj.mode.clicked;
              obj.SpecMode.Def = 0;
              if obj.SpecMode.HzPpm == 0
                  obj.SpecMode.SigRegion(2) = coord(1);
              elseif obj.SpecMode.HzPpm == 1
                  obj.SpecMode.SigRegion(2) = coord(1)*obj.datapars.frequency;
              end
              fprintf('Signal region set to %f to %f\n',obj.SpecMode.SigRegion(1),obj.SpecMode.SigRegion(2))
           elseif obj.SpecMode.Def == 4 % define noise region (not implemented!)
              coord = obj.mode.clicked;
              obj.SpecMode.Def = 14;
              if obj.SpecMode.HzPpm == 0
                  obj.SpecMode.NoiseRegion(1) = coord(1);
              elseif obj.SpecMode.HzPpm == 1
                  obj.SpecMode.NoiseRegion(1) = coord(1)*obj.datapars.frequency;
              end
              obj.DefineXPoint(1);
           elseif obj.SpecMode.Def == 14    % noise region, second point
              coord = obj.mode.clicked;
              obj.SpecMode.Def = 0;
              if obj.SpecMode.HzPpm == 0
                  obj.SpecMode.NoiseRegion(2) = coord(1);
              elseif obj.SpecMode.HzPpm == 1
                  obj.SpecMode.NoiseRegion(2) = coord(1)*obj.datapars.frequency;
              end
              fprintf('Noise region set to %f to %f\n',obj.SpecMode.NoiseRegion(1),obj.SpecMode.NoiseRegion(2))

          end
    case 'ChooseFidSpec'
        if strcmp(evt.NewValue.String, 'Spectrum') == 1    % now switched to Spec
            obj.SpecMode.FidSpec = 1;
            obj.CalcCurrentSpec;
            if obj.SpecMode.Compl <= 1  || obj.SpecMode.Compl == 4        % real or imaginary or both
                set(obj.Plot.PlotAxes,'YLim',[-max(abs(obj.spectrum(:))),max(abs(obj.spectrum(:)))]);
            elseif obj.SpecMode.Compl == 3     % abs
                set(obj.Plot.PlotAxes,'YLim',[0,max(abs(obj.spectrum))]);
            elseif obj.SpecMode.Compl == 2
                set(obj.Plot.PlotAxes,'YLim',[-pi,pi]);
            end
            obj.DrawSpec;
            set(widgets.specmenu1,'Visible',1);
            set(widgets.SetZero,'Visible',1);
            set(widgets.Phase1Ref,'Visible',1);
            set(widgets.SetSig,'Visible',1);
            set(widgets.SetNoise,'Visible',1);
        elseif strcmp(evt.NewValue.String, 'FID') == 1
            obj.SpecMode.FidSpec = 0;
            set(obj.Plot.PlotAxes,'YLimMode','auto');
            obj.DrawSpec;
            set(widgets.specmenu1,'Visible',0);
            set(widgets.SetZero,'Visible',0);
            set(widgets.Phase1Ref,'Visible',0);
        end
    case 'ChooseComplex'
        if strcmp(evt.NewValue.String, 'Real') == 1    % now switched to Real
            obj.SpecMode.Compl = 0;
            obj.Plot.showreal;
            set(obj.Plot.PlotAxes,'YLim',[-max(abs(obj.spectrum(:))),max(abs(obj.spectrum(:)))]);
            obj.DrawSpec;
        elseif strcmp(evt.NewValue.String, 'Imaginary') == 1
            obj.SpecMode.Compl = 1;
            obj.Plot.showimag;
            set(obj.Plot.PlotAxes,'YLim',[-max(abs(obj.spectrum(:))),max(abs(obj.spectrum(:)))]);
            obj.DrawSpec;
        elseif strcmp(evt.NewValue.String, 'Abs') == 1
            obj.SpecMode.Compl = 3;
            obj.Plot.showabs;
            set(obj.Plot.PlotAxes,'YLim',[0,max(abs(obj.spectrum(:)))]);
            obj.DrawSpec;
        elseif strcmp(evt.NewValue.String, 'Phase') == 1
            obj.SpecMode.Compl = 2;
            obj.Plot.showphase;
            set(obj.Plot.PlotAxes,'YLim',[-pi,pi]);
            obj.DrawSpec;
        elseif strcmp(evt.NewValue.String, 'Both') == 1
            obj.SpecMode.Compl = 4;
            obj.Plot.showboth;
            set(obj.Plot.PlotAxes,'YLim',[-max(abs(obj.spectrum)),max(abs(obj.spectrum))]);
            obj.DrawSpec;
        end
    case 'ChooseHzPpm'
        if strcmp(evt.NewValue.String, 'Hz') == 1    % now switched to Hz
            obj.SpecMode.HzPpm = 0;
            if obj.SpecMode.FidSpec == 1;
                obj.DrawSpec;            
            end
        elseif strcmp(evt.NewValue.String, 'ppm') == 1
            obj.SpecMode.HzPpm = 1;
            if obj.SpecMode.FidSpec == 1;
                obj.DrawSpec;
            end
        end
    case 'Phase'
        obj.SpecMode.SliderPhase = 0;
        set(widgets.phase1slider,'Visible',1);
        set(widgets.phase1label,'Visible',1);
        set(widgets.phase1value,'Visible',1);
        set(widgets.phase0slider,'Min',-180,'Max',180,'SliderStep',[1/360,5/360],'Value',obj.SpecMode.Phase0*180/pi);
        set(widgets.phase0label,'String','Phase 0');
        set(widgets.phase0value,'String',[num2str(obj.SpecMode.Phase0*180/pi),'°']);
    case 'Phase0'
        val = get(widgets.phase0slider,'Value');
        if obj.SpecMode.SliderPhase == 0    % The slider controls the phase
            if val == 180
                val = -180;
                set(widgets.phase0slider,'Value',val);
            elseif val == -180
                val = 180;
                set(widgets.phase0slider,'Value',val);
            end            
            set(widgets.phase0value,'String',[num2str(val),'°']);
            obj.SpecMode.Phase0 = val/180*pi;
            obj.CalcCurrentSpec;
            obj.DrawSpec;
                
            
        elseif obj.SpecMode.SliderPhase == 1   % The slider controls the filter
            set(widgets.phase0value,'String',[num2str(val),' Hz']);
            obj.SpecMode.Filter = val;
            obj.CalcCurrentSpec;
            obj.DrawSpec;   
            if obj.SpecMode.FidSpec == 0
                obj.widgets.HelperPlot.XData = linspace(0,obj.proFIDSize-1,obj.proFIDSize)/obj.datapars.bandwidth*1000;
                obj.widgets.HelperPlot.YData = exp(-obj.widgets.HelperPlot.XData *obj.SpecMode.Filter/1000)*mean(abs(obj.proFID(1:4)));
                obj.widgets.HelperPlot.Visible = 1;
            end
        end
    case 'Phase1'
        val = get(widgets.phase1slider,'Value');
            set(widgets.phase1value,'String',[num2str(val),'°/ppm']);
            obj.SpecMode.Phase1 = val/180*pi;
            obj.CalcCurrentSpec;
            obj.DrawSpec;
    case 'Filter'
        obj.SpecMode.SliderPhase = 1;
        set(widgets.phase1slider,'Visible',0);
        set(widgets.phase1label,'Visible',0);
        set(widgets.phase1value,'Visible',0);
        set(widgets.phase0slider,'Min',0,'Max',100,'SliderStep',[1/100,5/100],'Value',obj.SpecMode.Filter);
        set(widgets.phase0label,'String','Filter');
        set(widgets.phase0value,'String',[num2str(obj.SpecMode.Filter),' Hz']);
    case 'Zerofill'
        zf = 0;
        while zf < obj.rawFIDSize || zf > 20*obj.rawFIDSize
            def = {num2str(max([obj.rawFIDSize, obj.SpecMode.Zerofill]))};
            zfstring = inputdlg('Zerofill to','Zerofilling',1, def);
            if numel(zfstring) == 0
                return
            end
            zf = str2num(zfstring{1});
            
        end
        obj.SpecMode.Zerofill = zf;
        obj.CalcCurrentSpec;
        obj.DrawSpec;
    case 'DCOffset'
        def = round((obj.rawFIDSize)/32);
        npoints =  inputdlg('Number of final points for correction:','DC offset',1, {num2str(def)});
        if numel(npoints) == 0
            return;
        end
        npoints = str2double(npoints);
        obj.rawFID = obj.rawFID - mean(obj.rawFID(end-npoints,:));
        obj.CalcCurrentSpec;
        obj.DrawSpec;
    case 'Normalize'
        Hz = linspace(-obj.datapars.bandwidth/2,obj.datapars.bandwidth/2,obj.SpecMode.Zerofill);
        lims = sort([find(abs(Hz - obj.SpecMode.SigRegion(1)) == min(abs(Hz - obj.SpecMode.SigRegion(1)))),find(abs(Hz - obj.SpecMode.SigRegion(2)) == min(abs(Hz - obj.SpecMode.SigRegion(2))))]);
        obj.SpecMode.Scaling = max(abs(obj.spectrum(lims(1):lims(2),:))).*obj.SpecMode.Scaling;
        obj.CalcCurrentSpec;
        obj.DrawSpec;
    case 'Linewidth'
        if obj.SpecMode.Compl == 0 || obj.SpecMode.Compl == 4    % Real or both
            s = real(obj.spectrum);
        elseif obj.SpecMode.Compl == 1    % Imaginary
            s = imag(obj.spectrum);
        elseif obj.SpecMode.Compl == 3  % absolute
            s = abs(obj.spectrum);
        elseif obj.SpecMode.Compl == 2  % phase: no linewidth calculation possible
            sprintf('No linewidth calculation possible for phase!\n');
            return
        end
        peak = max(s);
        peakpos = find(s == peak);
        hp = peak/2;
        % Now start at the maximum and search until the first value below half is found
        p = peakpos;
        while s(p) > hp
            p = p+1;    % p will be the first point below 1/2
            if p > numel(obj.spectrum)
                printf('Cannot determine linewidth\n');
                return
            end
        end
       x = obj.Plot.Dat.x;
       val1 = (s(p-1)-hp)/(s(p-1)-s(p))*(x(p)-x(p-1))+x(p-1);
       % and now backwards
        p = peakpos;
        while s(p) > hp
            p = p-1;    % p will be the first point below 1/2
            if p < 1
                printf('Cannot determine linewidth\n');
                return
            end
        end
       val2 = (hp-s(p))/(s(p+1)-s(p))*(x(p+1)-x(p))+x(p);
       lw = val1 - val2;
       fprintf('Linewidth = %f (from %f to %f, max value = %f, half value = %f).\n',lw,val1,val2,hp,peak);
       if numel(s) == 2048
            noise = std(s(512:712));
            fprintf('Noise= %f, SNR= %f\n',noise,peak/noise);
       end
    case 'SetZero'
        if obj.SpecMode.FidSpec == 0
            return;
        end
        obj.DefineXPoint(1);
        obj.SpecMode.Def = 1;        
    case 'Phase1Ref'
        if obj.SpecMode.FidSpec == 0
            return;
        end
        obj.DefineXPoint(1);
        obj.SpecMode.Def = 2;
    case 'SigReg'
        if obj.SpecMode.FidSpec == 0
            return;
        end
        obj.DefineXPoint(2);
        obj.SpecMode.Def = 3;   
    case 'NoiseReg'
        obj.SpecMode.Def = 4;    % not implemented
    case 'ExportjMRUI'
        [f,p] = uiputfile([{'*.txt', ...
        'jMRUI datafile';'*.*','All Files'}], ... 
        'Set jMRUI-filename for writing')
        if f == 0
           return;
        end
        path = strcat(p,f);
        file = fopen(path,'w');
        if file == -1
            return
        end
        % Write header
        fprintf(file,'jMRUI Data Textfile\n\nFilename: %s\n\n',f);
        fprintf(file,'PointsInDataset: %d\n',obj.rawFIDSize);
        fprintf(file,'DatasetsInFile: 1\n');
        fprintf(file,'SamplingInterval: %5.4f\n',1/obj.datapars.bandwidth*1000);
        fprintf(file,'ZeroOrderPhase: 0E0\n');
        fprintf(file,'BeginTime: 0E0\n');
        fprintf(file,'TransmitterFrequency: %f\n', obj.datapars.frequency*1e6);
        fprintf(file,'MagneticField: %6.4E\n',obj.datapars.field);
        fprintf(file,'TypeOfNucleus: 1E0\n');
        fprintf(file,'NameOfPatient: hb\n');
        fprintf(file,'DateOfExperiment: 20150702\n');
        %fprintf(file,'DateOfExperiment: %s\n',data.ExpDate);
        fprintf(file,'Spectrometer: Investigational_Device_94T\n');
        fprintf(file,'AdditionalInfo: \n');
        fprintf(file,'\nSignal and FFT\n');
        fprintf(file,'sig(real)	sig(imag)	fft(real)	fft(imag)\n');
        %Write data
        fprintf(file,'Signal 1 out of 1 in file\n');
        %spec = circshift(fft(obj.rawFID),-obj.rawFIDSize/2);
        for cnt = 1:obj.rawFIDSize
                 %fprintf(file,'%11.4E %11.4E %11.4E %11.4E\n',real(obj.rawFID(cnt)),imag(obj.rawFID(cnt)),real(spec(cnt)),imag(spec(cnt))) ;
                 fprintf(file,'%11.4E %11.4E\n',real(obj.rawFID(cnt)),imag(obj.rawFID(cnt))) ;
        end

        fclose(file);
    case 'FastExport'        
        %assignin('caller','SpecWinFID',obj.proFID);
        %assignin('caller','SpecWinSpectrum',obj.spectrum);
        %assignin('caller','SpecWinMode',obj.SpecMode);
        assignin('base','SpecWinFID',obj.proFID);
        assignin('base','SpecWinSpectrum',obj.spectrum);
        assignin('base','SpecWinMode',obj.SpecMode);
        display(' Stored "SpecWinFID", "SpecWinSpectrum","SpecWinMode".');
    case 'zoomOut'
        PlotWindowCallback(obj, 0, 0, 'zoomOut');
        if obj.SpecMode.FidSpec == 1
             if obj.SpecMode.HzPpm == 0
                obj.SpecMode.SpecRange = obj.Plot.xRange;
             else
                 obj.SpecMode.SpecRange = obj.Plot.xRange*obj.datapars.frequency;
             end
        elseif obj.SpecMode.FidSpec == 0
            obj.SpecMode.FidRange = obj.Plot.xRange;
        end
end
end

function obj = SetSpec(obj, fid)
    obj.rawFID = fid;
    obj.rawFIDSize = numel(fid);
    obj.proFIDSize = numel(fid);
    obj.CalcCurrentSpec;
    if obj.SpecMode.FidSpec == 1    % show spectrum
        if obj.SpecMode.Compl <= 1  || obj.SpecMode.Compl == 4        % real or imaginary or both
            set(obj.Plot.PlotAxes,'YLim',[-max(abs(obj.spectrum)),max(abs(obj.spectrum))]);
        elseif obj.SpecMode.Compl == 3     % abs
            set(obj.Plot.PlotAxes,'YLim',[0,max(abs(obj.spectrum))]);
        elseif obj.SpecMode.Compl == 2
            set(obj.Plot.PlotAxes,'YLim',[-pi,pi]);
        end
    else     % show FID
        set(obj.Plot.PlotAxes,'YLimMode','auto');
    end
    obj.DrawSpec;
end

function obj = CalcCurrentSpec(obj)
    obj.proFID = obj.rawFID*exp(i*(obj.SpecMode.Phase0)).*exp(-linspace(0,obj.proFIDSize-1,obj.proFIDSize)'/obj.datapars.bandwidth*obj.SpecMode.Filter)./obj.SpecMode.Scaling;
    obj.spectrum = circshift(fft(obj.proFID,obj.SpecMode.Zerofill),-round(obj.SpecMode.Zerofill/2)).*exp(i*obj.SpecMode.Phase1/obj.datapars.frequency*(linspace(-obj.datapars.bandwidth/2,obj.datapars.bandwidth/2,obj.SpecMode.Zerofill)'-obj.SpecMode.Phase1Ref));%+data.refOffset-Mode.Phase1Ref));
end

function obj = DrawSpec(obj)
    set(obj.widgets.HelperPlot,'Visible','off');
    if obj.SpecMode.FidSpec == 0  % Show the FID
        obj.Plot.xReverse = 0;
        obj.Plot.xRange = obj.SpecMode.FidRange;
        obj.Plot.NewPlot(linspace(0,obj.proFIDSize,obj.proFIDSize)/obj.datapars.bandwidth*1000,obj.proFID);
    else                      % Show the spectrum
        obj.Plot.xReverse = 1;
        if obj.SpecMode.HzPpm == 0
            obj.Plot.xRange = obj.SpecMode.SpecRange;
            obj.Plot.NewPlot(linspace(-obj.datapars.bandwidth/2,obj.datapars.bandwidth/2,obj.SpecMode.Zerofill)-obj.SpecMode.Zero,obj.spectrum);
        elseif obj.SpecMode.HzPpm == 1
            obj.Plot.xRange = obj.SpecMode.SpecRange/obj.datapars.frequency;
            obj.Plot.NewPlot((linspace(-obj.datapars.bandwidth/2,obj.datapars.bandwidth/2,obj.SpecMode.Zerofill)-obj.SpecMode.Zero)/obj.datapars.frequency,obj.spectrum);
        end
        %obj.Plot.yDat(obj.spectrum);
    end
    obj.DrawPlot;
end
end
end

