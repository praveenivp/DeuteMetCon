classdef PlotWindow < handle 
    % Window for watching plots
    %   Detailed explanation goes here
    
    properties
        Plot
        WindowSize
        CurrentPlot
        widgets
        mode
        files
        
        
    end
    
    methods
        function obj = PlotWindow(PlotObject, WindowSize, baseWidget, varargin)
            % Interpret varargin
            if nargin > 3
                for cnt = 1:numel(varargin)
                    if ischar(varargin{cnt})
                        if contains(varargin{cnt},'ExportPath','IgnoreCase',true)
                            if numel(varargin) > cnt
                                obj.files.ExportPath = varargin{cnt+1};
                            end
                        elseif contains(varargin{cnt},'ExportFileExcel','IgnoreCase',true)
                            if numel(varargin) > cnt
                                obj.files.ExcelFileName = varargin{cnt+1};
                            end
                        elseif contains(varargin{cnt},'ExportFileMat','IgnoreCase',true)
                            if numel(varargin) > cnt
                                obj.files.MatFileName = varargin{cnt+1};
                            end
                        elseif contains(varargin{cnt},'PlotNames','IgnoreCase',true)
                            if numel(varargin) > cnt && numel(varargin{cnt+1})>0
                                plotnames = varargin{cnt+1};
                            end
                        end
                    end
                end
            end
            obj.mode.Click = 0;
            obj.CurrentPlot = 1;
            obj.mode.define = 0;
            obj.mode.clicked = [];
            if nargin < 2
                WindowSize = [700,400];
            end
            
            if numel(WindowSize)==1
                WindowSize = [WindowSize,WindowSize];
            end
            if numel(WindowSize)<4
                Pos = [100,100,WindowSize(1), WindowSize(2)];
            else
                Pos = WindowSize;
            end
            obj.WindowSize = [Pos(3),Pos(4)];
            if nargin > 0
            if isobject(PlotObject) && strcmp(class(PlotObject),'PlotClass')
                obj.Plot = PlotObject;
            else
                obj.Plot = PlotClass(PlotObject);
                if isobject(obj.Plot) == 0
                    return
                end
            end
            end
            if nargin < 3
                obj.widgets.base = figure('Name','Evaluate Series','NumberTitle','off','DockControls','off','Pointer','crosshair',...
                    'Toolbar','none','Menubar','none','Unit','Pixels','Position',Pos,'Color',[0.8,0.8,0.8], 'WindowButtonUpFcn',{@obj.PlotWindowCallback,'ClickUp'},'NextPlot','new','Tag','rPlotWindow','UserData',obj);%,...
                wpos = [0,0,Pos(3),Pos(4)];
            else
                obj.widgets.base = baseWidget;
                set(obj.widgets.base,'Tag','rpPlotWindow_res');
                if numel(obj.widgets.base.UserData) == 0
                    set(obj.widgets.base,'UserData',obj);
                end
                wpos = Pos;
            end
            set(obj.widgets.base, 'WindowButtonUpFcn',{@obj.PlotWindowCallback,'ClickUp'});
            obj.Plot.FigWin = obj.widgets.base;
            obj.widgets.main = uipanel(obj.widgets.base,'Unit','Pixels','Position',wpos,'BackgroundColor',[0.8,0.8,0.8]);
            obj.widgets.TextData = uicontrol(obj.widgets.main,'Style','text','Unit','Pixels','Position',[240,obj.WindowSize(2)-20,50,15],'BackgroundColor',[0.8,0.8,0.8], 'String','Data:', 'HorizontalAlignment','left');
            obj.widgets.LabelData = uicontrol(obj.widgets.main,'Style','text','Unit','Pixels','Position',[275,obj.WindowSize(2)-20,150,15],'BackgroundColor',[0.8,0.8,0.8], 'String','               ', 'HorizontalAlignment','left');
            obj.widgets.TextCursor = uicontrol(obj.widgets.main,'Style','text','Unit','Pixels','Position',[10,obj.WindowSize(2)-20,45,15],'BackgroundColor',[0.8,0.8,0.8], 'String','Cursor:', 'HorizontalAlignment','left');
            obj.widgets.LabelCursor = uicontrol(obj.widgets.main,'Style','text','Unit','Pixels','Position',[55,obj.WindowSize(2)-20,180,15],'BackgroundColor',[0.8,0.8,0.8], 'String','              ', 'HorizontalAlignment','left');
            
            axespos = [50,30,obj.WindowSize(1)-90, obj.WindowSize(2)-70];
            obj.widgets.PlotAxes = axes('Parent',obj.widgets.main,'Unit','Pixels','Position',axespos,'Box','On','XTickLabelMode','auto','YTickLabelMode','auto');
            obj.Plot.PlotAxes = obj.widgets.PlotAxes;
            obj.DrawPlot;
            
            obj.widgets.MarkerDot = line('Parent',obj.Plot.PlotAxes,'XData',0,'YData',0,'Visible','off','Color','r','Marker','o');
            obj.widgets.MarkerLine = line('Parent',obj.Plot.PlotAxes,'XData',[0,0],'YData',[-1,1],'Visible','off','Color','g');
            obj.widgets.MarkerLine2 = line('Parent',obj.Plot.PlotAxes,'XData',[0,0],'YData',[-1,1],'Visible','off','Color','g');
            obj.widgets.HelperPlot = line('Parent',obj.Plot.PlotAxes,'XData',[0,0],'YData',[-1,1],'Visible','off','Color','r');
            names = string(linspace(1,numel(obj.Plot.Plots),numel(obj.Plot.Plots))');
            if exist("plotnames","var") 
                names(1:numel(plotnames)) = plotnames;
            end
            obj.widgets.CurrPlot = uicontrol(obj.widgets.main,'Style','popupmenu','Unit','Pixels','Position',[obj.WindowSize(1)-70,obj.WindowSize(2)-40,55,35],'BackgroundColor',[0.8,0.8,0.8], 'String',names, 'HorizontalAlignment','right','Value',obj.CurrentPlot, 'Callback',{@obj.PlotWindowCallback,'SelCurrPlot'});
            obj.widgets.cm = uicontextmenu;
            obj.widgets.cm0 = uimenu(obj.widgets.cm,'Label','Zoom out','Callback',{@obj.PlotWindowCallback,'zoomOut'});
            obj.widgets.cm1 = uimenu(obj.widgets.cm,'Label','Plot');
            uimenu(obj.widgets.cm1,'Label','Remove','Callback',{@obj.PlotWindowCallback,'PlotRemove'});
            obj.widgets.cm2 = uimenu(obj.widgets.cm,'Label','Range');
            uimenu(obj.widgets.cm2,'Label','Dynamic y-range','Callback',{@obj.PlotWindowCallback,'dynamicyRange'});
            uimenu(obj.widgets.cm2,'Label','Set y range','Callback',{@obj.PlotWindowCallback,'SetyRange'});
            obj.widgets.cm3 = uimenu(obj.widgets.cm,'Label','File');
            uimenu(obj.widgets.cm3,'Label','Export Data','Callback',{@obj.PlotWindowCallback,'PlotExport'});
            uimenu(obj.widgets.cm3,'Label','Export Image','Callback',{@obj.PlotWindowCallback,'ImageExport'});
            set(obj.widgets.main,'UIContextMenu',obj.widgets.cm);
            set(obj.widgets.PlotAxes,'UIContextMenu',obj.widgets.cm);
            set(obj.widgets.base, 'WindowButtonMotionFcn',{@obj.PlotWindowCallback,'Motion'},'CurrentAxes',obj.widgets.PlotAxes,'KeyPressFcn',{@obj.PlotWindowCallback,'Key'},'ResizeFcn',{@obj.PlotWindowCallback,'Resize'});
            set(obj.widgets.PlotAxes, 'ButtonDownFcn',{@obj.PlotWindowCallback,'Click'});
            set(obj.Plot.Plots(1),'ButtonDownFcn', {@obj.PlotWindowCallback,'SelectPlot'});
            % Set variable files
            obj.files.ExportPath = cd;
            obj.files.ExportDataType = 'var';
            % Set default filenames
            obj.files.ExcelFileName = '';
            obj.files.MatFileName = '';

        end
        
        function obj = PlotWindowCallback(obj, src, evt, action)
            
            switch action
                case 'Motion'
                     set(obj.widgets.base,'CurrentAxes',obj.widgets.PlotAxes);
                     p = obj.GetCurrentPoint;
                     set(obj.widgets.LabelCursor, 'String',[num2str(p(1),6),', ',num2str(p(2),6)]);
                     xl = get(obj.widgets.PlotAxes,'XLim');
                     if p(1)>xl(1) & p(1)<xl(2)           
                         xpoint = find(abs(obj.Plot.Dat(obj.CurrentPlot).x - p(1)) == min(abs(obj.Plot.Dat(obj.CurrentPlot).x - p(1))));
                         set(obj.widgets.LabelData, 'String',[num2str(obj.Plot.Dat(obj.CurrentPlot).x(xpoint(1)),'%6.4g'),', ',num2str(obj.Plot.Dat(obj.CurrentPlot).y(xpoint(1)),'%6.4g')]);                         
                         if obj.Plot.ComplexMode == 1    % imaginary
                             set(obj.widgets.MarkerDot,'XData',obj.Plot.Dat(obj.CurrentPlot).x(xpoint(1)),'YData',imag(obj.Plot.Dat(obj.CurrentPlot).y(xpoint(1))),'Visible','on');
                         elseif obj.Plot.ComplexMode == 2    % phase
                             set(obj.widgets.MarkerDot,'XData',obj.Plot.Dat(obj.CurrentPlot).x(xpoint(1)),'YData',angle(obj.Plot.Dat(obj.CurrentPlot).y(xpoint(1))),'Visible','on');                             
                         elseif obj.Plot.ComplexMode == 3    % abs
                             set(obj.widgets.MarkerDot,'XData',obj.Plot.Dat(obj.CurrentPlot).x(xpoint(1)),'YData',abs(obj.Plot.Dat(obj.CurrentPlot).y(xpoint(1))),'Visible','on');                             
                         else    % real part or both
                           set(obj.widgets.MarkerDot,'XData',obj.Plot.Dat(obj.CurrentPlot).x(xpoint(1)),'YData',real(obj.Plot.Dat(obj.CurrentPlot).y(xpoint(1))),'Visible','on');
                         end
                     else
                         set(obj.widgets.MarkerDot,'Visible','off');
                     end
                     if obj.mode.Click == 1  || obj.mode.define == 1
                        set(obj.widgets.MarkerLine,'XData',[p(1),p(1)]);
                     elseif obj.mode.Click == 3  || obj.mode.define == 2
                        set(obj.widgets.MarkerLine2,'XData',[p(1),p(1)]); 
                     end

                case 'Click'
                      cp = obj.GetCurrentPoint();
                      if evt.Button == 1          % left click: select a plot
                          qqq = 1;
                      elseif evt.Button == 2       % zoom
                        if obj.mode.Click == 0
                            obj.mode.Click = 1;
                            ydat = get(obj.Plot.PlotAxes.YAxis(obj.Plot.Dat(obj.CurrentPlot).Axis),'Limits');
                            set(obj.widgets.MarkerLine,'Parent',obj.Plot.PlotAxes,'Visible','on','XData',[cp(1),cp(1)],'YData',ydat);                            
                        elseif obj.mode.Click == 2
                            obj.mode.Click = 3;
                            ydat = get(obj.Plot.PlotAxes,'YLim');
                            set(obj.widgets.MarkerLine2,'Parent',obj.Plot.PlotAxes,'Visible','on','XData',[cp(1),cp(1)],'YData',ydat); 
                        end
                      end
                case 'ClickUp' 
                      if obj.mode.Click == 1
                          obj.mode.Click = 2;
                      elseif obj.mode.Click == 3
                            obj.mode.Click = 0;
                            %set(obj.widgets.MarkerLine,'Visible','off','Parent',[]);
                            %set(obj.widgets.MarkerLine2,'Visible','off','Parent',[]);
                            set(obj.widgets.MarkerLine,'Visible','off');
                            set(obj.widgets.MarkerLine2,'Visible','off');
                            p = obj.GetCurrentPoint;
                            pp = get(obj.widgets.MarkerLine,'XData');
                            if p(1) > pp(1)
                                ppp = p(1);
                                p = pp(1);
                                pp = ppp;
                            end
                            obj.Plot.xRange=[p(1),pp(1)];
                            obj.DrawPlot;
                      elseif obj.mode.define > 0    % define a point
                              p = obj.GetCurrentPoint;
                              obj.mode.define = obj.mode.define - 1;
                              obj.mode.clicked = [obj.mode.clicked,p];    % returns the selected x-point
                              if obj.mode.define == 0
                                  set(obj.widgets.MarkerLine,'Visible',0);
                                  set(obj.widgets.MarkerLine2,'Visible',0);
                              else
                                  set(obj.widgets.MarkerLine,'Visible',1);
                                  ydat = get(obj.Plot.PlotAxes.YAxis(obj.Plot.Dat(obj.CurrentPlot).Axis),'Limits');
                                  set(obj.widgets.MarkerLine,'Visible','on','XData',[p(1),p(1)],'YData',ydat);
                              end
                          
                      end
                case 'SelectPlot'
                    if evt.Button == 1     % Click on a plot to select it
                        obj.CurrentPlot = find(obj.Plot.Plots == evt.Source);
                        set(obj.widgets.CurrPlot,'Value',obj.CurrentPlot);
                        if obj.Plot.Dat(obj.CurrentPlot).Axis == 2
                            yyaxis right
                        else
                            yyaxis left
                        end
                        set(obj.widgets.MarkerDot,'Parent',[]);
                        set(obj.widgets.MarkerDot,'Parent',obj.Plot.PlotAxes);
                        PlotWindowCallback(obj, src, evt, 'Motion');
                        return
                    end
                    if evt.Button == 2     % Middle mouse button: zoom in
                        obj = PlotWindowCallback(obj, src, evt, 'Click');
                        return
                    end
                case 'SelCurrPlot'
                        obj.CurrentPlot = evt.Source.Value;
                        obj.Plot.CurrentPlot = obj.CurrentPlot;
                        if obj.Plot.Dat(obj.CurrentPlot).Axis == 2
                            yyaxis right
                        else
                            yyaxis left
                        end
                        % Move the marker dot to the correct axis
                        set(obj.widgets.MarkerDot,'Parent',[]);
                        set(obj.widgets.MarkerDot,'Parent',obj.Plot.PlotAxes);
                        % Plot the selected Plot at the top
                        set(obj.Plot.Plots(obj.CurrentPlot),'Parent',[]);
                        set(obj.Plot.Plots(obj.CurrentPlot),'Parent',obj.Plot.PlotAxes);
                        %set(obj.widgets.CurrPlot,'Value',obj.CurrentPlot);                        
                        obj.DrawPlot;
                        PlotWindowCallback(obj, src, evt, 'Motion');
                        return
                case 'Key'
                case 'Resize'
                    p = get(obj.widgets.base,'Position');
                    s = p(3:4);
                    obj.WindowSize(1) = s(1);
                    obj.WindowSize(2) = s(2);
                    set(obj.widgets.main, 'Position',[0,0, obj.WindowSize(1),  obj.WindowSize(2)]);
                    obj.PlotWindowCallback(src,evt,'PanelResize');
                case 'PanelResize'
                    p = get(obj.widgets.main,'Position');
                    s = p(3:4);
                    obj.WindowSize(1) = s(1);
                    obj.WindowSize(2) = s(2);
                    set(obj.widgets.main, 'Position',p);
                    set(obj.widgets.PlotAxes,'Position',[50,30,obj.WindowSize(1)-90, obj.WindowSize(2)-70]);
                    set(obj.widgets.LabelData,'Position',[45,obj.WindowSize(2)-20,180,15]);
                    set(obj.widgets.LabelCursor,'Position',[285,obj.WindowSize(2)-20,120,15]);
                    set(obj.widgets.TextData,'Position',[10,obj.WindowSize(2)-20,35,15]);
                    set(obj.widgets.TextCursor,'Position',[240,obj.WindowSize(2)-20,45,15]);
                    set(obj.widgets.CurrPlot, 'Position',[obj.WindowSize(1)-40,obj.WindowSize(2)-40,35,35]);

                case 'zoomOut'
                    obj.Plot.Zoom();
                    obj.DrawPlot;
                case 'PlotRemove'
                    if obj.Plot.NData == 1
                        return;
                    end
                    obj.Plot.RemovePlot(obj.CurrentPlot);
                    if obj.CurrentPlot > obj.Plot.NData
                        obj.CurrentPlot = obj.CurrentPlot - 1;
                    end
                    obj.DrawPlot;
                    obj.widgets.CurrPlot.String = num2str(linspace(1,obj.Plot.NData,obj.Plot.NData)');
                    obj.widgets.CurrPlot.Value = obj.CurrentPlot;

                case 'SetyRange'
                    if max([obj.Plot.Dat(:).Axis]) > 1   % there are two axes
                        %yyaxis left
                        leftlims = get(obj.Plot.PlotAxes.YAxis(1), 'Limits');
                        %yyaxis right
                        rightlims = get(obj.Plot.PlotAxes.YAxis(2), 'Limits');
                        pos = get(obj.widgets.base,'Position');
                        pos(1:2) = pos(1:2) + pos(3:4)/2;
                        pos(3) = 240;
                        pos(4) = 100;
                        rangefig = figure('Name','Set y-range','NumberTitle','off','DockControls','off','Pointer','crosshair',...
                        'Toolbar','none','Menubar','none','Unit','Pixels','Position',pos,'Color',[0.8,0.8,0.8]);
                        OKbutton = uicontrol('Parent',rangefig, 'Style', 'pushbutton','Unit','Pixels','Position',[10,10,60,20],'String','OK'); 
                        cancelbutton = uicontrol('Parent',rangefig, 'Style', 'pushbutton','Unit','Pixels','Position',[170,10,60,20],'String','Cancel');
                        commonzerobutton = uicontrol('Parent',rangefig, 'Style', 'pushbutton','Unit','Pixels','Position',[85,10,70,20],'String','Common Zero');
                        label = uicontrol('Parent',rangefig, 'Style', 'text','Unit','Pixels','Position',[10,30,30,20],'String','Min','Backgroundcolor',[0.8,0.8,0.8]);
                        label = uicontrol('Parent',rangefig, 'Style', 'text','Unit','Pixels','Position',[10,52,30,20],'String','Max','Backgroundcolor',[0.8,0.8,0.8]);
                        label = uicontrol('Parent',rangefig, 'Style', 'text','Unit','Pixels','Position',[20,75,50,20],'String','left axis:','Backgroundcolor',[0.8,0.8,0.8]);
                        textminleft = uicontrol('Parent',rangefig, 'Style', 'edit','Unit','Pixels','Position',[50,32,60,20],'String',num2str(leftlims(1)));
                        textmaxleft = uicontrol('Parent',rangefig, 'Style', 'edit','Unit','Pixels','Position',[50,54,60,20],'String',num2str(leftlims(2)));
                        label = uicontrol('Parent',rangefig, 'Style', 'text','Unit','Pixels','Position',[120,30,30,20],'String','Min','Backgroundcolor',[0.8,0.8,0.8]);
                        label = uicontrol('Parent',rangefig, 'Style', 'text','Unit','Pixels','Position',[120,52,30,20],'String','Max','Backgroundcolor',[0.8,0.8,0.8]);
                        label = uicontrol('Parent',rangefig, 'Style', 'text','Unit','Pixels','Position',[140,75,50,20],'String','right axis:','Backgroundcolor',[0.8,0.8,0.8]);
                        textminright = uicontrol('Parent',rangefig, 'Style', 'edit','Unit','Pixels','Position',[160,32,60,20],'String',num2str(rightlims(1)));
                        textmaxright = uicontrol('Parent',rangefig, 'Style', 'edit','Unit','Pixels','Position',[160,54,60,20],'String',num2str(rightlims(2)));
                        set(OKbutton, 'Callback',{@obj.PlotWindowCallback,'yRangeOkay'});
                        set(cancelbutton, 'Callback',{@obj.PlotWindowCallback,'yRangeCancel'});
                        set(commonzerobutton, 'Callback',{@obj.PlotWindowCallback,'yRangeCommon0'});
                        waitfor(OKbutton);
                        if isvalid(rangefig)    % OK  pressed
                            %yyaxis left
                            set(obj.Plot.PlotAxes.YAxis(1),'Limits',[str2num(get(textminleft,'String')),str2num(get(textmaxleft,'String'))]);
                            %yyaxis right
                            set(obj.Plot.PlotAxes.YAxis(2),'Limits',[str2num(get(textminright,'String')),str2num(get(textmaxright,'String'))]);
                            delete(rangefig);
                        else                    % Cancel pressed
                            return
                        end
                        
                    else                                 % there is only one axis
                        oldlims = get(obj.Plot.PlotAxes, 'YLim');
                        limitstr = inputdlg({'Min','Max'},'y-Axis limits',1,strtrim(cellstr(num2str(oldlims'))));
                        if numel(limitstr) == 0
                            return;
                        end
                        limits = cellfun(@str2num,limitstr);
                        set(obj.Plot.PlotAxes, 'YLim',limits);
                    end
                case 'PlotExport'
                    obj.files = ExportData(obj.Plot.Dat,obj.files.ExportPath, obj.files); 
                case 'ImageExport'
                    ExportImage(obj.Plot.PlotAxes) 
                case 'dynamicyRange'
                    if max([obj.Plot.Dat(:).Axis]) > 1   % there are two axes
                        yyaxis right
                        set(obj.Plot.PlotAxes, 'YLimMode','auto');
                        %yyaxis left
                    end
                    set(obj.Plot.PlotAxes, 'YLimMode','auto');
                case 'yRangeOkay'
                    delete(src);
                case 'yRangeCancel'
                    delete(get(src,'Parent'));
                case 'yRangeCommon0'
                    ws = get(get(src,'Parent'),'Children');
                    leftlims = [str2num(get(ws(7),'String')),str2num(get(ws(6),'String'))];
                    rightlims = [str2num(get(ws(2),'String')),str2num(get(ws(1),'String'))];
                    leftminmax = leftlims(1)/leftlims(2);
                    rightminmax = rightlims(1)/rightlims(2);
                    leftnewlims = [min(rightminmax*leftlims(2), leftlims(1)),max(leftlims(1)/rightminmax, leftlims(2))];
                    rightnewlims = [min(leftminmax*rightlims(2), rightlims(1)),max(rightlims(1)/leftminmax, rightlims(2))];
                    leftnewminmax = leftnewlims(1)/leftnewlims(2);
                    rightnewminmax = rightnewlims(1)/rightnewlims(2);
                    if rightnewminmax ~= leftminmax
                        leftlims = leftnewlims;
                    elseif leftnewminmax ~= rightminmax
                        rightlims = rightnewlims;    
                    else
                        rightlims = rightnewlims;
                    end
                    set(ws(7),'String',num2str(leftlims(1)));
                    set(ws(6),'String',num2str(leftlims(2)));
                    set(ws(2),'String',num2str(rightlims(1)));
                    set(ws(1),'String',num2str(rightlims(2)));
                    
            end
        end
        
        function obj = DefineXPoint(obj, NumPoints)
            cp = obj.GetCurrentPoint();
            ydat = get(obj.Plot.PlotAxes.YAxis(obj.Plot.Dat(obj.CurrentPlot).Axis),'Limits');
            if NumPoints == 1
                set(obj.widgets.MarkerLine,'Visible','on','XData',[cp(1),cp(1)],'YData',ydat);
            elseif NumPoints == 2
                set(obj.widgets.MarkerLine2,'Visible','on','XData',[cp(1),cp(1)],'YData',ydat);
            end
            obj.mode.define = NumPoints;
            obj.mode.clicked = [];
        end
        
        function p = GetCurrentPoint(obj)
             %XData = obj.Plot.Dat(obj.CurrentPlot).x;
             %YData = obj.Plot.Dat(obj.CurrentPlot).y;
%              if max([obj.Plot.Dat(:).Axis]) > 1   % there are two axes
%                  if obj.Plot.Dat(obj.CurrentPlot).Axis == 1
%                      yyaxis left;
%                  else
%                      yyaxis right;
%                  end
%              end
             p = get(obj.widgets.PlotAxes,'CurrentPoint');
             p = squeeze(p(1,1:2));        
            return
        end
        
        
        function AddPlot(obj,xData, yData, MakeCurrent)
            if nargin < 4
                MakeCurrent = 1;
            end
            obj.Plot.AddPlot(xData, yData, MakeCurrent);
            obj.DrawPlot;
            set(obj.widgets.CurrPlot,'String',num2str(linspace(1,obj.Plot.NData,obj.Plot.NData)'),'Value',obj.CurrentPlot);
        end
        
        function DrawPlot(obj)
            obj.Plot.plot();
            for cnt = 1:obj.Plot.NData
                set(obj.Plot.Plots(cnt),'ButtonDownFcn',{@obj.PlotWindowCallback,'SelectPlot'});
            end
%             obj.widgets.MarkerDot = line('Parent',obj.Plot.PlotAxes,'XData',0,'YData',0,'Visible','off','Color','r','Marker','o');
%             obj.widgets.MarkerLine = line('Parent',obj.Plot.PlotAxes,'XData',[0,0],'YData',[-1,1],'Visible','off','Color','g');
%             obj.widgets.MarkerLine2 = line('Parent',obj.Plot.PlotAxes,'XData',[0,0],'YData',[-1,1],'Visible','off','Color','g');
        end
        

        

    end
    
end

