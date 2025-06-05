classdef PlotClass < handle 
    % Class containing one or several 1D datasets
    % 
    
	
    properties
        NData                       % number of datasets
        Dat                        % data values on y axis as array of structures
        xRange                      % range of the plotted region in x direction
        yRange                      % range of the plotted region in y direction
        CurrentPlot                 % index of current plot
        yUnit                       % unit of y axis
        xUnit                       % unit of x axis
        FigWin                      % ID of Figure the plot is plotted on
        PlotAxes                    % ID of axes the plot is plotted on
        Title                       % Title of the entire plot
        Names                       % Names of the single plot elements
        ComplexMode                 % How to treat complex numbers: 0: real part, 1: imaginary, 2: phase, 3: abs, 4: real + imaginary 
        Plots                       % Line objects containing the actual plots
        xReverse                    % if 1 then x-axis is plotted in reverse direction
    end
    
    methods
        %% Konstruktor: Passes an array to form the plot
        function obj = PlotClass(xData,yData,name)   
            if nargin == 0      % No parameters passed: generate empty object
                obj.Dat(1).x = [];
                obj.Dat(1).y = [];
                obj.Dat(1).name = '';
                obj.Dat(1).Axis = 1;
                obj.yUnit = '';
                obj.xUnit = '';
                obj.NData = 0;
                obj.xRange = [0,1];
                obj.xReverse = 0;
            else
                if nargin == 1   % One parameters passed: y values: generated x values
                    yData = xData;
                    s = size(xData);
                    s = max(s(1:2));
                    xData = linspace(1,s,s);
                end
                obj = obj.NewPlot(xData,yData);
                obj.xRange = [obj.Dat(1).x(1),obj.Dat(1).x(end)];
            end
            if nargin < 3
                
            end
            obj.FigWin = 0;
            obj.PlotAxes = 0;
            obj.Title = 'Plot';
            obj.ComplexMode = 0;
            obj.Dat(1).Axis = 1;
            obj.xReverse = 0;
            return;
        end
        
        %% Method that sets the yDat
        function obj = yDat(obj, yData)
            if nargin <= 1
                return;
            end
            if isnumeric(yData)                    
                yData = squeeze(yData);
                s = size(yData);
                if numel(s) > 2 || (s(2)>1 && s(1)>1)
                    display('Only one graph, please!');
                    return;
                end
                if (s(1) == 1 && s(2) > 1)
                    yData = yData';
                    s = size(yData);
                end
                if s(1)~= numel(obj.Dat(obj.CurrentPlot).x)
                    display('The number of points mustn''t change!')
                    display(['There should be ',num2str(numel(obj.Dat(obj.CurrentPlot).x)),' points.']);
                    return
                end
                    obj.Dat(obj.CurrentPlot).y = yData;
            end
            return
        end 
 
        %% Method that sets the xDat
        function obj = xDat(obj, xData)
            if nargin <= 1
                return;
            end
            if isnumeric(xData)                    
                xData = squeeze(xData);
                s = size(xData);
                if numel(s) > 2 || (s(2)>1 && s(1)>1)
                    display('Only one graph, please!');
                    return;
                end
                if (s(1) == 1 && s(2) > 1)
                    xData = xData';
                    s = size(xData);
                end
                if s(1)~= numel(obj.Dat(obj.CurrentPlot).y)
                    display('The number of points mustn''t change!')
                    display(['There should be ',num2str(numel(obj.Dat(obj.CurrentPlot).y)),' points.']);
                    return
                end
                    obj.Dat(obj.CurrentPlot).x = xData;
            end
            return
        end 
        
        
        %% Method to generate a new plot
        function obj = NewPlot(obj,xData, yData)
            if nargin <= 1
                return;
            end
            if nargin == 2
                yData = xData;
                xData = [];
            end
            if isnumeric(yData)                    
                yData = squeeze(yData);
                s = size(yData);
                sx = size(xData);
                if numel(s) > 2 || (s(1) <3 && s(2) < 3) 
                    display('These are no 1D datasets!');
                    return;
                end
                if (s(1) == 1 && s(2) > 1)
                    yData = yData';
                    s = size(yData);
                end
                if (sx(1) == 1 && sx(2) > 1)
                    xData = xData';
                    sx = size(xData);
                end
                if numel(xData)~=0 && (sx(1)~=s(1) || sx(2) > s(2) || (sx(2) < s(2) && sx(2) ~= 1))
                    display('x and y dimensions don''t agree!')
                    return
                end
                if s(2) < obj.NData
                    obj.Dat = [];
                end
                for cnt = 1:s(2)
                    obj.Dat(cnt).y = yData(:,cnt);
                    obj.Dat(cnt).Axis = 1;
                    obj.Dat(cnt).name = '';
                    if numel(xData) == 0
                        obj.Dat(cnt).x = linspace(0,s(1),s(1));
                    else
                        if sx(2)==1
                            obj.Dat(cnt).x = xData;
                        else
                            obj.Dat(cnt).x = xData(:,cnt);
                        end
                    end
                end                
                obj.NData = s(2);
                obj.CurrentPlot = 1;
            end
        end 

        function obj = AddPlot(obj,xData, yData, MakeCurrent)
            if nargin <= 1
                return;
            end
            if nargin == 2
                yData = xData;
                xData = [];
            end
            if nargin < 4
                MakeCurrent = 1;
            end
            if isnumeric(yData)                    
                yData = squeeze(yData);
                s = size(yData);
                sx = size(xData);
                if numel(s) > 2 || (s(1) <3 && s(2) < 3) 
                    display('These are no 1D datasets!');
                    return;
                end
                if (s(1) == 1 && s(2) > 1)
                    yData = yData';
                    s = size(yData);
                end
                if (sx(1) == 1 && sx(2) > 1)
                    xData = xData';
                    sx = size(xData);
                end
                if numel(xData)~=0 && (sx(1)~=s(1) || sx(2) > s(2) || (sx(2) < s(2) && sx(2) ~= 1))
                    display('x and y dimensions don''t agree!')
                    return
                end
                for cnt = obj.NData+1:obj.NData+s(2)
                    obj.Dat(cnt).y = yData;
                    obj.Dat(cnt).Axis = 1;
                    if numel(xData) == 0
                        obj.Dat(cnt).x = linspace(0,s(1),s(1));
                    else
                        if sx(2)==1
                            obj.Dat(cnt).x = xData;
                        else
                            obj.Dat(cnt).x = xData(:,cnt);
                        end
                    end
                end
                obj.NData = obj.NData + s(2);
                if MakeCurrent == 1
                    obj.CurrentPlot = obj.NData;
                end
            end
            return
        end 

        function obj = RemovePlot(obj,plotnr)
            if nargin <= 1 || isnumeric(plotnr)~=1
                return;
            end
            if plotnr <= 0 || plotnr > obj.NData
                return;
            end
            if plotnr == 1
                obj.Dat = obj.Dat(2:end);
            elseif plotnr == obj.NData
                obj.Dat = obj.Dat(1:end-1);
            else
                obj.Dat = [obj.Dat(1:plotnr-1),obj.Dat(plotnr+1:end)];
            end
            obj.NData = obj.NData-1;
            if (obj.CurrentPlot == plotnr)
                if plotnr == 1 && obj.NData > 0
                    obj.CurrentPlot = obj.CurrentPlot+1;
                elseif plotnr > 1 
                    obj.CurrentPlot = obj.CurrentPlot - 1;
                elseif obj.NData == 0
                    obj.CurrentPlot = 0;
                end
            end
            return
        end 
        
        %% Move current plot to the second axis
        function obj = SetAxis(obj, ax)
            if isnumeric(ax)==0 || (ax~=1 && ax ~=2)
                display('The axis must be 1 or 2!');
            end
            obj.Dat(obj.CurrentPlot).Axis = ax;
            return
        end

  
        %% Method to set ComplexMode to real
        function obj = showreal(obj)
            obj.ComplexMode = 0;
        end
        %% Method to set ComplexMode to imaginary
        function obj = showimag(obj)
            obj.ComplexMode = 1;
        end
        %% Method to set ComplexMode to phase
        function obj = showphase(obj)
            obj.ComplexMode = 2;
        end
        %% Method to set ComplexMode to abs
        function obj = showabs(obj)
            obj.ComplexMode = 3;
        end
        %% Method to set ComplexMode to real + imaginary
        function obj = showboth(obj)
            obj.ComplexMode = 4;
        end
        
        
        %% Method to plot the plot in a new figure window
        function obj = plot(obj, plotnr)
            if obj.NData == 0
                return
            end
            if obj.PlotAxes == 0 | ishandle(obj.PlotAxes) == 0   % ??? isvalid(obj.PlotAxes) == 0
                obj.FigWin = figure('Name',obj.Title);
                obj.PlotAxes = axes;
            else
                figure(obj.FigWin);
            end            
            if obj.xReverse == 1
                rev = 'reverse';
            else
                rev = 'normal';
            end
            set(obj.PlotAxes,'XLim',obj.xRange, 'XDir',rev);
            if nargin < 2
                plotnr = linspace(1, obj.NData, obj.NData);
            end
            nplots = numel(plotnr);
            col = 'kbrgmc';
            while nplots > numel(col)
                col = [col,col];
            end
            % Destroy all old line objects
            for cnt = 1:numel(obj.Plots)             
                    delete(obj.Plots(cnt));                    
            end 

            obj.Plots = matlab.graphics.primitive.Line.empty;
            if strcmp(get(obj.PlotAxes,'YAxisLocation'),'right')
                yyaxis left
            end
            currax = 1;
            plotorder = linspace(1,nplots,nplots);
            plotorder(obj.CurrentPlot:end) = circshift(plotorder(obj.CurrentPlot:end),-1);
            for cnt = plotorder 
                if max([obj.Dat(:).Axis]) >1
                    if obj.Dat(cnt).Axis == 1
                        yyaxis left
                    else
                        yyaxis right
                    end
                end
                if isreal(obj.Dat(cnt).y)
                    obj.Plots(cnt) = line('Parent',obj.PlotAxes,'YData',obj.Dat(plotnr(cnt)).y,'XData',obj.Dat(plotnr(cnt)).x,'Color',col(cnt),'visible','on');
                else
                    switch obj.ComplexMode
                        case 0                           
                            obj.Plots(cnt) = line('Parent',obj.PlotAxes,'YData',real(obj.Dat(plotnr(cnt)).y),'XData',obj.Dat(plotnr(cnt)).x,'Color',col(cnt),'visible','on');
                        case 1
                            obj.Plots(cnt) = line('Parent',obj.PlotAxes,'YData',imag(obj.Dat(plotnr(cnt)).y),'XData',obj.Dat(plotnr(cnt)).x,'Color',col(cnt),'visible','on');
                        case 2
                            obj.Plots(cnt) = line('Parent',obj.PlotAxes,'YData',angle(obj.Dat(plotnr(cnt)).y),'XData',obj.Dat(plotnr(cnt)).x,'Color',col(cnt),'visible','on');
                        case 3
                            obj.Plots(cnt) = line('Parent',obj.PlotAxes,'YData',abs(obj.Dat(plotnr(cnt)).y),'XData',obj.Dat(plotnr(cnt)).x,'Color',col(cnt),'visible','on');
                        case 4
                            obj.Plots(cnt,1) = line('Parent',obj.PlotAxes,'YData',real(obj.Dat(plotnr(cnt)).y),'XData',obj.Dat(plotnr(cnt)).x,'Color',col(cnt),'visible','on');
                            obj.Plots(cnt,2) = line('Parent',obj.PlotAxes,'YData',imag(obj.Dat(plotnr(cnt)).y),'XData',obj.Dat(plotnr(cnt)).x,'Color',col(cnt),'Linestyle','--','visible','on');
                    end
                    hold on 
                end
            end
            hold off
        end
        
        %% Method to fouriertransform all or single plots of the object
        function obj = fft(obj, plotnr)
            if obj.NData == 0
                return
            end
           
            if nargin < 2   % Fouriertransform all plots
                plotnr = linspace(1,obj.NData, obj.nData);
            end
            for cnt = 1:numel(plotnr)
                obj.Dat(cnt).y = fftshift(fft(obj.Dat(cnt).y),1);
            end
            
            return
        end

        %% Method to zoom in 
        function Zoom(obj, limits)
            if obj.NData == 0
                return
            end
            if nargin < 2   % if there are no limits given: zoom out
                obj.xRange = [obj.Dat(1).x(1),obj.Dat(1).x(end)];
            else
                if numel(limits) == 1 % If there is only one number, it is assumed that it is the right end
                    if limits > obj.xRange(1)
                        obj.xRange(2) = min([limits,obj.Dat(1).x(end)]);
                    end
                else
                    if limits(2) > limits(1)
                    obj.xRange = [max([limits(1),obj.Dat(1).x(1)]),min([limits(2),obj.Dat(1).x(end)])];
                    end
                end
            end
            return
        end
        
        %% Method to reverse the x-axis direction
        function reverse(obj)
            if obj.xReverse == 0    % the axis is currently not reversed
                obj.xReverse = 1;
                set(obj.PlotAxes, 'XDir','reverse');
            else
                obj.xReverse = 0;
                set(obj.PlotAxes, 'XDir','normal');
            end
        end

        %% Method to list all methods of the class
        function help(obj)
            display('fft(obj, plotnr):      Fouriertransforms the plot number plotnr (if missing: all plots in the object)');
            display('plot(obj, plotnr):     Plots the plot number plotnr in the plot-figure. plotnr can be an array to plot several plots in the same window.');
            display('                       If plotnr is missing, all plots of the object are plotted.');
        end
    end     
end

