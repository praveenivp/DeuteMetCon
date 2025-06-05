function data = ImportData(modal, defaultName)
    % Imports data of different kinds.
    % Currently available data types:
    %   - Import Matlab variable from base workspace

    data = [];
    if nargin < 1
        modal = 0;
    end
    if modal == 1
        WinStyle = 'modal';
    else
        WinStyle = 'normal';
    end
    % Get variables in basic workspace
    basevars = evalin('base','who');

    %% Generate widget
    widgets.base = figure('Name','Import Data','NumberTitle','off','DockControls','off','WindowStyle',WinStyle,...
                        'Toolbar','none','Menubar','none','Unit','Pixels','Position',[100,100,600, 150],'Color',[0.8,0.8,0.8]);
    widgets.DataName = uicontrol(widgets.base, 'Style','edit', 'Unit','Pixels', 'Position',[5,120,580,20]);
    widgets.SelectType = uibuttongroup(widgets.base, 'Unit','Pixels', 'Position',[5,70,150,40],'Title','Data Type','BackgroundColor',[0.8,0.8,0.8]);
    widgets.SelectVariable(1) = uicontrol(widgets.SelectType, 'Style','radiobutton','String','MATLAB variable', 'Unit','Pixels','Position',[1,5,150,15],'BackgroundColor',[0.8,0.8,0.8]);
    %widgets.SelectVariable(2) = uicontrol(widgets.SelectType, 'Style','radiobutton','String','.mat file', 'Unit','Pixels','Position',[1,5,150,15],'BackgroundColor',[0.8,0.8,0.8]);
    widgets.ListFile = uicontrol(widgets.base, 'Style','popupmenu','String',basevars, 'Unit','Pixels','Position',[170,75,300,20]);
    widgets.OKButton = uicontrol(widgets.base,'Style','pushbutton', 'Unit','Pixels', 'Position',[400,30,100,30], 'String','OK');
    widgets.CancelButton = uicontrol(widgets.base,'Style','pushbutton', 'Unit','Pixels', 'Position',[100,30,100,30], 'String','Cancel');
    widgets.output = uicontrol(widgets.base,'Style','text', 'Unit','Pixels', 'Position',[10,1,580,20], 'String','     ','BackgroundColor',[0.8,0.8,0.8]);
    widgets.Over = uicontrol(widgets.base,'Style','text', 'Unit','Pixels', 'Position',[0,0,1,1], 'String','     ','BackgroundColor',[0.8,0.8,0.8]);
    set(widgets.SelectType,'SelectionChangedFcn',{ @ImportDataCallback, 'SelectDataType', widgets});
    set(widgets.OKButton,'Callback',{ @ImportDataCallback, 'OK', widgets});
    set(widgets.CancelButton,'Callback',{ @ImportDataCallback, 'Cancel', widgets},'UserData',1);
    set(widgets.ListFile,'Callback',{ @ImportDataCallback, 'SelectFile', widgets});
    DataFound = 0;
    % If defaultName set: search for the variable
    selected = 0;
    if nargin > 1
        for cnt = 1:numel(basevars)
            if strcmp(basevars{cnt},defaultName)
                selected = cnt;
                set(widgets.ListFile,'Value',cnt);
                set(widgets.DataName,'String',defaultName);
            end
        end
    end

    while DataFound == 0
    waitfor(widgets.CancelButton,'UserData',0);
    if isgraphics(widgets.output) == 0       % widget is deleted by clicking on cross
        return;
    end
    if get(widgets.output,'UserData') == 0   % Cancel
        delete(widgets.base);
        return;
    else                                     % OK: Data selected
        DataType = get(get(widgets.SelectType,'SelectedObject'),'String');
        switch DataType
            case 'MATLAB variable'
                try
                    data = evalin('base',[get(widgets.DataName,'String'),';']);
                    DataFound = 1;
                catch
                    set(widgets.output,'String',['Variable ',get(widgets.DataName,'String'),' does not exist!!']);
                    DataFound = 0;
                end
        end
        if DataFound == 0
            set(widgets.CancelButton, 'UserData',1);
        end
    end
    end
    delete(widgets.base);
    return
end


function res=ImportDataCallback(src, evt, action, widgets)
    %% Callback Routines for the ImportData widget
    res = 0;
    switch action
        case 'SelectDataType'
            DataType = get(get(widgets.SelectType,'SelectedObject'),'String');
        case 'SelectFile'
            u = get(src,'String');
            name = u{get(src,'Value')};
            set(widgets.DataName,'String',name);
        case 'OK'
            set(widgets.output,'UserData',1);
            set(widgets.CancelButton,'UserData',0);
            %delete(widgets.Over);
        case 'Cancel'
            set(widgets.output,'UserData',0);
            set(widgets.CancelButton,'UserData',0);
            %delete(widgets.Over);
    end
end


