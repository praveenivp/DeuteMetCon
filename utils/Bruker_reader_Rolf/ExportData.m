function newexportpars = ExportData(data, defaultPathName, exportpars, write_immediately, varname)
    % Exports 1D-data of different kind.
    % Currently available data types:
    %   - Export Matlab variable from base workspace
    %   - csv file
    %   - mat file
    %   -Excel file
    
    % Expected data structure
    % either: 
    % (not yet implemented) 1D array: simply exports this array. If x-axis is required, it will be 1, 2, ...
    % (not yet implemented) 2D array: as above. If necessary, it will export all (:,i) separately.
    % data structure array:
    %   data.y: y-values. Can be a 2D-array, then all (:,i) stored separately with the same x values and name
    %   data.x: x-values
    %   data.name
    % can be an array of arbitrary size

    path = [];
    newexportpars = exportpars;
    % Set default path name
    defaultpath = '';
    if nargin > 1 && ischar(defaultPathName) && exist(defaultPathName,'dir') == 7
        defaultpath = defaultPathName;
    end
    if nargin > 2 && isstruct(exportpars) && isfield(exportpars,'ExcelFileName') && ischar(exportpars.ExcelFileName)
        defaultfile.Excel = exportpars.ExcelFileName;
    end
    if nargin > 2 && isstruct(exportpars) && isfield(exportpars,'MatFileName') && ischar(exportpars.MatFileName)
        defaultfile.Mat = exportpars.MatFileName;
    end
    if nargin < 4
        write_immediately = 0;
    end
    if nargin < 5
        write_immediately = 0;
        varname = '';
    end
    
    % Get variables in basic workspace
    if numel(exportpars.ExportDataType) < 3
        exportpars.ExportDataType = 'var';
    end
    if strcmp(exportpars.ExportDataType,'var')
        basevars = evalin('base','who');
    elseif strcmp(exportpars.ExportDataType,'mat')
        b = dir([defaultPathName,filesep, '*.mat']);
        basevars = {b.name};
    elseif strcmp(exportpars.ExportDataType,'csv')
        b = dir([defaultPathName,filesep, '*.csv']);
        basevars = {b.name};
    elseif strcmp(exportpars.ExportDataType,'exc')
        b = dir([defaultPathName,filesep, '*.xls*']);
        basevars = {b.name};
    end
    if numel(basevars) == 0
        basevars = ' ';
    end

    noGUI = 0;
    if write_immediately == 1
        if strcmp(exportpars.ExportDataType,'mat') && numel(defaultfile.Mat) > 0
            noGUI = 1;
        elseif strcmp(exportpars.ExportDataType,'exc') && numel(defaultfile.Excel) == 0
            noGUI = 1;
        end
    end
    
    % Determine whether this is 1D (graph) or 2D (image)
    s = size(data);
    % if the data is to be stored without widget, this is done here
    if noGUI == 1
        if strcmp(exportpars.ExportDataType,'mat')
            v.(varname) = data;
            save([defaultPathName,filesep,defaultfile.Mat],'-struct','v','-append');
        elseif strcmp(exportpars.ExportDataType,'exc')
        end
        return
    end
    
    %% Generate widget
    wsize = [600,200];
    widgets.base = figure('Name','Export Data','NumberTitle','off','DockControls','off',...
                        'Toolbar','none','Menubar','none','Unit','Pixels','Position',[100,100,wsize(1), wsize(2)],'Color',[0.8,0.8,0.8]);
    widgets.DataName = uicontrol(widgets.base, 'Style','edit', 'Unit','Pixels', 'Position',[5,wsize(2)-30,580,20], 'HorizontalAlignment','left', 'String',defaultpath, 'UserData',defaultfile);
    widgets.SelectType = uibuttongroup(widgets.base, 'Unit','Pixels', 'Position',[5,70,150,wsize(2)-110],'Title','Data Type','BackgroundColor',[0.8,0.8,0.8]);
    widgets.SelectVariable(1) = uicontrol(widgets.SelectType, 'Style','radiobutton','String','MATLAB variable', 'Unit','Pixels','Position',[1,wsize(2)-145,150,15],'BackgroundColor',[0.8,0.8,0.8]);
    widgets.SelectVariable(2) = uicontrol(widgets.SelectType, 'Style','radiobutton','String','mat file', 'Unit','Pixels','Position',[1,wsize(2)-160,150,15],'BackgroundColor',[0.8,0.8,0.8]);
    widgets.SelectVariable(3) = uicontrol(widgets.SelectType, 'Style','radiobutton','String','csv file', 'Unit','Pixels','Position',[1,wsize(2)-175,150,15],'BackgroundColor',[0.8,0.8,0.8]);
    widgets.SelectVariable(4) = uicontrol(widgets.SelectType, 'Style','radiobutton','String','Excel file', 'Unit','Pixels','Position',[1,wsize(2)-190,150,15],'BackgroundColor',[0.8,0.8,0.8]);
    widgets.DirName = uicontrol(widgets.base, 'Style','text', 'Unit','Pixels', 'Position',[170,wsize(2)-65,300,20], 'HorizontalAlignment','left', 'String',defaultpath,'BackgroundColor',[0.8,0.8,0.8]);
    widgets.ListFile = uicontrol(widgets.base, 'Style','popupmenu','String',basevars, 'Unit','Pixels','Position',[170,wsize(2)-90,300,20]);
    widgets.OKButton = uicontrol(widgets.base,'Style','pushbutton', 'Unit','Pixels', 'Position',[400,30,100,30], 'String','OK');
    widgets.CancelButton = uicontrol(widgets.base,'Style','pushbutton', 'Unit','Pixels', 'Position',[100,30,100,30], 'String','Cancel');
    widgets.output = uicontrol(widgets.base,'Style','text', 'Unit','Pixels', 'Position',[10,1,580,20], 'String','     ','BackgroundColor',[0.8,0.8,0.8]);
    widgets.Over = uicontrol(widgets.base,'Style','text', 'Unit','Pixels', 'Position',[0,0,1,1], 'String','     ','BackgroundColor',[0.8,0.8,0.8]);
    if strcmp(exportpars.ExportDataType,'var')
        widgets.SelectType.SelectedObject = widgets.SelectVariable(1);
        widgets.DirName.Visible = 'off';
    elseif strcmp(exportpars.ExportDataType,'mat')
        widgets.SelectType.SelectedObject = widgets.SelectVariable(2);
        if numel(exportpars.MatFileName) > 0
            ind = find(strcmp(widgets.ListFile.String, exportpars.MatFileName) == 1);
            if ind > 0
                widgets.ListFile.Value = ind;
            end
            widgets.DataName.String = [defaultpath,filesep,exportpars.MatFileName];
        end
    elseif strcmp(exportpars.ExportDataType,'csv')
        widgets.SelectType.SelectedObject = widgets.SelectVariable(3);
    elseif strcmp(exportpars.ExportDataType,'exc')
        widgets.SelectType.SelectedObject = widgets.SelectVariable(4);
    end
    set(widgets.SelectType,'SelectionChangedFcn',{ @ExportDataCallback, 'SelectDataType', widgets});
    set(widgets.OKButton,'Callback',{ @ExportDataCallback, 'OK', widgets});
    set(widgets.CancelButton,'Callback',{ @ExportDataCallback, 'Cancel', widgets},'UserData',1);
    set(widgets.ListFile,'Callback',{ @ExportDataCallback, 'SelectFile', widgets});
    ExportDone = 0;
    while ExportDone == 0
        waitfor(widgets.CancelButton,'UserData',0);
        if isgraphics(widgets.output) == 0       % widget is deleted by clicking on cross
            return;
        end
        if get(widgets.output,'UserData') == 0   % Cancel
            delete(widgets.base);
            return;
        else                                     % OK: Data selected
            p = strtrim(get(widgets.DirName,'String'));
            if numel(p)>0
                path = p;
            end
            
            DataType = get(get(widgets.SelectType,'SelectedObject'),'String');
            switch DataType
                case 'MATLAB variable'
                    name = strtrim(get(widgets.DataName,'String'));
                    doesexist = strcmp(name,basevars);
                    if max(doesexist)>0                % The variable name exists already in the base workspace
                        b = questdlg('Variable already exists!','Copy variable to workspace','Overwrite','Append','Cancel','Cancel');
                        if strcmp(b, 'Overwrite')
                            assignin('base',name,data);
                            ExportDone = 1;
                        elseif strcmp(b, 'Append')
                            val = evalin('base',[name,';']);
                            if (isstruct(val))     %Variable is a structure -> Check if it has the right members
                                if isfield(val,'x') && isfield(val, 'y') && isfield(val,'name')
                                    val = [val,data];
                                    assignin('base',name,val);
                                    ExportDone = 1;
                                else
                                    errordlg('Wrong structure');
                                    ExportDone = 0;
                                end
                            else                  % Variable is not a structure: must be array with correct dimensionality
                                s = size(val);
                                ds = size(data.y);
                                if s(1) == ds(1) && numel(s)<=2
                                    assignin('base',name,[val,data.y]);                                    
                                else
                                    errordlg('Wrong array dimensions');
                                    ExportDone = 0;
                                end
                            end
                        elseif strcmp(b, 'Cancel')
                            ExportDone = 0;
                        end
                    else                               % The variable name does not yet exist                       
                        assignin('base',name,data);
                        ExportDone = 1;
                    end
                case 'mat file'
                    name = strtrim(get(widgets.DataName,'String'));
                    if exist(name,'dir') == 7        % The file name already exists and is a folder
                        errordlg([name, ' is a directory!']);
                        ExportDone = 0;
                    elseif exist(name, 'file') == 2  % The file name already exists and is a file                        
                        b = questdlg('File exists!','Store as mat file','Append','Overwrite','Cancel','Append');
                        if strcmp(b, 'Cancel')
                            ExportDone = 0;
                        else
                            if exist('varname','var') == 1
                                newname{1} = varname;
                            else
                                newname = inputdlg('Enter variable name:','Save in mat file');
                                if numel(newname) == 0
                                    ExportDone = 0;
                                    return
                                end
                            end
                            v.(newname{1}) = data;
                            if strcmp(b, 'Overwrite')
                                save(name,'-struct','v');
                                ExportDone = 1;
                            end
                            if strcmp(b, 'Append')
                                cont = whos('-file',name);
                                filevarnames = {cont(:).name};
                                while max(strcmp(filevarnames,newname)) > 0
                                    newname = inputdlg('Variable already exists, enter new name:','Save in mat file');
                                end
                                if numel(newname) == 0
                                    ExportDone = 0;
                                    return
                                else
                                    v.(newname{1}) = data;
                                    if exist(name,'file') == 2
                                        save(name,'-struct','v','-append');
                                    else
                                        save(name,'-struct','v');
                                    end
                                    ExportDone = 1;
                                end
                            end
                            clipboard('copy',newname{1});
                        end
                    else
                        if exist('varname','var') == 1
                            newname{1} = varname;
                        else
                            newname = inputdlg('Enter variable name:','Save in mat file');
                            if numel(newname) == 0
                                return;
                            end
                        end
                        v.(newname{1}) = data;
                        save(name,'-struct','v');
                        ExportDone = 1;     
                        clipboard('copy',newname{1});
                    end
                    if ExportDone == 1
                        newexportpars.ExportDataType = 'mat';
                        [p,f,e] = fileparts(name);
                        newexportpars.MatFileName = [f,e];
                        newexportpars.ExportPath = p;
                    end
                case 'csv file'
                    name = strtrim(get(widgets.DataName,'String'));
                    if numel(name)<=4 || strcmp(name(end-3:end),'.csv')==0
                        name = [name,'.csv'];
                    end
                    if exist(name,'dir') == 7        % The file name already exists and is a folder
                        errordlg([name, ' is a directory!']);
                        ExportDone = 0;
                        return;
                    elseif exist(name, 'file') == 2  % The file name already exists and is a file
                        b = questdlg('File exists!','Store as csv file','Overwrite','Cancel','Cancel');
                        if strcmp(b, 'Cancel')
                            ExportDone = 0;
                        end
                    end
                    [file, errormsg] = fopen(name, 'w');
                    if file == -1
                        errordlg(['Cannot open file: ', errormsg]);
                        return
                    end
                    datsize = numel(data);
                    for cnt = 1:datsize
                        siz(cnt) = numel(data(cnt).y);
                        if numel(data(cnt).name) > 0
                            fprintf(file,['x;',data(cnt).name]);
                            if cnt < datsize
                                fprintf(file,';');
                            end
                        else
                            if isreal(data(cnt).y) 
                                fprintf(file,'x; data%d',cnt);
                            else
                                fprintf(file,'x; real(data%d);imag(data%d)',cnt,cnt);
                            end
                            if cnt < datsize
                                fprintf(file,';');
                            end
                        end
                    end
                    fprintf(file,'\n');
                    for cnt = 1:max(siz)
                        for cnt2 = 1:datsize
                            if cnt <= siz(cnt2)
                                if isreal(data(cnt2).y) 
                                    fprintf(file,'%f; %f',data(cnt2).x(cnt),data(cnt2).y(cnt));
                                else
                                    fprintf(file,'%f; %f; %f',data(cnt2).x(cnt),real(data(cnt2).y(cnt)),imag(data(cnt2).y(cnt)));
                                end
                            else
                                fprintf(file,';');
                            end
                            if cnt2 < datsize
                                fprintf(file,';');
                            end

                        end
                        fprintf(file,'\n');
                    end
                    fclose(file);
                    ExportDone = 1;
                case 'Excel file'
                    name = strtrim(get(widgets.DataName,'String'));
                    datsize = numel(data);
                    occupiedcols = 0;
                    sheet = 'data';
                    if contains(name, '.xls')==0
                        name = [name,'.xlsx'];
                    end
                    if exist(name,'dir') == 7        % The file name already exists and is a folder
                        errordlg([name, ' is a directory!']);
                        ExportDone = 0;
                        return;
                    elseif exist(name, 'file') == 2  % The file name already exists and is a file
                        b = questdlg('File exists!','Store as csv file','Overwrite','Append','Cancel','Cancel');
                        if strcmp(b, 'Cancel')
                            ExportDone = 0;
                            return
                        elseif strcmp(b, 'Append')
                            [status,sheets,xlFormat] = xlsfinfo(name);
                            if numel(status)==0
                                errordlg('Not a valid Excel file!');
                                return
                            end
                            sheet = sheets{1};
                            d = xlsread(name, sheet, 'a1:zz5');
                            occupiedcols = size(d);
                            occupiedcols = occupiedcols(2);
                            %col = char([floor(num/27),mod(num-1,26)+1]+64);
                        end
                    end
                    for cnt = 1:datsize
                        siz(cnt) = numel(data(cnt).y);
                    end

                    for cnt = 1:datsize
                        fieldnames{2*cnt-1} = 'x';
                        fieldnames{2*cnt} = data(cnt).name;
                    end
                    startcol = occupiedcols + 1;
                    sc = char([floor(startcol/27),mod(startcol-1,26)+1]+64);
                    if startcol <=26
                        sc = sc(2);
                    end
                    endcol = occupiedcols + datsize*2;
                    ec = char([floor(endcol/27),mod(endcol-1,26)+1]+64);
                    if endcol <=26
                        ec = ec(2);
                    end
                    [succ, errormsg] = xlswrite(name,{['imported by rpwin on ',date]}, sheet, [sc,'1:',sc,'1']);
                    if succ == 0
                        errormsg(['Could not write file: ',errormsg.message]);
                    end
                    [succ, errormsg] = xlswrite(name,fieldnames, sheet,[sc,'2:',ec,'2']);
                    datmat = zeros(max(siz),datsize*2);
                    for cnt = 1:datsize
                        datmat(1:numel(data(cnt).x),cnt*2-1) = data(cnt).x;
                        datmat(1:numel(data(cnt).y),cnt*2) = data(cnt).y;
                    end
                    [succ, errormsg] = xlswrite(name,datmat, sheet,[sc,'3:',ec,num2str(max(siz)+2)]);
                    ExportDone = 1;

            end
            if ExportDone == 0
                set(widgets.CancelButton, 'UserData',1);
            end
        end
    end
    delete(widgets.base);
    return
end


function res=ExportDataCallback(src, evt, action, widgets)
    %% Callback Routines for the ImportData widget
    res = 0;
    switch action
        case 'SelectDataType'
            DataType = get(get(widgets.SelectType,'SelectedObject'),'String');
            switch DataType
                case 'MATLAB variable'
                    basevars = evalin('base','who');
                    set(widgets.ListFile, 'String', basevars);
                    set(widgets.DirName, 'String',' ');
                    widgets.DirName.Visible = 'off';
                case 'mat file'
                    widgets.DirName.Visible = 'on';
                    newdir = uigetdir;
                    def = get(widgets.DataName,'UserData');
                    if isfield(def,'Mat') 
                        filename = def.Mat;
                    else
                        filename = '';
                    end
                    set(widgets.DataName,'String',[newdir, filesep,filename]);
                    listing = dir([newdir,filesep, '*.mat']);
                    isdir = [listing.isdir];
                    listing = listing(~isdir);
                    set(widgets.DirName, 'String',newdir);
                    if numel(listing) == 0
                        set(widgets.ListFile,'String',' ');
                        return
                    end
                    set(widgets.ListFile,'String',{listing.name});
                case 'csv file'
                    widgets.DirName.Visible = 'on';
                    newdir = uigetdir;
                    set(widgets.DataName,'String',[newdir, filesep]);
                    listing = dir([newdir,filesep, '*.csv']);
                    isdir = [listing.isdir];
                    listing = listing(~isdir);
                    set(widgets.DirName, 'String',newdir);
                    if numel(listing) == 0
                        set(widgets.ListFile,'String',' ');
                        return
                    end
                    set(widgets.ListFile,'String',{listing.name});
                case 'Excel file'
                    widgets.DirName.Visible = 'on';
                    def = get(widgets.DataName,'UserData');
                    if isfield(def,'Excel') 
                        filename = def.Excel;
                    else
                        filename = '';
                    end
                    direc = strtrim(get(widgets.DataName,'String'));
                    if numel(direc)>0
                        if (exist(direc,'dir')==7)
                            newdir = uigetdir(direc);
                        else
                            [direc,name,ext] = fileparts(direc);
                            if (exist(direc,'dir')==7)
                                newdir = uigetdir(direc);
                            end
                        end
                    else
                        newdir = uigetdir;
                    end
                    if ischar(newdir)== 0 || exist(newdir,'dir')==0
                        return
                    end
                    set(widgets.DataName,'String',[newdir, filesep,filename]);
                    listing = dir([newdir,filesep, '*.xls*']);
                    isdir = [listing.isdir];
                    listing = listing(~isdir);
                    set(widgets.DirName, 'String',newdir);
                    if numel(listing) == 0
                        set(widgets.ListFile,'String',' ');
                        return
                    end
                    set(widgets.ListFile,'String',{listing.name});
                    ind = find(strcmp({listing.name},filename)==1);
                    if numel(ind) == 0
                        ind(1) = 1;
                    end
                    set(widgets.ListFile,'Value',ind(1));
            end
        case 'SelectFile'
            u = get(src,'String');
            name = u{get(src,'Value')};
            DataType = get(get(widgets.SelectType,'SelectedObject'),'String');
            def = get(widgets.DataName,'UserData');
            switch DataType
                case 'MATLAB variable'
                case 'mat file'
                    DirName = get(widgets.DirName,'String');
                    name = [DirName,filesep,name];
                case 'csv file'
                    DirName = get(widgets.DirName,'String');
                    name = [DirName,filesep,name];
                case 'Excel file'
                    DirName = get(widgets.DirName,'String');
                    def.Excel = name;
                    name = [DirName,filesep,name];
            end
            set(widgets.DataName,'String',name);
            set(widgets.DataName,'UserData',def);
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


