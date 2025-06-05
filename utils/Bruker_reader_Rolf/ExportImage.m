function out = ExportImage(imaxes)
    % Exports images
    % Currently available data types:
    %   - fig-file
    %   - graphics files in different formats

    


    out = 0;
    % Get variables in basic workspace
% F = getframe(handles.axes1);
% Image = frame2im(F);
% imwrite(Image, 'Image.jpg')
% oder
% f_new = figure;
% ax_new = copyobj(ax_old,f_new)
% set(ax_new,'Position','default')
    newdir = cd;
    extensions = {'*.jpg','*.png','*.gif','*.tif','*.pdf','*.eps','*.svg'};
    listing = [];
    for cnt = 1:numel(extensions)
        list = dir([newdir, filesep, extensions{cnt}]);
        if numel(list)>0
            listing = [listing; list];
        end
    end
    if numel(listing) == 0
        listing.name = '';
    end

    %% Generate widget
    wsize = [600,200];
    widgets.base = figure('Name','Export Data','NumberTitle','off','DockControls','off',...
                        'Toolbar','none','Menubar','none','Unit','Pixels','Position',[100,100,wsize(1), wsize(2)],'Color',[0.8,0.8,0.8]);
    widgets.DataName = uicontrol(widgets.base, 'Style','edit', 'Unit','Pixels', 'Position',[5,wsize(2)-30,580,20], 'HorizontalAlignment','left','String',[newdir, filesep]);
    widgets.SelectType = uibuttongroup(widgets.base, 'Unit','Pixels', 'Position',[5,70,150,wsize(2)-110],'Title','Data Type','BackgroundColor',[0.8,0.8,0.8]);
    widgets.SelectVariable(1) = uicontrol(widgets.SelectType, 'Style','radiobutton','String','graphics file', 'Unit','Pixels','Position',[1,wsize(2)-145,150,15],'BackgroundColor',[0.8,0.8,0.8]);
    widgets.SelectVariable(2) = uicontrol(widgets.SelectType, 'Style','radiobutton','String','fig file', 'Unit','Pixels','Position',[1,wsize(2)-160,150,15],'BackgroundColor',[0.8,0.8,0.8]);
    widgets.DirName = uicontrol(widgets.base, 'Style','text', 'Unit','Pixels', 'Position',[170,wsize(2)-65,300,20], 'HorizontalAlignment','left', 'String','','BackgroundColor',[0.8,0.8,0.8], 'String',newdir);
    widgets.ListFile = uicontrol(widgets.base, 'Style','popupmenu', 'Unit','Pixels','Position',[170,wsize(2)-90,300,20],'String',{listing.name});
    widgets.OKButton = uicontrol(widgets.base,'Style','pushbutton', 'Unit','Pixels', 'Position',[400,30,100,30], 'String','OK');
    widgets.CancelButton = uicontrol(widgets.base,'Style','pushbutton', 'Unit','Pixels', 'Position',[100,30,100,30], 'String','Cancel');
    widgets.output = uicontrol(widgets.base,'Style','text', 'Unit','Pixels', 'Position',[10,1,580,20], 'String','     ','BackgroundColor',[0.8,0.8,0.8]);
    widgets.Over = uicontrol(widgets.base,'Style','text', 'Unit','Pixels', 'Position',[0,0,1,1], 'String','     ','BackgroundColor',[0.8,0.8,0.8]);
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
            DataType = get(get(widgets.SelectType,'SelectedObject'),'String');
            switch DataType
                case 'graphics file'
                    name = strtrim(get(widgets.DataName,'String'));
                    if exist(name,'dir') == 7        % The file name already exists and is a folder
                        errordlg([name, ' is a directory!']);
                        ExportDone = 0;
                        return;
                    elseif exist(name, 'file') == 2  % The file name already exists and is a file
                        b = questdlg('File exists!','Store as graphics file','Overwrite','Cancel','Cancel');
                        if strcmp(b, 'Cancel')
                            ExportDone = 0;
                        end
                    end
                    exportgraphics(imaxes,name);
%                     F = getframe(imaxes);  % This solution doesn't copy the axes labels
%                     Image = frame2im(F);
%                     [~,~,filetype] = fileparts(name);
%                     if strcmp(filetype,'.gif')
%                         [a, map] = rgb2ind(Image,256);
%                         imwrite(a,map,name);
%                     else
%                         imwrite(Image, name);
%                     end
                    ExportDone = 1;
                case 'fig file'
                    f_new = figure('Position',imaxes.Parent.Position);
                    if numel(imaxes.YAxis) == 1
                        ax_new = copyobj(imaxes,f_new)
                        set(ax_new,'Position',imaxes.Position)
                    else
                        yyaxis(imaxes, 'left');
                        ax = gca;
                        set(ax,'XLim',imaxes.XLim,'YLim',imaxes.YAxis(1).Limits,'YScale',imaxes.YAxis(1).Scale);
                        copyobj(imaxes.Children,ax);
                        yyaxis(imaxes, 'right');
                        yyaxis(ax, 'right');
                        set(ax,'XLim',imaxes.XLim,'YLim',imaxes.YAxis(2).Limits,'YScale',imaxes.YAxis(1).Scale);
                        copyobj(imaxes.Children,ax);

                    end
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
                case 'graphics file'
                    newdir = uigetdir;
                    extensions = {'*.jpg','*.png','*.gif','*.tif','*.pdf','*.eps','*.svg'};
                    listing = [];
                    for cnt = 1:numel(extensions)
                        list = dir([newdir, filesep, extensions{cnt}]);
                        if numel(list)>0
                            listing = [listing; list];
                        end
                    end
                    if numel(listing) == 0
                        listing.name = '';
                    end
                    set(widgets.DataName,'String',[newdir, filesep]);
                    set(widgets.DirName, 'String',newdir);
                    set(widgets.ListFile,'String',{listing.name});
                case 'fig file'
%                     newdir = uigetdir;
%                     set(widgets.DataName,'String',[newdir, filesep]);
%                     listing = dir(newdir); %,filesep, {'*.jpg','*.png','*.gif','*.tif','*.pdf','*.eps','*.svg'});
%                     isdir = [listing.isdir];
%                     listing = listing(~isdir);
%                     set(widgets.DirName, 'String',newdir);
%                     if numel(listing) == 0
%                         set(widgets.ListFile,'String',' ');
%                         return
%                     end
%                     set(widgets.ListFile,'String',{listing.name});
            end
        case 'SelectFile'
            u = get(src,'String');
            name = u{get(src,'Value')};
            DataType = get(get(widgets.SelectType,'SelectedObject'),'String');
            switch DataType
                case 'graphics file'
                    DirName = get(widgets.DirName,'String');
                    name = [DirName,filesep,name];
                case 'fig file'
                    DirName = get(widgets.DirName,'String');
                    name = [DirName,filesep,name];
            end
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


