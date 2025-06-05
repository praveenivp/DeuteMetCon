function images = SiReadProData(path, startpath)
if nargin < 2
    startpath = '.';
end
if nargin < 1
    path = '.';
end

files = SiSelectData(path,1);
if numel(files)==0
    images = [];
    return
end
images = struct('Info',{},'Image',{});

%if max([files.NSlices]) > 1 || max([files.NContrasts]) > 1
    params.complex = 0;
    params.phase = 0;
    params.slice = 1;
    params.contrast = 1;
    params.phase = 0;
    params.set = 0;
    wsize = [300,250];
    w(1) = figure('Name','Set Parameters','NumberTitle','off','DockControls','off',...
              'Toolbar','none','Menubar','none','Unit','Pixels','Position',[50,50,wsize(1),wsize(2)],'Color',[0.8,0.8,0.8],'UserData',params);
    w(2) = uicontrol('Parent',w(1),'Style','pushbutton','Units','pixels','Position',[10,10,100,20],'String',' OK ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);
    w(3) = uicontrol('Parent',w(1),'Style','radiobutton','Units','pixels','Position',[10,wsize(2)-30,100,20],'String',' Abs ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);
    w(4) = uicontrol('Parent',w(1),'Style','radiobutton','Units','pixels','Position',[10,wsize(2)-50,100,20],'String',' Complex ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8],'Value',params.complex);
    w(5) = uicontrol('Parent',w(1),'Style','radiobutton','Units','pixels','Position',[10,wsize(2)-70,100,20],'String',' Phase ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8],'Value',params.phase);
    w(6) = uicontrol('Parent',w(1),'Style','radiobutton','Units','pixels','Position',[10,wsize(2)-150,120,20],'String',' All Slices ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);
    junk = uicontrol('Parent',w(1),'Style','text','Units','pixels','Position',[10,wsize(2)-175,110,20],'String',' Slices: ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);
    w(7) = uicontrol('Parent',w(1),'Style','popupmenu','Units','pixels','Position',[10,wsize(2)-190,100,20],'String',num2str(linspace(1,files(1).NSlices,files(1).NSlices)'),'FontName','FixedWidth','BackgroundColor',[1,1,1],'Value',params.slice);
    w(8) = uicontrol('Parent',w(1),'Style','radiobutton','Units','pixels','Position',[150,wsize(2)-150,120,20],'String',' All Contrasts ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);
    junk = uicontrol('Parent',w(1),'Style','text','Units','pixels','Position',[150,wsize(2)-175,120,20],'String',' Contrasts: ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);
    w(9) = uicontrol('Parent',w(1),'Style','popupmenu','Units','pixels','Position',[150,wsize(2)-190,100,20],'String',num2str(linspace(1,files(1).NContrasts,files(1).NContrasts)'),'FontName','FixedWidth','BackgroundColor',[1,1,1],'Value',params.slice);
    w(10) = uicontrol('Parent',w(1),'Style','text','Units','pixels','Position',[150,wsize(2)-30,120,20],'String','Phase:','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);
    w(11) = uicontrol('Parent',w(1),'Style','radiobutton','Units','pixels','Position',[150,wsize(2)-50,120,20],'String',' None ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);
    w(12) = uicontrol('Parent',w(1),'Style','radiobutton','Units','pixels','Position',[150,wsize(2)-70,120,20],'String',' Each ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);
    w(13) = uicontrol('Parent',w(1),'Style','radiobutton','Units','pixels','Position',[150,wsize(2)-90,120,20],'String',' First ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);
    w(14) = uicontrol('Parent',w(1),'Style','radiobutton','Units','pixels','Position',[150,wsize(2)-110,120,20],'String',' Last ','FontName','FixedWidth','BackgroundColor',[0.8,0.8,0.8]);

    if numel(files) == max([files.NSlices]);
        set(w(6),'Value',1);
        params.slice = -1;
        set(w(1),'UserData',params);
    else
        set(w(6),'Value',0);
        params.slice = 1;
    end
    if max([files.NSlices]) == 1
        set(w(6),'Visible','off');
        set(w(7),'Enable','off');
    end
        
    if params.complex == 0 & params.phase == 0
        set(w(3),'Value',1);
    end
    set(w(2), 'Callback',{@PCallback, 'OK',w});
    set(w(3), 'Callback',{@PCallback, 'ABS',w});
    set(w(4), 'Callback',{@PCallback, 'COMPL',w});
    set(w(5), 'Callback',{@PCallback, 'PHASE',w});
    set(w(6), 'Callback',{@PCallback, 'ALLSLICE',w});
    set(w(7), 'Callback',{@PCallback, 'SLICE',w});
    set(w(11), 'Callback',{@PCallback, 'PH_NONE',w});
    set(w(12), 'Callback',{@PCallback, 'PH_EACH',w});
    set(w(13), 'Callback',{@PCallback, 'PH_FIRST',w});
    set(w(14), 'Callback',{@PCallback, 'PH_LAST',w});
    waitfor(w(2));
    if ishandle(w(1))
        params = get(w(1),'UserData');
        params.set = 1;
        delete(gcf);
        if params.slice > -1
            slpos = files(params.slice).Position;
            ind = find([files.Position]>slpos-0.01 & [files.Position]<slpos+0.01);
            files = files(ind);
        end
    end



%end
if params.slice == -1
    NSlices = max([files.NSlices]);
else
    NSlices = numel(params.slice);
end
NContrasts = numel(files)/NSlices;%max([files.Contrast]);
[junk,order] = sort([files.Position]);
files = files(order);
[junk, order] = sort([files.Contrast]);
files = files(order);
files = reshape(files, [NContrasts, NSlices]);
m=0;
for cnt=1:NSlices
    for cnt1 = 1:NContrasts
    %im = dicomread(files(cnt).Name);
    display(['Reading ',files(cnt1,cnt).Name]);
    im = cast(dicomread(files(cnt1,cnt).Name),'double');
    images(cnt1,cnt) = struct('Info',files(cnt1,cnt),'Image',im);
    ma = max(max(max(im)));
    if ma > m
        m = ma;
    end
    end
end
%%Scale images
m = cast(m,'double');
scale = 256./m;
for cnt=1:NSlices
    for cnt1 = 1:NContrasts
    images(cnt1,cnt).Image = images(cnt1,cnt).Image*scale;
    end
end
allinfo = [images.Info];
return
end


function res = PCallback(src, evt, action, w)
switch action
    case 'OK'
        delete(src);
    case 'ABS'
        if get(src, 'Value') == 1
            set(w(4),'Value',0);
            set(w(5),'Value',0);
            p=get(w(1),'UserData');
            p.complex = 0;
            p.phase = 0;
            set(w(1),'UserData',p);
        else
            set(src, 'Value',1);
        end
    case 'COMPL'
        if get(src, 'Value') == 1
            set(w(3),'Value',0);
            set(w(5),'Value',0);
            p=get(w(1),'UserData');
            p.complex = 1;
            p.phase = 0;
            set(w(1),'UserData',p);
        else
            set(src, 'Value',1);
        end
    case 'PHASE'
        if get(src, 'Value') == 1
            set(w(3),'Value',0);
            set(w(4),'Value',0);
            p=get(w(1),'UserData');
            p.complex = 0;
            p.phase = 1;
            set(w(1),'UserData',p);
        else
            set(src, 'Value',1);
        end
    case 'ALLSLICE'
            if get(src, 'Value')==1
                p=get(w(1),'UserData');
                p.slice = -1;
                set(w(1),'UserData',p);  
                set(w(7),'Enable','off');
            else
                p=get(w(1),'UserData');
                p.slice = get(w(7),'Value');
                set(w(1),'UserData',p);        
                set(w(7),'Enable','on');
            end
    case 'SLICE'
            p=get(w(1),'UserData');
            get(src,'Value');
            p.slice = get(src,'Value');
            set(w(1),'UserData',p);
end
return
end