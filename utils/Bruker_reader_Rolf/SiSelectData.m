function files = SiSelectData(path, return_fileinfo, read_immediately, filter)

%% Check whether a path was passed and whether it's a directory name
correctpath = 1;
correctfile = 0;

if nargin < 4
    filter = 'IMA';
end
extension = ['\*.',filter];

if nargin<3
    read_immediately = 0;
end
if nargin<2
    return_fileinfo=0;
end
if nargin<1
    correctpath=0;
else
    if isfile(path)          % The path is directly an existing file. It is assumed that the fie type is alright.
        if return_fileinfo==0
            files = path;
            return;
        else
            correctpath=1;
            correctfile = 1;
            f = dir(path);
        end 
        files=0;
    elseif ~isfolder(path)      % Path is a directory that doesn't exist - could be a path with wildcard
        if endsWith(path,'.ima','IgnoreCase',true)
            d = dir(path);
        else
            d = dir([path,'.IMA']);
        end
        if numel(d) > 0    % yes, it is a wildcard
            correctpath = 1;
            read_immediately = 1;
            if ~endsWith(path,'.ima','IgnoreCase',true)
                 path = [path,'.IMA'];
            end
            [p,f,ext] = fileparts(path);
            path = p;
            extension = [filesep,f,ext];
        else
            correctpath = 0;   % No file found - ask for new path
        end
    else                                            % Path is an existing directory
        if numel(dir([path,filesep,'*.',filter]))>0 || numel(dir([path,filesep,'*.',lower(filter)]))>0 || numel(dir([path,filesep,'*.',upper(filter)]))>0         % The directoy contains DICOM files
            correctpath = 1;
            read_immediately = 1;
        else
            d = dir(path);
            if isnumeric(str2double({d.name})) == 1
                correctpath = 1;
                read_immediately = 1;
            end
        end
    end
end
% If no valid directory was passed, ask for it
if correctpath == 0 | read_immediately == 0
    if correctpath == 0
        path = '.';
    end
    path = uigetdir(path, 'Select data directory');
    if path == 0
        files = 0;
        return 
    end
end

%% Find all the data files. Searches current directory and one level up.
if correctfile == 0                          % Only, if the passed path isn't already the full file name
    disp(['Reading directory ',path,'.']);
    f = dir([path,extension]);
    if numel(f) == 0
        filter = lower(filter);
        extension = ['\*.',filter];
        f = dir([path,extension]);
        if numel(f) == 0
            filter = upper(filter);
            extension = ['\*.',filter];
            disp(['Reading directory ',path,'.']);
            f = dir([path,extension]);
        end
    end
end

counter = 1;
if numel(f)==0
    dirs = dir([path, '\*']);
    for cnt = 1:numel(dirs)
        
        if dirs(cnt).isdir == 1 && strcmp(dirs(cnt).name,'.')==0 && strcmp(dirs(cnt).name,'..')==0
            f = dir([path,'\',dirs(cnt).name,extension]);
            for cnt2 = 1:numel(f)
                fnames{counter} = [path,'\',dirs(cnt).name,'\',f(cnt2).name];
                counter = counter+1;
            end
        end
    end
else
    for cnt = 1:numel(f)
        fnames{counter} = [path,'\',f(cnt).name];
        counter = counter+1;        
    end
end
display(['Number of datafiles found: ', int2str(numel(fnames))]); 
if numel(fnames)==0
    return
end

%% Read the file infos
listwidth = zeros(8,1);

fileinfo = struct('filename',fnames,'path',path,'Number',0,'Date',0,'Dim',0,'Protocol','', ...
               'Sequence','','Slice',0.0,'TR',0.0,'TE',0.0,'TI',0.0,'Averages',0, 'MTOffset',0, ...
               'Nucleus','','ImageSize',0,'ACQSize',0,'FOV',0,'Position',0,'Coord',[],'Orientation',[],'PatientPos','','NSlices',0,'NContrasts',0,'Contrast',0);

ind = 1;
for cnt = 1: numel(fnames)
try
    fileinfo(ind).filename;
    info = dicominfo(fnames{cnt});
    fileinfo(ind).Name = fnames{cnt};
    fileinfo(ind).Number = info.SeriesNumber;
    fileinfo(ind).Date = info.AcquisitionDate;
    fileinfo(ind).Protocol = info.SeriesDescription;
    fileinfo(ind).Sequence = info.SequenceName;
    fileinfo(ind).Slice = info.SliceThickness;
    fileinfo(ind).TR = info.RepetitionTime;
    fileinfo(ind).TE = info.EchoTime;
    if isfield(info, 'InversionTime')==1
        fileinfo(ind).TI = info.InversionTime;
    end
    fileinfo(ind).Averages = info.NumberOfAverages;
    fileinfo(ind).Nucleus = info.ImagedNucleus;
    fileinfo(ind).ImageSize = [info.Rows, info.Columns];
    fileinfo(ind).Position = info.SliceLocation;
    fileinfo(ind).Coord = info.ImagePositionPatient;
    fileinfo(ind).Orientation = info.ImageOrientationPatient;
    fileinfo(ind).PatientPos = info.PatientPosition;
    fileinfo(ind).FOV = [info.PixelSpacing(1)*double(info.Rows), info.PixelSpacing(2)*double(info.Columns)];
    if isfield(info, 'EchoNumber')==1
        fileinfo(ind).Contrast = info.EchoNumber;
    elseif isfield(info, 'EchoNumbers')==1
        fileinfo(ind).Contrast = info.EchoNumbers;
    end
    ind = ind+1;
catch
    disp('No DICOM file!');
    fileinfo = [fileinfo(1:ind-1),fileinfo(ind+1:end)];
end
end
[a,s] = sort([fileinfo.Number]);
fileinfo = fileinfo(s);
num = 0;
val = 0;
indices = [];
slpos=0;
for cnt=1:ind-1
    if a(cnt) ~= val
        val = a(cnt);
        num = num+1;
        indices=[indices,cnt];
        if cnt>1
            NSlices(num-1) = numel(unique([fileinfo(indices(num-1):indices(num)-1).Position]));
            NContrasts(num-1) = numel(unique([fileinfo(indices(num-1):indices(num)-1).Contrast]));
            for cnt2 = indices(num-1):indices(num)-1
                fileinfo(cnt2).NSlices =NSlices(num-1);
                fileinfo(cnt2).NContrasts =NContrasts(num-1);
            end
        end
        NSlices(num) = 1;
        NContrasts(num) = 1;
        slpos = fileinfo(cnt).Position;
    else
        if fileinfo(cnt).Position~=slpos
            NSlices(num) = NSlices(num)+1;
            slpos = fileinfo(cnt).Position;
        else
            NContrasts(num) = NContrasts(num)+1;
        end
    end
end
sp = [fileinfo(indices(num):end).Coord];
for cnt2 = indices(num):numel(fileinfo)
    fileinfo(cnt2).NSlices =NSlices(num);
    fileinfo(cnt2).NContrasts =NContrasts(num);
end
display([num2str(num), ' datasets']);
%%Create the widget
labelstrings = {' ','Protocol','TR','TE','TI','Ave','Slices', 'Contrasts'};
%num = numel(fileinfo);
maxlistwidth = 200;
items = num2str([fileinfo(indices).Number]','%3d');
s = size(items);
listwidth(1) = s(2)+2;
for cnt = 1:num
    for cnt2 = 1:numel(listwidth)
        switch cnt2
            case 1
                var = num2str(fileinfo(indices(cnt)).Number,3);
            case 2
                var = fileinfo(indices(cnt)).Protocol;
            case 3
                var = num2str(fileinfo(indices(cnt)).TR,5);
            case 4
                var = num2str(fileinfo(indices(cnt)).TE,5);
            case 5
                var = num2str(fileinfo(indices(cnt)).TI,5);
            case 6
                var = num2str(fileinfo(indices(cnt)).Averages,3);
            case 7
                var = num2str(NSlices(cnt),5);
            case 8
                var = num2str(NContrasts(cnt),5);
        end
        if numel(var)+2>listwidth(cnt2)
            listwidth(cnt2) = numel(var)+2;
        end
    end
end
for cnt = 1:numel(listwidth)
    if listwidth(cnt)>0
        if listwidth(cnt)<numel(labelstrings{cnt})+2
            listwidth(cnt)=numel(labelstrings{cnt})+2;
        end
    end
end

listheight = 17;
labelheight = 17;
buttonheight = 25;
%listwidth(3) = 7;
%listwidth(4) = 7;
%listwidth(5) = 7;
%listwidth(6) = 7;
%listwidth(7) = 7;
totallistwidth = sum(listwidth);
listwidthpix = 8*listwidth;
totallistwidthpix = sum(listwidthpix)+16;
titles = blanks(totallistwidth);
entries = repmat(' ',[num,totallistwidth]);
currpos = 0;
for cnt = 1:numel(listwidth)
    if listwidth(cnt)>0
        titles(currpos + 1:currpos + numel(labelstrings{cnt}))=labelstrings{cnt};
        if cnt < numel(listwidth)
            titles(currpos + listwidth(cnt)) = '|';
        end
        
        for cnt2 = 1:num
            switch cnt
                case 1
                    var = num2str(fileinfo(indices(cnt2)).Number,3);
                    typ = 'num';
                case 2
                    var = fileinfo(indices(cnt2)).Protocol;
                    typ = 'str';
                case 3
                    var = num2str(fileinfo(indices(cnt2)).TR,5);
                    typ = 'num';
                case 4
                    var = num2str(fileinfo(indices(cnt2)).TE,5);
                    typ = 'num';
                case 5
                    var = num2str(fileinfo(indices(cnt2)).TI,5);
                    typ = 'num';
                case 6
                    var = num2str(fileinfo(indices(cnt2)).Averages,3);
                    typ = 'num';
                case 7
                    var = num2str(NSlices(cnt2),5);
                    typ = 'num';
                case 8
                    var = num2str(NContrasts(cnt2),5);
                    typ = 'num';
            end
            if typ == 'str'
                entries(cnt2,currpos + 1:currpos + numel(var)) = var;
            else
                entries(cnt2,currpos + listwidth(cnt)- numel(var)-1:currpos+ listwidth(cnt) - 2) = var;
            end
            if cnt < numel(listwidth)
                entries(cnt2,currpos + listwidth(cnt)) = '|';
            end
        end
        
        
        currpos = currpos + listwidth(cnt);
    end
end

s = size(entries);
if (s(1) == 1)
    entries(2,1) = ' ';
end

lpos = 0;
base = figure('Name','Select Dataset','NumberTitle','off','DockControls','off',...
              'Toolbar','none','Menubar','none','Unit','Pixels','Position',[50,50,totallistwidthpix,min(num*listheight+buttonheight + labelheight,30*listheight+buttonheight+labelheight)],...
              'CloseRequestFcn',{@buttonCallback,'Figure'});

s = get(base, 'Position');
labelbase = uipanel('Parent',base, 'BorderType','line','BorderWidth',1,'Unit','Pixels','Position',[0,s(4)-labelheight,s(3),labelheight]);
label(2) = uicontrol('Parent',labelbase,'Style','text','Unit','Pixels','Position',[0,0,totallistwidthpix,labelheight], ...
                   'String',titles,'FontSize',10,'FontName','FixedWidth', ...
                   'HorizontalAlignment','Left', 'Callback',{@buttonCallback,'Label_Protocol',0});               

               
list = uicontrol('Parent',base,'Style','listbox','Unit','Pixels','Position',[lpos,buttonheight,totallistwidthpix,s(4)-buttonheight-labelheight], ...
                   'String',entries,'FontSize',10,'FontName','FixedWidth', ...
                   'HorizontalAlignment','Right', 'Callback',{@buttonCallback,'List',0}, 'Max',10);
buttons = uipanel('Parent',base, 'BorderType','line','BorderWidth',2,'Unit','Pixels','Position',[0,0,s(3),buttonheight]);
cancelbutton = uicontrol('Parent',buttons,'Style','pushbutton','String','Cancel','Unit','Pixels',...
                         'Position',[s(3)-70,1,50,20], 'Callback',{@buttonCallback,'Cancel', 0});
okbutton = uicontrol('Parent',buttons,'Style','pushbutton','String','OK','Unit','Pixels','Position',[20,1,50,20], ...
                     'Callback',{@buttonCallback,'OK',cancelbutton});
handles.list= list;
handles.cancelbutton = cancelbutton;
%w = uipanel('Parent',base,'Unit','Pixels','Position',[200,40,300,100]);
%a = axes('Parent',w,'Unit','Pixels','Position',[0,0,300,100],'Layer','top','DrawMode','fast');  
%q = line('XData',[0,1],'YData',[0,1]);

guidata(base, handles);
                     
waitfor(cancelbutton);
index = get(list,'Value');
delete(gcf);
files = [];
indices = [indices,ind];

for cnt = 1:numel(index)
    if return_fileinfo==1
        files = [files,fileinfo(indices(index(cnt)):indices(index(cnt)+1)-1)];
    else
        files = {files,fileinfo(indices(index(cnt)):indices(index(cnt)+1)-1).Name};
    end
end

return
end

function res=buttonCallback(src, evt, buttontype, cancelbutton)
handle = guidata(gcf);

switch buttontype
    case 'Cancel'
        set(handle.list,'Value',[]);
        delete(src);
        %close(gcf);
    case 'OK'
        %result.r = get(lists(1),'Value');
        delete(handle.cancelbutton);
    case 'List'
        index = get(src,'Value');
        top = get(src,'ListboxTop');
        cnt = 1;

        %while handle.lists(cnt)>0
            for cnt = 1:numel(handle.list)
                if handle.list(cnt)>0
            set(handle.list(cnt),'Value',index);
            set(handle.list(cnt),'ListboxTop',top);
            disp([cnt,index, top]);
            %cnt=cnt+1;
                end
        end
    case 'Figure'
        
        set(handle.list(1),'Value',[]);
        delete(handle.cancelbutton);
    case 'ButtonUp'
        disp('In button up');
    otherwise
end

return
end
