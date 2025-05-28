
function [res, fileinfo] = BrPickDataset(path, return_fileinfo, requiredStatus, title, fileinfo)
%Reads Bruker raw data with one to four dimensions and reconstructs the
%images.
%The parameter path can be either a directory, which serves as starting
%point for the search, or a raw data file. If path is completely missing, a
%file search dialogue is used to find the raw data file.
%requiredStatus can be either 'all', 'fid' or 'image'. With 'all', all scan
%directories are listed, with 'fid' only those where an fid-file exists,
%and with 'image' only those where an 2dseq-file exists. Default = fid

% return_fileinfo can be 0, which means that the first return parameter res contains only the filename of the file to read. 
% The return parameter fileinfo holds the parameters anyway.
% if return_filename is 1, res contains a structure with the main parameters.

res = '';
if nargin < 2
    return_fileinfo = 0;
end
if nargin < 3 | (strcmp(requiredStatus, 'all') && ~strcmp(requiredStatus, 'image'))
    requiredStatus = 'fid';
end
if nargin == 0
    path = '.';
    isdir=1;
else
    isfile = exist(path,'file');
    switch isfile 
        case 0
            path = '.';
            isdir = 1;
        case 2
            isdir = 0;
        case 7
            isdir = 1;
            path = strcat(path,filesep);
            %It may be already the desired scan directory
            if exist([path,'method'],'file') == 2
                if exist([path,'fid'],'file') == 2
                    res = [path,'fid'];
                elseif exist([path,'rawdata.job0'],'file') == 2
                    res = [path,'rawdata.job0'];
                end
                return
            end
        otherwise
            disp('Sorry, could not identify the file!');
            path = '.';
            isdir = 1;
    end
end
if isdir == 0
    if return_fileinfo == 0
        res = path;
    else
        [p,n,e] = fileparts(path);
        path = p;
        if strcmp(n,'2dseq')
            p = fileparts(fileparts(p));
        end
        acqpfile = [p,filesep,'acqp'];
        methodfile = [p,filesep,'method'];
        fileinfo = ReadParams(acqpfile, methodfile);
        fileinfo.filename = n;
        fileinfo.path = path;
        res = fileinfo;
    end
    return;
end

nStudies = 0;
nScans = 0;
if exist('fileinfo') == 1 & isstruct(fileinfo)
    nScans = numel(fileinfo);
    IsScans = 1;
    ReadFileinfo = 0;
    ScanOrder = zeros(nScans,1);
    for cnt = 1:nScans
        ScanOrder(cnt) = str2num(fileinfo(cnt).filename);
        ScanInfo(cnt,1:numel(fileinfo(cnt).filename)) = fileinfo(cnt).filename;
        ScanInfo(cnt,4) = '|';
        if fileinfo(cnt).version ~= 360
            filename{cnt} = [path,filesep,fileinfo(cnt).filename,filesep,'fid'];
        else
            filename{cnt} = [path,filesep,fileinfo(cnt).filename,filesep,'rawdata.job0'];
        end
    end
else
    ReadFileinfo = 1;
%     fileinfo = struct('filename','','path',path,'Number',0,'Date',0,'Dim',0,'Protocol','', ...
%            'Sequence','','NSlices',0.0,'TR',0.0,'TE',0.0,'TI',0.0,'Averages',0, ...
%            'Nucleus','','ImageSize',0,'ACQSize',0,'FOV',0,'Position',0,'TE2',0.0, 'NContrasts',0);
%Check whether the current directory contains studies or scans
numbers = ['1','2','3','4','5','6','7','8','9','0'];
while nScans == 0
%while nStudies == 0 & nScans == 0
    direc = dir(path);
    for cnt = 1:numel(direc)
        s = size(direc(cnt).name);
        if direc(cnt).isdir==1
            if max(s) > 3 & ismember(direc(cnt).name(s(1),s(2)),numbers) & direc(cnt).name(s(1),s(2)-3)=='.'
                nStudies = nStudies+1;
            elseif min(ismember(direc(cnt).name,numbers))==1 & s(2) <=4
                nScans = nScans+1;
            end
        end
    end
    if nScans == 0
    %if nStudies == 0 & nScans == 0
         path = uigetdir(path, 'Select study directory')
        if path == 0 
            return
        end
    end
end
if nScans == 0 %& nStudies == 0 
    return;
end
if nScans == 0 & nStudies > 0 
    IsScans = 0;
end
%if nScans > 0 & nStudies == 0  
if nScans > 0  
    IsScans = 1;
end
%if nScans > 0 & nStudies > 0  
%    IsScans = 0;
%end
end
%% If the directory contains studies: launch study selection widget
%if IsScans == 0
%    display('Studies');
%else
if IsScans > 0

%% If the directory contains scans: launch scan selection widget
    display('Scans');
    index = 1;
    displength = 200;
    protlength = 1;
    psize = zeros(11,1);
    ppos = 1;
if ReadFileinfo == 1
    for cnt = 1:numel(direc) 
        if direc(cnt).isdir & min(ismember(direc(cnt).name,numbers))==1 
            if requiredStatus(1) == 'f'
                if exist([path,filesep,direc(cnt).name,filesep,'fid'])~=2 && exist([path,filesep,direc(cnt).name,filesep,'rawdata.job0'])~=2 && exist([path,filesep,direc(cnt).name,filesep,'ser'])~=2
                    NoScan = 1;
                else                    
                    NoScan = 0;
                end
            end
            if requiredStatus(1) == 'i'
                if exist([path,filesep,direc(cnt).name,filesep,'pdata',filesep,'1',filesep,'2dseq'])~=2
                    NoScan = 1;
                else
                    NoScan = 0;
                end
            end
            if NoScan==0
                %Scan number
                ScanOrder(index) = str2num(direc(cnt).name);
                ScanInfo(index,1:numel(direc(cnt).name)) = direc(cnt).name;
                ScanInfo(index,4) = '|';
                %Protocol
                acqpfile = [path,filesep,direc(cnt).name,filesep,'acqp'];
                methodfile = [path,filesep,direc(cnt).name,filesep,'method'];
                fi = ReadParams(acqpfile, methodfile);
                acqparam = fi.Protocol;
                protlength = max([numel(acqparam),protlength]);
                if numel(acqparam)>0
                    ScanInfo(index,6:numel(acqparam)+5) = acqparam;
                end
                if index == 1
                    fileinfo = fi;
                else
                fileinfo(index) = fi;
                end
                filename{index} = [path,filesep,direc(cnt).name,filesep,'fid'];
                if exist(filename{index},'file')~=2
                    filename{index} = [path,filesep,direc(cnt).name,filesep,'rawdata.job0'];
                    if exist(filename{index},'file')~=2
                        filename{index} = [path,filesep,direc(cnt).name,filesep,'ser'];
                        if exist(filename{index},'file')~=2
                            filename{index}='';
                        end
                    end
                end
                fileinfo(index).filename = direc(cnt).name;
                [a,b,c] = fileparts(filename{index});
                fileinfo(index).datfile = [b,c];
                display(['Scanning dataset ',fileinfo(index).filename]);
                fileinfo(index).path = path;
                fileinfo(index).number = cnt;
            
                index = index+1;
            end

        end
    end
end
    display('Scans');
    position = protlength+7;
    % Create list entry
    for cnt1 = 1:11
        for index = 1:numel(fileinfo)
                switch cnt1
                    case 1
                        if index ==1
                            pmax = max(ScanOrder);
                            pstring = num2str(pmax,'%5.0f');
                            plength = numel(pstring);
                            psize(cnt1) = plength;
                            ParamNames(ppos:ppos+plength) = ' ';
                        end
                        plist(index,ppos:ppos+plength) = ' ';
                        if numel(ScanOrder(index))>0
                            pstring = num2str(ScanOrder(index),'%5.0f');
                            plist(index,ppos+plength-numel(pstring)+1:ppos+plength) = pstring;
                        end
                        
                    case 2
                        if index == 1
                            plength = 0;
                            for cnt2 = 1:numel(fileinfo)
                                if plength<numel(fileinfo(cnt2).Protocol)
                                plength = numel(fileinfo(cnt2).Protocol);
                                end
                            end
                            if plength>0
                                if plength<numel('Protocol')
                                    plength = numel('Protocol');
                                end
                                ParamNames(ppos-1) = '|';
                                ParamNames(ppos:ppos+plength) = ' ';
                                npos = ppos + floor((plength-numel('Protocol'))/2);
                                ParamNames(npos:npos+numel('Protocol')-1) = 'Protocol';
                            end
                            psize(cnt1) = plength;
                        end
                        if plength>0
                            plist(index,ppos-1) = '|';
                            plist(index,ppos:ppos+plength) = ' ';
                        
                            if numel(fileinfo(index).Protocol)>0
                                plist(index,ppos:ppos+numel(fileinfo(index).Protocol)-1) = fileinfo(index).Protocol;
                                if psize(cnt1) < numel(fileinfo(index).Protocol)
                                    psize(cnt1) = numel(fileinfo(index).Protocol);
                                end
                            end
                        end
                    case 3
                        if index == 1
                            plength = 0;
                            for cnt2 = 1:numel(fileinfo)
                                if plength<numel(fileinfo(cnt2).Sequence)
                                plength = numel(fileinfo(cnt2).Sequence);
                                end
                            end
                            if plength>0
                                if plength<numel('Method')
                                    plength = numel('Method');
                                end
                                ParamNames(ppos-1) = '|';
                                ParamNames(ppos:ppos+plength) = ' ';
                                npos = ppos + floor((plength-numel('Method'))/2);
                                ParamNames(npos:npos+numel('Method')-1) = 'Method';
                            end
                            psize(cnt1) = plength;
                        end
                        if plength>0
                            plist(index,ppos-1) = '|';                        
                            plist(index,ppos:ppos+plength) = ' ';
                            if numel(fileinfo(index).Sequence)>0
                                plist(index,ppos:ppos+numel(fileinfo(index).Sequence)-1) = fileinfo(index).Sequence;
                            end
                        end
                    case 4
                        if index ==1
                            pmax = max([fileinfo.TR]);
                            pstring = num2str(pmax,'%5.0f');
                            plength = numel(pstring);
                            if plength>0
                                s = numel('TR');
                                if plength<s
                                    plength = s;
                                end
                                ParamNames(ppos-1) = '|';
                                ParamNames(ppos:ppos+plength) = ' ';
                                npos = ppos + floor((plength-s)/2);
                                ParamNames(npos:npos+s-1) = 'TR';
                            end
                            psize(cnt1) = plength;
                        end
                        if plength>0
                            plist(index,ppos-1) = '|';
                            plist(index,ppos:ppos+plength) = ' ';
                            if numel(fileinfo(index).TR)>0
                                pstring = num2str(fileinfo(index).TR,'%5.0f');
                                plist(index,ppos+plength-numel(pstring)+1:ppos+plength) = pstring;
                            end
                        end
                    case 5
                        if index ==1
                            pmax = max([fileinfo.TE]);
                            pstring = num2str(pmax,'%5.1f');
                            plength = numel(pstring);
                            if plength>0
                                s = numel('TE');
                                if plength<s
                                    plength = s;
                                end
                                ParamNames(ppos-1) = '|';
                                ParamNames(ppos:ppos+plength) = ' ';
                                npos = ppos + floor((plength-s)/2);
                                ParamNames(npos:npos+s-1) = 'TE';
                            end
                            psize(cnt1) = plength;
                        end
                        if plength>0
                            plist(index,ppos-1) = '|';
                            plist(index,ppos:ppos+plength) = ' ';
                            if numel(fileinfo(index).TE)>0
                                pstring = num2str(fileinfo(index).TE,'%5.1f');
                                plist(index,ppos+plength-numel(pstring)+1:ppos+plength) = pstring;
                            end
                        end
                    case 6
                        if index ==1
                            plength = 0;
                            pmax = max([fileinfo.TI]);
                            if numel(pmax) ~= 0
                            pstring = num2str(pmax(1),'%5.1f');
                            plength = numel(pstring);
                            if plength>0
                                s = numel('TI');
                                if plength<s
                                    plength = s;
                                end
                                ParamNames(ppos-1) = '|';
                                ParamNames(ppos:ppos+plength) = ' ';
                                npos = ppos + floor((plength-s)/2);
                                ParamNames(npos:npos+s-1) = 'TI';
                            end
                            end
                            psize(cnt1) = plength;
                        end
                        if plength>0
                            plist(index,ppos-1) = '|';
                            plist(index,ppos:ppos+plength) = ' ';
                            if numel(fileinfo(index).TI)>0
                                pstring = num2str(fileinfo(index).TI(1),'%5.1f');
                                plist(index,ppos+plength-numel(pstring)+1:ppos+plength) = pstring;
                            end
                        end
                    case 7
                        if index ==1
                            pmax = max([fileinfo.Averages]);
                            pstring = num2str(pmax,'%5.0f');
                            plength = numel(pstring);
                            if plength>0
                                s = numel('Ave');
                                if plength<s
                                    plength = s;
                                end
                                ParamNames(ppos-1) = '|';
                                ParamNames(ppos:ppos+plength) = ' ';
                                npos = ppos + floor((plength-s)/2);
                                ParamNames(npos:npos+s-1) = 'Ave';
                            end
                            psize(cnt1) = plength;

                        end
                        if plength>0
                            plist(index,ppos-1) = '|';
                            plist(index,ppos:ppos+plength) = ' ';
                            if numel(fileinfo(index).Averages)>0
                                pstring = num2str(fileinfo(index).Averages,'%5.0f');
                                plist(index,ppos+plength-numel(pstring)+1:ppos+plength) = pstring;
                            end
                        end
                    case 8
                        if index ==1
                            pmax = max([fileinfo.NSlices]);
                            pstring = num2str(pmax,'%5.0f');
                            plength = numel(pstring);
                            if plength>0
                                s = numel('Slices');
                                if plength<s
                                    plength = s;
                                end
                                ParamNames(ppos-1) = '|';
                                ParamNames(ppos:ppos+plength) = ' ';
                                npos = ppos + floor((plength-s)/2);
                                ParamNames(npos:npos+s-1) = 'Slices';
                            end
                            psize(cnt1) = plength;
                        end
                        if plength>0
                            plist(index,ppos-1) = '|';
                            plist(index,ppos:ppos+plength) = ' ';
                            if numel(fileinfo(index).NSlices)>0
                                pstring = num2str(fileinfo(index).NSlices,'%5.0f');
                                plist(index,ppos+plength-numel(pstring)+1:ppos+plength) = pstring;
                            end
                        end
                    case 9
                        if index ==1
                            pmax = max([fileinfo.NContrasts]);
                            pstring = num2str(pmax,'%5.0f');
                            plength = numel(pstring);
                            if pmax <= 1
                                plength = 0;
                            end
                            if plength>0
                                s = numel('Contrasts');
                                if plength<s
                                    plength = s;
                                end
                                ParamNames(ppos-1) = '|';
                                ParamNames(ppos:ppos+plength) = ' ';
                                npos = ppos + floor((plength-s)/2);
                                ParamNames(npos:npos+s-1) = 'Contrasts';
                            end
                            psize(cnt1) = plength;
                        end
                        if plength>0
                            plist(index,ppos-1) = '|';
                            plist(index,ppos:ppos+plength) = ' ';
                            if numel(fileinfo(index).NContrasts)>0 && plength > 0
                                pstring = num2str(fileinfo(index).NContrasts,'%5.0f');
                                plist(index,ppos+plength-numel(pstring)+1:ppos+plength) = pstring;
                            end
                        end
                    case 10
                        if index ==1
                            pmax = max([fileinfo.Dim]);
                            pstring = num2str(pmax,'%5.0f');
                            plength = numel(pstring);
                            if pmax <= 2
                                plength = 0;
                            end
                            if plength>0
                                s = numel('Dims');
                                if plength<s
                                    plength = s;
                                end
                                ParamNames(ppos-1) = '|';
                                ParamNames(ppos:ppos+plength) = ' ';
                                npos = ppos + floor((plength-s)/2);
                                ParamNames(npos:npos+s-1) = 'Dims';
                            end
                            psize(cnt1) = plength;
                        end
                        if plength>0
                            plist(index,ppos-1) = '|';
                            plist(index,ppos:ppos+plength) = ' ';
                            if numel(fileinfo(index).Dim)>0 & plength > 0
                                pstring = num2str(fileinfo(index).Dim,'%5.0f');
                                plist(index,ppos+plength-numel(pstring)+1:ppos+plength) = pstring;
                            end
                        end
                    case 11
                end
        end
        if psize(cnt1)>0
            ppos = ppos + psize(cnt1) +2;
        end
                
    end
            

        
        
    [ScanOrder, order] = sort(ScanOrder);
    plist = plist(order,:);
    filename = {filename{order}};
    fileinfo = fileinfo(order);
    % Create the widgets
    % Create the widget structure
    %  1: The base Figure-object.
    %  2: The scan list
    widgets = zeros(1,15);
    xsize = sum(psize)*11;
    if exist('title') == 1
        name = title;
    else
        name = 'Select Scan';
    end
    widgets(1) = figure('Name',name,'NumberTitle','off','DockControls','off',...
                        'Toolbar','none','Menubar','none','Unit','Pixels','Position',[50,50,xsize,min([numel(fileinfo)*15+55,30*15+55])]);
    s = get(widgets(1), 'Position');

    w = uicontrol('Parent',widgets(1),'Style','text','Unit','Pixels','Position',[2,s(4)-20,xsize,20], ...
                            'String',ParamNames,'FontSize',10,'FontName','FixedWidth', ...
                            'HorizontalAlignment','Left');                    
    widgets(2) = uicontrol('Parent',widgets(1),'Style','listbox','Unit','Pixels','Position',[0,30,xsize,s(4)-50], ...
                            'String',plist,'FontSize',10,'FontName','FixedWidth', ...
                            'HorizontalAlignment','Center', 'Callback',{@buttonCallback,'List',0}, 'Max',10);
    buttons = uipanel('Parent',widgets(1), 'BorderType','line','BorderWidth',2,'Unit','Pixels','Position',[0,0,xsize,30]);
    cancelbutton = uicontrol('Parent',buttons,'Style','pushbutton','String','Cancel','Unit','Pixels',...
                             'Position',[100,5,50,20], 'Callback',{@buttonCallback,'Cancel', 0});
    okbutton = uicontrol('Parent',buttons,'Style','pushbutton','String','OK','Unit','Pixels','Position',[20,5,50,20], ...
                         'Callback',{@buttonCallback,'OK',cancelbutton});
end
waitfor(cancelbutton);
if ishandle(widgets(2))
    index = get(widgets(2),'Value');
    delete(gcf);
    if return_fileinfo == 0
        res = {filename{index}};
    else
        res = fileinfo(index);
    end
else
    res = '';
end


return;
end


function res=buttonCallback(src, evt, buttontype, cancelbutton)
switch buttontype
    case 'Cancel'
        %set(handle.lists(1),'Value',[]);
        delete(src);
        close(gcf);
    case 'OK'
        %result.r = get(lists(1),'Value');
        delete(cancelbutton);
    case 'List'
        %index = get(src,'Value');
        %top = get(src,'ListboxTop');

    otherwise
end
res = 0;
return
end

function f = ReadParams(acqpfile, methodfile)
pnames = {'Method','PVM_EchoTime','PVM_RepetitionTime','PVM_NAverages','PVM_NRepetitions', 'PVM_SelIrInvTime','PVM_MagTransOffset','PVM_MagTransOnOff','PVM_MagTransPower','IRTime'};
acqpnames = {'ACQ_protocol_name','NSLICES','ACQ_dim','ACQ_size','ACQ_jobs','ACQ_fov','NI','ACQ_echo_time','NR','ACQ_spatial_size_1','ACQ_spatial_size_2'};
[acqp, version] = BrReadParams(acqpnames,acqpfile);
[method, version] = BrReadParams(pnames, methodfile);
f.filename = '';
f.version = version;
f.path = '';
f.datfile = '';
f.number = 0;
f.Protocol = acqp.ACQ_protocol_name;
f.Sequence = method.Method;
f.TR = method.PVM_RepetitionTime;
f.TE = method.PVM_EchoTime;
f.Averages = method.PVM_NAverages;
f.TI = method.PVM_SelIrInvTime;
if numel(method.IRTime) > 0 
    f.TI = method.IRTime;
end
if version == 5 || version == 6
    f.ImageSize = acqp.ACQ_size;
elseif version == 360
    f.ImageSize = [acqp.ACQ_jobs(1).ScanSize,acqp.ACQ_spatial_size_1,acqp.ACQ_spatial_size_2];
end
f.ImageSize(1) = f.ImageSize(1)/2;
f.Dim = acqp.ACQ_dim;
f.FOV = acqp.ACQ_fov;
f.NSlices = acqp.NSLICES;
if strcmp(method.PVM_MagTransOnOff, 'On') && method.PVM_MagTransPower >=0.1
    f.MTOffset = method.PVM_MagTransOffset;
else
    f.MTOffset = 0;
end
if numel(acqp.NSLICES)>0
    f.NContrasts = acqp.NI/acqp.NSLICES*acqp.NR;
end
f.TE2 = acqp.ACQ_echo_time;

end
