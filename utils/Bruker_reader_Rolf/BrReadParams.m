function [arr, version] = BrReadParams(params, path)
% Reads specified parameters from a Bruker Parameter file (p.ex. acqp)
% Opens the file path or, if not given, opens file selection dialogue 
% to select a file, reads it and returns the values of the parameters 
% contained in the string array params.
% For the params-array, the following possibilities exist:
% - A single parameter name: In this case function returns the value of
% this parameter, with the correct type.
% - An arrays of parameters: The function returns a structure, the tag
% names of which are the parameter names.
% The parameter names in the params-array may also contain array indices ()
% or structure tag-numbers (p.ex. 'TPQQ.1').
% If a parameter is an array, the names of the returned array are named
% value1, ....

% Possibilities for variable types:
% String (enum):   ##$Method=<Bruker:FLASH>
% String (string): ##$PVM_TxCoilAdjStat1=( 1024 )
%                  <System (Coil)>
% String (??):     ##$PVM_IsotropicFovRes=Isot_None
% Number:          ##$PVM_GradCalConst=41856.5
% Array (numbers): ##$PVM_MaxMatrix=( 2 )
%                  2048 2048
% Array (strings): ##$ACQ_ReceiverSelect=( 11 )
%                  No Yes Yes No No No No No No No No
% Structure:       ##$ExcPulse1=(1.11111111111111, 5400, 45, Yes, 4, 6000, 0.163826239916836, 
%                  0.142428111523141, 0, 50, 8.70591848971339, <$ExcPulse1Shape>)
% Array (structs): ##$ACQ_jobs=( 1 )
%                  (1536, 70, 280, 71680, 101, 151515.151515152, 71680, 1, <job0>)







% Written by Rolf Pohmann, Max Planck Institute for Biological Cybernetics,
% Tübingen, Germany

if nargin < 1
    disp('Reads Bruker Parameter files and returns parameter values.');
    disp('Usage: result = BrReadParams(ParameterArray, Filename');
    arr = -1;
    return;
end
if ~isequal(class(params), 'char') &  ~isequal(class(params), 'cell')
    disp('Reads Bruker Parameter files and returns parameter values.');
    disp('Usage: result = BrReadParams(ParameterArray, Filename');
    arr = -1;
    return;
end
if isequal(class(params), 'char')
    params = {params};
end
if nargin < 2
    [f,p] = uigetfile([{'acqp;method;imnd;subject' ...
            'Bruker parameter files';'*.*','All Files'}], ... 
            'Select parameter-file for reading')
    if f == 0
        arr = -1;
        return;
    end
    path = strcat(p,f);
end
%% Open file 
file = fopen(path,'r');
if file == -1
    arr = -2;
    return
end
% The first line is used to ensure that the file has the right type and to determine the ParaVision version
line = strtrim(fgets(file));
if strncmp(line,'##TITLE=Parameter List',22)
    if numel(line) > 25
        vs = line(25:end);
        if strcmp(vs(1:10),'ParaVision')
            version = sscanf(vs(11:end),'%d');
        else
            version = 0;
        end
    else
        version = 5;
    end
else
    display('Not a Bruker parameter file');
end

line = 0;
paramfound = zeros(numel(params),1);
if numel(params)>1
    arr.version = version;
%     for cnt = 1: numel(params)
%         arr.(params{cnt})={};
%     end
end
%% Read file line by line
linenum = 0;
while line ~= -1
  line = fgets(file);
  linenum = linenum + 1;
  if line ~= -1
    line = strtrim(line);
    %Does the line contain a parameter name?
    if (strfind(line, '##$'))==1
        [pname, pvalue] = strtok(line(4:end), '=');
        %Does this parameter appear in the params-array?
        nfoundparams = 0;
          for cnt=1:numel(params)

            if strcmpi(pname,strtok(params{cnt},'([.'))   %Parameter found
               nfoundparams = nfoundparams + 1;
               %disp([pname,pvalue]);
               if nfoundparams > 1
                   res = saveres;
               else
                   res = [];
                   IsStruct = 0;
                   IsArray = 0;
                   %What type of parameter is it?
                   pvalue = pvalue(2:end);  %Remove =
                   pvalue = strtrim(pvalue);
                   
                   if pvalue(1) == '('          % Array or structure
                        % If the closing bracket is not on the same line, it's certainly a structure
                        % Otherwise, read the start of the next line. If it's the next parameter, then it's a structure
                        if pvalue(end) == ')'
                            l = ' ';
                            ar = [];
                            while l(1) ~= '#' && l(1) ~= '$' && ~feof(file)
                                pos = ftell(file);
                                l = strtrim(fgets(file));
                                if l(1) == '#' || l(1) == '$'
                                    fseek(file,pos,'bof');
                                else
                                    ar = [ar,' ', l];
                                end
                            end
                            if numel(ar) == 0         % Structure
                                IsStruct = 1;
                            else
                                IsArray = 1;
                            end
                        else
                            IsStruct = 1;
                            stopread = 0;
                            % The structure is longer than one line: Read till closing bracket
                            while pvalue(end) ~= ')' && stopread == 0 && ~feof(file)
                               pos = ftell(file);
                               l = strtrim(fgets(file));
                               if numel(l) == 0 || l(1) == '#' || l(1) == '$'
                                   stopread = 1;
                                   fseek(file,pos,'bof');
                               else
                                    pvalue = [pvalue,l];
                               end
                            end
                        end
                        if IsStruct == 1
                                pvalue = pvalue(2:end-1);
                                res = split(pvalue,',');
                                for elcnt = 1:numel(res)
                                    if isnan(str2double(res{elcnt}))
                                        res{elcnt} = strtrim(res{elcnt});
                                        if res{elcnt}(1) == '<'
                                            res{elcnt} = res{elcnt}(2:end-1);
                                        end
                                    else
                                        res{elcnt} = str2double(res{elcnt});
                                    end
                                end
                        else                      % array
                            %Determine the array dimensions
                            NumDim = str2num(pvalue(2:end-1));
                            if numel(NumDim) == 1
                                NumDim = [1,NumDim];
                            end
                            ar = strtrim(ar);
                            if ar(1) == '('      % array of structures
                                ar = ar(2:end-1);
                                ar = split(ar,') (');
                                nels = numel(ar);
                                for elscnt = 1:nels
                                    a = split(ar{elscnt},', ');
                                    for acnt = 1:numel(a)
                                        if isnan(str2double(a{acnt}))
                                            res{elscnt, acnt} = a{acnt};
                                        else
                                            res{elscnt, acnt} = str2double(a{acnt});
                                        end
                                    end
                                end
                                IsStruct = 1;
                            elseif ar(1) == '<'  % A string (character array) enclosed in < .. >
                                res = ar(2:end-1);
                            else
                                ar = split(ar);
                                numarr = str2double(ar);
                                nonum  = isnan(numarr);
                                if max(nonum) == 0    % only numbers
                                    res = numarr;
                                else
                                    nonumind = find(nonum > 0);
                                    indcnt = 1;
                                    offset = 0;
                                    while indcnt <= numel(nonumind)
                                        if regexp(ar{nonumind(indcnt)}, '@\d*\*(')
                                            q = textscan(ar{nonumind(indcnt)},'@%d*(%f)');     % Only for numbers
                                            if numel(numarr) == 1
                                                numarr = ones(q{1},1)*double(q{2});
                                            elseif nonumind(indcnt) == 1
                                                numarr = [ones(q{1},1)*double(q{2});numarr(nonumind(indcnt)+1:end)];
                                            elseif nonumind(indcnt) == numel(ar)
                                                numarr = [numarr(1:nonumind(indcnt)-1+offset);ones(q{1},1)*double(q{2})];
                                            else
                                                numarr = [numarr(1:nonumind(indcnt)-1+offset);ones(q{1},1)*double(q{2});numarr(nonumind(indcnt)+1+offset:end)];
                                            end
                                            offset = offset + q{1} - 1;
                                            indcnt = indcnt+1;
                                            res = numarr;
                                        else                             % array of strings
                                            res = ar;
                                            indcnt = numel(nonumind)+1;
                                        end
                                    end
                                end
                                res = reshape(res,NumDim);
                                res = squeeze(res);
                            end
                                
                       end
                   elseif pvalue(1) == '<'      % String
                       res = pvalue(2:end-1);
                   else                         % Number or single string
                       num = str2double(pvalue);
                       if ~isnan(num)        % Number
                           res = num;
                       else                     %String
                           res = pvalue;
                       end
                   end
                   
                   
                   
                   
                   
                   
%                    if pvalue(1)~='('    %It's not an array: number or String
%                        num = find((pvalue>57 | pvalue<43 | pvalue==44) ...
%                            & pvalue~=101 & pvalue~=69);
%                        if numel(num)==0    %It's a number
%                            %disp([pname,': ',pvalue]);
%                            res = str2double(pvalue);
%                        else %end of parameter is a single number
%                            res = pvalue;                       
%                        end  %end of parameter is a single string
%                    else     %end of parameter is a scalar
%                        %parameter is an array or a structure
%                        %Is it an array?
%                        if numel(find((pvalue(2:end-1)>57 | pvalue(2:end-1)<43) ...
%                                & pvalue(2:end-1)~=32 & pvalue(2:end-1)~=44)) == 0 ...
%                                & pvalue(end)==')'
%                            NumDim = str2num(pvalue(2:end-1));
%                            if numel(NumDim) == 1
%                                NumDim = [1,NumDim];
%                            end
%                            Ind = 0;
%                            r = '';
%                            res = '';
%                            while Ind < prod(NumDim)
%                                l = strtrim(fgets(file));
%                                %Determine type
%                                if numel(res)==0
%                                     if numel(find((l>57 | l<40) & l~=32 & l~=44 & l~=69 & l~=101 & l~=64)) == 0  %only numbers
%                                         IsNum = 1;
%                                         res = zeros(NumDim);
%                                     elseif l(1)=='<'    %A single string
%                                         IsNum = -1;
%                                         res = '';
%                                     elseif l(1)=='('    %An array of structures
%                                         IsNum = -2;
%                                         res = {prod(NumDim)};
%                                         IsStruct = 1;
%                                     else                %Any other kind of array
%                                         IsNum = 0;
%                                         res = cell(NumDim);
%                                     end
%                                     res = squeeze(res)';
%                                end
%                                while numel(l) > 0
%                                    if IsNum==0
%                                        [res{Ind+1}, l] = strtok(l);
%                                        Ind = Ind+1;
%                                    elseif IsNum == -1
%                                        if Ind==0 & l(end)=='>'
%                                            res = l(2:end-1)';
%                                            l='';Ind =prod(NumDim); 
%                                        elseif Ind==0 & l(end)~='>'
%                                            res = l(2:end)';
%                                        elseif Ind~=0 & l(end)=='>'
%                                            res = [res,l(1:end-1)];
%                                            l='';Ind =prod(NumDim);
%                                        elseif Ind~=0 & l(end)~='>'
%                                            res = [res,l(1:end)];
%                                        end
%                                        Ind = Ind+1;
%                                     elseif IsNum==-2
%                                         if numel(r)>0
%                                             l = [r,l];
%                                         end
%                                        [r, l] = strtok(l(strfind(l,'(')+1:end),')');
%                                        if numel(l)>0
%                                            if Ind == 0
%                                                res = BrConvertToStruc(r); 
%                                            else
%                                                res(Ind+1) = BrConvertToStruc(r);
%                                            end
%                                            r='';
%                                            Ind = Ind+1;
%                                        else
%                                            r = ['(',r];
%                                        end
%                                     elseif IsNum==1
%                                        [r, l] = strtok(l);
%                                        if r(1) == '@'
%                                           q = textscan(r,'@%d*(%d)');
%                                           res(Ind+1:Ind+q{1}) = q{2};
%                                           Ind = Ind+q{1};
%                                        else
%                                             res(Ind+1) = str2num(r);
%                                             Ind = Ind+1;
%                                        end                                       
%                                    end
% 
%                                end
%                         
%                            end
%                            res = res';
%                        else
%                            IsStruct=1;
%                            while (pvalue(end) ~= ')')
%                                l = strtrim(fgets(file));
%                                pvalue = [pvalue,l];
%                            end
%                        end  %end of parameter is an array else ...
%                    
%               end  %end of parameter is array or structure
            end
            saveres = res;
            nameappend = '';
            %If structure: Check whether the params{cnt} contains a .
            if IsStruct==1 
                    % Determine the names of the structure fields
                    sres = size(res);
                    if sres(2) > sres(1)
                        res = res';
                        sres = size(res);
                    end
                    if contains(params{cnt},'ACQ_jobs') && sres(1) == 9
                        fields = {'ScanSize','transactionBlocks','dummyScans','nTotalScans','receiverGain','swh','nStoredScans','ChanNum','title'};
                    elseif contains(params{cnt},'Pulse') && sres(1) == 12
                        fields = {'Length','Bandwidth','Flipangle','Calculated','Sharpness','Bwfac','Sint','Pint','Type','Rpfac','Pow','Shape'};
                    elseif contains(params{cnt},'Pulse') && contains(params{cnt},'Ampl') && sres(1) == 3
                        fields = {'ppow','pampl','patt'};
                    else
                        fields = strcat({'val'},num2str(linspace(1,sres(1),sres(1))','%02d'));
                    end
                    res = cell2struct(res, fields);
%                    pvalue = pvalue(2:end-1);
%                     ind = 1;
%                     res = struct;
%                     while numel(pvalue)>0
%                         [n,pvalue] = strtok(pvalue,',');
%                         [num, isnum] = str2num(n);
%                         if isnum == 0
%                             res.(['val',num2str(ind)]) = n;
%                         else
%                             res.(['val',num2str(ind)]) = num;
%                         end
%                         ind = ind+1;
%                     end
            end
            %Return result
            % If the params-Entry contains an index: Only return
            % this element
            [p,ind] = strtok(params{cnt},'([');
            if numel(ind) > 1
            ind = strtok(ind(2:end),'])');
            index = str2num(ind);
            if index<-1
                index=-1;
            end
            if index>numel(res)-1
                index=numel(res)-1;
            end
            if iscell(res)
                res = res{index+1};
            else
                res = res(index+1);
            end
            params{cnt} = [strtok(params{cnt},'(['),'_',num2str(index)];
            end
            
            if (numel(res)==1)
              if (numel(params)==1)
                  arr = res;
                  fclose(file);
                  return
              else
                  arr.(params{cnt})=res;
              end
            else
              %Add data to output structure of return value
              if (numel(params)==1)
                arr = res;
                fclose(file);
                return
              else
                arr.([strtok(params{cnt},'[(.'),nameappend]) = res;
              end
                
            end
            paramfound(cnt)=1;   
            end  %end of parameter found
           
          end      %end of loop through parameters
      end          %end of line contains parameter
  end  %end of next line read
end  %end of loop through file  
%Are there any not-found parameters?
if numel(params)==1
    arr = [];
else
    for cnt=1:numel(params)
        if paramfound(cnt)==0
            arr.(strtok(params{cnt},'[(.'))=[];
        end
    end
end
fclose(file);
return

function struc = BrConvertToStruc(String)

cnt = 1;
while numel(String)>0
    [item, String] = strtok(String,',');
    item = strtrim(item);
    if item(1) == '<'         %String
        item = item(2:end-1);
    elseif  numel(find((item>57 | item<43 | item==44) ...
                       & item~=101 & item~=69))==0;     %It's a number
        item = str2num(item);
    end
    struc.(['value',num2str(cnt)]) = item;
    cnt = cnt+1;
end    

return


