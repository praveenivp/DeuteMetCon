classdef ImageClass < handle
    %Class containing 2D or 3D images
    %   Detailed explanation goes here
    
    properties
        NDims                     % Number of Dimensions: 2 or 3
        MatrixSize                % Size in pixels. Has NDims Elements
        ImageData                 % Two or three-dimensional array containing the image data
        TrueColorData             % [1D, 2D, 3D, 3] array containing the true color values for the image (only if necessary)
        CurrentSlice              % Contains the current slice
        ColorMap                  % Current colormap of the image
        ImageParams               % An object of class rImParams containing some image parameters
        Coords                    % Image coordinates in real units
    end
    
    methods
        %% Konstruktor: Passes an array or a DataClass object to form the image
        function obj = ImageClass(ImData)   
            if nargin == 0      % No parameters passed: generate empty object
                obj.NDims = 0;
                obj.MatrixSize = [];
                obj.CurrentSlice = 0;
            else
                if isobject(ImData) 
                    if isa(ImData,'DataClass')
                        obj.NDims = numel(find(ImData.RecoParameters.MatrixSize > 1));
                        obj.MatrixSize = ones(3,1);
                        obj.MatrixSize(1:min(numel(ImData.RecoParameters.MatrixSize),3)) = ImData.RecoParameters.MatrixSize(1:min(numel(ImData.RecoParameters.MatrixSize),3));
                        obj.CurrentSlice = round(obj.MatrixSize(3)/2);
                        obj = obj.Data(ImData.ProData);
                    else
                    end
                else
                    obj = obj.Data(ImData);
                    if numel(obj.NDims) == 0
                        return;
                    end
                end
            end
            obj.ColorMap.cm = gray(256);
            obj.ColorMap.clim = [0,0];
            return;
        end
        
        
        %% Method to set the ImageData
        function obj = Data(obj,ImData)
            if nargin == 0
                return;
            end
            if isnumeric(ImData)                    
                ImData = squeeze(ImData);
                s = size(ImData);
                if numel(s) < 2 | numel(s) > 3 | s(1) <3 | s(2) < 4
                    display('This isn''t an image!');
                    return;
                end
                obj.ImageData = ImData;
                obj.NDims = numel(s);
                obj.MatrixSize = ones(1,3); %ones(1,3, 'int16');
                obj.MatrixSize(1:obj.NDims) = s;
                if obj.CurrentSlice < 1
                    obj.CurrentSlice = 1;
                elseif obj.CurrentSlice > obj.MatrixSize(end)
                    obj.CurrentSlice = obj.MatrixSize(end);
                end
                obj.Coords.x = linspace(0,s(2)-1,s(2))-(s(2)-1)/2;
                obj.Coords.y = linspace(0,s(1)-1,s(1))-(s(1)-1)/2;
                if obj.MatrixSize(3) <= 1
                    obj.Coords.z = 0;
                else
                    obj.Coords.z = linspace(0,s(3)-1,s(3))-(s(3)-1)/2;
                end
                obj.Coords.unit = 'px';
            else                                     % The argument is non-numeric: For the moment, generate an empty array
                obj.NDims = 0;
                obj.MatrixSize = [];
                obj.CurrentSlice = 0;                    
            end
        end

        % Method to set the ImageParams property
        % Can be either not existing -> sets default values, derived from current image object
        %                empty array -> removes ImageParams property
        %                ImageParams object -> sets it
        function Params(obj, acqpars,recopars)
            if nargin < 2
                obj.ImageParams = rImParams(obj);
            else
                if numel(acqpars) == 0
                    obj.ImageParams = [];
                else
                    if isobject(acqpars)
                        if isa(acqpars,'rImParams')
                            obj.ImageParams = acqpars;
                        end
                    elseif isstruct(acqpars)
                        if nargin > 1 && isstruct(recopars)
                            obj.ImageParams = rImParams(acqpars,recopars);
                        else
                            obj.ImageParams = rImParams(acqpars);
                        end
                    end
                end
            end
            if isobject(obj.ImageParams)
			% Bruker version:
                %obj.Coords.y = (linspace(1,obj.ImageParams.MatrixSize(1),obj.ImageParams.MatrixSize(1))-(obj.ImageParams.MatrixSize(1)+1)/2)/obj.ImageParams.MatrixSize(1)*obj.ImageParams.RealSize(1) + obj.ImageParams.CenterCoords(1);  % Stimmt das dann novh für Siemens?
                %obj.Coords.x = (linspace(1,obj.ImageParams.MatrixSize(2),obj.ImageParams.MatrixSize(2))-(obj.ImageParams.MatrixSize(2)+1)/2)/obj.ImageParams.MatrixSize(2)*obj.ImageParams.RealSize(2) - obj.ImageParams.CenterCoords(2); 
				% Siemens Version
                % The first versions of the next lines are replaced by lines that implement the 1/2 voxel shift. Is that correct?
                %obj.Coords.y = (linspace(1,obj.ImageParams.MatrixSize(1),obj.ImageParams.MatrixSize(1))-(obj.ImageParams.MatrixSize(1)+1)/2)/obj.ImageParams.MatrixSize(1)*obj.ImageParams.RealSize(1) - obj.ImageParams.CenterCoords(2);
                obj.Coords.y = (linspace(1,obj.ImageParams.MatrixSize(1),obj.ImageParams.MatrixSize(1))-(obj.ImageParams.MatrixSize(1)+2)/2)/obj.ImageParams.MatrixSize(1)*obj.ImageParams.RealSize(1) - obj.ImageParams.CenterCoords(2);
                %obj.Coords.x = (linspace(1,obj.ImageParams.MatrixSize(2),obj.ImageParams.MatrixSize(2))-(obj.ImageParams.MatrixSize(2)+1)/2)/obj.ImageParams.MatrixSize(2)*obj.ImageParams.RealSize(2) + obj.ImageParams.CenterCoords(1);
                obj.Coords.x = (linspace(1,obj.ImageParams.MatrixSize(2),obj.ImageParams.MatrixSize(2))-(obj.ImageParams.MatrixSize(2)+2)/2)/obj.ImageParams.MatrixSize(2)*obj.ImageParams.RealSize(2) + obj.ImageParams.CenterCoords(1);
                if obj.ImageParams.NDims >= 3
                    %obj.Coords.z = (linspace(1,obj.ImageParams.MatrixSize(3),obj.ImageParams.MatrixSize(3))-(obj.ImageParams.MatrixSize(3)+1)/2)/obj.ImageParams.MatrixSize(3)*obj.ImageParams.RealSize(3) + obj.ImageParams.CenterCoords(3);
                    obj.Coords.z = (linspace(1,obj.ImageParams.MatrixSize(3),obj.ImageParams.MatrixSize(3))-(obj.ImageParams.MatrixSize(3)+2)/2)/obj.ImageParams.MatrixSize(3)*obj.ImageParams.RealSize(3) + obj.ImageParams.CenterCoords(3);
                else
                    obj.Coords.z = 0;
                end
                obj.Coords.unit = obj.ImageParams.Units;
            end
        end

        % Method to set the real coordinates
        % realcoords is a structure that can be
        % either realcoords.x, realcoords.y, realcoords.z, where the number of elements of each array has to be equal to the image size
        % of realcoords.dist and realcoords.center, where dist and center are 3-element arrays that contain the distance between two voxels for x,y,z direction, and the coordinates of the image center.
        % realcoords.center can also be missing, in which case the center is 0.
        % There can also be an entry called realcoords.unit, which contains the unit as a string. If missing, no unit will be shown.
        function SetRealCoords(obj, realcoords)
            if nargin == 1 || isstruct(realcoords) == 0
                fprintf('Coordinate information must be a structure with either n-element arrays *.x, *.y, *.z or 3-element arrays *.dist, *.center, plus a unit (*.unit)');
                return
            end
            if isfield(realcoords,'unit')
                obj.Coords.unit = realcoords.unit;
            else
                obj.Coords.unit = '  ';
            end
            if isfield(realcoords,'x') && isfield(realcoords,'y') && numel(realcoords.x) == obj.MatrixSize(1) && numel(realcoords.y) == obj.MatrixSize(2)
                obj.Coords.x = realcoords.x;
                obj.Coords.y = realcoords.y;
                if (isfield(realcoords,'z') && numel(realcoords.z) == obj.MatrixSize(3))|| obj.NDims == 2
                    obj.Coords.z = realcoords.z;
                end  
                return
            end
            if isfield(realcoords,'dist')
                obj.Coords.y = (linspace(1,obj.MatrixSize(1),obj.MatrixSize(1))-(obj.MatrixSize(1)+1)/2)*realcoords.dist(1);
                obj.Coords.x = (linspace(1,obj.MatrixSize(2),obj.MatrixSize(2))-(obj.MatrixSize(2)+1)/2)*realcoords.dist(2);
                obj.Coords.z = (linspace(1,obj.MatrixSize(3),obj.MatrixSize(3))-(obj.MatrixSize(3)+1)/2)*realcoords.dist(3);
            end
            if isfield(realcoords,'center')
            end
        end

        %% Shows the image in a new figure 
        function show(obj, Slice)                        
            if nargin == 1
                Slice = 1;
            end
            if (Slice < 1)
                return
            end
            if obj.NDims == 2
                Slice = 1;
            else
                if Slice > obj.MatrixSize(3)
                    Slice = obj.MatrixSize(3);
                    display(['Showing slice number ',num2str(Slice),'.']);
                end
            end
            imagesc(obj.ImageData(:,:,Slice));colormap(obj.ColorMap.cm);
        end
        
        %% Method to show all slices in a figure window. Steps through slices with <Return>
        function move(obj)                        
            if obj.NDims == 2
                imagesc(obj.ImageData(:,:));
            else
                figure
                for cnt = 1:obj.MatrixSize(3)
                    imagesc(obj.ImageData(:,:,cnt));
                    input(num2str(cnt));                    
                end
            end
        end
        
       %% Method to set the color map as a [256,3] array
       function obj = SetColorMap(obj, cm)                        
           scm = size(cm);
           if scm(2) ~= 3 || scm(1) < 64 || scm(1) > 512
               display('Color map has invalid dimensions!');
               return
           end
           obj.ColorMap.cm = cm;
       end
       
       %% Method to flip the orientation of the image by 90° about one axis
       function obj = FlipOrient(obj, axis)
           scoords = obj.Coords;
           switch abs(axis)
               case 1
                                % an axis through the current plane -> inplane-rotation
                   obj.ImageData = permute(obj.ImageData, [2,1,3]);
                   m = obj.MatrixSize(1:2);
                   obj.MatrixSize(1:2) = [m(2),m(1)];
               case 2        
                                % horizontal axis along the current plane
                   if obj.NDims > 2
                       m = obj.MatrixSize;
                       if axis == 2
                            obj.ImageData = permute(obj.ImageData, [3,1,2]);
                            obj.MatrixSize = [m(3),m(1),m(2)];
                            obj.Coords.z = scoords.x;
                            obj.Coords.x = scoords.y;
                            obj.Coords.y = scoords.z;
                       elseif axis == -2
                           obj.ImageData = permute(obj.ImageData, [2,3,1]);
                            m = obj.MatrixSize;
                            obj.MatrixSize = [m(2),m(3),m(1)];
                            obj.Coords.z = scoords.y;
                            obj.Coords.y = scoords.x;
                            obj.Coords.x = scoords.z;
                       end
                       obj.CurrentSlice = min([obj.CurrentSlice, obj.MatrixSize(3)]);
                   end
               case 3       
                                % vertical axis along current plane
                   if obj.NDims > 2
                       m = obj.MatrixSize;
                       if axis == 3
                            obj.ImageData = permute(obj.ImageData, [3,2,1]);
                       elseif axis == -3
                            obj.ImageData = permute(obj.ImageData, [3,2,1]);
                       end
                       obj.MatrixSize = [m(3),m(2), m(1)]; 
                       obj.CurrentSlice = min([obj.CurrentSlice, obj.MatrixSize(3)]);
                       obj.Coords.z = scoords.y;
                       obj.Coords.y = scoords.z;
                   end                   
               case 4       
                                % vertical axis along current plane
                   if obj.NDims > 2
                       obj.ImageData = permute(obj.ImageData, [2,3,1]);
                       m = obj.MatrixSize;
                       obj.MatrixSize = [m(2),m(3), m(1)];
                       obj.CurrentSlice = min([obj.CurrentSlice, obj.MatrixSize(3)]);
                       obj.Coords.y = scoords.x; 
                       obj.Coords.x = scoords.z;
                       obj.Coords.z = scoords.y;
                   end                   
               case 5       
                                % vertical axis along current plane
                   if obj.NDims > 2
                       m = obj.MatrixSize;
                       if axis == 5
                            obj.ImageData = permute(obj.ImageData, [1,3,2]);
                            obj.MatrixSize = [m(1),m(3), m(2)];
                       elseif axis == -5
                            %obj.ImageData = fliplr(permute(fliplr(obj.ImageData), [1,3,2]));
                            obj.ImageData = permute(obj.ImageData, [1,3,2]);
                            obj.MatrixSize = [m(1),m(3), m(2)];
                       end
                       obj.CurrentSlice = min([obj.CurrentSlice, obj.MatrixSize(3)]);
                       obj.Coords.z = scoords.x;
                       obj.Coords.x = scoords.z;
                       obj.Coords.y = scoords.y;
                   end                   
           end
       end
        
       % Get the real coordinates for a set of point idices
       function res = getrealcoords(obj, indices)
           if isstruct(obj.Coords)
               res = [obj.Coords.x(indices(1)), obj.Coords.y(indices(2)), obj.Coords.z(indices(3))];
           else
               res = [0,0,0];
           end
       end

       % Get the indices from real coordinates 
       function res = getindcoords(obj, realcoords)
           if isstruct(obj.Coords)
               diffx = find(abs(realcoords(1) - obj.Coords.x) == min(abs(realcoords(1) - obj.Coords.x)));
               diffy = find(abs(realcoords(2) - obj.Coords.y) == min(abs(realcoords(2) - obj.Coords.y)));
               diffz = find(abs(realcoords(3) - obj.Coords.z) == min(abs(realcoords(3) - obj.Coords.z)));
               res = [diffx(1),diffy(1),diffz(1)];
           else
               res = [0,0,0];
           end
       end

    end
    
end

