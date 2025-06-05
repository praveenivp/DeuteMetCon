classdef rImParams < handle 
    % Class that contains image parameters
    % and methods that use them
    
    properties
        NDims
        MatrixSize
        RealSize
        CenterCoords
        SliceOrient
        Units
        MRParams
    end
    
    methods
        function obj = rImParams(whatever, NSlicesOrRecopars)
            % whatever can be:
            %   nothing: sets default parameters
            %   an array consisting of MatrixSize. The the second argument NSlices is assumed to be just that
            %   a structure, containing some or all of the properties
            %   an image object, from which the parameters are determined
            obj.NDims = 2;
            obj.Units = '';
            obj.MRParams = [];
            if nargin < 2
                NSlices = 1;
                recopars.RearrangeDims = [1,2,3];
            elseif isstruct(NSlicesOrRecopars)
                recopars = NSlicesOrRecopars;
                NSlices = 1;
            else                
                recopars.RearrangeDims = [1,2,3];
                NSlices = NSlicesOrRecopars;
            end
            if nargin == 0                          % nothing given: default values
                obj.MatrixSize = [0,0,0];
                obj.RealSize = [0,0,0];
                obj.CenterCoords = [0,0,0];
                obj.SliceOrient = 'xy';
            elseif isnumeric(whatever) && numel(whatever) <=3   % an array given for MatrixSize
                if numel(whatever) == 1
                    obj.MatrixSize = [whatever, whatever, NSlices];
                    if NSlices == 1
                        obj.NDims = 2;
                    else
                        obj.NDims = 3;
                    end
                    obj.RealSize = [whatever, whatever, NSlices];
                elseif numel(whatever) == 2
                    obj.MatrixSize = [whatever(1), whatever(2), NSlices];
                    if NSlices < 2
                        obj.NDims = 2;
                    else
                        obj.NDims = 3;
                    end
                    obj.RealSize = [whatever(1), whatever(2), NSlices];
                elseif numel(whatever) == 3
                    obj.MatrixSize = whatever;
                    obj.NDims = 3;
                    obj.RealSize = whatever;
                end
            elseif isobject(whatever) && isa(whatever,'ImageClass')
                obj.MatrixSize = whatever.MatrixSize;
                obj.RealSize = whatever.MatrixSize;
                obj.NDims = whatever.NDims;
                obj.CenterCoords = round(obj.MatrixSize/2);
                obj.SliceOrient = 'xy';
                if isfield(whatever, 'MRParams') && isfield(whatever.MRParams, 'Units') && (isstring(whatever.MRParams.Units) || ischar(whatever.MRParams.Units))
                    obj.Units = whatever.Units;
                end
            elseif isstruct(whatever)
                if numel(recopars) > 0 && isstruct(recopars)
                    pars = recopars;
                else
                    pars = whatever;
                end
                if ~isfield(pars, 'MatrixSize')   %The field 'MatrixSize' is obligatory
                    obj.MatrixSize = [0,0,0];
                    obj.RealSize = [0,0,0];
                else
                    if isfield(pars, 'Position') && numel(pars.Position) == 3
                        obj.CenterCoords = pars.Position;
                    else
                        obj.CenterCoords = [0,0,0];
                    end
                    if sum(pars.MatrixSize>1) == 1
                        if isfield(whatever, 'NDims') && whatever.NDims <=3 && whatever.NDims >=2
                            obj.MatrixSize = repmat(pars.MatrixSize,[whatever.NDims,1]);
                            obj.NDims = whatever.NDims;
                            obj.RealSize = pars.FOV;
                        end
                    elseif sum(pars.MatrixSize>1) == 2
                        obj.MatrixSize = pars.MatrixSize([2,1]);
                        obj.NDims = 2;
                        obj.RealSize = pars.FOV([2,1]);
                        if numel(pars.Position(:,3)) > 1    % This is for multislice experiments: make it 3D
                            obj.NDims = 3;
                            obj.MatrixSize(3) = numel(pars.Position(:,3));
                        end
                    elseif sum(pars.MatrixSize>1) == 3
                        obj.MatrixSize = pars.MatrixSize([2,1,3]);
                        obj.NDims = 3;
                        obj.RealSize = pars.FOV([2,1,3]);
                    end
                    if isfield(whatever, 'RealSize')
                        obj.MatrixSize(1:min([numel(whatever.RealSize),obj.NDims])) = whatever.RealSize(1:min([numel(whatever.RealSize),obj.NDims]));
                    end
                    if isfield(whatever, 'Units') && (isstring(whatever.Units) || ischar(whatever.Units))
                        obj.Units = char(whatever.Units);
                    end
                end
            end
        end
        
    end
end

