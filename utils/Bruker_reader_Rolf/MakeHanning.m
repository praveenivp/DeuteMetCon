function mat = MakeHanning( dim,strength, elliptical, offset, Hamming, zeroedges )
% Calculates the Hanning filter in one to three dimensions
% Parameters:
% dim: Array with one to three elements, giving the matrix size in the
% dimensions
% strength: filter strength: 1 is a normal Hanning filter, higher numbers
% weaken the filter
% elliptical: if 1, the Hanning filter is calculated using the absolute
% distance to zero. If 0, the three dimensions are filtered independently.
% offset: offset of the generalized Hamming function. 0.5 is a Hanning function, 0.58 is optimum for 3D
% Default is 1.

% Generates a Hanning function matrix
if nargin < 6
    zeroedges = 0;
end
if nargin < 5
    Hamming = 0;
end
if nargin < 4
    offset = 0;
end
if nargin < 3
elliptical = 1;
end
ndims = numel(dim);
if nargin < 2
    strength = [1,1,1];
else
    if numel(strength)== 1
        strength = [strength, strength, strength];
    end
end
dim2 = 1;
dim3 = 1;
strength2 = 1;
strength3 = 1;
dim1 = dim(1);
strength1 = strength(1);
if ndims >= 2
    dim2 = dim(2);
    strength2 = strength(2);
end
if ndims==3
   dim3 = dim(3);
   strength3 = strength(3);
end

% Generate the matrix
mat = ones(dim1, dim2, dim3);

if elliptical == 1

    % Every element of the matrix should be the distance from the center
    for cnt1 = 1:dim1
        for cnt2 = 1:dim2
            for cnt3 = 1:dim3
                mat(cnt1, cnt2, cnt3) = sqrt((((dim1+1)/2-cnt1)/(strength1*dim1+2))^2 + (((dim2+1)/2-cnt2)/(strength2*dim2+2))^2 + (((dim3+1)/2-cnt3)/(strength3*dim3+2))^2);
            end
        end
    end
    if zeroedges == 1
        edges = ones(size(mat));
        edges(find(mat>0.5)) = 0;
    end
    mat(find(mat>0.5)) = 0.5;  % function is constant outside central sphere

    % Now make the Hanning function
    if Hamming == 0
        mat = 0.5*(1+cos(2*pi*mat));
    else
        if offset == 0
            mat = 0.54+0.46*cos(2*pi*mat);
        else
            mat = offset+(1-offset)*cos(2*pi*mat);
        end
        if zeroedges == 1
            mat = mat.*edges;
        end
    end
else
    %in that case we need three one-dimensional Hanning functions
    %first dimension:
    d = linspace(((dim1+1)/2-1)/(strength1*dim1+2),((dim1+1)/2-1-dim1)/(strength1*dim1+2),dim1);
    if Hamming == 0
        d = 0.5*(1+cos(2*pi*d));
    else
        %d = 0.54+0.46*cos(2*pi*d);
        if offset == 0
            d = 0.54+0.46*cos(2*pi*d);
        else
            d = offset+(1-offset)*cos(2*pi*d);
        end
    end
    for cnt2 = 1:dim2
        for cnt3 = 1:dim3
            mat(:,cnt2, cnt3) =  mat(:,cnt2, cnt3) .* d';
        end
    end
    if ndims>1
    %second dimension:
    d = linspace(((dim2+1)/2-1)/(strength2*dim2+2),((dim2+1)/2-1-dim2)/(strength2*dim2+2),dim2);
    if Hamming == 0
        d = 0.5*(1+cos(2*pi*d));
    else
        %d = 0.54+0.46*cos(2*pi*d);
        if offset == 0
            d = 0.54+0.46*cos(2*pi*d);
        else
            d = offset+(1-offset)*cos(2*pi*d);
        end
    end
    for cnt1 = 1:dim1
        for cnt3 = 1:dim3
            mat(cnt1,:, cnt3) =  mat(cnt1,:, cnt3) .* d;
        end
    end
    end
    if ndims > 2
    %third dimension:
    d = linspace(((dim3+1)/2-1)/(strength3*dim3+2),((dim3+1)/2-1-dim3)/(strength3*dim3+2),dim3);
    if Hamming == 0
        d = 0.5*(1+cos(2*pi*d));
    else
        %d = 0.54+0.46*cos(2*pi*d);
        if offset == 0
            d = 0.54+0.46*cos(2*pi*d);
        else
            d = offset+(1-offset)*cos(2*pi*d);
        end
    end
    for cnt1 = 1:dim1
        for cnt2 = 1:dim2
            mat(cnt1, cnt2,:) =  mat(cnt1, cnt2,:) .* reshape(d,[1,1,numel(d)]);
        end
    end
    end
    
end

if offset ~= 0 && Hamming == 0
    mat = mat*(1-offset)+offset;
end

end

