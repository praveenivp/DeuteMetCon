function Vout=MyNiftiRead(filename,orientation_desired)
%V=MyNiftiRead(filename,orientation_desired)
%Nifti reader which respects orientaion
%
% Inputs
% filename - 3D volume
% orientation_desired - 3 char orienation like 'RAS','LPS','ARI'
% R/L- right/left, A/P- anterior/posterior I/S-Inferioir/superior
%
% Author: praveen.ivp
%
%Examples: Vout=MyNiftiRead(uigetfile('*.nii.gz'),'RAS');
%
%test dataset can be found at https://github.com/rordenlab/NIfTIspace.git

Vin=niftiread(filename);
h=niftiinfo(filename);

TT=h.Transform.T;
res_t=sqrt(diag(tform2rotm(TT)*tform2rotm(TT)'));
rotm=tform2rotm(TT)./res_t;


curr_orient=printOrientation(rotm);
fprintf('Converting %s to %s orientation \n',curr_orient,orientation_desired);
perm_vec=[zeros(1,3),4:ndims(Vin)];
flip_vec=zeros(1,3);
for i=1:3
    flip_vec(i)=~any(curr_orient(i)==orientation_desired);
    switch(orientation_desired(i))
        case {'A','P'}
            perm_vec(i)=find(curr_orient=='A'|curr_orient=='P');
        case {'R','L'}
           perm_vec(i)=find(curr_orient=='R'|curr_orient=='L');
        case {'I','S'}
            perm_vec(i)=find(curr_orient=='I'|curr_orient=='S');
    end
end
flip_vec=find(flip_vec>0);
Vout=(permute(ndflip(Vin, flip_vec),perm_vec));

end

function str=printOrientation(rotm)
cv={'R','A','S','L','P','I'};
str=[];
for i=1:3
[~,idx]=max(abs(rotm(i,:)));
% if negative 'L','P','I'
idx=idx+(rotm(i,idx)<0)*3; 
str=[str,cv{idx}];
end

end

function ndMat=ndflip(ndMat,dim)
for i=dim
    ndMat=flip(ndMat,i);
end

end