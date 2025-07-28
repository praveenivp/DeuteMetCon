function [brainmask,brain,fname_mask]=nyx_synthstrip(fname_nii,varargin)
% [brainmask,brain,filename]=nyx_synthstrip(fname_nii,desired_orientation);
% perform synstrip on the host with appatainer
% Inputs
% filename - nifti filename
% orientation_desired - 3 char orienation like 'RAS','LPS','ARI' or empty ''
% R/L- right/left, A/P- anterior/posterior I/S-Inferioir/superior

[pn,fn]=fileparts(fname_nii);

[fn,ext]=strtok(fn,'.');

if(isempty(pn))
pn=pwd;
fname_nii=fullfile(pn,fname_nii);
end
fname_stripped=fullfile(pn,[fn,'_strip',ext]);
fname_mask=fullfile(pn,[fn,'_mask',ext]);

% here the synthstrip magic happens
tic
cmd = sprintf('ssh $(hostname) apptainer exec -B /ptmp /ptmp/neurocontainers/synthstrip_7.4.1_20240913.simg mri_synthstrip -i %s -o %s -m %s', fname_nii, fname_stripped, fname_mask);
fprintf('calling synthstrip:%s\n',cmd)
[status,out] = system(cmd);
t=toc;
if status==0
    fprintf('synthstrip successfull! (t=%.2fs)\n', t)
else
    disp(out);
end

brainmask=MyNiftiRead(fname_mask,varargin{:});
brain=MyNiftiRead(fname_nii,varargin{:});

end

function [Vout_all,TransformFucntion]=MyNiftiRead(filename,orientation_desired)
%V=MyNiftiRead(filename,orientation_desired)
%Nifti reader which respects orientaion
%
% Inputs
% filename - nifti filename/file pattern/
% orientation_desired - 3 char orienation like 'RAS','LPS','ARI' or empty ''
% R/L- right/left, A/P- anterior/posterior I/S-Inferioir/superior
%
% Author: praveen.ivp
%
%Examples: Vout=MyNiftiRead(uigetfile('*.nii.gz'),'RAS');
% Vout=MyNiftiRead(dir('*.nii.gz'),'ASR');
%
%test dataset can be found at https://github.com/rordenlab/NIfTIspace.git
if(~isstruct(filename))
    dirst_nii=dir(filename);
    assert(~isempty(dirst_nii),'No Matching with provided pattern');
else
    dirst_nii=filename;
end



Vout_all=cell(length(dirst_nii),1);
TransformFucntion=struct;
for fidx=1:length(dirst_nii)

    cfile=fullfile(dirst_nii(fidx).folder,dirst_nii(fidx).name);
    Vin=niftiread(cfile);
    h=niftiinfo(cfile);

    TT=h.Transform.T;
    res_t=sqrt(diag(tform2rotm(TT)*tform2rotm(TT)'));
    rotm=tform2rotm(TT)./res_t;


    curr_orient=printOrientation(rotm);
    if(~exist('orientation_desired','var')||isempty(orientation_desired))
        orientation_desired=curr_orient;
    else
        fprintf('Converting %s to %s orientation \n',curr_orient,orientation_desired);
    end
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
    Vout_all{fidx}=Vout;
    TransformFucntion(fidx).flip_vec=flip_vec;
    TransformFucntion(fidx).perm_vec=perm_vec;
end

%try to convery cell into matrix
try
    Vout_all=cat(ndims(Vout_all{1})+1,Vout_all{:});
catch
    warning('Volumes cannot be merged')
end

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