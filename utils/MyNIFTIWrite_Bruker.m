function [nifti_header]=MyNIFTIWrite_Bruker(Vol,data_struct,filename,Description)
%[nifti_header]=MyNIFTIWrite(Vol_PRS,twix,filename,Description)
%
%[INPUTS]:
% Vol : upto 7D Input volume
%       Please make sure the first three physical dim are in this order: PHASExREADxPARTITION
% data_struct : DataClass object from Rolf Pohmann
% filename : string (optinal)
% Description : string (optional)
%
%Example:
%   [nifti_header]=MyNIFTIWrite(Vol_PRS,twix_obj);
%   MyNIFTIWrite(Image_volume_PRS,twix_obj,'test.nii','Some Documentaion here');
%
%References:
% https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
% https://brainder.org/2012/09/23/the-nifti-file-format/
% https://github.com/aghaeifar/recoMRD
% praveen.ivp

%input paramter handling
if(nargin<3)
    filename=[data_struct.AcqParameters.Filename,'.nii'];
end
if(nargin<4)
    Description='MyNIFTIWrite';
end
if(isfield(data_struct.AcqParameters,'CSIPosition'))
    %permute onlt CSI data accoriding to Pohmann

    %     Vol=permute(Vol,[abs(data_struct.RecoParameters.RearrangeDims),4,5]);
    Vol=permute(Vol,[2,1,3,4,5]);
    %    Vol=ndflip(Vol,find(data_struct.RecoParameters.RearrangeDims<0));
    %  Vol=ndflip(Vol,[1,2,3]);
    FOV_RPS=data_struct.AcqParameters.FOV(abs(data_struct.RecoParameters.RearrangeDims)); %mm
    %   FOV_RPS=data_struct.AcqParameters.FOV(:); %mm
else
    %image you don't need to permute
    FOV_RPS=data_struct.RecoParameters.FOV(:); %mm
    FOV_RPS=FOV_RPS([2 1 3]);
end



MatSz=[size(Vol,1);size(Vol,2);size(Vol,3) ];
Res_RPS=FOV_RPS(:)./MatSz; %mm

AffineMat=getAffineNormal(data_struct,FOV_RPS,MatSz);
trans_vec=AffineMat(1:3,4);

%% try to write the nifti
deafult_header=images.internal.nifti.niftiImage.niftiDefaultHeader(Vol, true, 'NIfTI1');

%modify
deafult_header.pixdim=[1 Res_RPS' ones(1,ndims(Vol)-3)];
deafult_header.sform_code=1;
deafult_header.srow_x=AffineMat(1,:);
deafult_header.srow_y=AffineMat(2,:);
deafult_header.srow_z=AffineMat(3,:);
deafult_header.descrip=Description;
deafult_header.dim_info=bin2dec(strcat('10','01','11'));% (slice,phase,read) -> 3 2 1
deafult_header.xyzt_units=setSpaceTimeUnits('Millimeter', 'Second');
deafult_header.cal_max=max(Vol(:));
deafult_header.cal_max=min(Vol(:));

% qform not tested but should work
deafult_header.qform_code=1;
quat=rotm2quat(AffineMat(1:3,1:3));
deafult_header.quatern_b=quat(2);
deafult_header.quatern_c=quat(3);
deafult_header.quatern_d=quat(4);
deafult_header.qoffset_x=trans_vec(1); %mm
deafult_header.qoffset_y=trans_vec(2); %mm
deafult_header.qoffset_z=trans_vec(3); %mm

deafult_header.scl_slope=1;

%simplify header
NV = images.internal.nifti.niftiImage(deafult_header);
nifti_header=NV.simplifyStruct();
nifti_header.raw=deafult_header;
%
niftiwrite(Vol,filename,nifti_header,'Compressed',false)
end

function affine=getAffineNormal(data_struct,FOV_RPS,MatSize)

% get resolution and FOV
% FOV_RPS=csi_struct.AcqParameters.FOV; %mm
Res_RPS=FOV_RPS(:)./MatSize(1:3); %mm
scaling_affine=diag([Res_RPS(:); 1]);

% get normal

if(isfield(data_struct.AcqParameters,'CSIPosition'))
    Normal=data_struct.AcqParameters.Orientation;
    %Rolf did some permute already
    Normal=[0;0;1];
else
    %Rolf did some permute already
    Normal=[0;0;1];
end

[~,mainOrientation]=max(Normal);
switch(mainOrientation)
    case 1
        Rm1=[[0, 0, 1]; [1, 0, 0]; [0, 1, 0]];
    case 2
        Rm1= [[1, 0, 0]; [0, 0, 1]; [0,-1, 0]];
    case 3
        Rm1=[[0,-1, 0]; [1, 0, 0]; [0, 0, 1]];
end

Normal_main=zeros(3,1);
Normal_main(mainOrientation)=1;
%rodrigues formula
v = cross(Normal_main, Normal) ;
s = norm(v) ;
c = dot(Normal_main, Normal)  ;
if (s <= 1e-5)
    Rm2 = eye(3)*c;
else
    I=eye(3);
    v_x=eye(3);
    for i=1:3
        v_x(i,:)=cross(I(i,:),v);
    end
    Rm2 = eye(3) + v_x + (v_x*v_x)./(1+c);
end

Rotmat_normal=Rm2*Rm1;

%inplane rotation
theta=0;
rotmat_inplane =  [[cos(-theta), sin(+theta), 0];...
    [sin(-theta), cos(-theta), 0];...
    [0         , 0         , 1]];
rotation_affine=eye(4);
rotation_affine(1:3,1:3)=Rotmat_normal*rotmat_inplane;

%%
% translation (center of image)
corner_mm =[-1*FOV_RPS(:)/2 ; 1];
corner_mm(3)=corner_mm(3)+Res_RPS(3)/2;
T=rotation_affine;
% image is bruker reconstructed so need the image center
%csi is not centered so we don't need image center
if(isfield(data_struct.AcqParameters,'CSIPosition'))
    offcenter_SCT=data_struct.AcqParameters.CSIPosition*0;
else
    offcenter_SCT=data_struct.AcqParameters.Position*-1;
end
T(1:3,4)=offcenter_SCT;
offset = T * corner_mm;
translation_affine = eye(4);
translation_affine(:,4) = offset;
% LPS to RAS, Note LPS and PCS (patient coordinate system [Sag, Cor, Tra] ) are identical here (head first/supine).
LPS_to_RAS = eye(4);%diag([-1, -1, 1, 1]); % Flip mm coords in x and y directions

affine = LPS_to_RAS *translation_affine * rotation_affine * scaling_affine;

end


function xyztCode = setSpaceTimeUnits(spaceUnitText, timeUnitText)
%setSpaceTimeUnits: convert simplified space time units to standard
%form.
%    This is a helper function to convert simplified space and time
%    units to the raw format as specified in the NIfTI header.

spaceKey   = {0, 1, 2, 3};
spaceValue = {'Unknown', 'Meter', 'Millimeter', 'Micron'};

spaceMap = containers.Map(spaceValue, spaceKey);

if isempty(find(strcmp(spaceValue, spaceUnitText), 1))
    error(message('images:nifti:spaceUnitNotSupported'));
end

spaceUnits = spaceMap(spaceUnitText);

timeValue = {'None', 'Second', 'Millisecond', 'Microsecond', 'Hertz', 'PartsPerMillion', 'Radian'};
timeKey = {0, 8, 16, 24, 32, 40, 48};

timeMap = containers.Map(timeValue, timeKey);

if isempty(find(strcmp(timeValue, timeUnitText), 1))
    error(message('images:nifti:timeUnitNotSupported'));
end

timeUnits = timeMap(timeUnitText);

spaceUnitCode = bitand(uint8(spaceUnits),uint8(7));
timeUnitCode  = bitand(uint8(timeUnits),uint8(56)); % 0x38

xyztCode = bitor(spaceUnitCode, timeUnitCode);

end