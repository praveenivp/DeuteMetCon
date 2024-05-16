function [nifti_header]=MyNIFTIWrite_CSI(Vol,twix,filename,Description,shift_mm)
%[nifti_header]=MyNIFTIWrite(Vol_PRS,twix,filename,Description)
%
%[INPUTS]:
% Vol : upto 7D Input volume
%       Please make sure the first three physical dim are in this order: PHASExREADxPARTITION
% twix : twix object from MAPVBVD
% filename : string (optinal)
% Description : string (optional)
% shift_mm: SCT
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
    filename=[twix.hdr.Config.ProtocolName '.nii'];
end
if(nargin<4)
    Description='MyNIFTIWrite';
end

sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};

%   need to account for POCS and Elliptical scanning
% kspst=twix.hdr.Phoenix.sKSpace;
% volTR=kspst.lPhaseEncodingLines*kspst.lPartitions*twix.hdr.Phoenix.alTR{1}*1e-6/(R_3D*R_PE);

try
    FOV_PRS=[sa.dPhaseFOV sa.dReadoutFOV  sa.dThickness + sa.dThickness*kspst.dSliceOversamplingForDialog]; %mm
catch
    FOV_PRS=[ sa.dPhaseFOV sa.dReadoutFOV sa.dThickness]; %mm
end


MatSz=[size(Vol,1),size(Vol,2),size(Vol,3) ];
Res_PRS=FOV_PRS./MatSz; %mm

% res=sa.dReadoutFOV/kspst.lBaseResolution;
% Res_PRS=[res/kspst.dPhaseResolution   res FOV_PRS(3)/kspst.lPartitions];

AffineMat=getAffineMatrix(twix,[size(Vol),1],shift_mm);
trans_vec=AffineMat(1:3,4);

%% try to write the nifti
deafult_header=images.internal.nifti.niftiImage.niftiDefaultHeader(Vol, true, 'NIfTI1');

%modify
deafult_header.pixdim=[1 Res_PRS ones(1,ndims(Vol)-3)];
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
function affine=getAffineMatrix(twix,MatSize,shift_mm)
sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};
kspst=twix.hdr.Phoenix.sKSpace;

try
    FOV_PRS=[sa.dPhaseFOV sa.dReadoutFOV  sa.dThickness + sa.dThickness*kspst.dSliceOversamplingForDialog]; %mm
catch
    FOV_PRS=[ sa.dPhaseFOV sa.dReadoutFOV sa.dThickness]; %mm
end

if(exist('MatSize'))
    Res_PRS=FOV_PRS./MatSize(1:3);
else
    res=sa.dReadoutFOV/kspst.lBaseResolution;
    Res_PRS=[res/kspst.dPhaseResolution   res FOV_PRS(3)/kspst.lPartitions];
    MatSize=round(FOV_PRS./Res_PRS);
end

 offcenter_SCT=(twix.image.slicePos(1:3,1)+0*shift_mm(:));
% seems like siemens always measure CSI wrt isocenter
% have to check sequence whether they use NCO objects properly
% offcenter_SCT(1:3)=0;
Quat=twix.image.slicePos(4:end,1);

T=quat2tform(Quat(:)');
T(1:3,4)=offcenter_SCT;


% translated from recoMRD.py
% https://github.com/aghaeifar/recoMRD/blob/7e294802ad3ec56fdcf33047c9f92dfddad3f9dc/recoMRD/readMRD.py

R = T(:,1:2) * diag(Res_PRS(1:2));
x1 = [1,1,1,1];
x2 = [1,1,MatSize(3),1];


zmax = (FOV_PRS(3) - Res_PRS(3)) ./ 2;
y1_c = T * [0, 0, -zmax, 1]';
y2_c = T * [0, 0, +zmax, 1]';
% SBCS Position Vector points to slice center this must be recalculated for DICOM to point to the upper left corner.
y1 = y1_c - T(:,1) * FOV_PRS(1)/2 - T(:,2) * FOV_PRS(2)/2;
y2 = y2_c - T(:,1) * FOV_PRS(1)/2 - T(:,2) * FOV_PRS(2)/2;

DicomToPatient =   (cat(1,x1, x2, eye(2,4))\cat(2,y1, y2, R)')';
% Flip voxels in y
AnalyzeToDicom = cat(2,diag([1,-1,1]), [0, (MatSize(2)+1), 0]');
AnalyzeToDicom = cat(1,AnalyzeToDicom, [0,0,0,1]);
% Flip mm coords in x and y directions
PatientToTal   = diag([-1, -1, 1, 1]);
affine         = PatientToTal * DicomToPatient * AnalyzeToDicom;
affine         = affine * cat(2,eye(4,3), [1,1,1,1]'); % this part is implemented in SPM nifti.m

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