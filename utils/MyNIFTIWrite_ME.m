function [nifti_header]=MyNIFTIWrite_ME(Vol,twix,filename,Description)
%[nifti_header]=MyNIFTIWrite_ME(Vol_PRS,twix,filename,Description)
%
%[INPUTS]:
% Vol : upto 7D Input volume
%       Please make sure the first three physical dim are in this order: PHASExREADxPARTITION
% twix : twix object from MAPVBVD
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
    filename=[twix.hdr.Config.ProtocolName '.nii'];
end
if(nargin<4)
    Description='MyNIFTIWrite';
end

sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};

%   need to account for POCS and Elliptical scanning
 kspst=twix.hdr.Phoenix.sKSpace;
% volTR=kspst.lPhaseEncodingLines*kspst.lPartitions*twix.hdr.Phoenix.alTR{1}*1e-6/(R_3D*R_PE);

try
    FOV_PRS=[sa.dPhaseFOV  sa.dReadoutFOV sa.dThickness + sa.dThickness*kspst.dSliceOversamplingForDialog]; %mm
catch
    FOV_PRS=[ sa.dPhaseFOV  sa.dReadoutFOV sa.dThickness]; %mm
end


MatSz=[size(Vol,1),size(Vol,2),size(Vol,3) ];
Res_PRS=FOV_PRS./MatSz; %mm

% res=sa.dReadoutFOV/kspst.lBaseResolution;
% Res_PRS=[res/kspst.dPhaseResolution   res FOV_PRS(3)/kspst.lPartitions];



AffineMat=getAffineNormal(twix,[size(Vol),1]);
Vol=flip(Vol,2); % cartesian
% Vol=ndflip(Vol,[1 2 3]); % spiral


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
function affine=getAffineNormal(twix,MatSize)
% get affine matrix from nornal vector
% a little more cleaner code without dicom bullshit
% https://raw.githubusercontent.com/aghaeifar/RecoTwix/main/recotwix/transformation.py
% get resolution and FOV
sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};
try
    FOV_PRS=[sa.dPhaseFOV  sa.dReadoutFOV  sa.dThickness + sa.dThickness*kspst.dSliceOversamplingForDialog]; %mm
catch
    FOV_PRS=[ sa.dPhaseFOV  sa.dReadoutFOV sa.dThickness]; %mm
end

if(exist('MatSize','var'))
    Res_PRS=FOV_PRS./MatSize(1:3);
else
    res=sa.dReadoutFOV/kspst.lBaseResolution;
    Res_PRS=[res/kspst.dPhaseResolution   res FOV_PRS(3)/kspst.lPartitions];
    MatSize=round(FOV_PRS./Res_PRS);
end
scaling_affine=diag([Res_PRS(:); 1]);

% get normal
Normal=zeros(3,1);
if(isfield(sa.sNormal,'dSag')   ),Normal(1)=sa.sNormal.dSag;end
if(isfield(sa.sNormal,'dCor')   ),Normal(2)=sa.sNormal.dCor;end
if(isfield(sa.sNormal,'dTra')   ),Normal(3)=sa.sNormal.dTra;end
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
if(isfield(sa,'dInPlaneRot')), theta=sa.dInPlaneRot;end
  rotmat_inplane =  [[cos(-theta), sin(+theta), 0];...
                    [sin(-theta), cos(-theta), 0];...
                    [0         , 0         , 1]];
rotation_affine=eye(4);
rotation_affine(1:3,1:3)=Rotmat_normal*rotmat_inplane;

%%
   % translation (center of image)
    corner_mm =[-1*FOV_PRS(:)/2 ; 1];
    corner_mm(3)=corner_mm(3)+Res_PRS(3)/2;
    T=rotation_affine;
    offcenter_SCT=twix.image.slicePos(1:3,1);
    T(1:3,4)=offcenter_SCT;
    offset = T * corner_mm;
    translation_affine = eye(4);
    translation_affine(:,4) = offset;
    % LPS to RAS, Note LPS and PCS (patient coordinate system [Sag, Cor, Tra] ) are identical here (head first/supine).
    LPS_to_RAS = diag([-1, -1, 1, 1]); % Flip mm coords in x and y directions

    affine = LPS_to_RAS *translation_affine * rotation_affine * scaling_affine;

end


function affine=getAffineMatrix(twix,MatSize)
sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};
kspst=twix.hdr.Phoenix.sKSpace;

try
    FOV_PRS=[sa.dPhaseFOV  sa.dReadoutFOV  sa.dThickness + sa.dThickness*kspst.dSliceOversamplingForDialog]; %mm
catch
    FOV_PRS=[ sa.dPhaseFOV  sa.dReadoutFOV sa.dThickness]; %mm
end

if(exist('MatSize','var'))
    Res_PRS=FOV_PRS./MatSize(1:3);
else
    res=sa.dReadoutFOV/kspst.lBaseResolution;
    Res_PRS=[res/kspst.dPhaseResolution   res FOV_PRS(3)/kspst.lPartitions];
    MatSize=round(FOV_PRS./Res_PRS);
end

offcenter_SCT=twix.image.slicePos(1:3,1);
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
% https://github.com/robotology-legacy/mex-wholebodymodel/blob/master/mex-wholebodymodel/matlab/utilities/+WBM/+utilities/quat2rotm.m
% function TT = quat2tform(quat)
%     if (length(quat) ~= 4)
%         error('quat2rotm: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
%     end
%     dcm = zeros(3,3);
% 
%     qnorm = norm(quat);
%     if (qnorm > 1)
%         quat = quat./qnorm; % normalize
%     end
%     % scalar (real) part:
%     q_0 = quat(1);
%     % vector (imaginary) part:
%     q_1 = quat(2);
%     q_2 = quat(3);
%     q_3 = quat(4);
% 
%     %% Compute the Direction Cosine Matrix (DCM) by applying the Euler-Rodrigues Parameterization
%     % for efficent computing of the rotation matrix R(s,r) in SO(3) with,
%     %
%     %       R(s,r) = I + (2*s)*S_r + 2*(S_r*S_r),
%     %
%     % where s denotes the scalar part of the given quaternion, S_r is the skew-symmetric matrix of
%     % the vector part r of the quaternion, and I is a (3x3) identity matrix.
%     % For more details about the parametrization formula, please check:
%     %   [1] CONTRIBUTIONS TO THE AUTOMATIC CONTROL OF AERIAL VEHICLES, Minh Duc HUA, PhD-Thesis, 2009,
%     %       <https://www-sop.inria.fr/act_recherche/formulaire/uploads/phd-425.pdf>, p. 101, eq. (3.7) & (3.8).
%     %   [2] Lecture notes for 'Multibody System Dynamics', Prof. Bottasso, Aerospace Engineering Departement, Politecnico di Milano,
%     %       <https://home.aero.polimi.it/trainelli/downloads/Bottasso_ThreeDimensionalRotations.pdf>, p. 33, eq. (1.160).
%     %   [3] The Vectorial Parameterization of Rotation, Olivier A. Bauchau & Lorenzo Trainelli, Nonlinear Dynamics, 2003,
%     %       <http://soliton.ae.gatech.edu/people/obauchau/publications/Bauchau+Trainelli03.pdf>, p. 16, eq. (51).
%     %   [4] Theory of Applied Robotics: Kinematics, Dynamics, and Control, Reza N. Jazar, 2nd Edition, Springer, 2010, p. 103, eq. (3.82).
%     %
%     % Note: The direct assignment is computationally faster than using directly the formula above. Furthermore, if the above formula
%     %       will be used, the matrix multiplications are causing accumulative round-off errors in the main diagonal of the DCM-matrix
%     %       (around ??0.7*1e-03). To be more accurate, it is better to use the derived transformation matrix of the formula.
%     dcm(1,1) = q_0*q_0 + q_1*q_1 - q_2*q_2 - q_3*q_3;
%     dcm(1,2) = 2*(q_1*q_2 - q_0*q_3);
%     dcm(1,3) = 2*(q_0*q_2 + q_1*q_3);
% 
%     dcm(2,1) = 2*(q_0*q_3 + q_1*q_2);
%     dcm(2,2) = q_0*q_0 - q_1*q_1 + q_2*q_2 - q_3*q_3;
%     dcm(2,3) = 2*(q_2*q_3 - q_0*q_1);
% 
%     dcm(3,1) = 2*(q_1*q_3 - q_0*q_2);
%     dcm(3,2) = 2*(q_0*q_1 + q_2*q_3);
%     dcm(3,3) = q_0*q_0 - q_1*q_1 - q_2*q_2 + q_3*q_3;
% 
% 
%     TT=zeros(4,4);
%     TT(1:3,1:3)=dcm;
%     TT(4,4)=1;
% end
% % https://github.com/robotology-legacy/mex-wholebodymodel/blob/master/mex-wholebodymodel/matlab/utilities/%2BWBM/%2Butilities/rotm2quat.m
% function quat = rotm2quat(rotm)
%     if ( (size(rotm,1) ~= 3) || (size(rotm,2) ~= 3) )
%         error('rotm2quat: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
%     end
%     quat    = zeros(4,1);
%     epsilon = 1e-12; % min. value to treat a number as zero ...
% 
%     %% Compute the corresponding (unit) quaternion from a given rotation matrix R:
%     % The transformation uses the computational efficient algorithm of Stanley.
%     % To be numerically robust, the code determines the set with the maximum
%     % divisor for the calculation.
%     % For further details about the Stanley Algorithm, see:
%     %   [1] Optimal Spacecraft Rotational Maneuvers, John L. Junkins & James D. Turner, Elsevier, 1986, pp. 28-29, eq. (2.57)-(2.59).
%     %   [2] Theory of Applied Robotics: Kinematics, Dynamics, and Control, Reza N. Jazar, 2nd Edition, Springer, 2010, p. 110, eq. (3.149)-(3.152).
%     % Note: There exist also an optimized version of the Stanley method and is the fastest
%     %       possible computation method for Matlab, but it does not cover all special cases.
%     % Further details about the fast calculation can be found at:
%     %   [3] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
%     %       p. 36, formula (2.30).
%     tr = rotm(1,1) + rotm(2,2) + rotm(3,3);
%     if (tr > epsilon) % if tr > 0:
%             % scalar part:
%             quat(1,1) = 0.5*sqrt(tr + 1);
%             s_inv = 1/(quat(1,1)*4);
%             % vector part:
%             quat(2,1) = (rotm(3,2) - rotm(2,3))*s_inv;
%             quat(3,1) = (rotm(1,3) - rotm(3,1))*s_inv;
%             quat(4,1) = (rotm(2,1) - rotm(1,2))*s_inv;
%     else % if tr <= 0, find the greatest diagonal element for calculating
%          % the scale factor s and the vector part of the quaternion:
%         if ( (rotm(1,1) > rotm(2,2)) && (rotm(1,1) > rotm(3,3)) )
%             quat(2,1) = 0.5*sqrt(rotm(1,1) - rotm(2,2) - rotm(3,3) + 1);
%             s_inv = 1/(quat(2,1)*4);
% 
%             quat(1,1) = (rotm(3,2) + rotm(2,3))*s_inv;
%             quat(3,1) = (rotm(2,1) + rotm(1,2))*s_inv;
%             quat(4,1) = (rotm(3,1) + rotm(1,3))*s_inv;
%         elseif (rotm(2,2) > rotm(3,3))
%             quat(3,1) = 0.5*sqrt(rotm(2,2) - rotm(3,3) - rotm(1,1) + 1);
%             s_inv = 1/(quat(3,1)*4);
% 
%             quat(1,1) = (rotm(1,3) - rotm(3,1))*s_inv;
%             quat(2,1) = (rotm(2,1) + rotm(1,2))*s_inv;
%             quat(4,1) = (rotm(3,2) + rotm(2,3))*s_inv;
%         else
%             quat(4,1) = 0.5*sqrt(rotm(3,3) - rotm(1,1) - rotm(2,2) + 1);
%             s_inv = 1/(quat(4,1)*4);
% 
%             quat(1,1) = (rotm(2,1) - rotm(1,2))*s_inv;
%             quat(2,1) = (rotm(3,1) + rotm(1,3))*s_inv;
%             quat(3,1) = (rotm(3,2) + rotm(2,3))*s_inv;
%         end
%     end
% end
