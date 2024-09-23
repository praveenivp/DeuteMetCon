function [sig_shifted,pos_PRS]= performFOVShift(sig,twix)
% [sig_shifted,pos_PRS]= performFOVShift(sig,twix)
% Fucntion for correction off-isocenter volume positions if it is not corrected in during acquistion 

sliceInfo=getSlicePosition(twix);

imSize=size(sig,2:4);
res_PRS=sliceInfo{1}.FOV_PRS./imSize; %mm
kmax_PRS=2*pi*(0.5./(res_PRS*1e-3)); %rad/m

[RM,iRM]=getRotationMatrices(twix);
HFS=[1 0 0; 0 -1 0 ; 0 0 -1]; % head first-supine
pos_PRS=(HFS*RM*iRM)\(1e-3*[-1; 1; 1;].*sliceInfo{1}.Position(:));

%                 pos_PRS=GradientXYZ2PRS(1e-3*[-1 1 1].*sliceInfo{1}.Position,twix) %only work for head first-supine

[x,y,z]=ndgrid(linspace(-0.5,0.5-1/imSize(1),imSize(1)), ...
               linspace(-0.5,0.5-1/imSize(2),imSize(2)), ...
               linspace(-0.5,0.5-1/imSize(3),imSize(3)) ...
    );

B0= exp(2i*(kmax_PRS(1)*x.*pos_PRS(1) ...
           +kmax_PRS(2)*y.*pos_PRS(2) ...
           +kmax_PRS(3)*z*pos_PRS(3)));

sig_shifted=sig.*permute(B0,[4 1 2 3]);

end



function [Rotmat_normal,rotmat_inplane,FOV_PRS]=getRotationMatrices(twix)
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

end

%%

function sliceInfo=getSlicePosition(twix_obj)

p=twix_obj.hdr.Phoenix.sSliceArray;
sliceInfo=cell(1,size(p.asSlice,2));
for i=1:size(p.asSlice,2)
    if(isfield(p.asSlice{i},'sPosition'))
        sliceInfo{i}.PosFieldName=['.dSag';'.dCor';'.dTra';];
        sliceInfo{i}.Position=mygetfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition);
    else
        sliceInfo{i}.Position=[0;0;0];
        fprintf('Position: Isocenter\n');
    end
    if(isfield(p.asSlice{i},'sNormal'))
        sliceInfo{i}.Normal=mygetfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{i}.sNormal);
    end
    if(isfield(twix_obj.hdr.Meas,'SliceThickness'))
        sliceInfo{i}.thickness=twix_obj.hdr.Meas.SliceThickness; %mm
    end
    if(isfield(p.asSlice{i},'dThickness') &&isfield(p.asSlice{i},'dReadoutFOV'))
        sliceInfo{i}.FOV_PRS=[p.asSlice{i}.dReadoutFOV p.asSlice{i}.dPhaseFOV  p.asSlice{i}.dThickness]; %mm
    end

end
end

function field=mygetfield(Struct,Fieldname)

if(nargin<2)
    field=zeros(1,3);
    if(isfield(Struct,'dCor'))
        field(2)=Struct.('dCor');
    end
    if(isfield(Struct,'dSag'))
        field(1)=Struct.('dSag');
    end
    if(isfield(Struct,'dTra'))
        field(3)=Struct.('dTra');
    end
else
    field=0;
    if(isfield(Struct,Fieldname))
        field=Struct.(Fieldname);
    end
end
end