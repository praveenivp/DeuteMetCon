%% Demo Script 
pathFolder='S:\Deuterium\20250311_RatData\LaKu_Flex1_20_MRI2_001_DMIRat39_t4_LaKu_Fle_1_1_20250311_144001';

% check the measurement ID
[picked_data]=BrPickDataset(pathFolder);
%  picked_data = regexp(picked_data{1}, '(\d+)\\fid','tokens');
%  picked_data=picked_data{1}{1};

addpath('DeuteMetCon')
addpath('OXSA')

 %% read and parse bruker data
CSIpath = fullfile(pathFolder,'15');
imagepath = fullfile(pathFolder,'4');

imObj = DataClass;
imObj = imObj.BrReadPro(imagepath);
csiObj = DataClass;
csiObj = csiObj.BrReadRaw(CSIpath);
csiObj = csiObj.CSISpatialFT(2*[size(csiObj.RawData,2),size(csiObj.RawData,3),size(csiObj.RawData,4)]);
rCSIView(csiObj,imObj);

%% 7T rat AMARES

% Input voxel index
voxel_idx={10,6,19}; % x and y has to be swapped
fid=csiObj.ProData(:,voxel_idx{:});
faxis=linspace(-0.5*csiObj.AcqParameters.Bandwidth,0.5*csiObj.AcqParameters.Bandwidth,length(fid));
ppmAxis=faxis./csiObj.AcqParameters.Frequency;
timeAxis=csiObj.AcqParameters.TE*1e-3+linspace(0,csiObj.AcqParameters.DwellTime*(length(fid)-1),length(fid));
nMet=3;
amares_struct=struct('chemShift',[1.3;3.7;4.8], ...in ppm
    'phase', zeros(nMet,1),...
    'amplitude', ones(nMet,1),...
    'linewidth', 10*ones(nMet,1), ...
    'imagingFrequency', csiObj.AcqParameters.Frequency,... in MHz
    'BW', csiObj.AcqParameters.Bandwidth,...
    'timeAxis', timeAxis(:), ... [s]
    'dwellTime', csiObj.AcqParameters.DwellTime,... [s]
    'ppmAxis',ppmAxis(:), ...
    'beginTime',-1*csiObj.AcqParameters.TE*1e-3, ... [s]
    'offset',0,...
    'samples',length(fid), ... [s]
    'peakName',string({'lac','glc','water'}));

pk=PriorKnowledge_DMI(amares_struct);
%relax the bounds for better fit 
pk.bounds(3).chemShift=4.8+[-1 1]*0.5;
pk.bounds(2).chemShift=3.7+[-1 1]*0.5;
pk.bounds(1).chemShift=1.3+[-1 1]*0.5;

% test the fit
[fitResults, fitStatus, figureHandle, CRBResults] = AMARES.amaresFit(double(fid(:)), amares_struct, pk, true,'quiet',false);

%% fit all voxels
nVoxels=numel(csiObj.ProData)/size(csiObj.ProData,1);
met_con=zeros(nVoxels,nMet);
for vxl=1:nVoxels
    fid=csiObj.ProData(:,vxl);
    [fitResults, fitStatus, figureHandle, CRBResults] = AMARES.amaresFit(double(fid(:)), amares_struct, pk, false,'quiet',true);
    met_con(vxl,:)=fitResults.amplitude;
end
met_con=reshape(met_con,size(csiObj.ProData,2),size(csiObj.ProData,3),size(csiObj.ProData,4),nMet);

%% try export NIFTI
cd(CSIpath)
MyNIFTIWrite_Bruker(imObj.ProData,imObj,'im.nii');
MyNIFTIWrite_Bruker(met_con,csiObj,'csi.nii');

%% OLD IDEAL code
% timeSel=7:size(csi_fid,4);
% offset=csi_para.PVM_FrqWork(1);
% met_struct=[struct("T1_s",0.4498212,"T2_s",0.270852,"freq_shift_Hz",-82,"name","acetate","T2star_s",0);...
%     struct("T1_s",0.4498212,"T2_s",0.270852,"freq_shift_Hz",71,"name","glucose","T2star_s",0);...
%     struct("T1_s",0.4498212,"T2_s",0.270852,"freq_shift_Hz",183,"name","water","T2star_s",0);];
% TE_s=0:dt:dt*(size(csi_fid,4)-1)+csi_para.PVM_EchoTime*1e-3;
% IDEAL_obj=IDEAL(met_struct,TE_s(timeSel),'solver','IDEAL','fm',[],'maxit',5,'SmoothFM',2);
% dat=myfft(sum(csi_fid,5),1:3);
% 
% met=IDEAL_obj'*dat(:,:,:,timeSel);
% as(met)

