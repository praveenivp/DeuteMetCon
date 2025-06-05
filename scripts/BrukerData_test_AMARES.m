%% plot simulated and measured bSSFP phase cycle profiles in a single ROI
% pn='S:\Deuterium\20250311_RatData\LaKu_Flex1_20_MRI2_001_DMIRat39_t3_LaKu_Fle_1_1_20250306_150629';
pathFolder='S:\Deuterium\20250311_RatData\LaKu_Flex1_20_MRI2_001_DMIRat39_t4_LaKu_Fle_1_1_20250311_144001';
% pn='S:\Deuterium\20220413_111704_Deuterium_SSFPtest_1_10';
cd('S:\Deuterium\20250311_RatData\proc')

% overview
[picked_data]=BrPickDataset(pathFolder);
 picked_data = regexp(picked_data{1}, '(\d+)\\fid','tokens');
 picked_data=picked_data{1}{1};

 %% Rolf portion
CSIpath = fullfile(pathFolder,'15');
imagepath = fullfile(pathFolder,'4');

im = DataClass;
im = im.BrReadPro(imagepath);
csi = DataClass;
csi = csi.BrReadRaw(CSIpath);
csi = csi.CSISpatialFT(2*[size(csi.RawData,2),size(csi.RawData,3),size(csi.RawData,4)]);
rCSIView(csi,im);

%% 7T rat AMARES

% Input voxel index
voxel_idx={10,6,19}; % x and y has to be swapped
fid=csi.ProData(:,voxel_idx{:});
faxis=linspace(-0.5*csi.AcqParameters.Bandwidth,0.5*csi.AcqParameters.Bandwidth,length(fid));
ppmAxis=faxis./csi.AcqParameters.Frequency;
timeAxis=csi.AcqParameters.TE*1e-3+linspace(0,csi.AcqParameters.DwellTime*(length(fid)-1),length(fid));
nMet=3;
amares_struct=struct('chemShift',[1.3;3.7;4.8], ...in ppm
    'phase', zeros(nMet,1),...
    'amplitude', ones(nMet,1),...
    'linewidth', 10*ones(nMet,1), ...
    'imagingFrequency', csi.AcqParameters.Frequency,... in MHz
    'BW', csi.AcqParameters.Bandwidth,...
    'timeAxis', timeAxis(:), ... [s]
    'dwellTime', csi.AcqParameters.DwellTime,... [s]
    'ppmAxis',ppmAxis(:), ...
    'beginTime',-1*csi.AcqParameters.TE*1e-3, ... [s]
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
nVoxels=numel(csi.ProData)/size(csi.ProData,1);
met_con=zeros(nVoxels,nMet);
for vxl=1:nVoxels
fid=csi.ProData(:,vxl);
[fitResults, fitStatus, figureHandle, CRBResults] = AMARES.amaresFit(double(fid(:)), amares_struct, pk, false,'quiet',true);
met_con(vxl,:)=fitResults.amplitude;
end
met_con=reshape(met_con,size(csi.ProData,2),size(csi.ProData,3),size(csi.ProData,4),nMet);

%% try export NIFTI
MyNIFTIWrite_Bruker(im.ProData,im,'im.nii');
MyNIFTIWrite_Bruker(met_con,csi,'csi.nii');

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

