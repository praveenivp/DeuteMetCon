
% Tumor Rat
%pathname = ['\\mrz10\MRZ14T\AGKS\rolf\Deuterium\7T\LaKu_Flex1_20_MRI2_001_DMIRatSurface_LaKu_F_1_1_20231124_155111\18'];

%Lennarts mouse data
% Day 3
% Rat 1
% pathname = ['\\mrz10\MRZ14T\AGKS\rolf\Deuterium\7T\Lennart\Messtag3\LeWi_Int1_23_MRI2_010_Glc_LeWi_Int1_23_MRI2_1_1_20240119_094142\20'];
% imagename = ['\\mrz10\MRZ14T\AGKS\rolf\Deuterium\7T\Lennart\Messtag3\LeWi_Int1_23_MRI2_010_Glc_LeWi_Int1_23_MRI2_1_1_20240119_094142\5'];

% Rat2
pathname = ['\\mrz10\MRZ14T\AGKS\rolf\Deuterium\7T\Lennart\Messtag3\LeWi_Int1_23_MRI2_012_Glc_LeWi_Int1_23_MRI2_1_1_20240119_130353\18'];
imagename = ['\\mrz10\MRZ14T\AGKS\rolf\Deuterium\7T\Lennart\Messtag3\LeWi_Int1_23_MRI2_012_Glc_LeWi_Int1_23_MRI2_1_1_20240119_130353\4'];

% Rat3
%pathname = ['\\mrz10\MRZ14T\AGKS\rolf\Deuterium\7T\Lennart\Messtag3\LeWi_Int1_23_MRI2_013_Glc_LeWi_Int1_23_MRI2_1_1_20240119_155544\19'];
% pathname = ['\\mrz10\MRZ14T\AGKS\rolf\Deuterium\7T\Lennart\Messtag3\LeWi_Int1_23_MRI2_013_Glc_LeWi_Int1_23_MRI2_1_1_20240119_155544\17'];
% imagename = ['\\mrz10\MRZ14T\AGKS\rolf\Deuterium\7T\Lennart\Messtag3\LeWi_Int1_23_MRI2_013_Glc_LeWi_Int1_23_MRI2_1_1_20240119_155544\5'];

im = DataClass;
im = im.BrReadPro(imagename);
csi = DataClass;
csi = csi.BrReadRaw(pathname);

%csi = csi.CoilCombine;
csi = csi.CSISpatialFT([128,128,32]);

rCSIView(csi,im);
%rim(csi.AcqParameters.kspace);
%csi = csi.CSICalcPSF;