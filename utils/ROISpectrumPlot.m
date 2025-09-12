function ROI_CSI=ROISpectrumPlot(fn_anat,mcobj,ROI_mask)
% ROI_CSI=ROISpectrumPlot(fn_anat,mcobj,ROI_mask)

%fn_anat nifti filename

%ROI_mask is boolean and same size as anatomy


% mask_wm=MyNiftiRead('rsegmented.nii','RAS')==2 | MyNiftiRead('rsegmented.nii','RAS')==41;

% mask_gm=mod(MyNiftiRead('rsegmented.nii','RAS'),39)==3 ;
% fn_ROI='/ptmp/pvalsala/deuterium/patients/HHu5-HMzF/anat/test.nii.gz';
% mcobj=mcobj_csi;

fn_metcon=mcobj{1}.WriteImages([],{'image'});

% get transformation matrices
ni_anat=niftiinfo(fn_anat); %RAS
ni_metcon=niftiinfo(fn_metcon);
anatpxl2dist=ni_anat.Transform.T';
dist2csipxl=pinv(ni_metcon.Transform.T');
%to match orientation in mcobj
ori_vec.flip_vec=[1]; %first dim is flipped
ori_vec.perm_vec=[1:4]; % does nothing


[x, y,z] = ind2sub(size(ROI_mask), find(ROI_mask));

ROI_CSI=uint8(mcobj{end}.getMask(80)); % try to get brain region

% Alternatively, you can create a cell array where each voxel is a separate element
voxels_anat = arrayfun(@(i) [x(i), y(i),z(i) ], 1:length(x), 'UniformOutput', false);
voxels_CSI=cell(size(voxels_anat));
for vxl=1:length(voxels_anat)
    pxlIdx_anat=voxels_anat{vxl}(:);
    anat_dist=anatpxl2dist*[pxlIdx_anat(:);1];
    anat_dist(ori_vec.flip_vec)=anat_dist(ori_vec.flip_vec)*-1;
    anat_dist=anat_dist(ori_vec.perm_vec); % does nothing
    pxlIdx_csi=dist2csipxl*anat_dist;
    voxels_CSI{vxl}=round(pxlIdx_csi(1:3));
    % voxels_CSI{vxl}=voxels_CSI{vxl}+[1,-1,0]; %shift from spatial maps script

    ROI_CSI(voxels_CSI{vxl}(1),voxels_CSI{vxl}(2),voxels_CSI{vxl}(3))=2; % we mark ROI voxels in the brain region
end

% this remove duplicates from voxels_CSI
[x1, y1,z1] = ind2sub(size(ROI_CSI), find(ROI_CSI==2));
voxels_CSI2 = arrayfun(@(i) [x1(i), y1(i),z1(i) ], 1:length(x1), 'UniformOutput', false);


%debug
displayOrthogonalSlices(ROI_CSI,  2+round(mean(cell2mat(voxels_CSI2(:)),1)));


for ds=1:length(mcobj)
    mcobj{ds}.demoFit(voxels_CSI2,50+ds);
end

f51=figure(51);
f52=figure(52);
plotTogther(f51,f52)
fprintf('centroid voxel is (%d,%d,%d)\n',round(mean(cell2mat(voxels_CSI),2)));

end

function plotTogther(f51,f52)

figure(22),clf
tt=tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
ax=nexttile();
copyobj(f51.Children(4).Children(2),gca)
copyobj(f52.Children(4).Children(2),gca)
ax.Children(2).LineStyle='--';
xlim([-12,0]),ax.XAxis.Direction='reverse';
set(ax.Children,'LineWidth',1.5)
legend('75 mins','124 mins')
title('Data')

ax=nexttile();
copyobj(f51.Children(4).Children(1),gca)
copyobj(f52.Children(4).Children(1),gca)
xlim([-12,0]),ax.XAxis.Direction='reverse';
ax.Children(2).LineStyle='--';
set(ax.Children,'LineWidth',1.5)
title('fit')



ax=nexttile();
copyobj(f51.Children(3).Children,gca)
copyobj(f52.Children(3).Children,gca)

ax.Children(5).LineStyle='--';
ax.Children(6).LineStyle='--';
ax.Children(7).LineStyle='--';
ax.Children(8).LineStyle='--';
set(ax.Children,'LineWidth',1.5)
xlim([-12,0]),ax.XAxis.Direction='reverse';
title('peaks')

end
function displayOrthogonalSlices(volume_data, voxel_coords)
% Function to display XY, YZ, and XZ slices of a 3D volume
% Inputs:
%   volume_data - 3D matrix representing the volume
%   voxel_coords - 1x3 vector specifying the voxel coordinates [x, y, z]

% Extract the slices for each plane
sagittal_slice = squeeze(volume_data(voxel_coords(1), :, :)); % XY plane
coronal_slice = squeeze(volume_data(:, voxel_coords(2), :));  % YZ plane
axial_slice = squeeze(volume_data(:, :, voxel_coords(3)));    % XZ plane

% Create a figure to display the slices
figure(28),clf

% Display XY Slice
subplot(1, 3, 1);
imshow(sagittal_slice, []);
title('XY Plane');

% Display YZ Slice
subplot(1, 3, 2);
imshow(coronal_slice, []);
title('YZ Plane');

% Display XZ Slice
subplot(1, 3, 3);
imshow(axial_slice, []);
title('XZ Plane');


% Adjust layout
sgtitle('Orthogonal Slices at Specified Voxel');
end