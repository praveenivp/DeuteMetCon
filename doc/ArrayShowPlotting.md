# plotting with Arrayshow
Arrayshow can link to user plotting fucntions that are useful for investigation and debuging. Here, three use case are presented. you need two things
1. appropriate `userCursorPosFcn()` in the path
2. set up `asObj.UserData` accordingly.

Use only one arrayShow instance when doing this!

## plot the AMARES fit of the voxel

```
function userCursorPosFcn(asObj, pos, plotDim)
% Entry function for user-defined cursor position callbacks.
% Use it to call your own cursor position dependent functions.
%
% This function is called from the asObj when pressing 'c', 'C'
% or when selecting the respective option via context menu.
try
    cursorPos=sprintf(regexprep(asObj.selection.getValue,':','%d'),asObj.cursor.getPosition);
    cursorPos=regexprep(cursorPos,',\d+$',',:'); %replace last dim with colon
    pxlIdx= str2double(strsplit(cursorPos,','));
    pxlIdx=pxlIdx(1:3);
    asObj.UserData.mcobj.demoFit(pxlIdx);

catch err
    error(err.message);
end
end
```

launch

## plot AMARES fit of the voxel corresponding to an antomical image
```
function userCursorPosFcn(asObj, pos, plotDim)
% Entry function for user-defined cursor position callbacks.
% Use it to call your own cursor position dependent functions.
%
% This function is called from the asObj when pressing 'c', 'C'
% or when selecting the respective option via context menu.
try
    % figure(420),clf,
    cursorPos=sprintf(regexprep(asObj.selection.getValue,':','%d'),asObj.cursor.getPosition);
    pxlIdx= str2double(strsplit(cursorPos,','));
    pxlIdx_anat=pxlIdx(1:3);
    anat_dist=asObj.UserData.anatpxl2dist*[pxlIdx_anat(:);1];
    anat_dist(asObj.UserData.ori_vec.flip_vec)=anat_dist(asObj.UserData.ori_vec.flip_vec)*-1; %you do one flip before 
    anat_dist=anat_dist(asObj.UserData.ori_vec.perm_vec);
    pxlIdx_csi=asObj.UserData.dist2csipxl*anat_dist;
    asObj.UserData.mcobj.demoFit(round(pxlIdx_csi(1:3)));   
catch err
    error(err.message);
end
end

```

The script I used to lauch it!
```
fn_metcon=mcobj_csi.WriteImages()
anat_vol=MyNiftiRead('../anat.nii','RAS'); %RAS should n't chnage
% may be becuse we flip before metcon nifti export
ori_vec.flip_vec=[1];

%% water voxel
ni_anat=niftiinfo('../anat.nii'); %RAS
ni=niftiinfo(fn_metcon);
pxl_anat=[56,148,164,1]';
% stratergy
% anat_pxl -> dist in RAS -> dist in PLI -> pxl_CSI
dist=ni_anat.Transform.T'*pxl_anat;

pxl_csi=round(pinv(ni.Transform.T')*dist);
%matrix to transform anat pxl to CSI pxl
T=(pinv(ni.Transform.T')*ni_anat.Transform.T');

%% both should be in same orintation (RAS)

as(anat_vol,'select',':,:,166')
pause(1)
asObjs.UserData.mcobj=mcobj_csi;
asObjs.UserData.anatpxl2dist=ni_anat.Transform.T';
asObjs.UserData.ori_vec=ori_vec;
asObjs.UserData.dist2csipxl=pinv(ni.Transform.T');
```


## plot AMARES fit of the average spectrum of an ROI in a different space


```
function userCursorPosFcn(asObj, pos, plotDim)
% Entry function for user-defined cursor position callbacks.
% Use it to call your own cursor position dependent functions.
%
% This function is called from the asObj when pressing 'c', 'C'
% or when selecting the respective option via context menu.
try
    % figure(420),clf,
    %get spectrum from cursor position(spectrum should be last dim)
    cursorPos=sprintf(regexprep(asObj.selection.getValue,':','%d'),asObj.cursor.getPosition);
    % cursorPos=regexprep(cursorPos,',\d+$',',:'); %replace last dim with colon

    pxlIdx= str2double(strsplit(cursorPos,','));

    ROI_mask=poly2mask(asObj.roi.objPoly.Position(:,1),asObj.roi.objPoly.Position(:,2),size(asObj.data.dat,1),size(asObj.data.dat,2));
    imSize=size(asObj.data.dat,1:2);

    [x, y] = ind2sub(size(ROI_mask), find(ROI_mask));

    % Alternatively, you can create a cell array where each voxel is a separate element
    voxels_anat = arrayfun(@(i) [x(i), y(i),pxlIdx(3) ], 1:length(x), 'UniformOutput', false);
    voxels_CSI=cell(size(voxels_anat));
    for vxl=1:length(voxels_anat)
        pxlIdx_anat=voxels_anat{vxl}(:);
        anat_dist=asObj.UserData.anatpxl2dist*[pxlIdx_anat(:);1];
        anat_dist(asObj.UserData.ori_vec.flip_vec)=anat_dist(asObj.UserData.ori_vec.flip_vec)*-1;
        anat_dist=anat_dist(asObj.UserData.ori_vec.perm_vec);
        pxlIdx_csi=asObj.UserData.dist2csipxl*anat_dist;
        voxels_CSI{vxl}=round(pxlIdx_csi(1:3));
    end


    asObj.UserData.mcobj.demoFit(voxels_CSI);

catch err
    error(err.message);
end
end

```