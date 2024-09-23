%% process and write images
%% input files
MeasPath='/ptmp/pvalsala/deuterium/HOSJ-D6P2';
sn=fullfile(MeasPath,'TWIX');

dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));
% dirst_csi=dirst_csi([1 3 4]);
pn=fullfile(MeasPath,sprintf('csi_GRE_%s',datetime('today','Format','yyyyMMMdd')));
mkdir(pn)
dirst_csi=dirst_csi([1 3 2 4 ]);
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))


metabolites=getMetaboliteStruct('invivo1');


% fn_fm='/ptmp/pvalsala/deuterium/DA77-F3UY/dep/fm_MeasUID692.nii';
% ref_space=mcobj_csi{1}.WriteImages();
% resliced_fm=myspm_reslice(ref_space, fn_fm,'linear','r');

%% final setttings
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[],'Solver','IDEAL','fm',[]};


mcobj_csi=cell(length(dirst_csi),1);
for cf=2%1:length(dirst_csi)
       fn=fullfile(sn,dirst_csi(cf).name);
    mcobj_csi{cf}=MetCon_CSI(fn,CSI_setting{:});
end

%%
cd(pn)
intake_time=cellfun(@(x) x.getMinutesAfterIntake('08:05'),mcobj_csi,'UniformOutput',true);
cellfun(@(x) x.WriteImages,mcobj_csi,'UniformOutput',false);
%% reslice images to anatomy
dirst_nii=dir('Metcon_*.nii');
reslice_all_sag=myspm_reslice('anat_sag.nii',dirst_nii, 'linear','rs');
reslice_all_tra=myspm_reslice('anat_tra.nii',dirst_nii, 'linear','rt');


%%


ylabel_str=strsplit(sprintf('%.0f min\n',intake_time(:)),'\n');
%%
anat_nii=double(MyNiftiRead("anat_sag.nii",'IPR'));
[met_snr]=MyNiftiRead('rsMetcon_SNR*IDEAL*.nii','IPR');
[met_mm]=MyNiftiRead('rsMetcon_mm*IDEAL*.nii','IPR');
met_mm(:,:,:,1,:)=met_snr(:,:,:,1,:);
% close all
figure(11),clf
tt=tiledlayout(2,6,'TileSpacing','tight','Padding','compact');
transform=@(x) ndflip(permute(x(1:end,1:end,:,:,:),[1:5]),[]);
% title_str={'Glx','Glc','D20'};
title_str={'Water [SNR]','Glc [mM]','Glx [mM]'};
cax_met={[0 40],[0 3],[0 3]}

cb_all={};cax_all={};cnt=1;
slct=16;
slcs=18;
for i=1:3
    nexttile(1+(i-1)*2);
[cb_all{cnt},cax_all{cnt}]=overlayplot(anat_nii,met_mm,'MetIdx',i,'prctile',99.8,'SlcSel',slcs,'transform',transform,...
    'cax',cax_met{i});
if(i==1),yticklabels(ylabel_str(1:end-1)); else, yticks([]); end
title(title_str{i})
cnt=cnt+1;
end
%


anat_nii=double(MyNiftiRead("anat_tra.nii",'PRS'));
[met_snr]=MyNiftiRead('rtMetcon_SNR*IDEAL*.nii','PRS');
[met_mm]=MyNiftiRead('rtMetcon_mm*IDEAL*.nii','PRS');
 met_mm(:,:,:,1,:)=met_snr(:,:,:,1,:);
% transform=@(x) ndflip(permute(x(30:end-20,1:end,:,:,:),[1 2 3 4 5]),[]);
% title_str={'Glx','Glc','D20'};
% title_str={'Glx/Water_{14 min}','Glc/Water_{14 min}','Water/Water_{-9 min}'};
for i=1:3
    nexttile(2+(i-1)*2);
[cb_all{cnt},cax_all{cnt}]=overlayplot(anat_nii,met_mm,'MetIdx',i,'prctile',100,'SlcSel',slct,'transform',transform,...
  'cax',cax_all{i});
yticks([]);
title(title_str{i})
cb_all{cnt}=colorbar;
cnt=cnt+1;

end



for i=4:length(cb_all)
makeColorbar(cb_all{i},cax_all{i});
end

%% stats
% figure(13),clf
mask=mcobj_csi{end}.getMask(80);
mask=imerode(mask,strel('sphere',2));
col_snr=cellfun(@(x)x.getNormalized,mcobj_csi,'UniformOutput',false);
col_snr=cat(5,col_snr{:});
% as(col_snr.*mask)
col_snr=reshape(col_snr,[],size(col_snr,4),size(col_snr,5));
col_snr=abs(col_snr(mask(:),:,:));

col_mM=cellfun(@(x)x.getmM,mcobj_csi,'UniformOutput',false);
col_mM=cat(5,col_mM{:});
% as(col_mM.*mask)
col_mM=reshape(col_mM,[],size(col_mM,4),size(col_mM,5));
col_mM=abs(col_mM(mask(:),:,:));


for cMet=1:3
nexttile(tt,6+cMet),cla
violin(squeeze(col_snr(:,cMet,:)));
grid on
title([metabolites(cMet).name ,' [SNR]' ])
end
for cMet=1:3
nexttile(tt,9+cMet),cla
% col_m<col_mM(:,cMet,:)
violin(squeeze(col_mM(:,cMet,:)))
xticklabels(intake_time),

ylim([0,4]),grid on
title([metabolites(cMet).name ,' [mM]' ])
end



%% get mM
mc_mm=abs(mcobj_csi{3}.getmM);
mc_mm(mc_mm>5)=0;
mc_mm(:,:,:,3)=mc_mm(:,:,:,3)/(1-mean([41.5e-2,37.9e-2])); % 
as(mc_mm)
%%

function makeColorbar(cb_handle,cax)



cb_handle.Visible='off';
ax2 = axes('Position',cb_handle.Position);

imagesc(linspace(0,1,100)'),colormap(ax2,'jet')
set(ax2,'YAxisLocation','right','FontSize',10,'FontWeight','bold','YDir','normal')
yticks(linspace(0,100,4))
yticklabels(round(linspace(0,1,4)*cax(2),2))
xticks([])
end



