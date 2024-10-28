%% process and write images
%% input files
MeasPath='/ptmp/pvalsala/deuterium/HOSJ-D6P2';
sn=fullfile(MeasPath,'TWIX');

dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));
% dirst_csi=dirst_csi([1 3 4]);
pn=fullfile(MeasPath,sprintf('proc/csi_bSSFP_%s',datetime('today','Format','yyyyMMMdd')));
mkdir(pn)
dirst_csi=dirst_csi([ 4 ]);
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))


metabolites=getMetaboliteStruct('invivo3');

%% pinv mode
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[]};

       fn=fullfile(sn,dirst_csi(1).name);   
    mcobj_csi_ssfp=MetCon_CSI(fn,CSI_setting{:},'fm','IDEAL','Solver','pinv');


%%
 fn=fullfile(sn,dirst_csi(end).name);   
mcobj_us=cell(4,1);
for cPC=1:length(mcobj_us)      
    mcobj_us{cPC}=copy(mcobj_csi_ssfp);
mcobj_us{cPC}.flags.PCSel=cPC;
mcobj_us{cPC}.performMetCon();
end




%%


All_PC=cellfun(@(mcobj) mcobj.getNormalized(),mcobj_us,'UniformOutput',false);
All_PC=cat(5, All_PC{:});
All_PC=cat(5,mcobj_csi_ssfp.getNormalized(),2*All_PC);

% sum(All_PC,5)./(sqrt(size(All_PC,5)))
figure(4),clf
tt=tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

% imPlot=ndflip(squeeze(permute(All_PC(:,20,:,:,:),[3 2 1 5 4])),[1 ]);
imPlot=ndflip(squeeze(permute(All_PC(:,:,32,:,:),[1 2 3 5 4])),[ ]);
for i=1:4
    nexttile(tt)
imagesc(createImMontage(imPlot(:,:,:,i),size(imPlot,3)))
colorbar,axis image
title(mcobj_csi_ssfp.metabolites(i).name)
xticks((0.5:size(imPlot,3)+0.5)*size(imPlot,2)),xticklabels({'all PC','180 deg','270 deg','360 deg','90 deg'})
yticks([])
colormap('jet')
end
fontsize(gcf,"scale",1.5)

    %% write images
cd(pn)
intake_time=cellfun(@(x) x.getMinutesAfterIntake('08:36'),mcobj_csi_ssfp,'UniformOutput',true);
cellfun(@(x) x.WriteImages,mcobj_csi_ssfp,'UniformOutput',false);
%% reslice images to anatomy
dirst_nii=dir('Metcon_*.nii');
reslice_all_sag=myspm_reslice('../*anat_*sag*.nii',dirst_nii, 'nearest','rs');
reslice_all_tra=myspm_reslice('../*anat_*tra*.nii',dirst_nii, 'nearest','rt');


%%

ylabel_str=strsplit(sprintf('%.0f min\n',intake_time(1:2)),'\n');
%%
anat_nii=double(MyNiftiRead("../*anat*_sag*.nii",'IPR'));
[met_snr]=MyNiftiRead('rsMetcon_SNR*.nii','IPR');
[met_mm]=MyNiftiRead('rsMetcon_mM*.nii','IPR');
met_mm(:,:,:,1,:)=met_snr(:,:,:,1,:);
% close all
figure(13),clf
tt=tiledlayout(2,6,'TileSpacing','tight','Padding','compact');
transform=@(x) ndflip(permute(x(1:end,1:end,:,:,:),[1:5]),[]);
% title_str={'Glx','Glc','D20'};
title_str={'Water [SNR]','Glc [mM]','Glx [mM]'};
cax_met={[0 80],[0 3],[0 3]}

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


anat_nii=double(MyNiftiRead("../*anat*_tra*.nii",'PRS'));
[met_snr]=MyNiftiRead('rtMetcon_SNR*.nii','PRS');
[met_mm]=MyNiftiRead('rtMetcon_mM*.nii','PRS');
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
mask=mcobj_csi_ssfp{end}.getMask(90);
mask=imerode(mask,strel('sphere',2));
col_snr=cellfun(@(x)x.getNormalized,mcobj_csi_ssfp,'UniformOutput',false);
col_snr=cat(5,col_snr{:});
% as(col_snr.*mask)
col_snr=reshape(col_snr,[],size(col_snr,4),size(col_snr,5));
col_snr=abs(col_snr(mask(:),:,:));

col_mM=cellfun(@(x)x.getmM,mcobj_csi_ssfp,'UniformOutput',false);
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
mc_mm=abs(mcobj_csi_ssfp{3}.getmM);
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



