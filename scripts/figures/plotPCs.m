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

set(gcf,'Color','w','InvertHardcopy','off','Position',[157 52 1200 1200])