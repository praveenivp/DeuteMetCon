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


metabolites=getMetaboliteStruct('invivo');

%% pinv mode
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[]};

       fn=fullfile(sn,dirst_csi(1).name);   
    mcobj_csi_ssfp_pinv=MetCon_CSI(fn,CSI_setting{:},'fm','IDEAL','Solver','pinv');
     mcobj_csi_ssfp_IDEAL=MetCon_CSI(fn,CSI_setting{:},'fm','IDEAL','Solver','IDEAL');
     mcobj_csi_ssfp_modes=MetCon_CSI(fn,CSI_setting{:},'fm','IDEAL','Solver','IDEAL-modes');


%%
 fn=fullfile(sn,dirst_csi(end).name);   
mcobj_us=cell(4,1);
for cPC=1:length(mcobj_us)      
    mcobj_us{cPC}=copy(mcobj_csi_ssfp_pinv);
mcobj_us{cPC}.flags.PCSel=cPC;
mcobj_us{cPC}.performMetCon();
end




%%


All_PC=cellfun(@(mcobj) mcobj.getNormalized(),mcobj_us,'UniformOutput',false);
All_PC=cat(5, All_PC{:});
All_PC=cat(5,mcobj_csi_ssfp_pinv.getNormalized(),mcobj_csi_ssfp_modes.getNormalized(),mcobj_csi_ssfp_IDEAL.getNormalized(),sqrt(4)*All_PC);

% sum(All_PC,5)./(sqrt(size(All_PC,5)))
figure(4),clf
tt=tiledlayout(4,1,'TileSpacing','none','Padding','compact');

% imPlot=ndflip(squeeze(permute(All_PC(:,20,:,:,:),[3 2 1 5 4])),[1 ]);
imPlot=ndflip(squeeze(permute(All_PC(:,:,30,:,:),[1 2 3 5 4])),[ ]);
for i=1:4
    nexttile(tt)
imagesc(createImMontage(imPlot(:,:,:,i),size(imPlot,3)))
colorbar,axis image
ylabel(mcobj_csi_ssfp_pinv.metabolites(i).name)
if(i==4)
xticks((0.5:size(imPlot,3)+0.5)*size(imPlot,2)),
xticklabels({'Linear','IDEAL-modes','IDEAL',['180',char(176)],['270',char(176)],['360',char(176)],['90',char(176)]})
else
    xticks([])
end
if(i==1),title('metabolite maps in SNR units'),end
yticks([])
colormap('turbo')
set(gca,'clim',get(gca,'clim').*[1 0.9])

for jj=3:3+4
text(size(imPlot,1)*jj,5,'x2','FontSize',12,'Color','W')
end

end
fontsize(gcf,"scale",1.5)

set(gcf,'Color','w','InvertHardcopy','off','Position', [150 212 1471 824])