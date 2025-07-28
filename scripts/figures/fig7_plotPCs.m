%% process and write images
%% input files
MeasPath='/ptmp/pvalsala/deuterium/HOSJ-D6P2';
% MeasPath='/ptmp/pvalsala/deuterium/H4DK-64GH';
% MeasPath='/ptmp/pvalsala/deuterium/I3BL-CJ5O';
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
% imPlot=ndflip(squeeze(permute(All_PC(:,29,:,:,:),[3 1 2 5 4])),[ 1 ]); %sag
imPlot=ndflip(squeeze(permute(All_PC(:,:,30,:,:),[1 2 3 5 4])),[  ]); %tra

% sum(All_PC,5)./(sqrt(size(All_PC,5)))
figure(4),clf
tt=tiledlayout(5,3,'TileSpacing','tight','Padding','compact');


%create mask
% mask=All_PC(:,:,30,1)>10;
%figure,mask=CreateMask(imPlot(:,:,1,1));

nexttile(tt,3),
fm=mcobj_csi_ssfp_modes.Experimental.fm_est;
imagesc(imdilate(mask,strel("disk",5,0)).*fm(:,:,30)),cb=colorbar;, axis image,ylabel('fieldmap [Hz]')
xticks([]),yticks([]),colormap(gca,"turbo")
set(gca,'FontWeight','bold'),clim([-30 30])

%calc stats
stats=reshape(imPlot,[],7,4);
stats=stats(mask(:),:,:);
stats_95p=squeeze(prctile(stats,95,1));
stats_std=squeeze(std(stats,[],1));
stats_mean=squeeze(mean(stats,1));
SNR_gain_p95=max(stats_95p(4:7,:),[],1)./stats_95p(1:3,:);


for i=1:4
    nexttile(tt,1+3*i,[1 3])
imagesc(createImMontage(imPlot(:,:,:,i),size(imPlot,3)))
colorbar,axis image

% if(i<4)
% text(1:51:51*8,ones(8,1)*48,strsplit(sprintf('%.1f ',stats_median(1:7,i)),' '),'color','w','FontSize',12)
% text(35:51:51*8,ones(8,1)*48,strsplit(sprintf('%.1f ',stats_95p(1:7,i)),' '),'color','w','FontSize',12)
% end




ylabel(mcobj_csi_ssfp_pinv.metabolites(i).name)
if(i==4)
xticks((0.5:size(imPlot,3)+0.5)*size(imPlot,2)),
xticklabels({'Linear','IDEAL-modes','IDEAL',['180',char(176)],['270',char(176)],['360',char(176)],['90',char(176)]})
else
    xticks([])
end

if(i<4)
xticks(25:51:51*7)
set(gca,'XAxisLocation',"top")
lbl={};
for ll=1:7
end
xticklabels(strsplit(sprintf('%.1f\\pm%.0f/%.1f ',[stats_mean(1:7,i),stats_std(1:7,i),stats_95p(1:7,i)]'),' '))
end

if(i==1),title('CSI-PC-bSSFP metabolite maps [SNR]'),end
yticks([])
colormap(gca,'turbo')
set(gca,'clim',get(gca,'clim').*[1 0.9])

for jj=3:3+4
text(size(imPlot,1)*jj,5,'x2','FontSize',12,'Color','W')
end
 set(gca,'FontWeight','bold')
end
fontsize(gcf,"scale",1.4)

set(gcf,'Color','w','InvertHardcopy','off','Position', [150 52 1500 1100])