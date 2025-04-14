MeasPath='/ptmp/pvalsala/deuterium/HOSJ-D6P2';
sn=fullfile(MeasPath,'TWIX');

dirst_csi=dir(fullfile(sn,"*rpcsi*.dat"));
% dirst_csi=dirst_csi([1 3 4]);
pn=fullfile(MeasPath,sprintf('proc/PlotPC_ME_bSSFP_%s',datetime('today','Format','yyyyMMMdd')));
mkdir(pn);
dirst_me=dir(fullfile(sn,"*rh_trufi*.dat"));
dirst_me=dirst_me([1 3 4]);
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))

metabolites=getMetaboliteStruct('invivo');


%% 
% dirst_fm=dir(fullfile(sn,'*piv_gre_B0mapping*'));
%  fmobj=B0map(fullfile(sn,dirst_fm(1).name),'UnwrapMode','SpatialUnwrap','doRegularization',false);
% fmobj.exportNIFTI()
%% process noise anf get fieldmap
fn_noise=dir(fullfile(sn,"*rh_trufi*noise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

fn_fm=(fullfile(sn,dirst_me(3).name));
ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',true,'doZeropad',[1 1 1]*0.5,'parfor',true};
% mcobj_ideal=MetCon_ME(fn_fm,'Solver','IDEAL',ME_setting{:},'mask',80);
% fm_meas_Hz=mcobj_ideal.Experimental.fm_est*(-2*pi)/(6.536 /42.567);
%% process all ME- data

% ME_setting=[ME_setting,{'fm',fm_meas_Hz,'Solver','pinv'}];
mcobj_me=cell(length(dirst_me),1);
for cf=1:length(dirst_me)

       fn=fullfile(sn,dirst_me(cf).name);   
    mcobj_me{cf}=MetCon_ME(fn,ME_setting{:},'fm','IDEAL','Solver','pinv');
%     mcobj_me{cf-1}=MetCon_ME(fn,ME_setting{:},'fm','M01169_piv_gre_B0mapping_5Echoes_fmap.nii','Solver','pinv');
end

%%
 mcobj_me_pinv=MetCon_ME(fn,ME_setting{:},'fm','IDEAL','Solver','pinv');
  mcobj_me_modes=MetCon_ME(fn,ME_setting{:},'fm','IDEAL','Solver','IDEAL-modes');
    mcobj_me_IDEAL=MetCon_ME(fn,ME_setting{:},'fm','IDEAL','Solver','IDEAL');
%%
 fn=fullfile(sn,dirst_me(end).name);   
mcobj_us=cell(length(mcobj_me{3}.DMIPara.PhaseCycles),1);
im_all=mean(cell2mat(cellfun(@(x)x.img,mcobj_me,'UniformOutput',false)),1);
 mcobj_temp=copy(mcobj_me{end});
 mcobj_temp.flags.Solver='pinv';
% mcobj_temp.img=im_all;
 mcobj_temp.FieldMap='IDEAL';
% mcobj_temp.getFieldmap();
% mcobj_temp.performMetCon();
for cPC=1:length(mcobj_us)      
mcobj_us{cPC}=copy(mcobj_temp);
mcobj_us{cPC}.flags.PCSel=cPC;
mcobj_us{cPC}.getFieldmap();
mcobj_us{cPC}.performMetCon();
mcobj_us{cPC}.performPxlShiftCorrection();
end

%%


All_PC=cellfun(@(mcobj) mcobj.getNormalized(),mcobj_us,'UniformOutput',false);
All_PC=cat(5, All_PC{:});
All_PC=cat(5,mcobj_me_pinv.getNormalized(),mcobj_me_modes.getNormalized(),mcobj_me_IDEAL.getNormalized(),sqrt(18)*All_PC);
%%
% sum(All_PC,5)./(sqrt(size(All_PC,5)))
figure(5),clf
tt=tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

% imPlot=ndflip(squeeze(permute(All_PC(:,20,:,:,:),[3 2 1 5 4])),[1 ]);
% imPlot=ndflip(squeeze(permute(All_PC(10:end-7,20:end-20,15,:,:),[1 2 3 5 4])),[ ]); % sag
imPlot=ndflip(squeeze(permute(All_PC(10:end-7,36,:,:,:),[1 2 3 5 4])),[ ]); % sag
for i=1:4
    nexttile(tt)
imagesc(createImMontage(imPlot(:,:,:,i),size(imPlot,3)))
colorbar,axis image
title(mcobj_me_modes.metabolites(i).name)
if(i==4)
xticks((0.5:size(imPlot,3)+0.5)*size(imPlot,2)),
xticklabels([{'Linear',' IDEAL-\newlinemodes','IDEAL'},strsplit(num2str(mcobj_me_modes.DMIPara.PC_deg),' ')])
else
    xticks([])
end
yticks([])
colormap('jet')
end
fontsize(gcf,"scale",1.5)

set(gcf,'Color','w','InvertHardcopy','off','Position',[0 417 1907 588])


%%
% figure
% mask_me=CreateMask(imPlot(:,:,1,1,1));

stats_me=reshape(imPlot,[],18+3,4);
stats_me(:,4:end,:)=movmean(stats_me(:,4:end,:),1,2)*sqrt(1);
stats_me=stats_me(mask_me(:),:,:);
stats_95p_me=squeeze(prctile(stats_me,95,1));
stats_std_me=squeeze(std(stats_me,[],1));
stats_mean_me=squeeze(mean(stats_me,1));
SNR_gain_p95=max(stats_95p_me(4:end,:),[],1)./stats_95p_me(1:3,:)
SNR_gain_mean=max(stats_mean_me(4:end,:),[],1)./stats_mean_me(1:3,:);