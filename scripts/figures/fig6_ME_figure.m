%% process and write images
%% input files
sn='/ptmp/pvalsala/deuterium/DDRM-T4EU/TWIX';
sn='/ptmp/pvalsala/deuterium/20240920_newbSSFPseq/New folder';
dirst_me=dir(fullfile(sn,"*rh_trufi*.dat"));
dirst_me=dirst_me([1 2]);
pn='/ptmp/pvalsala/deuterium/DDRM-T4EU/proc/bSSFP2';
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))

metabolites=getMetaboliteStruct('phantom');


%% process noise anf get fieldmap
fn_noise=dir(fullfile(sn,"*rh_trufi*noise*.dat"));
twix_noise=mapVBVD(fullfile(sn,fn_noise(1).name),'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

fn_fm=(fullfile(sn,dirst_me(1).name));
ME_setting={'NoiseDecorr',D_image,'mask',[],'metabolites',metabolites,...
            'doPhaseCorr',false,'doZeropad',[1 1 1]*0.5,'parfor',true};
mcobj_ideal=MetCon_bSSFP(fn_fm,'Solver','IDEAL',ME_setting{:},'mask',80);
fm_meas_Hz=mcobj_ideal.Experimental.fm_est*(-2*pi)/(6.536 /42.567);
%% process all ME- data
ME_setting=[ME_setting,{'fm',fm_meas_Hz,'Solver','pinv'}];

mcobj_me=cell(length(dirst_me),1);
for cf=1%:length(dirst_me)
       fn=fullfile(sn,dirst_me(cf).name);
    mcobj_me{cf}=MetCon_bSSFP(fn,ME_setting{:});
end


%%
cd(pn)
intake_time=cellfun(@(x) x.getMinutesAfterIntake('08:05'),mcobj_me,'UniformOutput',true);
cellfun(@(x) x.WriteImages,mcobj_me,'UniformOutput',false);
%% reslice images to anatomy
dirst_nii=dir('Metcon_*.nii');
reslice_all_sag=myspm_reslice('anat_sag.nii',dirst_nii, 'linear','rs');
reslice_all_tra=myspm_reslice('anat_tra.nii',dirst_nii, 'linear','rt');


%%


ylabel_str=strsplit(sprintf('%.0f min\n',intake_time(:)),'\n');
%%
anat_nii=double(MyNiftiRead("anat_sag.nii",'IPR'));
[met_snr]=MyNiftiRead('rsMetcon_SNR*.nii','IPR');
[met_mm]=MyNiftiRead('rsMetcon_mm*.nii','IPR');
met_mm(:,:,:,1,:)=met_snr(:,:,:,1,:);
% close all
figure(16),clf
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
[met_snr]=MyNiftiRead('rtMetcon_SNR*.nii','PRS');
[met_mm]=MyNiftiRead('rtMetcon_mm*.nii','PRS');
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
mask=mcobj_me{1}.getMask(95);
mask=imerode(mask,strel('sphere',2));
col_snr=cellfun(@(x)x.getNormalized,mcobj_me,'UniformOutput',false);
col_snr=cat(5,col_snr{:});
% as(col_snr.*mask)
col_snr=reshape(col_snr,[],size(col_snr,4),size(col_snr,5));
col_snr=abs(col_snr(mask(:),:,:));

col_mM=cellfun(@(x)x.getmM,mcobj_me,'UniformOutput',false);
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



