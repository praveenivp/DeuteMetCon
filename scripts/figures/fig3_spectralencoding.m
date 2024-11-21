% To demonstrate spectral encoding of phase cycling dimension
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/MATLAB'))

%% all input files
sn='/ptmp/pvalsala/deuterium/20240813_spectral';
fn_noise=fullfile(sn,'meas_MID00849_FID14520_pvrh_trufi_Noise.dat');
fn_1PC=fullfile(sn,'meas_MID00845_FID14516_pvrh_trufi_5E_1PC_12P5mm_FA50_s4_r180.dat');
fn_18PC=fullfile(sn,'meas_MID00844_FID14515_pvrh_trufi_5E_18PC_12P5mm_FA50_s4_r180.dat');

% fn_1PC=fullfile(sn,'meas_MID00847_FID14518_pvrh_trufi_6E_1PC_12P5mm_FA50_s4_r180.dat');
% fn_18PC=fullfile(sn,'meas_MID00848_FID14519_pvrh_trufi_6E_18PC_12P5mm_FA50_s4_r180.dat');
%%
metabolites=getMetaboliteStruct('phantom',0);
common_settings={'metabolites',metabolites,'doPhaseCorr',false,'doZeroPad',[1 1 1],'mask',80};

%% get dependencies fieldmap and noise decorr
twix_noise=mapVBVD(fn_noise,'rmos');
[D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix_noise);

 mcobj_fm=MetCon_bSSFP(fn_18PC,'NoiseDecorr',D_image,'Solver','IDEAL',common_settings{:});
%field map in rad/s in proton scale
fm_ideal=imgaussfilt3(mcobj_fm.Experimental.fm_est*(-2*pi)/(6.536 /42.567),0.5);

%get field map from field maps
%  r=B0map('meas_MID00858_FID14529_piv_gre_B0mapping_5Echoes.dat');
%   mask_im=abs(squeeze(r.reco_obj.img(1,55,10:end-5,:,1)));
%   fm_1H=squeeze(r.Fmap(55,10:end-5,:)).*(mask_im>0.5);


common_settings_2=[common_settings,{'fm',fm_ideal,'NoiseDecorr',D_image','Solver','pinv'}];


%%
DMIPara=getDMIPara(mapVBVD(fn_1PC)); 
Nechoes=length(DMIPara.TE);
 mcobj_18PC=MetCon_bSSFP(fn_18PC,'EchoSel',1:Nechoes,common_settings_2{:});
 mcobj_1PC=MetCon_bSSFP(fn_1PC,'EchoSel',1:Nechoes,common_settings_2{:});


% process all retrpective echo undersampling
allmet=[];
cond_all=[];

for i=1:Nechoes
    mcobj_18PC.flags.EchoSel=1:i;
    mcobj_1PC.flags.EchoSel=1:i;

    mcobj_18PC.performMetCon();
    mcobj_1PC.performMetCon();
    allmet=cat(5,allmet, ...
        cat(6,mcobj_18PC.getNormalized(),mcobj_1PC.getNormalized()) );

        cond_all=cat(5,cond_all, ...
        cat(6,mcobj_18PC.Experimental.condition,mcobj_1PC.Experimental.condition) );
end

%% process condition number
cond_all1=reshape(cond_all,[],Nechoes,2);
mask80=mcobj_18PC.getMask(80);
 cond_all1=cond_all1(mask80(:),:,:);
cond_mean=squeeze(mean(cond_all1,1,"omitnan"))';
cond_std=squeeze(std(cond_all1,[],1,"omitnan"))';

cond_std(2,1:3)=nan;
cond_mean(2,1:3)=nan;
dataset_label={'18PC','1PC'};
cd(sn)
% save('latest_fig3 data.mat',"allmet","cond*","Nechoes","dataset_label")

%% cretate nice mask
  figure
mask=CreateMask(squeeze(im_plot(:,:,1,1,1)));
%%

disp_echo=2:Nechoes;
figure(4),clf
im_plot=permute(squeeze(mean(allmet(15:end-10,48-1:48+1,:,:,disp_echo,:),2)),[1 2 3 5 4]);

%adjust SNR for missing echoes
im_plot=im_plot./permute(sqrt(disp_echo),[1 3 4 5 2]);

im_plot(:,:,:,2,1:2)=nan;

im_plot=im_plot.*mask;
im_plot(im_plot==0)=nan;

tt=tiledlayout(5,3*4,"Padding","compact",'TileSpacing','tight');
imsize=size(im_plot,1:2);

phantom_im=imread('/ptmp/pvalsala/Packages/DeuteMetCon/doc/images/phantom.png');
nexttile(1,[2 4]),imshow(phantom_im),axis image
title('phantom')


fm=squeeze(fm_ideal(15:end-10,round(size(allmet,2)/2),:,:,:,:));
fm=fm.*mask;
fm(isnan(fm))=0;
fm=fm./((-2*pi)/(6.536 /42.567));
% fm=4.1521+fm_1H./((-2*pi)/(6.536 /42.567)); % 4.1521 Hz offset between field maps
nexttile(5,[2 4]),imagesc(fm),axis image,colorbar
title('fielmap [Hz]')
xticks([]),yticks([]);
clim([-15 20])

nexttile(9,[2 4]),
 errorbar(disp_echo,cond_mean(2,disp_echo)',cond_std(2,disp_echo)','linewidth',1.5,'CapSize',10),hold on
 errorbar(disp_echo,cond_mean(1,disp_echo)',cond_std(1,disp_echo)','linewidth',1.5,'CapSize',10)
%  boxchart(rmoutliers(cond_all1(:,disp_echo,1)),'JitterOutliers','off','MarkerStyle','.'),hold on
%   boxchart(rmoutliers(cond_all1(:,disp_echo,2)),'JitterOutliers','off','MarkerStyle','.')
%   xticklabels(disp_echo);
legend(dataset_label{2},dataset_label{1}),axis square,grid on
ylabel('condition number'),xlabel('number of echoes')


for i=1:4
    nexttile(tt,12*2+(i-1)*3+1,[3 3])
    im_montage=createImMontage(reshape(im_plot(:,:,i,:,:),imsize(1),imsize(2),[]),2);
    imagesc(im_montage,'AlphaData',~isnan(im_montage))
    set(gca,'color',[0,0,0])
    title(metabolites(i).name)
    yticks([imsize(1)/2:imsize(1):imsize(1)*length(disp_echo)]);
    yticklabels(disp_echo)
    xticks([imsize(2)/2:imsize(2):imsize(2)*size(im_plot,4)])
      xticklabels({'18PC','1PC'})
      if( i==1),ylabel('number of echoes'); else, yticklabels([]); end
      colormap turbo,axis image
end
fontsize(gcf,'scale',1.5)
set(gcf,'color','w','Position', [435 81 1307 947],'InvertHardcopy','off')
savefig('fig3_PhaseCycling')
print(gcf,'fig3_phaseCycling','-dpng','-r300')