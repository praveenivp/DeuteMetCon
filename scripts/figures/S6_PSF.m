%%
%% input files
MeasPath='/ptmp/pvalsala/deuterium/HOSJ-D6P2';
sn=fullfile(MeasPath,'TWIX');

dirst_csi=dir(fullfile(sn,"*rpcsi_fid*.dat"));
dirst_csi=dirst_csi(2);
dirst_csi_ssfp=dir(fullfile(sn,"*rpcsi_ssfp*.dat"));
dirst_me=dir(fullfile(sn,"*pvrh_trufi_5E_*.dat"));
dirst_me=dirst_me(1);

pn=fullfile(MeasPath,sprintf('proc/csi_GRE_%s',datetime('today','Format','yyyyMMMdd')));


addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))


%% read al data

fn=fullfile(sn,dirst_csi(1).name);
twix_csi=mapVBVD(fn);

fn=fullfile(sn,dirst_csi_ssfp(1).name);
twix_csi_ssfp=mapVBVD(fn);

fn=fullfile(sn,dirst_me(1).name);
twix_me=mapVBVD(fn);
%%
 [fw64_csifisp,PSF_csifisp,W_csifisp]=getPSF_CSI(twix_csi,true);
 [fw64_csibssfp,PSF_csibssfp,W_csibssfp]=getPSF_CSI(twix_csi_ssfp,true);
 [fw64_mebssfp,PSF_mebssfp,W_mebssfp]=getPSF_CSI(twix_me,true);
%%
figure(5),clf
cm=lines(3);
tt=tiledlayout(2,3,"TileSpacing","compact",'Padding','compact');
w=hamming(1024);
nexttile()
stem(round(PSF_csifisp.kaxis_E1*12),PSF_csifisp.W_E1,'LineWidth',1.5)
xticks(-12:3:12),xlabel('PE_{1-3}'),ylabel('averages')
hold on, plot(linspace(-12,12,1024),w*12,'LineWidth',1.5)
title('CSI weighting'),axis square,grid on
ylim([0 12])
nexttile()
stem(round(PSF_csifisp.kaxis_E1*12),PSF_csibssfp.W_E1./4,'LineWidth',1.5)
hold on, plot(linspace(-12,12,1024),w*6,'LineWidth',1.5)
xticks(-12:4:12),xlabel('PE_{1-3}'),ylabel('averages')
title('CSI-PC-bSSFP weighting')
axis square,grid on

nexttile()
imagesc(-8:7,-12:11,squeeze(W_mebssfp(1,:,:)))
axis image
xlabel('PE_1'),ylabel('PE_2'),xticks(-8:4:7),yticks(-12:4:12)
box on
title('ME-PC-bSSFP: elliptical scanning')

nexttile()
plot(PSF_csifisp.FOV_E1,PSF_csifisp.PSF_E1,'LineWidth',1.5)
xlim([-1 1]*50),axis square
grid minor,grid on
title('CSI: PSF')
xlabel('distance [mm]'),ylabel('amplitude [a.u]')
box on
 text(0,0.64,' 15.9\newline mm','HorizontalAlignment','center','FontSize',8,'Color',cm(1,:))
%15.93
nexttile()
plot(PSF_csibssfp.FOV_E1,PSF_csibssfp.PSF_E1,'LineWidth',1.5)
xlim([-1 1]*50),axis square
grid minor,grid on
title('CSI-PC-bSSFP: PSF')
xlabel('distance [mm]'),ylabel('amplitude [a.u]')
box on
text(0,0.64,' 16.2\newline mm','HorizontalAlignment','center','FontSize',8,'Color',cm(1,:))

nexttile()
hold on
plot(PSF_mebssfp.FOV_E1,PSF_mebssfp.PSF_E1,'LineWidth',1.5)
plot(PSF_mebssfp.FOV_E2,PSF_mebssfp.PSF_E2,'LineWidth',1.5)
plot(PSF_mebssfp.FOV_E3,PSF_mebssfp.PSF_E3,'LineWidth',1.5)
xlim([-1 1]*50),axis square,ylim([-0.25 1])
grid minor,grid on
legend('read','PE1','PE2')
title('ME-PC-bSSFP: PSF')

text(-25,0.5,' 12.5\newline mm','HorizontalAlignment','center','FontSize',8,'Color',cm(1,:))
text(0,0.5,' 14.4\newline mm','HorizontalAlignment','center','FontSize',8,'Color',cm(2,:))
text(25,0.5,' 14.8\newline mm','HorizontalAlignment','center','FontSize',8,'Color',cm(3,:))

xlabel('distance [mm]'),ylabel('amplitude [a.u]')
set(gcf,'Color','w','Position',[459 273 1079 671])
fontsize(gcf,'scale',1.3)
box on

cd('/ptmp/pvalsala/deuterium/paper/')
print(gcf,'S8_PSF.png','-dpng','-r300')

%% calculate signal energy lost in negative lobes of PSF
 load('/ptmp/pvalsala/deuterium/paper/weights.mat')
% W2=ndCircShift(W_csibssfp,floor(size(W_csibssfp)/2)+1,[1 2 3]);
W2=padarray(W_fispcsi,[1,1,1]*512,0,'both');
PSF_csifisp=myfft(W2,[1 2 3],[1 2 3]);
W2=padarray(W_csibssfp,[1,1,1]*512,0,'both');
PSF_csibssfp=real(myfft(W2,[1 2 3],[1 2 3]));
W2=padarray(W_mebssfp,[1,1,1]*512,0,'both');
PSF_me=real(myfft(W2,[1 2 3],[1 2 3]));

clear W2;

Int_csifisp=sum(real(PSF_csifisp)./max(real(PSF_csifisp(:))),'all');
Int_csibssfp=sum(real(PSF_csibssfp)./max(real(PSF_csibssfp(:))),'all');
Int_me=sum(real(PSF_me)./max(real(PSF_me(:))),'all');

FOV_me = [400, 300,200]; % Field of view in mm [x,y,z]
FOV_csi = [208, 208, 208]; % Field of view in mm [x,y,z]
FOV_fac=prod(FOV_csi)/prod(FOV_me);

res_me=[12.44 14.42 14.84 ]; %mm
new_vol=[4.0425    4.2831    2.6670]; %mL
res_fac=new_vol(1:2)/new_vol(3);
% 
%2.66 mL vs 4.3 mL(16.2579^3 ) vs  4.08 mLvs 
SNR_fac_csifisp=(Int_csifisp/Int_me)*(FOV_fac*res_fac(1))
SNR_fac_csibssfp=(Int_csibssfp/Int_me)*(FOV_fac*res_fac(2))




%after : 4.0425./[4.0425    4.2831    2.6670]  1.0000    0.9438    1.5157
%before :2.85./[2.85 2.98 2.18]  1.000    0.9564    1.3073

