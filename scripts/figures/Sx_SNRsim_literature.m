%% Peters et al 2021
% 0.1974           0.2401
%  0.80797            1.7652
%   0.28553        0.79157

metabolites=getMetaboliteStruct('Peters12T');
nMet=3;
metabolites=metabolites(1:3);

B0=0;%randn([1000,1])*2*6.5360*9.4; % 2ppm
PC=180;
TR=95e-3;
FA=90;
TE=2.1e-3;
DC=0.84;

[Msig_all,dc_fac]=MetSignalModel   (metabolites,TE,0, ...
    TR,B0,deg2rad(FA),'GRE',DC);
SE_csi=abs(Msig_all).*sqrt(dc_fac)./sqrt(TR);

TR=12.2e-3;
FA=80;
DC=0.70;
[Msig_all,dc_fac]=MetSignalModel(metabolites,TE,deg2rad(PC), ...
    TR,B0,deg2rad(FA),'bSSFP',DC);
SE_me=abs(Msig_all).*sqrt(dc_fac)./sqrt(TR);

figure(21),
tt=tiledlayout(3,1);


metnames=reordercats(categorical({metabolites.name}),{metabolites.name});
clo=flip(lines(3),1);
ax=nexttile(tt);
bh=barh(metnames,cat(2,SE_me,SE_csi),'FaceAlpha',0.75,'BarWidth',0.8);
set(gca,'ColorOrder',clo([1 3],:))%[[ 0.8500,0.3250,0.0980];[ 0,0.4470,0.7410];[0.4940,0.1840,0.5560];])
lh=legend('ME-bSSFP','CSI','Location','southeast');
lh.Direction='reverse';
grid minor
set(gca,'FontSize',12)
title('signal efficiencies (Peters et al 2022)')
xlim(get(gca,'xlim')+[0 0.7])

for jj=1:2
    xtips = bh(jj).XEndPoints;
    ytips = bh(jj).YEndPoints;
    switch(3-jj)
        case 1
            labels = string(strsplit(sprintf('%0.1fx ',(round([1 1 1])))));
        case 2
            labels = string(strsplit(sprintf('%0.1fx ',round(SE_me./SE_csi,2))));
    end
    text(ytips+0.01,xtips,labels(1:3),'HorizontalAlignment','left',...
        'VerticalAlignment','middle','FontSize',12)
end 
set(gca,'FontSize',12)

%% Montrazi et al 2023
metabolites=getMetaboliteStruct('Peters12T');
nMet=3;
metabolites=metabolites(1:3);

B0=0;%randn([1000,1])*2*6.5360*9.4; % 2ppm
PC=180;
TR=95e-3;
FA=90;
TE=2.1e-3;
DC=0.84;

[Msig_all,dc_fac]=MetSignalModel   (metabolites,TE,0, ...
    TR,B0,deg2rad(FA),'GRE',DC);
SE_csi=abs(Msig_all).*sqrt(dc_fac)./sqrt(TR);

TR=11.48e-3;
FA=60;
DC=0.80;
[Msig_all,dc_fac]=MetSignalModel(metabolites,TE,deg2rad(PC), ...
    TR,B0,deg2rad(FA),'bSSFP',DC);
SE_bcsi=abs(Msig_all).*sqrt(dc_fac)./sqrt(TR);

TR=11.48e-3;
FA=60;
DC=0.70;
[Msig_all,dc_fac]=MetSignalModel(metabolites,TE,deg2rad(PC), ...
    TR,B0,deg2rad(FA),'bSSFP',DC);
SE_me=abs(Msig_all).*sqrt(dc_fac)./sqrt(TR);



metnames=reordercats(categorical({metabolites.name}),{metabolites.name});

ax=nexttile(tt);
bh=barh(metnames,cat(2,SE_me,SE_bcsi,SE_csi),'FaceAlpha',0.75,'BarWidth',0.8);
set(gca,'ColorOrder',flip(lines(3),1))%[[ 0.8500,0.3250,0.0980];[ 0,0.4470,0.7410];[0.4940,0.1840,0.5560];])
lh=legend('ME-bSSFP','CSI-bSSFP','CSI','Location','southeast');
lh.Direction='reverse';
grid minor
set(gca,'FontSize',12)
title('signal efficiencies (Montrazi et al 2023)')
xlim(get(gca,'xlim')+[0 0.7])

for jj=1:3
    xtips = bh(jj).XEndPoints;
    ytips = bh(jj).YEndPoints;
    switch(4-jj)
        case 1
            labels = string(strsplit(sprintf('%0.1fx ',(round([1 1 1])))));
        case 2
            labels = string(strsplit(sprintf('%0.1fx ',round(SE_bcsi./SE_csi,2))));
        case 3
            labels = string(strsplit(sprintf('%0.1fx ',round(SE_me./SE_csi,2))));
    end
    text(ytips+0.01,xtips,labels(1:3),'HorizontalAlignment','left',...
        'VerticalAlignment','middle','FontSize',12)
end 
set(gca,'FontSize',12)

%%
%% Peters et al 2021
% 0.1974           0.2401
%  0.80797            1.7652
%   0.28553        0.79157

metabolites=getMetaboliteStruct('Roig7T');
nMet=3;
metabolites=metabolites([2 3 4]);

B0=0;%randn([1000,1])*2*6.5360*9.4; % 2ppm
PC=180;
TR=290e-3;
FA=86;
TE=2.1e-3;
DC=0.869;

[Msig_all,dc_fac]=MetSignalModel   (metabolites,TE,0, ...
    TR,B0,deg2rad(FA),'GRE',DC);
SE_csi=abs(Msig_all).*sqrt(dc_fac)./sqrt(TR);

TR=23e-3;
FA=50;
DC=0.76;
[Msig_all,dc_fac]=MetSignalModel(metabolites,TE,deg2rad(PC), ...
    TR,B0,deg2rad(FA),'bSSFP',DC);
SE_me=abs(Msig_all).*sqrt(dc_fac)./sqrt(TR);



metnames=reordercats(categorical({metabolites.name}),{metabolites.name});
clo=flip(lines(3),1);
ax=nexttile(tt);
bh=barh(metnames,cat(2,SE_me,SE_csi),'FaceAlpha',0.75,'BarWidth',0.8);
set(gca,'ColorOrder',clo([1 3],:))%[[ 0.8500,0.3250,0.0980];[ 0,0.4470,0.7410];[0.4940,0.1840,0.5560];])
lh=legend('CRT-bSSFP','CSI','Location','southeast');
lh.Direction='reverse';
grid minor
set(gca,'FontSize',12)
title('signal efficiencies (Frese et al 2025)')
xlim(get(gca,'xlim')+[0 0.7])

for jj=1:2
    xtips = bh(jj).XEndPoints;
    ytips = bh(jj).YEndPoints;
    switch(3-jj)
        case 1
            labels = string(strsplit(sprintf('%0.1fx ',(round([1 1 1])))));
        case 2
            labels = string(strsplit(sprintf('%0.1fx ',round(SE_me./SE_csi,2))));
    end
    text(ytips+0.01,xtips,labels(1:3),'HorizontalAlignment','left',...
        'VerticalAlignment','middle','FontSize',12)
end 
set(gca,'FontSize',12)