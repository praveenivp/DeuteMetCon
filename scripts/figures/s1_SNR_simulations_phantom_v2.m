
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))

refVoltage=550; % upto 550-610 V
kFactor=0.83;
RFfac=447/refVoltage; % Flip angle scale factor

metabolites=getMetaboliteStruct('phantom',0);
nMet=4;
metabolites=metabolites(1:nMet);

 pc_range_csi=linspace(0,360-360/4,4)+180;
 pc_range_me=linspace(0,360-360/18,18)+180;
TR_all_bssfp=linspace(15e-3,25e-3,50);
TR_all_gre=linspace(15e-3,100e-3,80);
FA_all=linspace(10,90,60);
TE=2.3e-3;
apply_T2star=true;
calc_signaleff=true;


B0=randn([1e3,1])/3*(2*6.5360*9.4); % 2ppm

%% SSFP signal signal

sig_all_bssfp=zeros(length(TR_all_bssfp),length(FA_all),length(metabolites));
sig_all_me=zeros(length(TR_all_bssfp),length(FA_all),length(metabolites));
sig_all_bssfp_noPC=zeros(length(TR_all_bssfp),length(FA_all),length(metabolites));
sig_all_me_noPC=zeros(length(TR_all_bssfp),length(FA_all),length(metabolites));
dc_fac_bssfp=zeros(length(TR_all_bssfp),1,length(metabolites));
dc_fac_me=zeros(length(TR_all_bssfp),1,length(metabolites));
maxFA_bssfp=zeros([length(TR_all_bssfp) 3]);
DC_ssfp=@(TR_s) ((TR_s-4.2e-3)/TR_s).*double(TR_s>4.1e-3); % 4.2 ms non-encoing time csi-bSSFP
DC_me=@(TR_s) ((TR_s-5.5e-3)/TR_s).*double(TR_s>5.5e-3); % 5.5 ms non-encoing time ME-bSSFP
% figure,tiledlayout("flow")
for cTR=1:size(sig_all_bssfp,1)
    [~,maxFA_bssfp(cTR,1)]=SimpleSARModel(1,500e-6,TR_all_bssfp(cTR),refVoltage,kFactor);
    [~,maxFA_bssfp(cTR,2)]=SimpleSARModel(1,1400e-6,TR_all_bssfp(cTR),refVoltage,kFactor);
    [~,maxFA_bssfp(cTR,3)]=SimpleSARModel(1,2000e-6,TR_all_bssfp(cTR),refVoltage,kFactor);
    for CM=1:size(sig_all_bssfp,3)
        %off-resonance is set to -1*chemical shift otherwise increase phase-cyles
        [Msig_all,dc_fac_bssfp(cTR,1,CM)]=MetSignalModel(metabolites(CM),TE,deg2rad(pc_range_csi), ...
            TR_all_bssfp(cTR),B0,deg2rad(FA_all),'bSSFP',DC_ssfp);
        sig_all_bssfp(cTR,:,CM)=abs(mean(abs(Msig_all),[3,5]));
                [Msig_all,dc_fac_me(cTR,1,CM)]=MetSignalModel(metabolites(CM),TE,deg2rad(pc_range_me), ...
            TR_all_bssfp(cTR),B0,deg2rad(FA_all),'bSSFP',DC_me);
        sig_all_me(cTR,:,CM)=abs(mean(abs(Msig_all),[3,5]));

   %no phase cycling
        [Msig_all,dc_fac_csi(cTR,1,CM)]=MetSignalModel(metabolites(CM),TE,deg2rad(180), ...
            TR_all_bssfp(cTR),0,deg2rad(FA_all),'bSSFP',DC_ssfp);
        sig_all_bssfp_noPC(cTR,:,CM)=abs(mean(abs(Msig_all),[3,5]));
                [Msig_all,dc_fac_me(cTR,1,CM)]=MetSignalModel(metabolites(CM),TE,deg2rad(180), ...
            TR_all_bssfp(cTR),0,deg2rad(FA_all),'bSSFP',DC_me);
        sig_all_me_noPC(cTR,:,CM)=abs(mean(abs(Msig_all),[3,5]));
    end
end

%% GRE signal efficiency
sig_all_gre=zeros(length(TR_all_gre),length(FA_all),length(metabolites));
sig_all_FISP=zeros(length(TR_all_gre),length(FA_all),length(metabolites));

dc_fac_fisp=zeros(length(TR_all_gre),1,length(metabolites));
dc_fac_gre=zeros(length(TR_all_gre),1,length(metabolites));

DC=@(TR_s) ((TR_s-7.56e-3)/TR_s).*double(TR_s>7.56e-3); % 7.56 ms non-encoing time

maxFA_gre=zeros([length(TR_all_gre) 3]);
for cTR=1:size(sig_all_gre,1)
    [~,maxFA_gre(cTR,1)]=SimpleSARModel(1,500e-6,TR_all_gre(cTR),refVoltage,kFactor);
    [~,maxFA_gre(cTR,2)]=SimpleSARModel(1,1400e-6,TR_all_gre(cTR),refVoltage,kFactor);
    [~,maxFA_gre(cTR,3)]=SimpleSARModel(1,2000e-6,TR_all_gre(cTR),refVoltage,kFactor);
    for CM=1:size(sig_all_gre,3)
        
        [Msig_all,dc_fac_fisp(cTR,1,CM)]=MetSignalModel   (metabolites(CM),TE,0, ...
            TR_all_gre(cTR),0,deg2rad(FA_all),'FISP',DC);
        sig_all_FISP(cTR,:,CM)=mean(abs(Msig_all),3);

        [Msig_all,dc_fac_gre(cTR,1,CM)]=MetSignalModel   (metabolites(CM),TE,0, ...
            TR_all_gre(cTR),0,deg2rad(FA_all),'FLASH',DC);
        sig_all_gre(cTR,:,CM)=mean(abs(Msig_all),3);

    end
end

% apply T2* and duty cycle penalty
if(apply_T2star)

    sig_all_gre=sig_all_gre.*sqrt(dc_fac_gre);
    sig_all_FISP=sig_all_FISP.*sqrt(dc_fac_fisp);
    sig_all_bssfp=sig_all_bssfp.*sqrt(dc_fac_bssfp);
     sig_all_me=sig_all_me.*sqrt(dc_fac_me);
       sig_all_bssfp_noPC=sig_all_bssfp_noPC.*sqrt(dc_fac_csi)*fac;
     sig_all_me_noPC=sig_all_me_noPC.*sqrt(dc_fac_me)*fac;
end

%calcualte signal efficiency
if(calc_signaleff)
    sig_all_gre=sig_all_gre./sqrt(TR_all_gre(:));
    sig_all_FISP=sig_all_FISP./sqrt(TR_all_gre(:));
    sig_all_bssfp=sig_all_bssfp./sqrt(TR_all_bssfp(:));
    sig_all_me=sig_all_me./sqrt(TR_all_bssfp(:));
    sig_all_bssfp_noPC=sig_all_bssfp_noPC./sqrt(TR_all_bssfp(:));
    sig_all_me_noPC=sig_all_me_noPC./sqrt(TR_all_bssfp(:));
end



%some plot handles
plotProtTR74= @ () plot(74,57*RFfac,'r*','MarkerSize',10);
plotProtTR36= @ () plot(36,41*RFfac,'r*','MarkerSize',10);
plotProtbssfp= @ () plot(19,50*RFfac,'r*','MarkerSize',10);

%% Plot everything
fh=figure(35); clf
set(fh,'OuterPosition', [60 49 1162 1011],'color','w'),clf

tt=tiledlayout(9,nMet,"TileSpacing","compact",'Padding','tight');
all_clim={[0 2.5],[0,2.5],[0 2.5],[0 2.5]};
all_clim2={[0 1.5],[0,1.5],[0 1.5],[0 1.5]};
for i=1:nMet


    % plot FISP signal efficiency
    ax=nexttile(tt,i,[2 1]);
    imagesc((TR_all_gre*1e3),FA_all,sig_all_FISP(:,:,i)'),colorbar
    hold on
    plot(TR_all_gre*1e3,maxFA_gre,'LineWidth',2);
    mycolors = [1 0 0; 0 1 0; 0 0 1];
    ax = gca;
    ax.ColorOrder = mycolors;
    %     plotMax((TR_all*1e3),FA_all,sig_all_gre(:,:,i),ax)
    % imagesc((TR_all*1e3),FA_all,sig_all_bssfp(:,:,i)'),colorbar
    title(sprintf('%s | T1/T2* =%.0f/%.f ms',metabolites(i).name,metabolites(i).T1_s*1e3,metabolites(i).T2star_s*1e3))
    ax=gca;
    ax.XAxis.Direction="normal";
    ax.YAxis.Direction="normal";
    axis square
   set(gca,'clim',round(get(gca,'clim').*[0 1],1),'FontSize',10)
   clim(all_clim2{i})
    % xlabel('TR [ms]'),
    if(i==1),ylabel('FA [deg]'),end
    if(exist("plotProtTR36",'var')),hold on,plotProtTR36(), end %plotProtGre2()
    % legend('800 us','1400 us','2000 us','Location','northwest')



    ax=nexttile(tt,2*nMet +i,[2 1]);
    imagesc((TR_all_bssfp*1e3),FA_all,sig_all_bssfp(:,:,i)'),colorbar
    hold on
    plot(TR_all_bssfp*1e3,maxFA_bssfp,'LineWidth',2);
    mycolors = [1 0 0; 0 1 0; 0 0 1];
    ax = gca;
    ax.ColorOrder = mycolors;
    % plotMax((TR_all*1e3),FA_all,sig_all(:,:,i),ax)
    % imagesc((TR_all*1e3),FA_all,sig_all_bssfp(:,:,i)'),colorbar
    title(sprintf('%s | T1/T2 =%.0f/%.0f ms',metabolites(i).name,metabolites(i).T1_s*1e3,metabolites(i).T2_s*1e3))
    ax=gca;
    ax.XAxis.Direction="normal";
    ax.YAxis.Direction="normal";
    axis square
    set(gca,'clim',round(get(gca,'clim').*[0 1],1),'FontSize',10)
    clim(all_clim{i})
     % xlabel('TR [ms]'),
    if(i==1),ylabel('FA [deg]'),end
    if(exist("plotProtbssfp",'var')),hold on,plotProtbssfp(),end



    ax=nexttile(tt,2*2*nMet+i,[2 1]);
    imagesc((TR_all_bssfp*1e3),FA_all,sig_all_me(:,:,i)'),colorbar
    hold on
    plot(TR_all_bssfp*1e3,maxFA_bssfp,'LineWidth',2);
    mycolors = [1 0 0; 0 1 0; 0 0 1];
    ax = gca;
    ax.ColorOrder = mycolors;
    % plotMax((TR_all*1e3),FA_all,sig_all(:,:,i),ax)
    % imagesc((TR_all*1e3),FA_all,sig_all_bssfp(:,:,i)'),colorbar
    title(sprintf('%s | T1/T2 =%.0f/%.0f ms',metabolites(i).name,metabolites(i).T1_s*1e3,metabolites(i).T2_s*1e3))
    ax=gca;
    ax.XAxis.Direction="normal";
    ax.YAxis.Direction="normal";
    axis square
    set(gca,'clim',round(get(gca,'clim').*[0 1],1),'FontSize',10)
    clim(all_clim{i})
    xlabel('TR [ms]')
    if(i==1),ylabel('FA [deg]'),end
    if(exist("plotProtbssfp",'var')),hold on,plotProtbssfp(),end

end
% 
ax=nexttile(tt,nMet,[2 1]);
legend('0.5 ms','1.4 ms','2 ms','Location',[0.852487456415608 0.793195026049927 0.0791738371584378 0.0523012538272966])
sgtitle('signal efficiency [1/\surds]','Interpreter','tex','fontsize',24)

ax=nexttile(tt,1,[2 1]);
annotation(gcf,'textbox',...
    [0.01,ax.Position(2)+ax.Position(3)/2,0,0],...
    'String',{'CSI'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',15,'HorizontalAlignment','center');
ax=nexttile(tt,1+2*nMet,[2 1]);
annotation(gcf,'textbox',...
    [0.01,ax.Position(2)+ax.Position(3)/2,0,0],...
    'String',{'CSI-PC-bSSFP'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',15,'HorizontalAlignment','center');

ax=nexttile(tt,2*2*nMet+1,[2 1]);
annotation(gcf,'textbox',...
    [0.01,ax.Position(2)+ax.Position(3)/2,0,0],...
    'String',{'ME-PC-bSSFP'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',15,'HorizontalAlignment','center');
% plot protocol
TR_ssfp=19e-3; %s

[~,idxTR]=min(abs(TR_all_bssfp-TR_ssfp));
FA_bSSFP=maxFA_bssfp(idxTR,2);
[~,idxFA]=min(abs(FA_all-FA_bSSFP)); %use max flip angle
sig_bSSFP=squeeze(abs(sig_all_bssfp(idxTR,idxFA,1:nMet)));
sig_me=squeeze(abs(sig_all_me(idxTR,idxFA,1:nMet)));

sig_bSSFP_noPC=squeeze(abs(sig_all_bssfp_noPC(idxTR,idxFA,1:nMet)));
sig_me_noPC=squeeze(abs(sig_all_me_noPC(idxTR,idxFA,1:nMet)));

%FISP
TR_gre=36e-3; %s
FA_gre=41*RFfac; %deg
[~,idxTR]=min(abs(TR_all_gre-TR_gre));
[~,idxFA]=min(abs(FA_all-FA_gre)); %use max flip angle
sig_fisp=squeeze(abs(sig_all_FISP(idxTR,idxFA,1:nMet)));

TR_gre=36e-3; %s
FA_gre=41*RFfac; %deg
[~,idxTR]=min(abs(TR_all_gre-TR_gre));
[~,idxFA]=min(abs(FA_all-FA_gre)); %use max flip angle
sig_flash=squeeze(abs(sig_all_gre(idxTR,idxFA,1:nMet    )));


metnames=reordercats(categorical({metabolites.name}),{metabolites.name});
 nexttile(tt,nMet*2*3+1,[3 nMet])
if(true) % plot only phase-cycled
bh=barh(metnames,cat(2,sig_me,sig_bSSFP,sig_fisp),'FaceAlpha',0.75);
set(gca,'ColorOrder',flip(lines(3),1));%,[[ 0.8500,0.3250,0.0980];[ 0,0.4470,0.7410];[0.4940,0.1840,0.5560];])
lh=legend('ME-PC-bSSFP','CSI-PC-bSSFP','CSI','Location','southeast');
lh.Direction='reverse';
grid minor
set(gca,'FontSize',12)
title('signal efficiencies of the study protocols')
xlim(get(gca,'xlim')+[0 0.2])
lab_csi=string(strsplit(sprintf('%0.0fx ',(round([1 1 1 1])))));
text(bh(3).YEndPoints,bh(3).XEndPoints,lab_csi(1:nMet),"FontSize",12)

lab_csi=string(strsplit(sprintf('%0.2fx ',(round(sig_bSSFP./sig_fisp,2)))));
text(bh(2).YEndPoints,bh(2).XEndPoints,lab_csi(1:nMet   ),"FontSize",12)

lab_csi=string(strsplit(sprintf('%0.2fx ',(round(sig_me./sig_fisp,2)))));
text(bh(1).YEndPoints,bh(1).XEndPoints,lab_csi(1:nMet),"FontSize",12)


 else
hold on
dat=cat(1,[sig_flash,0*sig_flash],[sig_bSSFP,sig_bSSFP_noPC-sig_bSSFP],[sig_me,sig_me_noPC-sig_me]);
xindx=[(1:4)+0.25,1:4,(1:4)-0.25];
co=flip(lines(3),1);
clear bh;
bh{1}=barh (xindx(mod(xindx*10,10)==7.5),dat(mod(xindx*10,10)==7.5,:),    'stacked','BarWidth',0.2,'FaceColor',co(1,:));
bh{2}=barh (xindx(mod(xindx*10,10)==0),dat(mod(xindx*10,10)==0,:),    'stacked','BarWidth',0.2,'FaceColor',co(2,:));
bh{3}=barh (xindx(mod(xindx*10,10)==2.5),dat(mod(xindx*10,10)==2.5,:),    'stacked','BarWidth',0.2,'FaceColor',co(3,:));
bh{1}(2).FaceAlpha=0.2;
bh{2}(2).FaceAlpha=0.2;
bh{3}(2).FaceAlpha=0.2;

bh{1}(1).FaceAlpha=0.7;
bh{2}(1).FaceAlpha=0.7;
bh{3}(1).FaceAlpha=0.7;
% barh([1 2 3]*2+0.2,[sig_bSSFP,sig_bSSFP_noPC],'stacked','LineWidth',0.2)
lh=legend([bh{1}(1),bh{2}(1),bh{3}(1)],'ME-PC-bSSFP','CSI-PC-bSSFP','CSI','Location','southeast');
lh.Direction='reverse';
grid minor,grid on
set(gca,'FontSize',12)
title('signal efficiencies of the study protocols')
 xlim(get(gca,'xlim')+[-1 1]),ylim([0.4 4.6])
yticks(1:nMet),yticklabels({metabolites.name})

FtSize=12;
lab_csi=compose('%0.0fx',(round([1 1 1 1])));
text(bh{3}(1).YEndPoints,bh{3}(1).XEndPoints,lab_csi,"FontSize",FtSize)

lab_csi=compose('%0.2fx ',(round(sig_bSSFP./sig_fisp,2)));
text(bh{2}(1).YEndPoints,bh{2}(1).XEndPoints,lab_csi,"FontSize",FtSize)

lab_csi=compose('%0.2fx ',(round(sig_me./sig_fisp,2)));
text(bh{1}(1).YEndPoints,bh{1}(1).XEndPoints,lab_csi,"FontSize",FtSize)

lab_csi=compose('%+.0f%% ',100*(round(sig_bSSFP_noPC./sig_fisp-sig_bSSFP./sig_fisp,2)));
text(bh{2}(2).YEndPoints,bh{2}(2).XEndPoints,lab_csi,"FontSize",FtSize,'Color',[1,1,1]*0.6)

lab_csi=compose('%+0.0f%% ',100*(round(sig_me_noPC./sig_fisp-sig_me./sig_fisp,2)));
text(bh{1}(2).YEndPoints,bh{1}(2).XEndPoints,lab_csi,"FontSize",FtSize,'Color',[1,1,1]*0.6)

tab=table(sig_bSSFP,sig_me,sig_fisp,sig_flash,sig_bSSFP./sig_fisp,sig_me./sig_fisp, ...
    'RowNames',{metabolites.name},'VariableNames',{'CSI-PC-bSSFP','ME-PC-BSSFP','FISP','gre','CSI-bSSFP/FISP','ME-bSSFP/FISP'})
xlim([-0.5 3]),box on

end



%%

function plotMax(x,y,sig,ax)


[~,idx]=max(sig(:));
[xi,yi]=ind2sub([length(x) length(y)],idx);
% [~,yi]=max(sig,[],2);

plot(x(xi),y(yi),'r+','MarkerSize',10)
% scatter(gca,x(xi)*1e3,y(yi),50)


end
