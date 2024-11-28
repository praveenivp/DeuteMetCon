
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))

refVoltage=480; % upto 480-520 V
kFactor=0.83;
RFfac=447/refVoltage; % Flip angle scale factor

metabolites=getMetaboliteStruct('invivo');
nMet=3;
metabolites=metabolites(1:nMet);

pc_range=linspace(0,360-360/4,4)+180;
TR_all_bssfp=linspace(15e-3,25e-3,50);
TR_all_gre=linspace(15e-3,100e-3,80);
FA_all=linspace(10,90,60);
TE=2.3e-3;
apply_T2star=true;
calc_signaleff=true;


B0=randn([1000,1])*2*6.5360*9.4; % 2ppm

%% SSFP signal signal

sig_all_bssfp=zeros(length(TR_all_bssfp),length(FA_all),length(metabolites));
dc_fac_ssfp=zeros(length(TR_all_bssfp),1,length(metabolites));
maxFA_bssfp=zeros([length(TR_all_bssfp) 3]);
 DC=@(TR_s) ((TR_s-4.2e-3)/TR_s).*double(TR_s>4.1e-3); % 4.2 ms non-encoing time csi-bSSFP
%  DC=@(TR_s) ((TR_s-5.5e-3)/TR_s).*double(TR_s>5.5e-3); % 5.5 ms non-encoing time ME-bSSFP
% figure,tiledlayout("flow")
for cTR=1:size(sig_all_bssfp,1)
    [~,maxFA_bssfp(cTR,1)]=SimpleSARModel(1,500e-6,TR_all_bssfp(cTR),refVoltage,kFactor);
    [~,maxFA_bssfp(cTR,2)]=SimpleSARModel(1,1400e-6,TR_all_bssfp(cTR),refVoltage,kFactor);
    [~,maxFA_bssfp(cTR,3)]=SimpleSARModel(1,2000e-6,TR_all_bssfp(cTR),refVoltage,kFactor);
    for CM=1:size(sig_all_bssfp,3)
        %off-resonance is set to -1*chemical shift otherwise increase phase-cyles
        [Msig_all,dc_fac_ssfp(cTR,1,CM)]=MetSignalModel(metabolites(CM),TE,deg2rad(pc_range), ...
            TR_all_bssfp(cTR),B0,deg2rad(FA_all),'bSSFP',DC);
        sig_all_bssfp(cTR,:,CM)=abs(mean(abs(Msig_all),[3,5]));
        %            nexttile(), plot(pc_range,abs(Msig_all(:)))
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
    sig_all_bssfp=sig_all_bssfp.*sqrt(dc_fac_ssfp);
end

%calcualte signal efficiency
if(calc_signaleff)
    sig_all_gre=sig_all_gre./sqrt(TR_all_gre(:));
    sig_all_FISP=sig_all_FISP./sqrt(TR_all_gre(:));
    sig_all_bssfp=sig_all_bssfp./sqrt(TR_all_bssfp(:));

end



%some plot handles
plotProtTR74= @ () plot(74,57*RFfac,'r*','MarkerSize',10);
plotProtTR36= @ () plot(36,41*RFfac,'r*','MarkerSize',10);
plotProtbssfp= @ () plot(19,50*RFfac,'r*','MarkerSize',10);

%% Plot everything
fh=figure(3);
set(fh,'OuterPosition', [72 42 1e3 1e3]),clf

tt=tiledlayout(4,nMet*2,"TileSpacing","compact",'Padding','compact');

for i=1:nMet
    % plot GRE signal efficiency
    ax=nexttile(tt,2*(i-1)+1,[1 2]);
    imagesc((TR_all_gre*1e3),FA_all,sig_all_gre(:,:,i)'),colorbar
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
    xlabel('TR [ms]'),ylabel('FA [deg]')
     if(exist("plotProtTR36",'var')),hold on,plotProtTR36(), end %plotProtGre2()
    % legend('800 us','1400 us','2000 us','Location','northwest')

    % plot GRE signal efficiency
    ax=nexttile(tt,2*(i-1)+1+2*nMet,[1 2]);
    imagesc((TR_all_gre*1e3),FA_all,sig_all_FISP(:,:,i)'),colorbar
    hold on
    plot(TR_all_gre*1e3,maxFA_gre,'LineWidth',2);
    mycolors = [1 0 0; 0 1 0; 0 0 1];
    ax = gca;
    ax.ColorOrder = mycolors;
    %     plotMax((TR_all*1e3),FA_all,sig_all_gre(:,:,i),ax)
    % imagesc((TR_all*1e3),FA_all,sig_all_bssfp(:,:,i)'),colorbar
    title(sprintf('%s | T1/T2 =%.0f/%.0f ms',metabolites(i).name,metabolites(i).T1_s*1e3,metabolites(i).T2_s*1e3))
    ax=gca;
    ax.XAxis.Direction="normal";
    ax.YAxis.Direction="normal";
    axis square
   set(gca,'clim',round(get(gca,'clim').*[0 1],1),'FontSize',10)
    xlabel('TR [ms]'),ylabel('FA [deg]')
    if(exist("plotProtTR36",'var')),hold on,plotProtTR36(), end %plotProtGre2()
    % legend('800 us','1400 us','2000 us','Location','northwest')



    ax=nexttile(tt,2*(i-1)+1+2*2*nMet,[1 2]);
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
    xlabel('TR [ms]'),ylabel('FA [deg]')
    if(exist("plotProtbssfp",'var')),hold on,plotProtbssfp(),end


end
ax=nexttile(tt,2*(i-1)+1+2*nMet,[1 2]);
legend('700 us','1400 us','2000 us','Location','southeast')
sgtitle('signal efficiency [1/\surds]','Interpreter','tex','fontsize',24)

ax=nexttile(tt,1);
annotation(gcf,'textbox',...
    [0.01,ax.Position(2)+ax.Position(3)/2-0.04,0,0],...
    'String',{'FLASH'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',20);


ax=nexttile(tt,1+nMet*2);
annotation(gcf,'textbox',...
    [0.01,ax.Position(2)+ax.Position(3)/2-0.02,0,0],...
    'String',{'FISP'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',20);

ax=nexttile(tt,1+nMet*2*2);
annotation(gcf,'textbox',...
    [0.01,ax.Position(2)+ax.Position(3)/2-0.04,0,0],...
    'String',{'bSSFP'},...
    'Rotation',90,...
    'FontWeight','bold',...
    'FontSize',20);
% plot protocol
TR_ssfp=19e-3; %s

[~,idxTR]=min(abs(TR_all_bssfp-TR_ssfp));
FA_bSSFP=maxFA_bssfp(idxTR,2);
[~,idxFA]=min(abs(FA_all-FA_bSSFP)); %use max flip angle
sig_prot1=squeeze(abs(sig_all_bssfp(idxTR,idxFA,1:nMet)));

%FISP
TR_gre=36e-3; %s
FA_gre=41*RFfac; %deg
[~,idxTR]=min(abs(TR_all_gre-TR_gre));
[~,idxFA]=min(abs(FA_all-FA_gre)); %use max flip angle
sig_prot2=squeeze(abs(sig_all_FISP(idxTR,idxFA,1:nMet)));

TR_gre=36e-3; %s
FA_gre=41*RFfac; %deg
[~,idxTR]=min(abs(TR_all_gre-TR_gre));
[~,idxFA]=min(abs(FA_all-FA_gre)); %use max flip angle
sig_prot3=squeeze(abs(sig_all_gre(idxTR,idxFA,1:nMet    )));


metnames=reordercats(categorical({metabolites.name}),{metabolites.name});
nexttile(tt,[1 3])
barh(metnames,cat(2,sig_prot1,sig_prot2,sig_prot3))
legend('bSSFP','FISP','FLASH','Location','southeast')
grid minor
set(gca,'FontSize',12)
title('Signal efficiencies for the investigated protocols')
xlim(get(gca,'xlim')+[0 0.7])

nexttile(tt,[1 3])
barh(metnames,cat(2,sig_prot1./sig_prot2,sig_prot2./sig_prot3))
legend('bSSFP/FISP','FISP/FLASH','Location','southeast')
grid minor
title('ratio')
set(gcf,'color','w')
xlim(get(gca,'xlim')+[1 -0.2])
set(gca,'FontSize',12)

tab=table(sig_prot1,sig_prot2,sig_prot3,sig_prot1./sig_prot2,sig_prot1./sig_prot3,'RowNames',{metabolites.name},'VariableNames',{'bSSFP','FISP','gre','bSSFP/FISP','bSSFP/GRE'})


% print(gcf,'fig2_SNRsiminvivo3','-dpng','-r300')

%%

function plotMax(x,y,sig,ax)


[~,idx]=max(sig(:));
[xi,yi]=ind2sub([length(x) length(y)],idx);
% [~,yi]=max(sig,[],2);

plot(x(xi),y(yi),'r+','MarkerSize',10)
% scatter(gca,x(xi)*1e3,y(yi),50)


end
