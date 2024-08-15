
refVoltage=520; % upto 520
kFactor=0.83;
RFfac=447/refVoltage; % Flip angle scale factor

metabolites=getMetaboliteStruct('phantom');
pc_range=linspace(0,359,18)+180;
TR_all_bssfp=linspace(15e-3,25e-3,50);
TR_all_gre=linspace(15e-3,100e-3,80);
FA_all=linspace(10,80,60);
TE=2e-3;

%% SSFP signal signal

sig_all_bssfp=zeros(length(TR_all_bssfp),length(FA_all),length(metabolites));
maxFA_bssfp=zeros([length(TR_all_bssfp) 3]);
% figure,tiledlayout("flow")
for cTR=1:size(sig_all_bssfp,1)
    [~,maxFA_bssfp(cTR,1)]=SimpleSARModel(1,800e-6,TR_all_bssfp(cTR),refVoltage,kFactor);
    [~,maxFA_bssfp(cTR,2)]=SimpleSARModel(1,1400e-6,TR_all_bssfp(cTR),refVoltage,kFactor);
    [~,maxFA_bssfp(cTR,3)]=SimpleSARModel(1,2000e-6,TR_all_bssfp(cTR),refVoltage,kFactor);
    for CM=1:size(sig_all_bssfp,3)
    %off-resonance is set to -1*chemical shift otherwise increase phase-cyles
        [Msig_all]=MetSignalModel   (metabolites(CM),TE,deg2rad(pc_range),TR_all_bssfp(cTR),-1*metabolites(CM).freq_shift_Hz,deg2rad(FA_all),'bSSFP-peters');
        sig_all_bssfp(cTR,:,CM)=mean(abs(Msig_all),3);
        %            nexttile(), plot(pc_range,abs(Msig_all(:)))
    end
end
DC=0.79;
%calcualte signal efficiency
sig_all_bssfp=sig_all_bssfp./sqrt(TR_all_bssfp(:))*sqrt(DC);
%% GRE signal efficiency
sig_all_gre=zeros(length(TR_all_gre),length(FA_all),length(metabolites));
maxFA_gre=zeros([length(TR_all_gre) 3]);
for cTR=1:size(sig_all_gre,1)
    [~,maxFA_gre(cTR,1)]=SimpleSARModel(1,800e-6,TR_all_gre(cTR),refVoltage,kFactor);
    [~,maxFA_gre(cTR,2)]=SimpleSARModel(1,1400e-6,TR_all_gre(cTR),refVoltage,kFactor);
    [~,maxFA_gre(cTR,3)]=SimpleSARModel(1,2000e-6,TR_all_gre(cTR),refVoltage,kFactor);
    for CM=1:size(sig_all_gre,3)
        DC=(TR_all_gre(cTR)-7e-3)/TR_all_gre(cTR);
        [Msig_all]=MetSignalModel   (metabolites(CM),TE,0,TR_all_gre(cTR),-1*metabolites(CM).freq_shift_Hz,deg2rad(FA_all),'GRE-peters',DC);
        sig_all_gre(cTR,:,CM)=mean(abs(Msig_all),3);

    end
end
%calcualte signal efficiency
sig_all_gre=sig_all_gre./sqrt(TR_all_gre(:));


fh=figure(3);
set(fh,'OuterPosition',[1 41 1920 1082]),clf
tt=tiledlayout(3,4,"TileSpacing","compact",'Padding','compact');
for i=1:length(metabolites)
    

    % plot GRE signal efficiency
    ax=nexttile(tt,i);
    imagesc((TR_all_gre*1e3),FA_all,sig_all_gre(:,:,i)'),colorbar
    hold on
    plot(TR_all_gre*1e3,maxFA_gre,'LineWidth',2);
    mycolors = [1 0 0; 0 1 0; 0 0 1];
    ax = gca;
    ax.ColorOrder = mycolors;
%     plotMax((TR_all*1e3),FA_all,sig_all_gre(:,:,i),ax)
    % imagesc((TR_all*1e3),FA_all,sig_all_bssfp(:,:,i)'),colorbar
    title(sprintf('%s | T1/T2* =%.2f/%.2f ms',metabolites(i).name,metabolites(i).T1_s*1e3,metabolites(i).T2star_s*1e3))
    ax=gca;
    ax.XAxis.Direction="normal";
    ax.YAxis.Direction="normal";
    axis square
    % clim([0 0.5])
    xlabel('TR [ms]'),ylabel('FA [deg]')
    if(exist("plotProtGre2",'var')),hold on,plotProtGre1(),plotProtGre2(), end
    % legend('800 us','1400 us','2000 us','Location','northwest')



    ax=nexttile(tt,i+4);
    imagesc((TR_all_bssfp*1e3),FA_all,sig_all_bssfp(:,:,i)'),colorbar
    hold on
    plot(TR_all_bssfp*1e3,maxFA_bssfp,'LineWidth',2);
    mycolors = [1 0 0; 0 1 0; 0 0 1];
    ax = gca;
    ax.ColorOrder = mycolors;
    % plotMax((TR_all*1e3),FA_all,sig_all(:,:,i),ax)
    % imagesc((TR_all*1e3),FA_all,sig_all_bssfp(:,:,i)'),colorbar
    title(sprintf('%s | T1/T2 =%.2f/%.2f ms',metabolites(i).name,metabolites(i).T1_s*1e3,metabolites(i).T2_s*1e3))
    ax=gca;
    ax.XAxis.Direction="normal";
    ax.YAxis.Direction="normal";
    axis square
    % clim([0 0.5])
    xlabel('TR [ms]'),ylabel('FA [deg]')
    if(exist("plotProtbssfp",'var')),hold on,plotProtbssfp(),end


end
legend('800 us','1400 us','2000 us','Location','southeast')
sgtitle('signal efficiency [1/\surds]','Interpreter','tex')


%% plot protocol
TR_ssfp=19e-3; %s

[~,idxTR]=min(abs(TR_all_bssfp-TR_ssfp));
FA_bSSFP=maxFA_bssfp(idxTR,2);
[~,idxFA]=min(abs(FA_all-FA_bSSFP)); %use max flip angle
sig_prot1=squeeze(abs(sig_all_bssfp(idxTR,idxFA,:)));

%1.26 T2 star readout
TR_gre=36e-3; %s
FA_gre=41*RFfac; %deg
[~,idxTR]=min(abs(TR_all_gre-TR_gre));
[FA_gre,idxFA]=min(abs(FA_all-FA_gre)); %use max flip angle
sig_prot2=squeeze(abs(sig_all_gre(idxTR,idxFA,:)));

TR_gre=74e-3; %s
FA_gre=57*RFfac; %deg
[~,idxTR]=min(abs(TR_all_gre-TR_gre));
[FA_gre,idxFA]=min(abs(FA_all-FA_gre)); %use max flip angle
sig_prot3=squeeze(abs(sig_all_gre(idxTR,idxFA,:)));


metnames=reordercats(categorical({metabolites.name}),{metabolites.name});
nexttile(tt,9,[1 2])
barh(metnames,cat(2,sig_prot1,sig_prot2,sig_prot3))
legend('bSSFP','GRE-TR36','GRE-TR74','Location','southeast')
grid minor
title('Signal efficiencies for the investigated protocols')

nexttile(tt,11,[1 2])
barh(metnames,cat(2,sig_prot1./sig_prot2,sig_prot2./sig_prot3))
legend('bSSFP/GRE-TR36','GRE-TR36/TR74','Location','southwest')
grid minor
title('ratio')

tab=table(sig_prot1,sig_prot2,sig_prot3,'RowNames',{metabolites.name},'VariableNames',{'bssfp19ms','gre36ms','gre74ms'})

%some plot handles: needs second execition
plotProtGre1= @ () plot(74,57*RFfac,'r*','MarkerSize',10);
plotProtGre2= @ () plot(36,41*RFfac,'r*','MarkerSize',10);
plotProtbssfp= @ () plot(19,50*RFfac,'r*','MarkerSize',10);


%%

function plotMax(x,y,sig,ax)


[~,idx]=max(sig(:));
[xi,yi]=ind2sub([length(x) length(y)],idx);
% [~,yi]=max(sig,[],2);

plot(x(xi),y(yi),'r+','MarkerSize',10)
% scatter(gca,x(xi)*1e3,y(yi),50)


end
