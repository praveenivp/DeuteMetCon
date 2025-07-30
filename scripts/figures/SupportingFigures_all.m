%% why phase cycles decreases SNR
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'))

refVoltage=520; % upto 480-520 V
RFfac=447/refVoltage; % Flip angle scale factor

metabolites=getMetaboliteStruct('invivo');
metabolites=metabolites(1:3);

pc_range=linspace(0,360,40);
TR=19e-3;
FA=50;
TE=2e-3;
B0=0;

figure,tt=tiledlayout(1,2,"TileSpacing","tight",'Padding','compact')
nexttile()
DC=@(TR_s) ((TR_s-4.2e-3)/TR_s).*double(TR_s>4.1e-3); % 4.2 ms non-encoing time csi-bSSFP
[Msig_all,dc_fac_ssfp]=MetSignalModel(metabolites,TE,deg2rad(pc_range), ...
    TR,-1*[metabolites.freq_shift_Hz],deg2rad(FA*RFfac),'bSSFP',DC);
DC=@(TR_s) ((TR_s-7.56e-3)/TR_s).*double(TR_s>7.56e-3); % 4.2 ms non-encoing time csi-FISP
[Msig_all_FISP,dc_fac_ssfp]=MetSignalModel(metabolites,TE,0, ...
    36e-3,0,deg2rad(41*RFfac),'FISP',DC);
Msig_all=cat(2,squeeze(Msig_all(1,:,:,:,1)),squeeze(Msig_all(2,:,:,:,2)),squeeze(Msig_all(3,:,:,:,3)));
plot(pc_range,squeeze(abs(Msig_all)),'LineWidth',2)
hold on
set(gca,'ColorOrder',lines(4),'ColorOrderIndex',1)
plot(pc_range, mean(abs(Msig_all_FISP),10)'.*ones(size(Msig_all)),'--','LineWidth',2)
legend({metabolites.name},'Location','south')
title('invivo'),grid on,xlabel('RF phase increment [deg]'),ylabel('amplitude [n.u]')
set(gcf,'Color','w','Position',[322 427 982 442])

metabolites=getMetaboliteStruct('phantom');
refVoltage=600; % V
RFfac=447/refVoltage; % Flip angle scale factor
nexttile()

DC=@(TR_s) ((TR_s-4.2e-3)/TR_s).*double(TR_s>4.1e-3); % 4.2 ms non-encoing time csi-bSSFP
[Msig_all,dc_fac_ssfp]=MetSignalModel(metabolites,TE,deg2rad(pc_range), ...
    TR,-1*[metabolites.freq_shift_Hz],deg2rad(FA*RFfac),'bSSFP',DC);
DC=@(TR_s) ((TR_s-7e-3)/TR_s).*double(TR_s>4.1e-3); % 4.2 ms non-encoing time csi-FISP
[Msig_all_FISP,dc_fac_ssfp]=MetSignalModel(metabolites,TE,0, ...
    36e-3,0,deg2rad(41*RFfac),'FISP',DC);
Msig_all=cat(2,squeeze(Msig_all(1,:,:,:,1)),squeeze(Msig_all(2,:,:,:,2)),squeeze(Msig_all(3,:,:,:,3)),squeeze(Msig_all(4,:,:,:,4)));
plot(pc_range,squeeze(abs(Msig_all)),'LineWidth',2)
hold on
set(gca,'ColorOrder',lines(4),'ColorOrderIndex',1)
plot(pc_range, mean(abs(Msig_all_FISP),10)'.*ones(size(Msig_all)),'--','LineWidth',2)
legend({metabolites.name},'Location','south')
title('Phantom'),grid on,xlabel('RF phase increment [deg]')
fontsize(gcf,"scale",1.2)
print(gcf,'sx_phase_cycle_SNR','-dpng','-r300')
%% plot simulated mode and Eigen vector
figure(7),clf

refVoltage=500; % upto 480-520 V
RFfac=447/refVoltage; % Flip angle scale factor

metabolites=getMetaboliteStruct('invivo');
metabolites=metabolites(1:3);
 metabolites(3).T2_s=40e-3;
pc_range=linspace(0,360,40);
TR=19e-3;
FA=50;
TE=4e-3;
B0=0;

figure,tt=tiledlayout(1,2,"TileSpacing","tight",'Padding','compact')
nexttile()
DC=@(TR_s) ((TR_s-4.2e-3)/TR_s).*double(TR_s>4.1e-3); % 4.2 ms non-encoing time csi-bSSFP
[Msig_all,dc_fac_ssfp]=MetSignalModel(metabolites,TE,deg2rad(pc_range), ...
    TR,B0,deg2rad(FA*RFfac),'bSSFP',DC);
plot(-10:10,abs(squeeze(calc_Fn2(Msig_all,deg2rad(pc_range),10)))','-x','LineWidth',2)
legend({metabolites.name})
 xlim([-4,4]),
xticks(-10:10)
xticklabels(strsplit(sprintf('F_{%d} ',-10:10),' '))
title('In vivo modes decay : simulated')
xlabel('SSFP Configuration order'),ylabel('amp [a.u]')
% hold on,yyaxis right
% set(gca,'ColorOrder',lines(4))
% plot(-10:10,unwrap(angle(squeeze(calc_Fn2(Msig_all,deg2rad(pc_range),10)))'),'--','LineWidth',2)


nexttile()
% H4 ME-bSSFP mcobj
SF=cell2mat(mcobj_me{3}.Experimental.sing_val(1:3)').^2;
plot(-8:8,abs(cell2mat(mcobj_me{3}.Experimental.sing_vec(1:3)'))./SF(1,:),'-x','LineWidth',2)
legend({metabolites.name})
 xlim([-4,4]),
xticks(-10:10)
xticklabels(strsplit(sprintf('F_{%d} ',-10:10),' '))
title('In vivo : scaled eigen vectors')
xlabel('SSFP Configuration order'),ylabel('amp [a.u]')
fontsize(gcf,"scale",1.5)
set(gcf,'color','w','Position',[511 571 1157 420])
print(gcf,'sx_modes_eigen_vec','-dpng','-r300')


%% plot conditioning/NSA
metabolites=getMetaboliteStruct('invivo');


% % 3T sim
% for i=1:4
% metabolites(i).freq_shift_Hz=metabolites(i).freq_shift_Hz.*(3/9.4);
% end
% deltaT=(0.01:0.05:1)*1e-3;
% N=64*4;
% deltaTE=(0.1:0.01:20)*1e-3;
% NE=5;


deltaT=(0.01:0.01:0.5)*1e-3;
N=64;
deltaTE=(0.1:0.01:6)*1e-3;
NE=5;

NSA_CSI=zeros(length(metabolites),length(deltaT));
for cDT=1:length(deltaT)
    TE=2e-3+(0:N-1)*deltaT(cDT);
    IDEALobj_pinv=IDEAL(metabolites,TE,'solver','phaseonly','PhaseCorr',true);
    A=(IDEALobj_pinv.getA);
    NSA_CSI(:,cDT)=abs(diag(inv(A'*A)));
end

NSA_ME=zeros(length(metabolites),length(deltaTE));
for cDT=1:length(deltaTE)
    TE=2e-3+(0:NE-1)*deltaTE(cDT);
    IDEALobj_pinv=IDEAL(metabolites,TE,'solver','phaseonly','PhaseCorr',true);
    A=(IDEALobj_pinv.getA);
    NSA_ME(:,cDT)=abs(diag(inv(A'*A)));
end


figure(56),
tt=tiledlayout(1,2,"TileSpacing","compact",'Padding','compact');

nexttile()
plot(deltaTE*1e3,1./(NE*NSA_ME),'LineWidth',2)
hold on
line([1 1]*3.7,[0 1],'linewidth',2,'color','magenta')
xlabel('echo spacing [ms]'),ylabel('Normalized NSA')
title('multi-echo | N=5'),grid on
legend([{metabolites.name},'ME-bSSFP'],'Location','northwest')

nexttile()
plot(deltaT*1e3*N,1./(N*NSA_CSI),'LineWidth',2)
title('CSI '),grid on
xlabel('readout length [ms]'),ylabel('Normalized NSA')
hold on
line([1 1]*14.9,[0 1],'linewidth',2,'color','magenta' ),line([1 1]*28.44,[0 1],'linewidth',2,'color','cyan')
legend([{metabolites.name},'CSI-bSSFP','CSI-FISP'],'Location','northwest')
fontsize(gcf,'scale',1.5)
set(gcf,'color','w','Position', [335 381 1289 556])
print(gcf,'sx_NSA','-dpng','-r300')
%% optimal flip angles differs for bSSFP with phase-cycling
offset=0;
pc_range=180;%linspace(0,360-360/18,18)+offset;
refVoltage=500; % upto 480-520 V
RFfac=447/refVoltage; % Flip angle scale factor

TR=19e-3;
FA=1:120*RFfac;
TE=TR/2;
B0=0;
% SSFP signal signal

 DC=@(TR_s) ((TR_s-4.2e-3)/TR_s).*double(TR_s>4.1e-3); % 4.2 ms non-encoing time csi-bSSFP
        %off-resonance is set to -1*chemical shift otherwise increase phase-cyles
        [Msig_all,dc_fac_ssfp]=MetSignalModel(metabolites(1:3),TE,deg2rad(pc_range), ...
            TR,B0,deg2rad(FA),'bSSFP',DC);

         pc_range_4pc=linspace(0,360-360/4,4)+offset;
       [Msig_all_4pc,dc_fac_ssfp]=MetSignalModel(metabolites(1:3),TE,deg2rad(pc_range_4pc), ...
            TR,B0,deg2rad(FA),'bSSFP',DC);
        figure(1),clf
        tiledlayout(1,2)
        nexttile()
        plot(FA./RFfac,squeeze(mean(abs(Msig_all_4pc),3))')
        xlabel('FA [deg]'),grid on
        title('4 PC')

             nexttile()
        plot(FA./RFfac,squeeze(mean(abs(Msig_all),3))')
        hold on
        plot(FA./RFfac,squeeze(mean(abs(Msig_all_4pc),3))','--')
        xlabel('FA [deg]'),grid on
        title('1 PC')

%%  FISP and FLASH signal for Glx at optimal flip angle vs TR
metabolites=getMetaboliteStruct('invivo');
cMet=3;
metabolites(cMet).T2star_s=23e-3;
TR_all=linspace(15e-3,metabolites(cMet).T1_s*3,500);
FA_all=acos(exp(-TR_all./metabolites(3).T1_s));
DC=@(TR_s) ((TR_s-7.56e-3)/TR_s).*double(TR_s>7.56e-3); % 7.56 ms non-encoing time
TE=2.3e-3;

xaxis=TR_all*1e3;
xaxis=TR_all./metabolites(cMet).T1_s;

[Msig_all,dc_fac]=MetSignalModel(metabolites(cMet),0,0, ...
    TR_all,0,FA_all,'FLASH',DC);
figure(131),clf%
 plot(xaxis,squeeze(max(Msig_all,[],6))./sqrt(TR_all'),'LineWidth',2)
hold on
% plot(xaxis,sqrt(DC(TR_all))'.*squeeze(max(Msig_all,[],6))./sqrt(TR_all'),'LineWidth',2)
  plot(xaxis,10*(dc_fac(:)).*squeeze(max(Msig_all,[],6)),'LineWidth',2) % T2star

[Msig_all,dc_fac]=MetSignalModel(metabolites(cMet),0,0, ...
    TR_all,0,FA_all,'FISP',DC);


plot(xaxis,squeeze(max(Msig_all,[],6))./sqrt(TR_all'),'LineWidth',2)
 % plot(xaxis,sqrt(DC(TR_all))'.*squeeze(max(Msig_all,[],6))./sqrt(TR_all'),'LineWidth',2)
   plot(xaxis,10*(dc_fac(:)).*squeeze(max(Msig_all,[],6)),'LineWidth',2) %T2star



xlabel('TR/$T_1$','Interpreter','latex'),ylabel('signal eff, [1/$\sqrt{TR}$]','Interpreter','latex')
grid on,ylim([0.4 2.1])

% stem((36e-3)/metabolites(cMet).T1_s ,2)

stem((metabolites(cMet).T2star_s*1.26+7.56e-3) /metabolites(cMet).T1_s ,2)
legend('FLASH','FLASH-T_2* effects','FISP','FISP-T_2* effects','TR=36 ms')
fontsize('scale',1.2)

%% show 1.26*T2* is optimum readout length

figure(131),clf,plot(TR_all*1e3,dc_fac.*sqrt(TR_all)) % SNR metric
% plot(TR_all*1e3,dc_fac) % SNR metric
hold on
 stem(metabolites(3).T2star_s*1e3*1.26+7.56,0.1)
