%%

pulseDur=[500e-6,1400e-6,2000e-6]; %[s]
sys_freq_offset=0; %[Hz]
met_st=getMetaboliteStruct('invivo');

figure(37),clf
for ii=1:length(pulseDur)
dt=1000e-9; %[s]


t=0:dt:2*pulseDur(ii);
rf=zeros(size(t));
rf(1:round(pulseDur(ii)/dt))=1;

rf_fft=fftshift(fft(rf,2^20));
faxis=linspace(-0.5/dt,0.5/dt,length(rf_fft));

rf_fft=rf_fft./max(abs(rf_fft));
plot(faxis,abs(rf_fft),'LineWidth',1.5)
hold on
scale_fac=interp1(faxis,abs(rf_fft),[met_st.freq_shift_Hz]+sys_freq_offset)
end



% plot simulated spectrum

% build n-peak Lorentzian model function
model = @(cf,a,g)a./(1 + ((faxis - cf)/(g/2)).^2) ; %  Lorentzian
sig_sim=zeros(size(faxis));
amp=[1 0.3 0.5 0.2];

for i=1:length(met_st)
    fwhm=0.1/(pi*met_st(i).T2star_s);
    sig_sim=sig_sim+model(met_st(i).freq_shift_Hz+sys_freq_offset,amp(i),fwhm);
end

xlim([-1 1]*300)
plot(faxis,abs(sig_sim),'LineWidth',1.5);
xlabel('frequency [Hz]'),ylabel('normalized amplitude [a.u]')
title (sprintf('pulse profile'))
grid on
legend('0.5 ms' ,'1.4 ms','2 ms','2H spectrum','Location','southeast')
fontsize(gcf,'scale',1.5)
set(gcf,'color','w')
% scalefactor of flip angle?

