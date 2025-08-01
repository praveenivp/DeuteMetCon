%% load sequence file
clearvars,%clc
MeasPath='/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/sub-01';
metabolites=getMetaboliteStruct('invivo');
flip=false;

% addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'));
addpath(genpath('/ptmp/pvalsala/Packages/pulseq'));
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))


%% sim
figure(13),clf

tt=tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
for jj=1:3
% simulate
SNR=2*jj; %20-30
noise=randn(512,15);
T2=100e-3;
T2star_s=23e-3;
taxis=(0:st.dwell_s:st.dwell_s*511)';
a=exp(-st.TE_array/100e-3);
model_time = @(cf,a,t2)a.*exp(-taxis/t2).*exp(2i*pi*cf*taxis); %  Lorentzian

fid=model_time(154,a',T2star_s);
spec1=specFft(fid);

% Calculate signal power (peak amplitude in spectrum)
signal_amplitude = a(1)*T2star_s*(1-exp(128e-3/T2star_s));%max(abs(spec1(:))); % Peak amplitude in frequency domain

% Generate noise
noise = randn(size(fid)) + 1i*randn(size(fid)); % Complex Gaussian noise
noise_std = ( signal_amplitude/ SNR); % Scale noise to achieve desired SNR
noise = noise * noise_std / std(noise(:)); % Normalize noise to desired std


data_Combined=model_time(154,a',T2star_s)+noise;
spec1=specFft(data_Combined);
SNR_measured=max(abs(spec1(:)))/std(real(spec1(350:450,end)));



nexttile(tt)
faxis=calcFreqAxis(st.dwell_s,size(spec1,1))./(9.938*6.56);
plot(faxis,real(spec1)),xlim([0 8]),
ax=gca;
ax.XAxis.Direction='reverse';
text(ax.XLim(2)-2,ax.YLim(2)*0.8,sprintf('SNR=%d',SNR))

%  amares fit

st.AcqDelay_s=0;
[expParams,pk]=getAMARES_structs_T1inv(twix,metabolites(3),st);

spec1=specFft(padarray(data_Combined,[512*4 0],'post'));

faxis=calcFreqAxis(st.dwell_s,size(spec1,1));
pk.bounds(1).chemShift=[1.5 3.5];
pk.bounds(1).phase=[0 2*pi];
pk.bounds(1).amplitude=[0 Inf];
pk.bounds(1).linewidth=[ 10 20 ];
for i=1:size(data_Combined,2)
    [fitResults{i}, fitStatus{i},~,CRBResults{i}] = AMARES.amaresFit(double(data_Combined(:,i)), expParams, pk,0,'quiet',true);
    amp_all(i)=[fitResults{i}.amplitude];
    % pause(1)
end

%  plot overview
nexttile(tt)
hold on
weights=1./cell2mat(cellfun(@(x)x.amplitude,CRBResults,'UniformOutput',false)');
lcolour='r'; typ='AMARES';

ft = fittype( 'a*(exp(-x/b))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares'  );
opts.Display = 'Off';
opts.Robust = 'on';
opts.StartPoint = [1 100];
opts.Lower = [-10 10];
opts.Upper = [Inf 300];

% Fit model to data.
[fitresult, gof] = fit( col(st.TE_array*1e3),col(amp_all) , ft, opts );

plot(st.TE_array*1e3,col(amp_all),[lcolour,'x'],'LineWidth',1.5)
hold on
plot(st.TE_array*1e3,fitresult(st.TE_array*1e3),lcolour,'LineWidth',1.5),
lege{1,1}=sprintf('%s , T2= %.2f ms',typ,fitresult.b);
lege{2,1}=sprintf('%.2f*exp(-t/%.0f)',fitresult.a,fitresult.b);


xlabel('TE [ms]')
title(peakName{3})
legend(lege{:},'Location','best')
grid on

end
fontsize(gcf,"scale",1.3)

% set(gcf,'color','w','Position',[249 195 1661 727])
% savefig(fullfile(MeasPath,'T2_sim.fig'))
