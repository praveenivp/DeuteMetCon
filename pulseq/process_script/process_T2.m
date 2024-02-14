%% load sequence file
sn='/ptmp/pvalsala/deuterium/EAZW-GUMK/TWIX/';

addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'));
addpath(genpath('/ptmp/pvalsala/Packages/pulseq'));
%%


system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
seq.read(fullfile(sn,'../EXPDATA/pulseq/T2_nonsel_2H_10min.seq'))

st.dwell_s=seq.getDefinition('dwell');
st.TE_array=seq.getDefinition('TE_array');
st.rf_dur=seq.getDefinition('rf_dur');
st.averages=seq.getDefinition('averages');
st.repetitions=seq.getDefinition('repetitions');
%%
dirst=dir(fullfile(sn,'*pulseq_T2_10mins*.dat'));
st.filename=dirst(end).name; %load last file
twix=mapVBVD(fullfile(sn,st.filename));
if(iscell(twix)), twix=twix{end}; end
st.Cfreq=twix.hdr.Dicom.lFrequency; %Hz
st.hdr=twix.hdr;
st.RefVoltage= twix.hdr.Spice.TransmitterReferenceAmplitude;

   %% coil combination

data=twix.image{''};
data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions);

% noise decorr
dirst_noise=dir(fullfile(sn,'*noise*.dat'));
 twix_noise=mapVBVD(fullfile(dirst_noise(1).folder,dirst_noise(1).name));
 if(iscell(twix_noise)), twix_noise=twix_noise{end}; end
 [D_noise,D_image]=CalcNoiseDecorrMat(twix_noise);

data=permute(data,[2 1 4 3]);
data_whiten=reshape(D_image*data(:,:),size(data));
%  data_whiten=mean(data_whiten,4);
data_whiten=permute(data_whiten,[2 1 3 4]);
% % 
%                  Acqdelay=st.rf_dur/2+330e-6; %s
%                  ext_size=round(Acqdelay/st.dwell_s);
%                  [data_whiten] = fidExtrp(data_whiten,ext_size);


% Combine coil data
data_whiten=mean(data_whiten,4);
DataSize=size(data_whiten);
%we already noise decorrelated data!
wsvdOption.noiseCov         =0.5*eye(DataSize(2));
data_Combined=zeros([DataSize(1) DataSize(3:end)]);
% CoilWeights = zeros([DataSize(2) DataSize(3)]);
for rep=1:size(data_whiten,3)
    rawSpectra=squeeze(data_whiten(:,:,rep));  % fid x Coil
    [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(rawSpectra, [], wsvdOption);
%     CoilWeights(:,rep) = wsvdWeights;
    data_Combined(:,rep)=wsvdCombination;
    wsvdQuality_all(rep)=wsvdQuality;
end

%% phase correction
Acqdelay=seq.getBlock(5).blockDuration+seq.getBlock(4).rf.shape_dur/2+seq.getBlock(4).rf.ringdownTime; %s

ext_size=3;%round(Acqdelay/st.dwell_s)+1;
                 [data_Combined2] = fidExtrp(data_Combined(1:end,:),ext_size);




 spec1=squeeze(fftshift(fft(data_Combined2,[],1),1));

 faxis=linspace(-0.5/st.dwell_s,0.5/st.dwell_s,size(spec1,1));

  [phi0]=CalcZerothPhase(faxis,spec1,1000);
% spec2=spec1*exp(1i*median(phi0));
spec2=1*spec1.*exp(1i*(median(phi0)+faxis(:).*0.e-1));
disp(median(phi0))

figure(1),clf,plot(faxis,real(spec1))


%  Nlorentz fit
freq=[3.9 -50  -140 -199];
fitf=cell(size(spec2,2),1);
gof=cell(size(spec2,2),1);
figure(58),clf,hold on
amp_all=zeros(length(st.TE_array),length(freq));
FWHM_all=zeros(length(st.TE_array),length(freq));
for crep=1:15
[fitf{crep},gof{crep},foptions]=NLorentzFit(faxis,real(spec2(:,crep)),freq);

 amp_all(crep,:)=[fitf{crep}.A1 fitf{crep}.A2 fitf{crep}.A3 fitf{crep}.A4];
  FWHM_all(crep,:)=[fitf{crep}.gamma1 fitf{crep}.gamma2 fitf{crep}.gamma3 fitf{crep}.gamma4];
 plot(faxis,real(spec2(:,crep)))
 xlim([-500 500])
 plot(faxis,fitf{crep}(faxis))
end


save(fullfile(sn,'../proc/T1T2/T2_fit.mat'),'fitf','gof','foptions','st')
%%

figure(12),clf

tt=tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
nexttile(tt,1,[2 2])
hold on
for crep=1:length(fitf)
plot(faxis,fitf{crep}(faxis))
end
plot(faxis,real(spec2),'.','LineWidth',1.5)
 xlim([-300 200]),xlabel('frequency [Hz]')
 legend(strsplit(sprintf('TE= %.1f ms\n',st.TE_array*1e3),'\n'))
   set(gca,'ColorOrder',jet(size(spec2,2)));
peakName={'water','Glc','Glx','lactate/lipids'};
title('Global T2 : spin echo')
grid on

plt_id=[3 4, 7 8];
for i=1:length(freq)
    nexttile(tt,plt_id(i));


    lcolour='r'; typ='Lorentz';
    % Set up fittype and options.
    ft = fittype( 'a*(exp(-x/b))+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'on';
    opts.StartPoint = [1 0.1e3 0];
    opts.Lower = [-Inf 0 -Inf];
    opts.Upper = [Inf 2e3 Inf];
    % Fit model to data.
    [fitresult, gof] = fit( st.TE_array*1e3,col(amp_all(:,i)) , ft, opts );

    plot(st.TE_array*1e3,amp_all(:,i),[lcolour,'.'])
    hold on
    plot(st.TE_array*1e3,fitresult(st.TE_array*1e3),lcolour),
    lege{1,1}=sprintf('%s , T2= %.2f ms',typ,fitresult.b);
    lege{2,1}=sprintf('%.1f*(1-2*exp(-x/%.2f))+%.1f',fitresult.a,fitresult.b,fitresult.c);
    disp(fitresult)

    xlabel('TE [ms]')
    title(peakName{i})
    legend(lege{:},'Location','northeast')
       grid on

end
fontsize(gcf,"scale",1.3)

set(gcf,'color','w')
savefig(fullfile(sn,'../proc/T1T2/T2_final.fig'))

%% picked
spectrum=abs(fftshift(fft(data_Combined,2048,1)));
faxis=linspace(-0.5/st.dwell_s,0.5/st.dwell_s,length(spectrum));
freq_idx_all=[257 248 237 229 ];
figure
for i=1:4
[~,freq_idx]=min(abs(faxis-freq(i)))
% freq_idx=freq_idx_all(i);
subplot(2,2,i),plot(st.TI_array, mean(spectrum(freq_idx+[-1:1],:),1))
end


