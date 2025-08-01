%% load sequence file
clearvars,%clc
MeasPath='/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/phantom';
metabolites=getMetaboliteStruct('phantom');
flip=false;

% addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'));
addpath(genpath('/ptmp/pvalsala/Packages/pulseq'));
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))
%%


system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
% seq.read(fullfile(sn,'../EXPDATA/pulseq/T1_nonsel_FOCI_2H_10min.seq'))
seq.read(fullfile(MeasPath,'../pulseq/','T2_nonsel_2H_10min.seq'))

st.dwell_s=seq.getDefinition('dwell');
st.TE_array=seq.getDefinition('TE_array');
st.rf_dur=seq.getDefinition('rf_dur');
st.averages=seq.getDefinition('averages');
st.repetitions=seq.getDefinition('repetitions');
%%
% dirst=dir(fullfile(sn,'*pulseq_T1_invFOCI_10mins*.dat'));
dirst=dir(fullfile(MeasPath,'*_T2*.dat'));
st.filename=dirst(end).name; %load last file
twix=mapVBVD(fullfile(MeasPath,st.filename));

% read noise data and calc Noise decorrelation
dirst_noise=dir(fullfile(MeasPath,'*oise*.dat'));
twix_noise=mapVBVD(fullfile(dirst_noise(1).folder,dirst_noise(1).name));
if(iscell(twix_noise)), twix_noise=twix_noise{end}; end
[D_noise,D_image]=CalcNoiseDecorrMat(twix_noise);

if(iscell(twix)), twix=twix{end}; end
data=twix.image{''};
data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions);
data=padarray(data,[2*size(data,1),0,0,0],0,'post');
Spectrum=squeeze(mean(fftshift(fft(data,[],1),1),[2 3]));
st.Cfreq=twix.hdr.Dicom.lFrequency; %Hz
st.hdr=twix.hdr;
st.RefVoltage= twix.hdr.Spice.TransmitterReferenceAmplitude;
st.VectorSize=512;


%% noise decorr and coil combination

data=twix.image{''};
data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions);
data=permute(data,[2 1 4 3]);
data_whiten=reshape(D_image*data(:,:),size(data));
data_whiten=permute(data_whiten,[2 1 3 4]);

% Combine coil data with wsvd
data_whiten=mean(data_whiten,4);
DataSize=size(data_whiten);
%we already noise decorrelated data!
wsvdOption.noiseCov         =0.5*eye(DataSize(2));
data_Combined=zeros([DataSize(1) DataSize(3:end)]);
for rep=1:size(data_whiten,3)
    rawSpectra=squeeze(data_whiten(:,:,rep));  % fid x Coil
    [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(rawSpectra, [], wsvdOption);
    data_Combined(:,rep)=-1*wsvdCombination;
    wsvdQuality_all(rep)=wsvdQuality;
end
%% Test amares
st.AcqDelay_s=seq.getBlock(5).blockDuration+seq.getBlock(6).blockDuration+seq.getBlock(7).adc.delay;
filter_delay=(2.9023e-6+st.dwell_s*1.4486);
st.AcqDelay_s=0;%(st.AcqDelay_s+filter_delay);

%    metabolites(2).freq_shift_Hz=-68;
[expParams,pk]=getAMARES_structs_T1inv(twix,metabolites,st);

spec1=specFft(padarray(data_Combined,[512*4 0],'post'));
% n=n+1;
%  faxis=linspace(-0.5/st.dwell_s,(0.5/st.dwell_s)-(1/(st.dwell_s*size(spec1,1))),size(spec1,1));
%   faxis=linspace(-0.5/st.dwell_s,(0.5/st.dwell_s)-1*(1/(st.dwell_s*size(spec1,1))),size(spec1,1));
faxis=calcFreqAxis(st.dwell_s,size(spec1,1));

%fix zeroth order phase offset so that amplitude doesn't flip

[phi0]=CalcZerothPhase(faxis,specFft(padarray(data_Combined,[512*4 0],'post')),300);

for cm=1:4
    pk.bounds(cm).phase=50*[-1 +1]-rad2deg(median(phi0))  -180*flip;
    pk.initialValues(cm).phase=-1*rad2deg(median(phi0)) -180*flip;
end
pk.bounds(1).linewidth=[3,25];
pk.bounds(3).linewidth=[3,25];
pk.bounds(2).linewidth=[10 40];
pk.bounds(4).linewidth=[10 30];
%pick the best acq_delay
clear relNorm;
acq_delay_arr=linspace(st.AcqDelay_s,1000e-6,50);
for jj=1:length(acq_delay_arr)
    expParams.beginTime=acq_delay_arr(jj);
    [~, fitStatus,~,~] = AMARES.amaresFit(double(data_Combined(:,end)), expParams, pk,0,'quiet',true);
    relNorm(jj)=fitStatus.relativeNorm;
end
[minRelNorm,idx]=min(relNorm);
st.AcqDelay_s=acq_delay_arr(idx);
expParams.beginTime=acq_delay_arr(idx);

clear fitResults fitStatus CRBResults
for i=1:size(data_Combined,2)
    [fitResults{i}, fitStatus{i},~,CRBResults{i}] = AMARES.amaresFit(double(data_Combined(:,i)), expParams, pk,0,'quiet',true);
    amp_all(i,:)=[fitResults{i}.amplitude];
    % pause(1)
end

%  figure, plot(cell2mat(cellfun(@(x) x.linewidth,fitResults,'UniformOutput',false)'))

%repeat with restricted linewidth
all_lw=cell2mat(cellfun(@(x) x.linewidth,fitResults,'UniformOutput',false)');
med_lw=median(all_lw(1:5,:),1);
std_lw=std(all_lw(1:5,:),[],1);

for cMet=1:length(metabolites)
    pk.bounds(cMet).linewidth=[-1,1]*0.1*std_lw(cMet)+med_lw(cMet);
end
clear fitResults fitStatus CRBResults
for i=1:size(data_Combined,2)
    [fitResults{i}, fitStatus{i},~,CRBResults{i}] = AMARES.amaresFit(double(data_Combined(:,i)), expParams, pk,0,'quiet',true);
    amp_all(i,:)=[fitResults{i}.amplitude];
    % pause(1)
end

%
%  plot overview
figure(13),clf

tt=tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
nexttile(tt,1,[2 2])
hold on

taxis=0:st.dwell_s:st.dwell_s*(length(faxis)-1);
taxis=taxis+st.AcqDelay_s;

cmap_lines=jet(15);
clear all_lines;
for crep=1:15
    %build n-peak Lorentzian model function
    % model = @(cf,a,g)a./(1 + ((faxis - cf)/(g/2)).^2) ; %  Lorentzian
    model_time = @(cf,a,t2)a.*exp(-taxis/t2).*exp(2i*pi*cf*taxis); %  Lorentzian
    spec_sim=zeros(size(faxis));
    amp=[1 0.5 0.50 0.5];
    fid=zeros(size(taxis));
    for i=1:length(metabolites)
        T2star_s=1/(pi*fitResults{crep}.linewidth(i)); %s
        freq_shift_Hz=fitResults{crep}.chemShift(i)*expParams.imagingFrequency;
        fid=fid+model_time(freq_shift_Hz,fitResults{crep}.amplitude(i),T2star_s);
    end
    phi0_=deg2rad(fitResults{crep}.phase(1));
    zero_order=exp(-1i*(phi0_));
    first_order_only=exp(-1i*(2*pi*faxis*st.AcqDelay_s));

    hold on,
    all_lines(crep)=plot(faxis,real(specFft(fid(:)).*first_order_only),'LineWidth',1.2,'Color',cmap_lines(crep,:),'DisplayName',sprintf('TE= %.1f ms',st.TE_array(crep)*1e3));
    plot(faxis,real(spec1(:,crep).*zero_order.*first_order_only),'.','MarkerSize',7,'LineWidth',1.2,'Color',cmap_lines(crep,:))

end


xlim([-250 100]),xlabel('frequency [Hz]'),%ylim([-0.5 0.5])
legend(all_lines,'location','northwest')
set(gca,'ColorOrder',jet(size(Spectrum,2)));
peakName={metabolites.name };
title('Global T2 : spin-echo ')
grid on

%     xlim([-75.6206 -38.8002]),ylim([-0.0377 0.0571]) %glucose

weights=1./cell2mat(cellfun(@(x)x.amplitude,CRBResults,'UniformOutput',false)');
plt_id=[3 4, 7 8];
for i=1:length(metabolites)
    nexttile(tt,plt_id(i));


    lcolour='r'; typ='AMARES';
    if(i==1 && contains(MeasPath,'sub'))
        % Set up fittype and options.
        ft = fittype( 'f*(a*(exp(-x/b)) +d*(exp(-x/e))) +c +1e5*(a+d-1)', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares','Weights',weights(:,i) );
        opts.Display = 'Off';
        opts.Robust = 'on';
        opts.StartPoint = [1 0.1e3 0 0.1 500 1];
        opts.Lower = [0.7 1 -Inf 0 200 -Inf];
        opts.Upper = [1 50 Inf 0.25  900 Inf];

        % Fit model to data.
        [fitresult, gof] = fit( st.TE_array*1e3,col(amp_all(:,i)) , ft, opts );

        plot(st.TE_array*1e3,amp_all(:,i),[lcolour,'x'],'LineWidth',1.5)
        hold on
        plot(st.TE_array*1e3,fitresult(st.TE_array*1e3),lcolour,'LineWidth',1.5),
        lege{1,1}=sprintf('%s , T2= %.2f ms',typ,fitresult.b);
        lege{2,1}=sprintf('%.2f*exp(-t/%.0f)+%.2f*exp(-t/%.0f)',fitresult.a,fitresult.b,fitresult.d,fitresult.e);


    else
        ft = fittype( 'a*(exp(-x/b))+c', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares','Weights',weights(:,i)  );
        opts.Display = 'Off';
        opts.Robust = 'on';
        opts.StartPoint = [1 0.1e3 0];
        opts.Lower = [-Inf 0 -Inf];
        opts.Upper = [Inf 2e3 Inf];

        % Fit model to data.
        [fitresult, gof] = fit( st.TE_array*1e3,col(amp_all(:,i)) , ft, opts );

        plot(st.TE_array*1e3,amp_all(:,i),[lcolour,'x'],'LineWidth',1.5)
        hold on
        plot(st.TE_array*1e3,fitresult(st.TE_array*1e3),lcolour,'LineWidth',1.5),
        lege{1,1}=sprintf('%s , T2= %.2f ms',typ,fitresult.b);
        lege{2,1}=sprintf('%.2f*exp(-t/%.0f)',fitresult.a,fitresult.b);
    end

    xlabel('TE [ms]')
    title(peakName{i})
    legend(lege{:},'Location','best')
    grid on
    fitresult_all{i}=fitresult;
end
fontsize(gcf,"scale",1.3)

set(gcf,'color','w','Position',[249 195 1661 727])
savefig(fullfile(MeasPath,'T2_final_amares.fig'))

% make table
fprintf('\n%% %s',MeasPath)
fprintf('\nT2=[%.4f,%.4f,%.4f,%.4f]*1e-3; %%s',cellfun(@(x)x.b ,fitresult_all))
CI=cellfun(@(x)confint(x),(fitresult_all),'UniformOutput',false);
fprintf('\nT2_CI=[%.4f,%.4f,%.4f,%.4f]*1e-3; %%s diff(CI95)/2',cellfun(@(x) diff(x(:,2))/2,CI))
if(isfield(fitresult_all{1},'d'))
    fprintf('\n%%T2(%.2f%%) =%.4f and T2(%.2f%%) =%.4f',fitresult_all{1}.a*100,fitresult_all{1}.b,fitresult_all{1}.d*100,fitresult_all{1}.e)
end

T2star_ms=median(1e3./(pi*cell2mat(cellfun(@(x) x.linewidth,fitResults,'UniformOutput',false)')),1);
fprintf('\nT2_star_s=[%.4f,%.4f,%.4f,%.4f]*1e-3; %%s median \n',T2star_ms)
fprintf('Acq_delay=%.4f; %%us\n',st.AcqDelay_s*1e6)
fprintf('ImagingFreq=%.4f; %%Hz\n',st.Cfreq)

chem_ppm=median(cell2mat(cellfun(@(x) x.chemShift,fitResults,'UniformOutput',false)'),1);
fprintf('chemShift_Hz=[%.4f,%.4f,%.4f,%.4f]; %%Hz median \n',chem_ppm*st.Cfreq*1e-6)
fprintf('AMARES_rel_norm=%.4f; %%median not mean here frac\n',median(cellfun(@(x)x.relativeNorm,fitStatus)))

