%% Script to
addpath(genpath('/ptmp/pvalsala/Packages/mapVBVD'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'));
addpath(genpath('/ptmp/pvalsala/Packages/pulseq'));
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))
%
% metabolites=getMetaboliteStruct('phantom');
% MeasPath: '/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/phantom';

% metabolites=getMetaboliteStruct('invivo');
% MeasPath='/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/sub-01';

%% read pulseq sequence
system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
seq.read(fullfile(MeasPath,'../pulseq','T1_nonsel_FOCI_2H_10min.seq'))


st.TI_array=seq.getDefinition('TI_array');
st.rf_dur=seq.getDefinition('rf_dur');
st.averages=seq.getDefinition('averages');
st.repetitions=seq.getDefinition('repetitions');
%%
% MeasPath='/ptmp/pvalsala/deuterium/phantoms/20250717_relaxometryphantom';
% MeasPath='/ptmp/pvalsala/deuterium/phantoms/20240925_Newsequence/TWIX';
rmos_flag={'rmos'};
% dirst=dir(fullfile(sn,'*pulseq_T1_invFOCI_10mins*.dat'));
dirst=dir(fullfile(MeasPath,'*_T1_*.dat'));
st.filename=dirst(end).name; %load last file
twix=mapVBVD(fullfile(MeasPath,st.filename),rmos_flag{:});

% read noise data and calc Noise decorrelation
dirst_noise=dir(fullfile(MeasPath,'*oise*.dat'));
twix_noise=mapVBVD(fullfile(dirst_noise(1).folder,dirst_noise(1).name),rmos_flag{:});
if(iscell(twix_noise)), twix_noise=twix_noise{end}; end
[D_noise,D_image]=CalcNoiseDecorrMat(twix_noise);

if(iscell(twix)), twix=twix{end}; end
% data=twix.image{''};
% data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions);
% data=padarray(data,[2*size(data,1),0,0,0],0,'post');
% Spectrum=squeeze(mean(fftshift(fft(data,[],1),1),[2 3]));
st.Cfreq=twix.hdr.Dicom.lFrequency; %Hz
st.hdr=twix.hdr;
st.RefVoltage= twix.hdr.Spice.TransmitterReferenceAmplitude;
st.VectorSize=twix.image.NCol/(1+~isempty(rmos_flag));
st.dwell_s=seq.getDefinition('dwell')*(1+~isempty(rmos_flag));

%% noise decorr and coil combination

data=twix.image{''};
data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions,[]);
data=permute(data,[2 1 4 3 5]);
data_whiten=reshape(D_image*data(:,:),size(data));
data_whiten=permute(data_whiten,[2 1 3 4 5]);

% Combine coil data with wsvd
data_whiten=mean(data_whiten,[4 5]);
DataSize=size(data_whiten);
%we already noise decorrelated data!
wsvdOption.noiseCov         =0.5*eye(DataSize(2));
data_Combined=zeros([DataSize(1) DataSize(3:end)]);
for rep=1:prod(DataSize(3:end))
    rawSpectra=squeeze(data_whiten(:,:,rep));  % fid x Coil
    [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(rawSpectra, [], wsvdOption);
    % [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvdApod(rawSpectra, [], 1/20e-3,timeAxis(:),wsvdOption);
    data_Combined(:,rep)=wsvdCombination;
    wsvdQuality_all(rep)=wsvdQuality;
end

%normalize data: some data may ne acquired with low receiver gain
data_Combined=data_Combined./max(abs(specFft(data_Combined(:,end))));

% if you want to correct phase manually
spec1=specFft(padarray(data_Combined,[512 0],0,'post'));
 faxis=linspace(-0.5/(st.dwell_s),0.5/st.dwell_s,length(spec1));
InteractivePhaseCorr(faxis,spec1)
%% full 2d fitting
st.AcqDelay_s=seq.getBlock(5).blockDuration+seq.getBlock(4).rf.ringdownTime+seq.getBlock(6).adc.delay;
filter_delay=(2.9023e-6+st.dwell_s*1.4486);
st.AcqDelay_s=(st.AcqDelay_s+filter_delay);




% Optimization options
tol=1e-9;
options = optimoptions(@fmincon, ...
    'Display', 'final-detailed', ...
    'Algorithm', 'interior-point', ...
    'UseParallel', true,'StepTol',tol,'OptimalityTol',tol,'FunctionTol',tol,...
    'MaxFunction',inf,'Maxiter',2000,'Hessian','bfgs');

taxis=0:st.dwell_s:st.dwell_s*(length(data_Combined)-1);
taxis=taxis+st.AcqDelay_s;




if(contains(MeasPath,'phantom'))

    %get additional delay=100e-6 from manual phase correction!
    % inputs
    x0 = [1 0.1 0.1 0.1, ...    % Amplitudes (a_i)
        [metabolites.freq_shift_Hz],...       % Frequencies (cf_i, Hz)
        [35,17.62,28,22], ...  % T2* times (s)
        [metabolites.T1_s], ...  % T1 times (s)
        pi,...       % Zero-order phase (rad)
        -50e-6,ones(1,4)*1.8];         %additional delay [s] and fit

    % Parameter bounds: [lb] and [ub]
    lb = [0,0,0,0, ...           % a_i >= 0
        [metabolites.freq_shift_Hz]-10, ... % Frequencies
        [15e-3, 10e-3, 15e-3, 15e-3], ...  % T2* times (s)
        200e-3,40e-3,100e-3,200e-3, ...           % T1 times [s]
        pi, ...                % phi0 in [-π, π]
        x0(18),ones(1,4)*1];
    ub = [Inf,Inf,Inf,Inf, ...    % a_i upper bound
        [metabolites.freq_shift_Hz]+10, ...    % Frequencies
        [55e-3, 20e-3, 50e-3, 80e-3],... %T2* times [s]
        500e-3,100e-3,300e-3,400e-3, ...    % T1 times [s]
        2*pi, ...                 % phi0
        x0(18),ones(1,4)*2];
    %glucose peak decreases in the last TIs : dynamic range problem?
    obj_func=@(x)objective(x,st.TI_array(1:1:end-3),taxis,data_Combined(:,1:1:end-3));
else
    % inputs
    x0 = [0.3664 0.4598 0.0835 0.0703, ...    % Amplitudes (a_i)
        [metabolites.freq_shift_Hz],...       % Frequencies (cf_i, Hz)
        [metabolites.T2star_s], ...  % T2* times (s)
        [metabolites.T1_s], ...  % T1 times (s)
        1.32,...       % Zero-order phase (rad
        0,ones(1,4)*1.8];         %additional delay [s] and fit

    % Parameter bounds: [lb] and [ub]
    lb = [0,0,0,0, ...           % a_i >= 0
        [metabolites.freq_shift_Hz]-10, ... % Frequencies
        [5e-3, 5e-3, 5e-3, 3e-3], ...  % T2* times (s)
        200e-3,50e-3,100e-3,10e-3, ...           % T1 times [s]
        0, ...                % phi0 in [-π, π]
        -2e-4,ones(1,4)*1];
    ub = [Inf,Inf,Inf,Inf, ...    % a_i upper bound
        [metabolites.freq_shift_Hz]+10, ...    % Frequencies
        [25e-3, 25e-3, 25e-3, 15e-3],... %T2* times [s]
        450e-3,200e-3,300e-3,200e-3, ...    % T1 times [s]
        2*pi, ...                 % phi0
        2e-4,ones(1,4)*2];

    obj_func=@(x)objective(x,st.TI_array(1:1:end),taxis,data_Combined(:,1:1:end));
end


inp_func=@(x) simulate_fid(x,st.TI_array,taxis);
s0=inp_func(x0);
[x_opt,fval,exitflag,output,lambda,grad,hessian] = fmincon(obj_func,x0,[],[],[],[],lb,ub,[],options);
st.AcqDelay_s=st.AcqDelay_s+x_opt(18);



% calculate effcetive sample size
N=size(data_Combined,1);
acf_r=autocorr(real(data_Combined(:,1)));
acf_i=autocorr(imag(data_Combined(:,1)));
Neff_real = N / (1 + 2 * acf_r(2));
Neff_imag = N / (1 + 2 * acf_i(2));
Neff = 15* (Neff_real + Neff_imag);

P=length(x0);
% 1. Calculate Mean Squared Error
mse = fval / (Neff - P);

% 2. Calculate parameter covariance matrix
covariance_matrix = mse * pinv(hessian);

% 3. Extract the standard deviations
std_devs = sqrt(diag(covariance_matrix));
assert(any(imag(std_devs)<eps)) % fit parameter at bounds
clear fitResults_T1
fitResults_T1.T1_std=std_devs([13:16])'*1e3;
fitResults_T1.amplitude = x_opt(1:4);
fitResults_T1.freq_shift_Hz = x_opt(5:8); %Hz
fitResults_T1.T2star_ms = x_opt(9:12)*1e3;
fitResults_T1.T1_ms = x_opt(13:16)*1e3;
fitResults_T1.phi0 = (x_opt(17));
fitResults_T1.additional_delay=x_opt(18);
fitResults_T1.b=x_opt(19:22);
fitResults_T1.resnorm=fval;
fitResults_T1.MeasPath=MeasPath;
fitResults_T1.ub=ub;
fitResults_T1.lb=lb;
fitResults_T1.x0=x0;
fitResults_T1.x_opt  =x_opt;
disp(fitResults_T1)
% fitResults_T1{str2double(MeasPath(end-1:end))}=fitResults_2D;

% plot
figure(4),clf,
tt2=tiledlayout(2,6,'TileSpacing','tight','Padding','tight');

taxis=0:st.dwell_s:st.dwell_s*(length(data_Combined)-1);
taxis=taxis+st.AcqDelay_s;
taxis2=linspace(min(taxis),max(taxis)*4,length(taxis)*4);
S_opt_plot= simulate_fid(x_opt,st.TI_array,taxis2);
data_Combined2=padarray(data_Combined,[length(taxis2)-size(data_Combined,1),0],0,'post');
spec_data=specFft(data_Combined2);
faxis=linspace(-0.5/(st.dwell_s),0.5/st.dwell_s,length(taxis2));
first_order_only=exp(-1i*(2*pi*faxis*st.AcqDelay_s)).';
zero_order=exp(1i*x_opt(17));

spec_data=real(spec_data.*first_order_only.*zero_order);
spec_fit=real(specFft(S_opt_plot).*first_order_only.*zero_order);


cax=[-1 1]*0.1;
nexttile([1 2])
imagesc(st.TI_array*1e3,faxis, spec_data)
axis square,ylim([-300 100]),colorbar,cb=colorbar;title('data'),xlabel('TI [ms]'),ylabel('frequency [Hz]'),set(gca,'FontWeight','bold')
clim(cax)
nexttile([1 2])
imagesc(st.TI_array*1e3,faxis,   spec_fit)
axis square,ylim([-300 100]),colorbar,clim(cb.Limits);title('fit'),xlabel('TI [ms]'),ylabel('frequency [Hz]'),set(gca,'FontWeight','bold')
clim(cax)
nexttile([1 2])
imagesc(st.TI_array*1e3,faxis,    abs(spec_data- spec_fit))
axis square,ylim([-300 100]),colorbar,clim(cb.Limits),title('residual'),xlabel('TI [ms]'),ylabel('frequency [Hz]'),set(gca,'FontWeight','bold')
colormap('turbo')
clim(cax)
cmap_lines=jet(15);
nexttile([1 3])
for crep=1:15
    hold on,
    lines_all{crep}=plot(faxis,spec_fit(:,crep),'LineWidth',1.5,'Color',cmap_lines(crep,:),'DisplayName',sprintf('TI=%.0f ms',st.TI_array(crep)*1e3));
    plot(faxis,real(spec_data(:,crep)),'.','MarkerSize',7,'LineWidth',1.5,'Color',cmap_lines(crep,:))
end
xlim([-300 100])
title('fit with data'),xlabel('frequency [Hz]'),ylabel('amplitude [a.u]')
l1=plot(nan,nan,'-black','LineWidth',2,'DisplayName','fit');
l2=plot(nan,nan,'.black','LineWidth',2,'DisplayName','data');
legend([l1,l2,lines_all{:}],'Location','northwest','numColumns',2),set(gca,'FontWeight','bold'),box on,grid minor

nexttile([1 3])
for crep=1:15
    hold on,
    lines_all{crep}=plot(faxis,real(specFft(S_opt_plot(:,crep)).*zero_order.*first_order_only),'LineWidth',1.5,'Color',cmap_lines(crep,:),'DisplayName',sprintf('TI=%.0f ms',st.TI_array(crep)*1e3));
    plot(faxis,-0.2+real(specFft(S_opt_plot(:,crep)).*zero_order.*first_order_only)-real(spec_data(:,crep)),'--','MarkerSize',7,'LineWidth',1.5,'Color',cmap_lines(crep,:))
end
xlim([-300 100]),title('fit with residuals'),xlabel('frequency [Hz]'),ylabel('amplitude [a.u]')
l1=plot(nan,nan,'-black','LineWidth',2,'DisplayName','fit');
l2=plot(nan,nan,'--black','LineWidth',2,'DisplayName','residual');
legend([l1,l2,lines_all{:}],'Location','northwest','numColumns',2),set(gca,'FontWeight','bold')

sgtitle(['Inversion recovery T1 measurements | ', MeasPath(regexp(MeasPath,'[^/]*$'):end)],'fontweight','bold','fontsize',16)
fontsize(gcf,'scale',1.1),box on,grid minor
set(gcf,'Position',[296 348 1200 730],'color','w')


% debug plots
if(false)

    crep=[1];
    figure(2),clf,plot(faxis,spec_data(:,crep)-spec_fit(:,crep),'-+'),xlim([-300 100]),hold on,plot(faxis,spec_data(:,crep),'-*'),plot(faxis,spec_fit(:,crep))

    figure(6),clf
    tt=tiledlayout(2,2,"TileSpacing","compact","Padding","compact");
    % nexttile()
    % plot(abs(specFft(s0)))
    faxis=linspace(-0.5/(st.dwell_s),0.5/st.dwell_s,length(data_Combined));
    nexttile()
    plot(faxis,real(specFft(inp_func(x_opt))))

    nexttile()
    plot(faxis,real(specFft(data_Combined)))

    nexttile(),
    s_opt=inp_func(x_opt);
    plot(real(data_Combined(:,1)))
    hold on
    plot(real(s_opt(:,1)))
    title('fid')
    nexttile()
    s_opt=inp_func(x_opt);

    taxis=0:st.dwell_s:st.dwell_s*(length(data_Combined)-1);
    first_order_only=exp(-1i*(2*pi*faxis*st.AcqDelay_s)).';
    zero_order=exp(1i*x_opt(17));

    for plotRep=[1,5,10,15]
        l=plot(Hz2ppm(faxis),real(zero_order.* first_order_only.*(specFft(data_Combined(:,plotRep)))),'-*');
        hold on
        plot(Hz2ppm(faxis),real(zero_order.* first_order_only.*specFft(s_opt(:,plotRep))),'LineWidth',2,'Color',l.Color)
        % plot(abs(specFft(s0(:,1))),'--b')
    end
    xlim([0 6]),%ylim([-0.04 0.14])
end
%% FID simulation with phase correction
function s = simulate_fid(params, TI, taxis)
% Unpack parameters
a = params(1:4);
freq_shift_Hz = params(5:8); %Hz
T2star_s = params(9:12);
T1_s = params(13:16);
phi0 = params(17);
taxis=taxis+params(18);

b=params(19:22);
% taxis=taxis+beginTime;
% Initialize FID matrix [N_TE x N_time]
s = zeros(length(taxis),length(TI));
model_time = @(cf,a,t2)a.*exp(-taxis/t2).*exp(2i*pi*cf*taxis); %  Lorentzian

for j=1:length(TI)
    fid=zeros(size(taxis));
    for i=1:length(a)
        fid=fid+model_time(freq_shift_Hz(i),a(i).*(1-b(i)*exp(-TI(j)/T1_s(i))),T2star_s(i));
    end
    % fid=fid+model_time(freq_shift_Hz(1),0.2*a(1).*(1-b*exp(-TI(j)/(1.75*T1_s(1)))),2.5*T2star_s(i));
    s(:,j)=fid(:);
end

% s(1,:)=s(1,:)*2; %Do we need this for anti-alias filter?

% Apply zero-order phase correction
s = s * exp(-1i * phi0);

end

%% Objective function to minimize
function obj = objective(params, TE, taxis, FID_data)
model_fid = simulate_fid(params, TE, taxis);
residual = model_fid - FID_data;

obj = sum(abs(residual(:)).^2);  % Sum of squared residuals
end
