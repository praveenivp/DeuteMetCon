%% Script to calculate T2 with custom non-linear constrained optimization
addpath(genpath('/ptmp/pvalsala/Packages/mapVBVD'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'));
addpath(genpath('/ptmp/pvalsala/Packages/pulseq'));
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))

% metabolites=getMetaboliteStruct('phantom');
% MeasPath='/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/phantom';
%
% metabolites=getMetaboliteStruct('invivo');
% MeasPath='/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/sub-01';
%%


system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
% seq.read(fullfile(sn,'../EXPDATA/pulseq/T1_nonsel_FOCI_2H_10min.seq'))
seq.read(fullfile(MeasPath,'../pulseq/','T2_nonsel_2H_10min.seq'))

st.dwell_s=seq.getDefinition('dwell');
st.TE_array=seq.getDefinition('TE_array')-0.46e-3;
st.rf_dur=seq.getDefinition('rf_dur');
st.averages=seq.getDefinition('averages');
st.repetitions=seq.getDefinition('repetitions');
%%
% MeasPath='/ptmp/pvalsala/deuterium/phantoms/20240925_Newsequence/TWIX';
% MeasPath='/ptmp/pvalsala/deuterium/phantoms/20250709_phantomSNR_4rep/TWIX';
% MeasPath='/ptmp/pvalsala/deuterium/phantoms/20250717_relaxometryphantom';
dirst=dir(fullfile(MeasPath,'*_T2*.dat'));
st.filename=dirst(end).name; %load last file
twix=mapVBVD(fullfile(MeasPath,st.filename));
if(iscell(twix)), twix=twix{end}; end

% read noise data and calc Noise decorrelation
dirst_noise=dir(fullfile(MeasPath,'*oise*.dat'));
twix_noise=mapVBVD(fullfile(dirst_noise(1).folder,dirst_noise(1).name));
if(iscell(twix_noise)), twix_noise=twix_noise{end}; end
[D_noise,D_image]=CalcNoiseDecorrMat(twix_noise);


st.Cfreq=twix.hdr.Dicom.lFrequency; %Hz
st.hdr=twix.hdr;
st.RefVoltage= twix.hdr.Spice.TransmitterReferenceAmplitude;
st.VectorSize=512;
st.AcqDelay_s=seq.getBlock(1).blockDuration/2;
filter_delay=(2.9023e-6+st.dwell_s*1.4486);
st.AcqDelay_s=(st.AcqDelay_s+filter_delay);

% noise decorr and coil combination
data=twix.image{''};
data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions,[]);
data=permute(data,[2 1 4 3 5]);
data_whiten=reshape(D_image*data(:,:),size(data));
data_whiten=permute(data_whiten,[2 1 3 4 5]);
timeAxis=0:st.dwell_s:st.dwell_s*(size(data_whiten,1)-1);

% Combine coil data with wsvd
data_whiten=mean(data_whiten,[4 5]);
DataSize=size(data_whiten);
%we already noise decorrelated data!
wsvdOption.noiseCov         =0.5*eye(DataSize(2));
data_Combined=zeros([DataSize(1) (DataSize(3:end))]);
for rep=1:prod(DataSize(3:end))
    rawSpectra=squeeze(data_whiten(:,:,rep));  % fid x Coil
    if(contains(MeasPath,'phantom'))
        [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(rawSpectra, [], wsvdOption);
    else
        [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvdApod(rawSpectra, [], 1/20e-3,timeAxis(:),wsvdOption);
    end
    data_Combined(:,rep)=wsvdCombination;
    wsvdQuality_all(rep)=wsvdQuality;
end

%normalize data: some data may ne acquired with low receiver gain
data_Combined=data_Combined./max(abs(specFft(data_Combined(:,1))));

%B0 correction for sub-03
if(contains(MeasPath,'sub-03'))
    %do zeroth order phase correction between rep
    spec1=specFft(padarray(data_Combined,[512 0],0,'post'));
    faxis=linspace(-0.5/(st.dwell_s),0.5/st.dwell_s,length(spec1));
    [phi0]=CalcZerothPhase(faxis,spec1,300);
    spec1=spec1.*exp(1i*phi0');

    first_order_only=exp(-1i*(2*pi*faxis*st.AcqDelay_s)).';
    zero_order=1;%exp(1i*1.7044);
    [~,idx]=max(abs(real(spec1.*zero_order.*first_order_only)));
    idx=idx-idx(1);
    idx(abs(idx)>10)=0;
    for crep=1:length(idx)
        spec1(:,crep)=ndCircShift(spec1(:,crep),-2*idx(crep),1);
    end
    data_Combined=specInvFft(spec1);
    data_Combined=data_Combined(1:512,:);
end


% if you want to correct phase manually
% spec1=specFft(padarray(data_Combined,[512 0],0,'post'));
% faxis=linspace(-0.5/(st.dwell_s),0.5/st.dwell_s,length(spec1));
% InteractivePhaseCorr(faxis,spec1)
%% full 2d fitting
st.AcqDelay_s=seq.getBlock(1).blockDuration/2;
filter_delay=(2.9023e-6+st.dwell_s*1.4486);
st.AcqDelay_s=(st.AcqDelay_s+filter_delay);

% some Initial condition/bounds are manually tuned to minimize the uncertainity!
if(contains(MeasPath,'phantom'))
    % inputs
    x0 = [0,0,0,0 ...    % Amplitudes (a_i)
        [metabolites.freq_shift_Hz]+10,...       % Frequencies (cf_i, Hz)
        30e-3, 10e-3, 25e-3, 10e-3, ...  % T2* times (s)
        0.3, 0.05, 0.1, 0.1, ...  % T2 times (s)
        pi,-150e-6,...       % Zero-order phase (rad), additional delay
        ];         %water 2
    % Parameter bounds: [lb] and [ub]
    lb = [0,0,0,0, ...           % a_i >= 0
        [metabolites.freq_shift_Hz]-10, ... % Frequencies
        [20e-3, 15e-3, 15e-3, 20e-3], ...  % T2* times (s)
        100e-3,10e-3,30e-3,100e-3, ...           % T2_i > 0
        pi,x0(18),...                % phi0 in [-π, π]
        ];
    ub = [Inf,Inf,Inf,Inf, ...    % a_i upper bound
        [metabolites.freq_shift_Hz]+10, ...    % Frequencies
        [60e-3, 25e-3, 40e-3, 80e-3],... %T2* times
        350e-3,100e-3,200e-3,300e-3, ...    % T2 times
        2*pi,x0(18), ...                 % phi0
        ];
else
    % inputs
    x0 = [0.3664 0.4598 0.0835 0.0703, ...    % Amplitudes (a_i)
        [metabolites.freq_shift_Hz],...       % Frequencies (cf_i, Hz)
        30e-3, 10e-3, 25e-3, 10e-3, ...  % T2* times (s)
        0.3, 0.05, 0.2, 0.3, ...  % T2 times (s)
        2.1933,0,...       % Zero-order phase (rad), additional delay
        400e-3,0.1];         %water 2
    % Parameter bounds: [lb] and [ub]
    lb = [0,0,0,0, ...           % a_i >= 0
        [metabolites.freq_shift_Hz]-10, ... % Frequencies
        [15e-3, 5e-3, 15e-3, 11e-3]-10e-3, ...  % T2* times (s)
        10e-3,10e-3,60e-3,10e-3, ...           % T2_i > 0
        0,-1e-3,...                % phi0 in [-π, π]
        200e-3,0];
    ub = [Inf,Inf,Inf,Inf, ...    % a_i upper bound
        [metabolites.freq_shift_Hz]+10, ...    % Frequencies
        [15e-3, 15e-3, 15e-3, 15e-3]+15e-3,... %T2* times
        50e-3,100e-3,200e-3,200e-3, ...    % T2 times
        2*pi,1e-3, ...                 % phi0
        500e-3,Inf];
end

% Optimization options
tol=1e-9;
options = optimoptions(@fmincon, ...
    'Display', 'final', ...
    'Algorithm', 'interior-point', ...
    'UseParallel', true,'StepTol',tol,'OptimalityTol',tol,'FunctionTol',tol,...
    'MaxFunction',inf,'Maxiter',200);

taxis=0:st.dwell_s:st.dwell_s*(length(data_Combined)-1);
taxis=taxis+st.AcqDelay_s;

inp_func=@(x) simulate_fid(x,st.TE_array,taxis);
s0=inp_func(x0);

obj_func=@(x)objective(x,st.TE_array(:),taxis,data_Combined(:,:));
[x_opt,fval,exitflag,output,lambda,grad,hessian] = fmincon(obj_func,x0,[],[],[],[],lb,ub,[],options);

clear fitResults_T2
fitResults_T2.amplitude = x_opt(1:4);
fitResults_T2.freq_shift_Hz = x_opt(5:8); %Hz
fitResults_T2.T2star_ms = x_opt(9:12)*1e3;
fitResults_T2.T2_ms = x_opt(13:16)*1e3;
fitResults_T2.phi0 = x_opt(17);
fitResults_T2.additional_delay=x_opt(18);
fitResults_T2.resnorm=fval;
fitResults_T2.MeasPath=MeasPath;
fitResults_T2.ub=ub;
fitResults_T2.lb=lb;
fitResults_T2.x0=x0;
fitResults_T2.x_opt  =x_opt;
st.AcqDelay_s=st.AcqDelay_s+x_opt(18);


% very suspicious uncertainity calculation!
% calculate effcetive sample size
N=size(data_Combined,1);
acf_r=autocorr(real(data_Combined(:,1)));
acf_i=autocorr(imag(data_Combined(:,1)));
Neff_real = N / (1 + 2 * acf_r(2));
Neff_imag = N / (1 + 2 * acf_i(2));
Neff = 15* (Neff_real + Neff_imag) / 2;

P=length(x0);
% 1. Calculate Mean Squared Error
mse = fval / (Neff - P);

% 2. Calculate parameter covariance matrix
covariance_matrix = mse * pinv(hessian);

% 3. Extract the standard deviations (standard errors)
std_devs = sqrt(diag(covariance_matrix));
% [std_devs([13:16 19])*1e3; std_devs(20)]
fitResults_T2.T2_std=std_devs([13:16])'*1e3;


if(~contains(MeasPath,'phantom'))
    fitResults_T2.water2_T2_ms=x_opt(19)*1e3;

    fitResults_T2.water2_fac=x_opt(20)./(x_opt(1)+x_opt(20));
    fitResults_T2.water2_T2_std_ms=std_devs(19)*1e3;
    % https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    f=x_opt(20)./(x_opt(1)+x_opt(20));A=x_opt(20);B=x_opt(1);
    fitResults_T2.water2_fac_std=abs(f/(A+B))*sqrt((B^2/A^2)*std_devs(20)^2+std_devs(1)^2);
    fitResults_T2_all2{str2double(MeasPath(end-1:end))}=fitResults_T2;

end
disp(fitResults_T2)
%  %%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(41),clf,
tt2=tiledlayout(2,6,'TileSpacing','tight','Padding','tight');

taxis=0:st.dwell_s:st.dwell_s*(length(data_Combined)-1);
taxis=taxis+st.AcqDelay_s;
taxis2=linspace(min(taxis),max(taxis)*2,length(taxis)*2);
S_opt_plot= simulate_fid(x_opt,st.TE_array,taxis2);
data_Combined2=padarray(data_Combined,[length(taxis),0],0,'post');
spec1=specFft(data_Combined2);
faxis=linspace(-0.5/(st.dwell_s),0.5/st.dwell_s,length(S_opt_plot));
first_order_only=exp(-1i*(2*pi*faxis*st.AcqDelay_s)).';
zero_order=exp(1i*x_opt(17));

cax=[-0.1 0.1+eps];
nexttile([1 2])
imagesc(st.TE_array*1e3,faxis,    real(spec1.*zero_order.*first_order_only))
axis square,ylim([-300 100]),title('data'),xlabel('TE [ms]'),ylabel('frequency [Hz]'),set(gca,'FontWeight','bold')
clim(cax)
nexttile([1 2])
imagesc(st.TE_array*1e3,faxis,    real(specFft(S_opt_plot).*zero_order.*first_order_only))
axis square,ylim([-300 100]),clim(cax);title('fit'),xlabel('TE [ms]'),ylabel('frequency [Hz]'),set(gca,'FontWeight','bold')
yticklabels([]),ylabel([])

nexttile([1 2])
imagesc(st.TE_array*1e3,faxis,    abs(real(spec1.*zero_order.*first_order_only)- real(specFft(S_opt_plot).*zero_order.*first_order_only)))
axis square,ylim([-300 100]),colorbar,clim(cax),title('residual'),xlabel('TE [ms]'),ylabel('frequency [Hz]'),set(gca,'FontWeight','bold')
colormap('turbo'), yticklabels([]),ylabel([])

cmap_lines=jet(15);
nexttile([1 3])
for crep=1:15
    hold on,
    lines_all{crep}=plot(faxis,real(specFft(S_opt_plot(:,crep)).*zero_order.*first_order_only),'LineWidth',1.5,'Color',cmap_lines(crep,:),'DisplayName',sprintf('TE= %.1f ms',st.TE_array(crep)*1e3));
    plot(faxis,real(spec1(:,crep).*zero_order.*first_order_only),'.','MarkerSize',7,'LineWidth',1.5,'Color',cmap_lines(crep,:))
end
xlim([-300 100])
title('fit with data'),xlabel('frequency [Hz]'),ylabel('amplitude [a.u]')
l1=plot(nan,nan,'-black','LineWidth',2,'DisplayName','fit');
l2=plot(nan,nan,'.black','LineWidth',2,'DisplayName','data');
legend([l1,l2,lines_all{:}],'Location','northwest','numColumns',2),set(gca,'FontWeight','bold'),box on,grid minor

nexttile([1 3])
for crep=1:15
    hold on,
    lines_all{crep}=plot(faxis,real(specFft(S_opt_plot(:,crep)).*zero_order.*first_order_only),'LineWidth',1.5,'Color',cmap_lines(crep,:),'DisplayName',sprintf('TE= %.1f ms',st.TE_array(crep)*1e3));
    plot(faxis,-0.1+real(specFft(S_opt_plot(:,crep)).*zero_order.*first_order_only)-real(spec1(:,crep).*zero_order.*first_order_only),'--','MarkerSize',7,'LineWidth',1.5,'Color',cmap_lines(crep,:))
end
xlim([-300 100]),title('fit with residuals'),xlabel('frequency [Hz]'),ylabel('amplitude [a.u]')
l1=plot(nan,nan,'-black','LineWidth',2,'DisplayName','fit');
l2=plot(nan,nan,'--black','LineWidth',2,'DisplayName','residual');
legend([l1,l2,lines_all{:}],'Location','northwest','numColumns',2),set(gca,'FontWeight','bold')

sgtitle(['Spin-echo T2 measurements | ', MeasPath(regexp(MeasPath,'[^/]*$'):end)],'fontweight','bold','fontsize',16)
fontsize(gcf,'scale',1.1),box on,grid minor
set(gcf,'Position',[296 348 1200 720],'color','w')

% debug plot
if(false)
    figure(6),clf
    tt=tiledlayout(2,2,"TileSpacing","compact","Padding","compact");
    nexttile()
    plot(real(specFft(inp_func(x_opt))))

    nexttile()
    plot(real(specFft(data_Combined)))

    nexttile(),
    s_opt=inp_func(x_opt);
    plot(real(data_Combined(:,1)))
    hold on
    plot(real(s_opt(:,1)))

    nexttile(),

    s_opt=inp_func(x_opt);

    for plotRep=[1 5 15]
        l=plot(abs((specFft(data_Combined(:,plotRep)))),'-*');
        hold on
        plot(abs(specFft(s_opt(:,plotRep))),'LineWidth',2,'Color',l.Color)
        % plot(abs(specFft(s0(:,1))),'--b')
    end
    xlim([200 300])

end
%% FID simulation with phase correction
function s = simulate_fid(params, TE, taxis)
% Unpack parameters
a = params(1:4);
freq_shift_Hz = params(5:8); %Hz
T2star_s = params(9:12);
T2_s = params(13:16);
phi0 = params(17);
taxis=taxis+params(18);



% taxis=taxis+beginTime;
% Initialize FID matrix [N_TE x N_time]
s = zeros(length(taxis),length(TE));
model_time = @(cf,a,t2)a.*exp(-taxis/t2).*exp(2i*pi*cf*taxis); %  Lorentzian

for j=1:length(TE)
    fid=zeros(size(taxis));
    for i=1:length(a)
        if(length(params)>18 && i==1)
            water2_T2=params(19);
            water2_fac=params(20);
            fid=fid+model_time(freq_shift_Hz(i),a(i).*exp(-TE(j)/T2_s(i)),T2star_s(i))+...
                model_time(freq_shift_Hz(i),water2_fac*exp(-TE(j)/water2_T2),T2star_s(i));
        else
            fid=fid+model_time(freq_shift_Hz(i),a(i).*exp(-TE(j)/T2_s(i)),T2star_s(i));

        end
    end
    s(:,j)=fid(:);
end

s(1,:)=s(1,:)*0.5;

% Apply zero-order phase correction
s = s * exp(-1i * phi0);
end

%% Objective function to minimize
function obj = objective(params, TE, taxis, FID_data)
model_fid = simulate_fid(params, TE, taxis);
residual = model_fid - FID_data;
% residual=residual.*exp(-TE'/100e-3);

% if(length(params)>18)
%     reg=(1*(params(20)/params(1)-0.1))^2; %maximize T2 values difference between water and free water
% else
% reg=0;
% end

obj = sum(abs(residual(:)).^2);%+reg;  % Sum of squared residuals
end
