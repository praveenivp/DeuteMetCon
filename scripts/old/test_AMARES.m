%% simulate spectrum
addpath(genpath('/ptmp/pvalsala/Packages/OXSA/main'))
addpath(genpath('/ptmp/pvalsala/Packages/OXSA/utils'))
dt=40e-6; %[s]
met_st=getMetaboliteStruct('phantom');
beginTime=1e-3;
taxis=(0:dt:100e-3)+beginTime;
sys_freq_offset=0;
faxis=linspace(-0.5/dt,0.5/dt,length(taxis));

% build n-peak Lorentzian model function
model = @(cf,a,g)a./(1 + ((faxis - cf)/(g/2)).^2) ; %  Lorentzian
model_time = @(cf,a,t2)a.*exp(-taxis/t2).*exp(2i*pi*cf*taxis); %  Lorentzian
spec_sim=zeros(size(faxis));
amp=[1 0.5 0.50 0.5];
fid=zeros(size(taxis));
for i=1:length(met_st)
    fwhm=2/(pi*met_st(i).T2star_s);
    spec_sim=spec_sim+model(met_st(i).freq_shift_Hz+sys_freq_offset,amp(i),fwhm);
    fid=fid+model_time(met_st(i).freq_shift_Hz+sys_freq_offset,amp(i),met_st(i).T2star_s);
end

figure, plot(faxis,abs(fftshift(fft(fid(:),length(faxis),1)))),xlim([-500 500])
%% amaresFit
DW=dt;
BW=1/dt; %Hz
samples=length(fid);
timeAxis=0:DW:(samples-1)*DW;
imagingFrequency=9.38*6.5e6*1e-6;
offset=0;
ppmAxis=linspace(-BW/2,BW/2-BW/samples,samples)/imagingFrequency;

freq=[met_st.freq_shift_Hz];
nMet=length(met_st);
chemShift=freq./imagingFrequency;
phase=ones(1,nMet).*0;
amplitude=ones(1,nMet);
linewidth=ones(1,nMet).*12;

amares_struct=struct('chemShift',chemShift, ...
    'phase', phase,...
    'amplitude', amplitude,...
    'linewidth',linewidth, ...
    'imagingFrequency', imagingFrequency,...
    'BW', BW,...
    'timeAxis', timeAxis(:), ...
    'dwellTime', DW,...
    'ppmAxis',ppmAxis(:), ...
    'beginTime',beginTime, ...
    'offset',offset,...
    'samples',samples, ...
    'peakName',string({met_st.name}));

pk=PriorKnowledge_DMI(amares_struct);
pk.bounds(2).linewidth=[5 30];
pk.bounds(3).linewidth=[5 20];


%             output format (4th dim) : use xFit order [chemicalshift x N] [linewidth xN] [amplitude xN] [phase x1] fitStatus.relativeNorm fitStatus.resNormSq]

                fid2= fidExtrp(fid(:),beginTime/dt);

                fid3 = AMARES.makeInitialValuesModelFid(pk, amares_struct);

[fitResults, fitStatus, figureHandle, CRBResults] = AMARES.amaresFit(double(fid(:)), amares_struct, pk, true,'quiet',false);
%                         metabol_con(i,:)= [fitStatus.xFit fitStatus.relativeNorm fitStatus.resNormSq]';
fitted_T2star=1e3./(fitResults.linewidth*pi)

real_T2star=[met_st.T2star_s]*1e3

%% AMARES with builit-in phasecorr

fid_struct=struct('signals',cell(1),'dwellTime',dt,'timeAxis',0:dt:dt*(length(fid)-1),'samples',samples, ...
        'imagingFrequency', imagingFrequency,...
        'ppmAxis',ppmAxis(:));
pk=PriorKnowledge_DMI(amares_struct);
pk.bounds(2).linewidth=[5 30];
pk.bounds(3).linewidth=[5 20];
fid_struct.signals{1}=double(fid(:));
% Results = AMARES.amares(spec, instanceNum ,voxelNum, beginTime, expOffset, pk, showPlot, varargin)
results    =AMARES.amares(fid_struct,1,1,              beginTime,0,pk,true);


results=call_amares(fid(:),fid_struct,1,1,              beginTime,0,pk,true);