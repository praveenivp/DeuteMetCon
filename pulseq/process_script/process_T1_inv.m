%% load sequence file
system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
seq.read('T1_sel_FOCI_2H_5min.seq')

st.dwell_s=seq.getDefinition('dwell');
st.TI_array=seq.getDefinition('TI_array');
st.rf_dur=seq.getDefinition('rf_dur');
st.averages=seq.getDefinition('averages');
st.repetitions=seq.getDefinition('repetitions');
%%
dirst=dir('meas_MID00325_FID04998_pulseq_T1_5min0phase.dat');
st.filename=dirst(end).name; %load last file
twix=mapVBVD(st.filename);
twix=twix{end};
data=twix.image{''};
  data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions);
  data=padarray(data,[2*size(data,1),0,0,0],0,'post');
Spectrum=squeeze(mean(fftshift(fft(data,[],1),1),[2 3]));
st.Cfreq=twix.hdr.Dicom.lFrequency; %Hz
st.hdr=twix.hdr;
st.RefVoltage= twix.hdr.Spice.TransmitterReferenceAmplitude;


faxis=linspace(-0.5/st.dwell_s,0.5/st.dwell_s,length(Spectrum));


figure
    subplot(121),plot(faxis,Spectrum),xlabel('Freq [Hz]'),title(st.filename,'Interpreter','none');
    xlim([-200 200])
    subplot(122)
%     plot(st.TR_array*1e3,max(Spectrum,[],1)),xlabel('TR [ms]')
% subplot(122),plot(rad2deg(st.fa_array),max(Spectrum,[],1))
%
lcolour='r'; typ='max';
% Set up fittype and options.
        ft = fittype( 'a*(1-exp(-x/b))+c', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Robust = 'LAR';
        opts.StartPoint = [1 0.4 0];
        opts.Lower = [-Inf 0 -Inf];
        opts.Upper = [Inf 2 Inf];
        % Fit model to data.
        [fitresult, gof] = fit( st.TI_array,col(max(Spectrum,[],1)) , ft, opts );
    
        plot(st.TI_array,max(Spectrum,[],1),[lcolour,'.'])
        hold on
        plot(st.TI_array,fitresult(st.TI_array),lcolour),
        lege{1,1}=sprintf('%s , T1= %.2f ms',typ,fitresult.b);
        lege{2,1}=sprintf('%.1d*(1-exp(-x/%.2f))+%.1d',fitresult.a,fitresult.b,fitresult.c);
        disp(fitresult)

           xlabel('TR(ms)')
    title('D20')
    legend(lege{:},'Location','southeast')

   %% AMARES

data=twix.image{''};
data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions);


% noise decorr
% addpath(genpath('C:\Users\pvalsala\Documents\Packages2\DeuteMetCon'));
% twix_noise=mapVBVD('meas_MID00059_FID03976_pvrh_Noise.dat');
% twix_noise=twix_noise{end};
% [D_noise,D_image]=CalcNoiseDecorrMat(twix_noise);

data=permute(data,[2 1 4 3]);
data_whiten=reshape(D_noise*data(:,:),size(data));
% data_whiten=mean(data_whiten(:,:,:,1),4);
data_whiten=permute(data_whiten,[2 1 3 4]).*exp(1i*-1.1);

                Acqdelay=st.rf_dur/2+120e-6; %s
                ext_size=round(Acqdelay/st.dwell_s);
%                  [data_whiten] = fidExtrp(data_whiten,ext_size);


%% Combine coil data
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
% Combine later?
%  CSI = sum(bsxfun (@times,CSI_Data_Filtered,permute(CoilWeights,[1 3 2])),3);

fids=data_Combined;
faxis=linspace(-0.5/st.dwell_s,0.5/st.dwell_s,length(fids));
%%

% do first order phase corr

%  Acqdelay=mcobj.twix.hdr.MeasYaps.alTE{1}*1e-6; %s
                Acqdelay=st.rf_dur/2+200e-6; %s
                ext_size=round(Acqdelay/st.dwell_s);

                [data_Combined] = fidExtrp(data_Combined,ext_size);
               phi0=1; %obj.flags.phaseoffset(1);
                data_Combined=data_Combined*exp(-1i*phi0);
Spectrum=squeeze(fftshift(fft(data_Combined,[],1),1));
%manual

% Spectrum=squeeze(fftshift(fft(data_Combined,[],1),1));
% faxis=linspace(-0.5/st.dwell_s,0.5/st.dwell_s,length(Spectrum));
% phi0=1; %obj.flags.phaseoffset(1);
% phi1=0.004637;
% Spectrum=Spectrum.*exp(-1i*(phi0+faxis(:).*phi1));


st.Cfreq=twix.hdr.Dicom.lFrequency; %Hz

figure(3),clf,plot(faxis,real(Spectrum),'LineWidth',1.5)
xlim([-300 300]),xlabel('frequency [Hz]')
legend(strsplit(sprintf('TE= %.1f ms\n',st.TI_array*1e3),'\n'))
set(gca,'ColorOrder',linspecer(size(Spectrum,2)));
fontsize(gcf,0.5,'centimeters')
fids=ifft(ifftshift(Spectrum,1),[],1);




%%


DW= st.dwell_s;
BW=1/DW; %Hz
samples=size(fids,1);
timeAxis=0:DW:(samples-1)*DW;
imagingFrequency=twix.hdr.Dicom.lFrequency*1e-6;
offset=0;
ppmAxis=linspace(-BW/2,BW/2-BW/samples,samples)/imagingFrequency;
nMet=4;
freq=[1 -70 -157 -215 ];
peakName={'water','Glu','Glx','Lactate'};
chemShift=freq./imagingFrequency;
phase=ones(nMet,1).*0;
amplitude=[4 1 1 1];
linewidth=[8 8  8 8];

amares_struct=struct('chemShift',chemShift, ...
    'phase', phase,...
    'amplitude', amplitude,...
    'linewidth',linewidth, ...
    'imagingFrequency', imagingFrequency,...
    'BW', BW,...
    'timeAxis', timeAxis(:), ...
    'dwellTime', DW,...
    'ppmAxis',ppmAxis(:), ...
    'beginTime',0, ...
    'offset',offset,...
    'samples',samples, ...
    'peakName',{peakName});

pk=PriorKnowledge_DMI(amares_struct);


metabol_con=zeros(size(fids,2),nMet*3+3);
%             output format (4th dim) : use xFit order [chemicalshift x N] [linewidth xN] [amplitude xN] [phase x1] fitStatus.relativeNorm fitStatus.resNormSq]
amp_all=[];ph_all=[];
for i=1:size(fids,2)
    %                 if(mod(i,1000)==0), fprintf('%.0f %d done\n',i/size(fids,1)*100); drawnow(); end
    %                     waitbar(i/size(fids,1),wbhandle,'Performing AMARES fitting');
    [fitResults, fitStatus, figureHandle, CRBResults] = AMARES.amaresFit((double(fids(:,i))), amares_struct, pk, 1);%,'quiet',true);
    metabol_con(i,:)= [fitStatus.xFit fitStatus.relativeNorm fitStatus.resNormSq]';
    amp_all=[amp_all; fitResults.amplitude];
    ph_all=[ph_all; fitResults.phase(1)];
     pause(1);
end


figure,plot(st.TI_array,amp_all(:,2:4))

%%

figure,


subplot(2,4,[1 2 5 6])
,plot(faxis,real(Spectrum),'LineWidth',1.5)
xlim([-300 300]),xlabel('frequency [Hz]')
legend(strsplit(sprintf('TE= %.1f ms\n',st.TI_array*1e3),'\n'))
set(gca,'ColorOrder',linspecer(size(Spectrum,2)));


plt_id=[3 4, 7 8];
for i=1:nMet
    subplot(2,4,plt_id(i));


    lcolour='r'; typ='AMARES';
    % Set up fittype and options.
    ft = fittype( 'a*(1-exp(-x/b))+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [1 0.1e3 0];
    opts.Lower = [-Inf 0 -Inf];
    opts.Upper = [Inf 2e3 Inf];
    % Fit model to data.
    [fitresult, gof] = fit( st.TI_array*1e3,col(amp_all(:,i)) , ft, opts );

    plot(st.TI_array*1e3,amp_all(:,i),[lcolour,'.'])
    hold on
    plot(st.TI_array*1e3,fitresult(st.TI_array*1e3),lcolour),
    lege{1,1}=sprintf('%s , T2= %.2f ms',typ,fitresult.b);
    lege{2,1}=sprintf('%.1d*(1-exp(-x/%.2f))+%.1d',fitresult.a,fitresult.b,fitresult.c);
    disp(fitresult)

    xlabel('TE [ms]')
    title(peakName{i})
    legend(lege{:},'Location','northeast')


end
fontsize(gcf,"scale",1.4)

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


