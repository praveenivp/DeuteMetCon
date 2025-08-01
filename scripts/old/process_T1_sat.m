%% load sequence file
system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
seq.read('T1_nonsel_FOCI_2H_10min.seq')

st.dwell_s=seq.getDefinition('dwell');
st.TR_array=seq.getDefinition('TR_array');
st.rf_dur=seq.getDefinition('rf_dur');
st.averages=seq.getDefinition('averages');
st.repetitions=seq.getDefinition('repetitions');
%%
dirst=dir('*T1*.dat');
st.filename=dirst(end).name %load last file
twix=mapVBVD(st.filename);
data=twix.image{''};
  data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions);
  data=padarray(data,[4*size(data,1),0,0,0],0,'post');
Spectrum=squeeze(sos(myfft(data,1),[2 3]));
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
        [fitresult, gof] = fit( st.TR_array,col(max(Spectrum,[],1)) , ft, opts );
    
        plot(st.TR_array,max(Spectrum,[],1),[lcolour,'.'])
        hold on
        plot(st.TR_array,fitresult(st.TR_array),lcolour),
        lege{1,1}=sprintf('%s , T1= %.2f ms',typ,fitresult.b);
        lege{2,1}=sprintf('%.1d*(1-exp(-x/%.2f))+%.1d',fitresult.a,fitresult.b,fitresult.c);
        disp(fitresult)

           xlabel('TR(ms)')
    title('D20')
    legend(lege{:},'Location','southeast')


    %% wsvd combine coils
    %% Combine coil data
    data=twix.image{''};
    data=reshape(data,size(data,1),size(data,2),st.averages,st.repetitions);
    data=permute(data,[2 1 3 4 5]);

    sz     = size(data);
    data   = mcobj.D*data(:,:);
    data    = reshape(data,sz(1),sz(2),[]);
    

    %we already noise decorrelated data!
    CSI_wsvdOption.noiseCov         =0.5*eye(10);
   
 
    CSI_Combined=zeros([sz(2) prod(sz(3:end))]);
    CoilWeights = zeros([prod(sz(3:end)) sz(1)]);
    for vx=1:size(data,3)
        rawSpectra=squeeze(data(:,:,vx)).';  % Spectrum x Coil
        [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(rawSpectra, [], CSI_wsvdOption);
%         CoilWeights(:,vx) = wsvdWeights;
        CSI_Combined(:,vx)=wsvdCombination;
    end
data_combined=reshape(CSI_Combined,sz(2:end));

% data_combined=sum(data_combined,3);

    %% do proper phase correction
    DW= 200e-6; %twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9;

%     Acqdelay=0.5e-3; %s %from pulseq
%     ext_size=round(Acqdelay/DW);
% 
%     [data_corr] = fidExtrp(permute(data_combined,[2 3 1]),ext_size);
% 
% data_corr=squeeze(data_corr).';


phi0=-pi;
phi1=-0.00;
% phi0=0;
% phi1=0.01; % phase evolution per Hz

BW=1/200e-6; %Hz
timeAxis=0:DW:(length(fid)-1)*DW;
spec1=fftshift(fft(data_combined,[],1),1);
faxis=linspace(-BW/2,BW/2-BW/length(spec1),length(spec1));

%  data_corr=spec1.*exp(1i*(phi0+faxis(:).*phi1));



 %%
 %rf spoiling not added to ADC object.
  rand_phase = mod(117*((9:24).^2 + (9:24) + 2), 360)*pi/180;
 spec_corr=spec1.*exp(-1i*rand_phase);
 spec_corr=sum(spec_corr,2);

faxis1=faxis;
% faxis1(1:225)=0;
% faxis1(270:end)=0;

  spec_corr=(spec_corr.*exp(-1i*(-2.414+faxis1(:).*0.001)));
 figure(7),clf,plot(faxis,real(spec_corr(:,:)))
 xlim([-300 100])
%%
data_corr=ifft(ifftshift(spec_corr,1),[],1);
 amp=[];
 clear fitResults_all;
 for tt=1:size(data_corr,3)
 fid =sum(squeeze(data_corr(:,:,tt)),2);
%  figure,plot(abs(fid))
  %
%   fid =(sum(squeeze(data_corr(:,:,30)),2));
%               DW= twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9;
            BW=1/DW; %Hz
            timeAxis=0:DW:(length(fid)-1)*DW;
            imagingFrequency=twix.hdr.Dicom.lFrequency*1e-6; %MHz
            offset=-5;
            ppmAxis=linspace(-BW/2,BW/2-BW/length(fid),length(fid))/imagingFrequency;
            nMet=3;
            freq=[-150 -48 7];
            chemShift=1*freq./imagingFrequency;
            phase=ones(nMet,1).*0;
            amplitude=[1 1 1];
            linewidth=[12 12  12];

            samples=length(fid);

            amares_struct=struct('chemShift',chemShift, ...
                'phase', phase,...
                'amplitude', amplitude,...
                'linewidth',linewidth, ...
                'imagingFrequency', imagingFrequency,...
                'BW', BW,...
                'timeAxis', timeAxis(:), ...
                'dwellTime', DW,...
                'ppmAxis',ppmAxis(:), ...
                'beginTime',0.e-3, ...
                'offset',offset,...
                'samples',samples);

            pk=PriorKnowledge_DMI(amares_struct);




            [fitResults, fitStatus, figureHandle, CRBResults] = AMARES.amaresFit(double(fid), amares_struct, pk, 11,'quiet',true);


           amp= [amp;fitResults.amplitude];
           fitResults_all{tt}=fitResults;
 end

%%

 figure(1),clf,hold on
 for ii=1:length(fitResults_all)
ft=fitResults_all{ii};

 ft.chemShift  =ft.chemShift;
 ft.phase=ft.phase*0;
 [modelFid, modelFids] = AMARES.makeModelFid(ft,timeAxis(:), imagingFrequency);
  plot(faxis,real(fftshift(fft(modelFids,[],1),1)))
 end

 %%
 figure
 met_names={'Glx','Glc','Water'};
 for i=1:3
     subplot(2,2,i)
%  figure,plot(st.TR_array, amp)
         ft = fittype( 'a*(1-exp(-x/b))+c', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';        
        opts.Robust = 'LAR';
%         opts.StartPoint = [0.89 0.96 ];
        opts.StartPoint = [577 0.25 0];
        opts.Lower = [-Inf 0 0];
        opts.Upper = [Inf 2 0];
        % Fit model to data.
        [fitresult, gof] = fit( st.TR_array(:),amp(:,i), ft, opts );
         plot(st.TR_array,amp(:,i),[lcolour,'.'])
        hold on
        plot(st.TR_array,fitresult(st.TR_array),lcolour),
        lege{1,1}=sprintf('%s , T1= %.2f ms','AMARES',fitresult.b*1e3);
        lege{2,1}=sprintf('%.1d*(1-exp(-x/%.2f))+%.1d',fitresult.a,fitresult.b,fitresult.c);
        disp(fitresult)

           xlabel('TR [ms]')
    title(met_names{i})
    legend(lege{:},'Location','southeast')
 end