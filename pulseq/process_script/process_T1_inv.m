%% load sequence file
system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
seq.read('T1_nonsel_FOCI_2H_10min.seq')

st.dwell_s=seq.getDefinition('dwell');
st.TI_array=seq.getDefinition('TI_array');
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
