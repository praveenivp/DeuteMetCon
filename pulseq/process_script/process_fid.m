
dirst=dir('allData#S94Tuebingen#F52182#M1025#D220523#T113251#rpfid-T1.dat');
st.filename=dirst(end).name %load last file
st.dwell_s=twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9;
twix=mapVBVD(st.filename);
data=twix.image{''};
  data=reshape(data,size(data,1),size(data,2),[]);
  data=padarray(data,[4*size(data,1),0,0,0],0,'post');
Spectrum=squeeze(sos(myfft(data,1),[2 3]));
st.Cfreq=twix.hdr.Dicom.lFrequency; %Hz
st.hdr=twix.hdr;
% st.RefVoltage= twix.hdr.Spice.TransmitterReferenceAmplitude;


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


%%
% load('/ptmp/pvalsala/deuterium/DDRM-T4EU/proc/NoiseMat.mat')
data=twix.image{''};
data=reshape(data,size(data,1),size(data,2),[]);
data=mean(permute(data,[2 1 4 3]),[3 4]);
data_whiten=reshape(D_noise*data(:,:),size(data));
% data_whiten=mean(data_whiten(:,:,:,1),4);
 data_whiten=permute(data_whiten,[2 1 3 4]).*exp(1i*0);


                Acqdelay= twix.hdr.Phoenix.sSpecPara.lAcquisitionDelay*1e-6; %s
                ext_size=round(Acqdelay/st.dwell_s);
%                  [data_whiten] = fidExtrp(data_whiten,ext_size);


    %%
    DataSize=[size(data_whiten) 1];
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

figure,plot(faxis(:),abs(myfft(fids(:))))
%%
