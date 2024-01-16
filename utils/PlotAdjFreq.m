function [Spectrum,faxis,header]=PlotAdjFreq(filename)
% [Spectrum,faxis,header]=PlotAdjFreq(filename)
% usage:
%PlotAdjFreq;
%PlotAdjFreq(TwixFolder);
% [Spectrum,faxis,header]=PlotAdjFreq(twixfile)

if(nargin==0)
    [fn,pn]=uigetfile('*Adj*Fre*.dat','select Adj Freq');
    filename=fullfile(pn,fn);
elseif(isfolder(filename))
    [fn,pn]=uigetfile(fullfile(filename,'*AdjFre*.dat'),'select Adj Freq');
    filename=fullfile(pn,fn);
end
% GammaH2=6.536e6 ; %Hz/T
sp_t=mapVBVD(filename);
sp_t=sp_t{end};
data=mean(sp_t.image{''},[ 3 4 5 6 7]);

data=padarray(data,[size(data,1)*2 0 0],0,'post'); %zeropadding

Spectrum=sos(fftshift(fft(data(:,:,1),[],1),1),2);

header.Centerfreq=sp_t.hdr.Dicom.lFrequency; %Hz
header.dwell=sp_t.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9; %s
header.Hz_pxl=round(1/(sp_t.hdr.Phoenix.sKSpace.lBaseResolution*header.dwell));
faxis=linspace(-0.5/header.dwell,0.5/header.dwell,length(Spectrum));
% figure,plot(faxis/(sp_t.hdr.Phoenix.sProtConsistencyInfo.flNominalB0*GammaH2-Cfreq),spectrum)%ppm
% ,plot(faxis,spectrum,'LineWidth',1.5),xlabel('Frequency(Hz)'),title('sos Spectrum');

header.FA=sp_t.hdr.Phoenix.adFlipAngleDegree{1};
try
header.pulseVoltage=sp_t.hdr.Phoenix.sTXSPEC.aRFPULSE{1}.flAmplitude;
end
%only work for rect pulse
%     header.ref_vol=sp_t.hdr.Spice.TransmitterReferenceAmplitude;
%     header.pulse_dur=(ref_vol/pulseVoltage)*(0.5*FA_All/90);
header.TR=sp_t.hdr.Phoenix.alTR{1}*1e-3; %ms
header.TE=sp_t.hdr.Phoenix.alTE{1}*1e-3; %ms
figure,plot(faxis,Spectrum),xlabel(sprintf('Frequency(Hz)+%d Hz',sp_t.hdr.Dicom.lFrequency)),title(sprintf('sos Spectrum M%d',sp_t.hdr.Config.MeasUID))
end

function res = sos(x ,dim, pnorm)
% res = sos(x [,dim, pnorm])
%
% function computes the square root of sum of squares along dimension dim.
% If dim is not specified, it computes it along the last dimension.
%
% (c) Michael Lustig 2009

if nargin < 2
    dim = size(size(x),2);
end

if nargin < 3
    pnorm = 2;
end


res = (sum(abs(x.^pnorm),dim)).^(1/pnorm);
end