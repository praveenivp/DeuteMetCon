%% MR sampling
FOV= 180e-3; %m
res=5.6e-3;% m;
DW=114.8e-6; %s
B0=0; %Hz
T2star=20e-3; %s

kmax=0.5/res;  %1/m
Tacq= DW*(FOV/res);%readoutime in s
Gread=1; %readout gradietnt mT/m

taxis=0:DW:Tacq-DW; 
%                       T2 star realax.         B0 dispersion
sig=ones(size(taxis)).*exp(-taxis/T2star).*exp(1i*2*pi*B0*taxis);
figure(4),clf,subplot(121)
plot(taxis*1e3,sig)
xlabel('time (ms)')
title(sprintf('T2*= %1.2f ms | B0 = %d Hz',T2star*1e3,B0))


im=ifftshift(ifft(sig,2^10));
im_axis=linspace(-FOV/2,FOV/2-FOV/length(im),length(im));
subplot(122),plot(im_axis*1e3,abs(im))
hold on
[fwhm,points]=getFWHM(im_axis,abs(im));
xlabel('distance (mm)')
plot(points.x*1e3,points.y,'LineWidth',1.4)
title(sprintf('PSF FWHM=%2.2f mm | Nominal=%1.2f mm',fwhm*1e3,res*1e3))

%% do for data
%   twix=mapVBVD('X:\mrdata\echtdata\studies\48\experiments\DDRM-T4EU\TWIX\allData#S94Tuebingen#F52142#M985#D220523#T100321#rprh_trufi_5mc_monopolar.dat');
%
sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};
kp=twix.hdr.Phoenix.sKSpace;

try
FOV=[sa.dReadoutFOV sa.dPhaseFOV sa.dThickness + sa.dThickness*kp.dSliceOversamplingForDialog]; %mm
catch
   FOV=[sa.dReadoutFOV sa.dPhaseFOV sa.dThickness]; %mm
end

MatSz=[kp.lBaseResolution  kp.lPhaseEncodingLines kp.lPartitions ];
res=FOV./MatSz; %mm

FOV= FOV(1)*1e-3; %m
res=FOV/MatSz(1);% m;
DW=2*twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9; %s
B0=0; %Hz
T2star=20e-3; %s

kmax=0.5/res;  %1/m
Tacq= DW*(FOV/res);%readoutime in s
Gread=1; %readout gradietnt mT/m

taxis=0:DW:Tacq-DW; 
%                       T2 star realax.         B0 dispersion
sig=ones(size(taxis)).*exp(-taxis/T2star).*exp(1i*2*pi*B0*taxis);
figure(4),clf,subplot(121)
plot(taxis*1e3,sig)
xlabel('time (ms)')
title(sprintf('T2*= %1.2f ms | B0 = %d Hz',T2star*1e3,B0))


im=ifftshift(ifft(sig,2^16));
im_axis=linspace(-FOV/2,FOV/2-FOV/length(im),length(im));
subplot(122),plot(im_axis*1e3,abs(im))
hold on
[fwhm,points]=getFWHM(im_axis,abs(im));
xlabel('distance (mm)')
plot(points.x*1e3,points.y,'LineWidth',1.4)
title(sprintf('PSF FW64%%=%2.2f mm | Nominal=%1.2f mm',fwhm*1e3,res*1e3))
%%
samp_pattern=zeros(twix.image.NLin,twix.image.NPar);
% Lin=twix.image.Lin(twix.image.Ave==1 & twix.image.Rep==1 & twix.image.Eco==1);
% Par=twix.image.Par(twix.image.Ave==1 & twix.image.Rep==1 &twix.image.Eco==1);
for i=1:length(Lin)
 samp_pattern(Lin(i),Par(i))=1;
end

figure,imagesc(samp_pattern),xlabel('PE1'),ylabel('PE2')
title('Elliptical scanning: DDRM dataset')
%% 
function [fwhm,points]=getFWHM(x,y)
x=x(:);
y=y(:);
halfMax = (min(y) + max(y)) *0.64;
index1 = find(y >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(y >= halfMax, 1, 'last');
fwhm = x(index2) - x(index1);

points.x=[x(index1); x(index2)];
points.y=[y(index1);y(index2)];
end