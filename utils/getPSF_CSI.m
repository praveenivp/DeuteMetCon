function [fw64,PSF,W]=getPSF_CSI(twix)
% [fw64,PSF,W]=getPSF_CSI(twix)
% function to plot the point spread function of CSI experiment.
% FW64 not FWHM because FW64(unweighted experiment)/nominal resolution ~=1


fw64=zeros(4,1);

%
sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};
kp=twix.hdr.Phoenix.sKSpace;

try
FOV=[sa.dReadoutFOV sa.dPhaseFOV*kp.dPhaseResolution sa.dThickness + sa.dThickness*kp.dSliceOversamplingForDialog]; %mm
catch
   FOV=[sa.dReadoutFOV sa.dPhaseFOV*kp.dPhaseResolution  sa.dThickness]; %mm
end

MatSz=[kp.lBaseResolution  kp.lPhaseEncodingLines kp.lPartitions ];
res=FOV./MatSz; %mm

    
isCSI= contains(twix.hdr.Phoenix.tSequenceFileName,'csi');
if(isCSI)

% twix=mapVBVD('X:\mrdata\echtdata\studies\48\experiments\KSRI-QYZ6\TWIX\allData#S94Tuebingen#F37802#M376#D170123#T103746#rpcsi_fid.dat');
% MatSz=[twix.image.NLin, twix.image.NPar, twix.image.NSeg];
LinIdx=twix.image.Lin-1;
ParIdx=twix.image.Par-1;
SegIdx=twix.image.Seg-1;

W=zeros(MatSz);
for i=1:length(LinIdx)
W(LinIdx(i),ParIdx(i),SegIdx(i))=W(LinIdx(i),ParIdx(i),SegIdx(i))+1;
end

% Spectral PSF
sp=twix.hdr.Phoenix.sSpecPara;
DW=twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9; %s  % readout oversampling x2?
B0=0; %Hz
T2star=20e-3; %s
fprintf('Using T2* =20 ms \n');

Tacq= DW*sp.lVectorSize;%readoutime in s
taxis=0:DW:Tacq-DW;
BW=0.5/DW;

sig=ones(size(taxis)).*exp(-taxis/T2star).*exp(1i*2*pi*B0*taxis);
Spec=ifftshift(ifft(sig,2^17));
spec_axis=linspace(-BW/2,BW/2-BW/length(Spec),length(Spec));

else % multi-echo
    % twix=mapVBVD('X:\mrdata\echtdata\studies\48\experiments\KSRI-QYZ6\TWIX\allData#S94Tuebingen#F37802#M376#D170123#T103746#rpcsi_fid.dat');
% MatSz=[twix.image.NLin, twix.image.NPar, twix.image.NSeg];
LinIdx=twix.image.Lin;
ParIdx=twix.image.Par;


W=zeros(MatSz);
for i=1:length(LinIdx)
W(:,LinIdx(i),ParIdx(i))=W(:,LinIdx(i),ParIdx(i))+1;
end

end
 W=ones(size(W));
PSF=fftshift(fftn(W,ones(1,ndims(W)).*2^9));



%
kaxis_e1=linspace(-1,1,MatSz(1));
kaxis_e2=linspace(-1,1,MatSz(2));
kaxis_e3=linspace(-1,1,MatSz(3));
kCenter=ceil(MatSz./2);

FOV_e1=linspace(-FOV(1)/2,FOV(1)/2-FOV(1)/size(PSF,1),size(PSF,1));
FOV_e2=linspace(-FOV(2)/2,FOV(2)/2-FOV(2)/size(PSF,2),size(PSF,2));
FOV_e3=linspace(-FOV(3)/2,FOV(3)/2-FOV(3)/size(PSF,3),size(PSF,3));
PSFc=round(size(PSF)/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(14),clf
subplot(2,4,1),
plot(kaxis_e1,W(:,kCenter(2),kCenter(3)))
xlabel('normalized kspace'),title('Acquisition weighting: E1')

subplot(2,4,5),
plot(FOV_e1,abs(PSF(:,PSFc(2),PSFc(3))))
[fw64(1),points]=getFW64(FOV_e1,abs(PSF(:,PSFc(2),PSFc(3))));
hold on
xlabel('distance (mm)')
plot(points.x,points.y,'LineWidth',1.4)
title(sprintf('PSF FW64%%=%2.2f mm | Nominal=%1.2f mm',fw64(1),res(1)))



subplot(2,4,2),
plot(kaxis_e2,W(kCenter(1),:,kCenter(3)))
xlabel('normalized kspace'),title('Acquisition weighting: E2')

subplot(2,4,6),
plot(FOV_e2,abs(PSF(PSFc(1),:,PSFc(3))))
[fw64(2),points]=getFW64(FOV_e2,abs(PSF(PSFc(1),:,PSFc(3))));
hold on
xlabel('distance (mm)')
plot(points.x,points.y,'LineWidth',1.4)
title(sprintf('PSF FW64%%=%2.2f mm | Nominal=%1.2f mm',fw64(2),res(2)))


subplot(2,4,3),
plot(kaxis_e3,squeeze(W(kCenter(1),kCenter(2),:)))
xlabel('normalized kspace'),title('Acquisition weighting: E3')

subplot(2,4,7),
plot(FOV_e3,squeeze(abs(PSF(PSFc(1),PSFc(2),:))))
[fw64(3),points]=getFW64(FOV_e3,abs(PSF(PSFc(1),PSFc(2),:)));
hold on
xlabel('distance (mm)')
plot(points.x,points.y,'LineWidth',1.4)
title(sprintf('PSF FW64%%=%2.2f mm | Nominal=%1.2f mm',fw64(3),res(3)))

if(isCSI)
subplot(2,4,4)
plot(taxis*1e3,sig)
xlabel('time (ms)')
title(sprintf('T2*= %1.2f ms | B0 = %d Hz',T2star*1e3,B0))

subplot(2,4,8),plot(spec_axis,abs(Spec))
hold on
[fw64(4),points]=getFW64(spec_axis,abs(Spec));
xlabel('frequency (Hz)')
plot(points.x,points.y,'LineWidth',1.4)
xlim([-1, 1].*200)
title(sprintf('PSF FW64%%=%2.2f Hz | Nominal=%1.2f Hz',fw64(4),(1/(2*Tacq))))
end

end


%% 
function [fw64,points]=getFW64(x,y)
x=x(:);
y=y(:);
halfMax = (max(y)-min(y)) *(0.64);
index1 = find(y >= halfMax, 1, 'first');
index2 = find(y >= halfMax, 1, 'last');
fw64 = x(index2) - x(index1);

points.x=[x(index1); x(index2)];
points.y=[y(index1);y(index2)];
end