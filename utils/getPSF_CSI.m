function [fw64,PSF_struct,W]=getPSF_CSI(twix,plotFlag)
% [fw64,PSF_struct,W]=getPSF_CSI(twix,plotFlag)
% function to plot the point spread function of CSI experiment.
% FW64 not FWHM because FW64(unweighted experiment)/nominal resolution ~=1

if(~exist('plotFlag','var'))
    plotFlag=false;
end

fw64=zeros(4,1);
nominal_Resolution=zeros(4,1);

%
sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};
kp=twix.hdr.Phoenix.sKSpace;

try
    FOV=[sa.dReadoutFOV sa.dPhaseFOV*kp.dPhaseResolution sa.dThickness + sa.dThickness*kp.dSliceOversamplingForDialog]; %mm
catch
    FOV=[sa.dReadoutFOV sa.dPhaseFOV*kp.dPhaseResolution  sa.dThickness]; %mm
end

MatSz=[kp.lBaseResolution  kp.lPhaseEncodingLines kp.lPartitions ];
nominal_Resolution(1:3)=FOV./MatSz; %mm


isCSI= contains(twix.hdr.Phoenix.tSequenceFileName,'csi');
if(isCSI)

    % twix=mapVBVD('X:\mrdata\echtdata\studies\48\experiments\KSRI-QYZ6\TWIX\allData#S94Tuebingen#F37802#M376#D170123#T103746#rpcsi_fid.dat');
    % MatSz=[twix.image.NLin, twix.image.NPar, twix.image.NSeg];
    LinIdx=twix.image.Lin;
    ParIdx=twix.image.Par;
    SegIdx=twix.image.Seg;

    W=zeros(MatSz);
    for i=1:length(LinIdx)
        W(LinIdx(i),SegIdx(i),ParIdx(i))=W(LinIdx(i),SegIdx(i),ParIdx(i))+1;
    end

    % Spectral PSF (spectral blurring along CSI time axis)
    sp=twix.hdr.Phoenix.sSpecPara;
    DW=2*twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9; %s  % readout oversampling x2?
    B0=0; %Hz
    T2star=20e-3; %s
    % fprintf('Using T2* =20 ms \n');

    Tacq= DW*sp.lVectorSize;%readoutime in s
    taxis=0:DW:Tacq-DW;
    BW=1/DW;

    SigDecay=ones(size(taxis)).*exp(-taxis/T2star).*exp(1i*2*pi*B0*taxis);
    PSF_time=real(fftshift(ifft(SigDecay)));
    PSF_time=PSF_time./max(PSF_time(:));
    faxis_orig=linspace(-BW/2,BW/2-BW/length(PSF_time),length(PSF_time));
    PSF_struct.faxis=linspace(-BW/2,BW/2-BW/length(PSF_time),20*length(PSF_time));%interp polation axis
    PSF_struct.PSF_time=interp1(faxis_orig,PSF_time,PSF_struct.faxis);
    PSF_struct.taxis=taxis;
    PSF_struct.SigDecay=SigDecay;
    PSF_struct.T2star=T2star;
    [fw64(4)] = calculate_FW64( PSF_struct.PSF_time, PSF_struct.faxis,0.5);
    nominal_Resolution(4)=1/(Tacq);

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




interp_factor = 10; % Interpolation factor

% do some zeropadding before
W2=padarray(W,[1,1,1]*256,0,'both');

% Calculate PSF
psf = real(fftshift(ifftn(ifftshift(W2)))); % 3D PSF
psf = psf / max(psf(:)); % Normalize to maximum value

% Get center indices
[Nx, Ny, Nz] = size(psf);
cx = floor(Nx/2)+ 1;
cy = floor(Ny/2)+ 1;
cz = floor(Nz/2)+ 1;
%alternative way
% [~, max_idx] = max(psf(:));
% [cx, cy, cz] = ind2sub(size(psf), max_idx);

% Create high-resolution axes for interpolation
FOV_E1 = linspace(-FOV(1)/2, FOV(1)/2-FOV(1)/Nx, Nx*interp_factor);
FOV_E2 = linspace(-FOV(2)/2, FOV(2)/2-FOV(2)/Ny, Ny*interp_factor);
FOV_E3 = linspace(-FOV(3)/2, FOV(3)/2-FOV(3)/Nz, Nz*interp_factor);

% Extract and interpolate profiles along each axis
x_profile = squeeze(psf(:, cy, cz));
y_profile = squeeze(psf(cx, :, cz));
z_profile = squeeze(psf(cx, cy, :));

PSF_E1 = interp1(linspace(-FOV(1)/2, FOV(1)/2-FOV(1)/Nx, Nx), x_profile, FOV_E1, 'linear');
PSF_E2 = interp1(linspace(-FOV(2)/2, FOV(2)/2-FOV(2)/Ny, Ny), y_profile, FOV_E2, 'linear');
PSF_E3 = interp1(linspace(-FOV(3)/2, FOV(3)/2-FOV(3)/Nz, Nz), z_profile, FOV_E3, 'linear');

% Calculate FW64% using high-resolution profiles
[fw64(1), x1, x2] = calculate_FW64(PSF_E1, FOV_E1);
[fw64(2), y1, y2] = calculate_FW64(PSF_E2, FOV_E2);
[fw64(3), z1, z2] = calculate_FW64(PSF_E3, FOV_E3);


kaxis_e1=linspace(-1,1,MatSz(1));
kaxis_e2=linspace(-1,1,MatSz(2));
kaxis_e3=linspace(-1,1,MatSz(3));
kCenter=round(MatSz./2);

W_E1=squeeze(W(:,kCenter(2),kCenter(3)));
W_E2=squeeze(W(kCenter(1),:,kCenter(3)));
W_E3=squeeze(W(kCenter(1),kCenter(2),:));

PSF_struct.W_E1=W_E1;
PSF_struct.W_E2=W_E2;
PSF_struct.W_E3=W_E3;
PSF_struct.kaxis_E1=kaxis_e1;
PSF_struct.kaxis_E2=kaxis_e2;
PSF_struct.kaxis_E3=kaxis_e3;

PSF_struct.PSF_E1=PSF_E1;
PSF_struct.PSF_E2=PSF_E2;
PSF_struct.PSF_E3=PSF_E3;
PSF_struct.FOV_E1=FOV_E1;
PSF_struct.FOV_E2=FOV_E2;
PSF_struct.FOV_E3=FOV_E3;
PSF_struct.fw64=fw64;
PSF_struct.FOV_mm=FOV;
PSF_struct.MatSize=MatSz;
PSF_struct.nominal_Resolution=nominal_Resolution;

if(plotFlag>0 ||nargout==0)
    plotPSF_(PSF_struct);
end

end

function fh=plotPSF_(PSF)
% Plot results
fh=figure('Position', [100, 100, 1200, 800]);


[~, x1, x2] = calculate_FW64(PSF.PSF_E1, PSF.FOV_E1);
[~, y1, y2] = calculate_FW64(PSF.PSF_E2, PSF.FOV_E2);
[~, z1, z2] = calculate_FW64(PSF.PSF_E3, PSF.FOV_E3);

for i=1:3
subplot(2,4,i);
switch (i)
    case 1
        plot(PSF.kaxis_E1,PSF.W_E1,'LineWidth', 2);
    case 2
        plot(PSF.kaxis_E2,PSF.W_E2,'LineWidth', 2);
    case 3
        plot(PSF.kaxis_E3,PSF.W_E3,'LineWidth', 2);
end

xlabel('time [ms]')
xlabel('normalized kspace');
ylabel('weight');
title(sprintf('weighting:E%d\n Nom. Res. = %.2f mm',i,PSF.nominal_Resolution(i)))
grid on;
axis tight;
end

% X-axis PSF
subplot(2,4,5);
plot(PSF.FOV_E1, PSF.PSF_E1, 'LineWidth', 2);
hold on;
plot([x1 x2], [0.64 0.64], 'r-', 'LineWidth', 1.5);
plot([x1 x1], [0 0.64], 'r:');
plot([x2 x2], [0 0.64], 'r:');
title(sprintf('X-axis PSF\nFW64%% = %.2f mm', PSF.fw64(1)));
xlabel('Position [mm]');
ylabel('Normalized Intensity');
grid on;
axis tight;

% Y-axis PSF
subplot(2,4,6);
plot(PSF.FOV_E2, PSF.PSF_E2, 'LineWidth', 2);
hold on;
plot([y1 y2], [0.64 0.64], 'r-', 'LineWidth', 1.5);
plot([y1 y1], [0 0.64], 'r:');
plot([y2 y2], [0 0.64], 'r:');
title(sprintf('Y-axis PSF\nFW64%% = %.2f mm', PSF.fw64(2)));
xlabel('Position [mm]');
ylabel('Normalized Intensity');
grid on;
axis tight;

% Z-axis PSF
subplot(2,4,7);
plot(PSF.FOV_E3, PSF.PSF_E3, 'LineWidth', 2);
hold on;
plot([z1 z2], [0.64 0.64], 'r-', 'LineWidth', 1.5);
plot([z1 z1], [0 0.64], 'r:');
plot([z2 z2], [0 0.64], 'r:');
title(sprintf('Z-axis PSF\nFW64%% = %.2f mm', PSF.fw64(3)));
xlabel('Position [mm]');
ylabel('Normalized Intensity');
grid on;
axis tight;

%time axis
if(isfield(PSF,'PSF_time'))

subplot(2,4,4);
plot(PSF.taxis*1e3,PSF.SigDecay,'LineWidth', 2);
xlabel('time (ms)')
title(sprintf('Signal Decay(T2*= %1.0f ms)\n Nom Res= %.2f Hz',PSF.T2star*1e3,PSF.nominal_Resolution(4)))
xlabel('time [ms]');
ylabel('Normalized Intensity');
grid on;
axis tight;


subplot(2,4,8);
[~, z1, z2] = calculate_FW64(PSF.PSF_time, PSF.faxis,0.5);
plot(PSF.faxis, real(PSF.PSF_time), 'LineWidth', 2);
hold on;
plot([z1 z2], [0.5 0.5], 'r-', 'LineWidth', 1.5);
plot([z1 z1], [0 0.5], 'r:');
plot([z2 z2], [0 0.5], 'r:');
title(sprintf('Time PSF\nFWHM = %.2f Hz', PSF.fw64(4)));
xlabel('Frequency [Hz]');
ylabel('Normalized Intensity');
grid on;
axis tight;
xlim([-1, 1].*200)
end

end

% user simpler fucntion as we do interpolation already!
function [fw64, x1, x2] = calculate_FW64(profile, axis,level)
if(~exist('level','var'))
    level=0.64;
end
% Copyright Dr. Rolf Pohmann;) 
fw64=numel(find(profile>=level))*diff(axis(1:2));
x1=axis(profile>=level);
x2=x1(end);x1=x1(1);
end


% FW64% calculation function (optimized for interpolated data)
function [fw64, x1, x2] = calculate_FW64_older(profile, axis,level)
    
if(~exist('level','var'))
    level=0.64;
end

    [peak_val, peak_idx] = max(profile);
    
    % Right side
    right_profile = profile(peak_idx:end);
    right_axis = axis(peak_idx:end);
    right_idx = find(right_profile <= level*peak_val, 1, 'first');
    
    % Left side
    left_profile = profile(1:peak_idx);
    left_axis = axis(1:peak_idx);
    left_idx = find(left_profile <=level*peak_val, 1, 'last');
    
    % More accurate interpolation
    if isempty(right_idx)
        x_right = axis(end);
    else
        if right_idx > 1
            x_right = interp1(right_profile(right_idx-1:right_idx), ...
                             right_axis(right_idx-1:right_idx), level*peak_val, 'spline');
        else
            x_right = right_axis(right_idx);
        end
    end
    
    if isempty(left_idx)
        x_left = axis(1);
    else
        if left_idx < length(left_profile)
            x_left = interp1(left_profile(left_idx:left_idx+1), ...
                            left_axis(left_idx:left_idx+1), level*peak_val, 'spline');
        else
            x_left = left_axis(left_idx);
        end
    end
    
    fw64 = abs(x_right - x_left);
    x1 = x_left;
    x2 = x_right;
end
