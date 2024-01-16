function ppm_axis= Hz2ppm(faxis,B0,gamma,center_freq)
%  ppm_axis= Hz2ppm(faxis,B0,gamma,center_freq)
% example
% for deuterium at 9.4T : ppm_axis= Hz2ppm(faxis_Hz)

if(~exist("center_freq",'var'))
    center_freq=0; % Hz
end

if(~exist("B0",'var'))
    B0=9.38; %[T]
end

if(~exist("gamma",'var'))
    gamma=6.536; %[MHz/T]
end

% 4.8 is water chemical shift should be same as center_freq
ppm_axis=(faxis-center_freq)./(B0*gamma)+4.8;


end