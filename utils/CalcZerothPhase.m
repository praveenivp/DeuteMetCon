function [phi0,spec2]=CalcZerothPhase(faxis,spec,freq_Hz)


% try permute to match the size
faxis=faxis(:);
if(size(spec,1)~=length(faxis))
    spec=spec.';
end

phi0=zeros(size(spec,2),1);
spec2=zeros(size(spec));


%get index of te frequency points
 [~,idx_freq(1)]= min(abs(faxis+freq_Hz));
 [~,idx_freq(2)]= min(abs(faxis-freq_Hz));


CalcSlope = @(x)  (imag(x(2)) -imag(x(1)))/(real(x(2)) -real(x(1)));
% line scan range and resolution
phi0_axis=linspace(-pi,pi,1024);

for crep=1:size(spec,2)

% line scan to maximize the slope of two points  
spec_freqPts=spec(idx_freq);
for iii=1:length(phi0_axis)
     slp(iii)=CalcSlope(spec_freqPts.*exp(1i*phi0_axis(iii)));
%      slp(iii)=tan(angle(diff(spec_freqPts.*exp(1i*phi0_axis(iii)))));
end

% figure,plot(phi0_axis,slp)

 [~,idx_max_slp]= max(abs(slp));
phi0(crep)=phi0_axis(idx_max_slp);
spec2(:,crep)=spec(:,crep)*exp(1i*phi0_axis(idx_max_slp));
end
%flip them if necesssary
% phi0=phi0+(max(real(-1*spec2)) >max(real(1*spec2)))'*pi;
% spec2=spec2.*exp(1i*(max(-1*spec2) >max(1*spec2))*pi);
end