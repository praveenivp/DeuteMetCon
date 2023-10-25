function bSSFP=bSSFP_profile_approx(TR,TE,flip,phi,off_resonance,lambda)
%function bSSFP=bSSFP_profile_approx(TR,TE,flip,phi,off_resonance,lambda)
%
%Calculation of the bSSFP profile (series of measurements with varying RF phase increment)
%Following Ganter, MRM, 2006: 
%Steady State of Gradient Echo Sequences with Radiofrequency Phase Cycling: Analytical Solution, Contrast Enhancement with Partial Spoiling
%Using the approximation: E1,2 ~ 1-TR/T1,2, which holds for TR << T2 < T1
%
%INPUT
%
%repetition time in ms (TR)
%echo time in ms (TE)
%flip angle in rad (flip)
%vector of RF phase increments in rad (phi)
%off-resonance-related shift in rad (off_resonance)
%relaxation time ratio, T1/T2 (lambda)
%
%OUTPUT
%
%complex bSSFP signal for the given range of RF phase increments phi (bSSFP)

M0 = 1;

theta = off_resonance - phi;

D = 1+2*lambda+(1-2*lambda).*cos(flip)-cos(theta)-cos(flip).*cos(theta);

bSSFP = -(complex(0,1)./D).*sin(flip).*(1-exp(-complex(0,1)*theta));

%% Read out signal at t=TE

bSSFP = M0*bSSFP.*exp(complex(0,1)*off_resonance*(TE/TR));
% bSSFP = M0*bSSFP;
