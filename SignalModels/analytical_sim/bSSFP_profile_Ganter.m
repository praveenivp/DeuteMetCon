function bSSFP=bSSFP_profile_Ganter(TR,TE,flip,phi,off_resonance,T1,T2)
%function bSSFP=bSSFP_profile_Ganter(TR,TE,flip,phi,T1,T2)
%
%Calculation of the bSSFP profile (series of measurements with varying RF phase increment)
%Following Ganter, MRM, 2006: 
%Steady State of Gradient Echo Sequences with Radiofrequency Phase Cycling: Analytical Solution, Contrast Enhancement with Partial Spoiling
%
%INPUT
%
%repetition time in ms (TR)
%echo time in ms (TE)
%flip angle in rad (flip)
%vector of RF phase increments in rad (phi)
%off-resonance-related shift in rad (off_resonance)
%spin-lattice relaxation time in ms (T1)
%spin-spin relaxation time in ms (T2)
%
%OUTPUT
%
%complex bSSFP signal for the given range of RF phase increments phi (bSSFP)
%
%
%Reference: eq 42-44
% Ganter, C. (2005). https://doi.org/10.1002/mrm.20736 
%
% Rahel Heule

M0 = 1;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);

theta = off_resonance - phi;

D = (1-E1.*cos(flip)).*(1-E2.*cos(theta))-(E1-cos(flip)).*(E2-cos(theta)).*E2;

bSSFP = -(1i./D).*(1-E1).*sin(flip).*(1-E2.*exp(-1i*theta));

%% Read out signal at t=TE

bSSFP = M0*bSSFP.*exp(-TE./T2).*exp(1i*off_resonance*(TE/TR));
end
