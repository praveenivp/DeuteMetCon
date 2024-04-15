function bSSFP=bSSFP_profile(TR,TE,flip,phi,shift,T1,T2)
%function bSSFP=bSSFP_profile(TR,TE,flip,phi,T1,T2)
%
%Calculation of the bSSFP profile (series of measurements with varying RF phase increment)
%
%INPUT
%
%repetition time in ms (TR)
%echo time in ms (TE)
%flip angle in rad (flip)
%vector of RF phase increments in rad (phi)
%off-resonance shift in rad (shift)
%spin-lattice relaxation time in ms (T1)
%spin-spin relaxation time in ms (T2)
%
%OUTPUT
%
%complex bSSFP signal for the given range of RF phase increments phi (bSSFP)

M0    = 1;
E1    = exp(-TR./T1);
E2    = exp(-TR./T2);

C     = E2.*(E1-1).*(1+cos(flip));
D     = (1-E1.*cos(flip))-(E1-cos(flip)).*E2.^2;
bSSFP = (1-E1).*sin(flip).*(1-E2.*exp(-complex(0,1)*(-phi+shift)))./(C.*cos(-phi+shift)+D);

%% Read out signal at t=TE

bSSFP = M0*bSSFP.*exp(-TE./T2).*exp(complex(0,-1)*phi*(TE/TR));
% bSSFP = M0*bSSFP.*exp(-TE./T2);
end