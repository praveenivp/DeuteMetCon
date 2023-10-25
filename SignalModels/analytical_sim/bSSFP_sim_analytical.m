function [Msig_all]=bSSFP_sim_analytical(Metabolites,TE,PhaseCyles,TR,dfreq,FA)
% [Msig_all]=bSSFP_sim_analytical(Metabolites,TE,PhaseCyles,TR,dfreq,FA)
% FA : flipangle in radians
% T1,T2,TE,TR are in seconds
% dfreq is frequency offset in Hz(meant for B0 map)
% Phasecyles in radians
% Species is a struct array with all relavent properties of your chemical compound
% Species(1)=struct('T1_s',300e-3,'T2_s',200e-3,'freq_shift_Hz',0,'name','water');



Msig_all=zeros(1,length(Metabolites),length(TE),length(PhaseCyles),length(TR),length(dfreq),length(FA));
for cMb=1:length(Metabolites)
%     if(cMb==3) scal=5; else scal=1;end
    for cTE=1:length(TE)
        for cPC=1:length(PhaseCyles)
            for cTR=1:length(TR)      
                for cFA=1:length(FA)
                    for cFreq=1:length(dfreq)
                        
                        freqoffset=Metabolites(cMb).freq_shift_Hz+dfreq(cFreq);
                        [Msig_all(:,cMb,cTE,cPC,cTR,cFreq,cFA)]=...
                            bSSFP_profile_Ganter2(FA(cFA),Metabolites(cMb).T1_s,Metabolites(cMb).T2_s,TE(cTE),TR(cTR),PhaseCyles(cPC),Freq2Phase(freqoffset,TR(cTR)));
                    end
                end
            end
        end
    end
end

end

function bSSFP=bSSFP_profile_Ganter2(flip,T1,T2,TE, TR,phi,off_resonance)
%function bSSFP=bSSFP_profile_Ganter(TR,TE,flip,phi,T1,T2)
%
%Calculation of the bSSFP profile (series of measurements with varying RF phase increment)
%Following Ganter, MRM, 2006: 
%Steady State of Gradient Echo Sequences with Radiofrequency Phase Cycling: Analytical Solution, Contrast Enhancement with Partial Spoiling
%
%INPUT
%
%repetition time in s (TR)
%echo time in s (TE)
%flip angle in rad (flip)
%vector of RF phase increments in rad (phi)
%off-resonance-related shift in rad (off_resonance)
%spin-lattice relaxation 
%spin-spin relaxation time in s (T2)
%
%OUTPUT
%
%complex bSSFP signal for the given range of RF phase increments phi (bSSFP)

M0 = 1;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);

theta = off_resonance - phi;
% Equation 44
D = (1-E1.*cos(flip)).*(1-E2.*cos(theta))-(E1-cos(flip)).*(E2-cos(theta)).*E2;
% equation 42
bSSFP = -(1i./D).*(1-E1).*sin(flip).*(1-E2.*exp(-1i*theta));

%% Read out signal at t=TE

bSSFP = M0*bSSFP.*exp(-TE./T2).*exp(1i*off_resonance*(TE/TR));
end


function bSSFP=bSSFP_profile2(flip,T1,T2,TE, TR,phi,off_resonance)
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
%off-resonance shift in rad (off_resonance)
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
bSSFP = (1-E1).*sin(flip).*(1-E2.*exp(-complex(0,1)*(-phi+off_resonance)))./(C.*cos(-phi+off_resonance)+D);

%% Read out signal at t=TE

bSSFP = M0*bSSFP.*exp(-TE./T2).*exp(complex(0,1)*off_resonance*(TE/TR));
% bSSFP = M0*bSSFP.*exp(-TE./T2);

end

%%

% ph_rad=freq2Phase(freq_hz,TR_s);

function phase_rad=Freq2Phase(freq_Hz,TR_s)
% ph_rad=freq2Phase(freq_Hz,TR_s)
% all these things are equivalent
% -1/TR to +1/TR  <====> -360 to 360 deg  <===> -2*pi to 2*pi

% if(freq_Hz<0)
%     freq_Hz= -1*mod(freq_Hz,1/TR_s);
% else
%     freq_Hz= mod(freq_Hz,1/TR_s);
% end
% 
% phase_rad=interp1(linspace(-1/TR_s,1/TR_s,1e3),linspace(-2*pi,2*pi,1e3),freq_Hz(:));
 phase_rad=freq_Hz*2*pi*TR_s;

end

function freq_Hz=Phase2Freq(phase_rad,TR_s)
% freq_Hz=Phase2Freq(phase_rad,TR_s)
% all these things are equivalent
% -1/TR to +1/TR  <====> -360 to 360 Hz  <===> -2*pi to 2*pi

if(phase_rad<0)
    phase_rad= -1*mod(phase_rad,2*pi);
else
    phase_rad= mod(phase_rad,2*pi);
end

freq_Hz=interp1(linspace(-2*pi,2*pi,1e3),linspace(-1/TR_s,1/TR_s,1e3),phase_rad(:));
end