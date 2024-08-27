function Msig_all=MetSignalModel(Metabolites,TE,PhaseCyles,TR,dfreq,FA,Type,DutyCycle)
% [Msig_all]=SignalModel(Metabolites,TE,PhaseCyles,TR,dfreq,FA,Type)
% FA : flipangle in radians
% T1,T2,TE,TR are in seconds
% dfreq is frequency offset in Hz(meant for B0 map)
% Phasecyles in radians
% Species is a struct array with all relavent properties of your chemical compound
% Species(1)=struct('T1_s',300e-3,'T2_s',200e-3,'freq_shift_Hz',0,'name','water');
% TYPE: {'bSSFP','GRE','FISP'}
% Dutycycle : 0-1 scalar
if(~exist('DutyCycle','var'))
    DutyCycle=1;
end

% first vectorize along dfreq which can get quite slow with 3D field map
Msig_all=zeros(length(Metabolites),length(TE),length(PhaseCyles),length(TR),length(dfreq),length(FA));
for cMb=1:length(Metabolites)
    freqOffset=Metabolites(cMb).freq_shift_Hz+dfreq;
    T1=Metabolites(cMb).T1_s;
    T2=Metabolites(cMb).T2_s;
    
    for cTE=1:length(TE)
        for cPC=1:length(PhaseCyles)
            for cTR=1:length(TR)
                for cFA=1:length(FA)

                    switch(Type)
                        case 'bSSFP'
                            phaseOffset=Freq2Phase(freqOffset,TR(cTR));
                            [Msig_all(cMb,cTE,cPC,cTR,:,cFA)]=...
                                bSSFP_profile_Ganter2(FA(cFA),T1,T2,TE(cTE),TR(cTR),PhaseCyles(cPC),phaseOffset);
                        case 'GRE'
                            T2star=Metabolites(cMb).T2star_s;
                            [Msig_all(cMb,cTE,cPC,cTR,:,cFA)]=...
                                GRESignal(FA(cFA),T1,T2star,TE(cTE),TR(cTR),freqOffset);
                        case 'bSSFP-peters'
                            [Msig_all(cMb,cTE,cPC,cTR,:,cFA)]=...
                                bSSFP_peters(FA(cFA),T1,T2,TE(cTE),TR(cTR));
                        case 'GRE-peters'
                            T2star=Metabolites(cMb).T2star_s;
                            [Msig_all(cMb,cTE,cPC,cTR,:,cFA)]=...
                                GRESignal_peters(FA(cFA),T1,T2star,TE(cTE),TR(cTR),freqOffset,DutyCycle);
                        otherwise
                            error('unknown mode: {''bSSFP'',''GRE'',''bSSFP-peters'',''GRE-peters''}')

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

bSSFP = M0*bSSFP.*exp(-TE./T2).*exp(1i*off_resonance*((TE)/TR));
end


function bSSFP=bSSFP_peters(flip,T1,T2,TE, TR)
% Equation 2 from  DOI: 10.1002/mrm.28906

M0 = 1;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);


num=M0*exp(-(TR/2)/T2)*(1-E1).*sin(flip);

bSSFP=num./(1- (E1-E2)*cos(flip)-E1*E2);
%phase evolution: not relavant for SNR
% bSSFP=bSSFP.*exp(-1i*2*pi*freqOffset*TE);
end


function Sflash=GRESignal(FA,T1,T2star,TE,TR,freqOffset)
%     sig=GREsignalmodel(flip,T1,T2star,TE,TR,freqOffset)

M0    = 1;
E1    = exp(-TR./T1);
% https://mriquestions.com/spoiled-gre-parameters.html
Sflash=M0*sin(FA)*(1-E1)*exp(-TE/T2star);
Sflash=Sflash./(1-cos(FA)*E1);
Sflash=Sflash.*exp(-1i*2*pi*freqOffset*TE);

end

function Sflash=GRESignal_peters(FA,T1,T2star,TE,TR,freqOffset,DutyCycle)
%     sig=GREsignalmodel(flip,T1,T2star,TE,TR,freqOffset)

M0    = 1;
E1    = exp(-TR./T1);
% https://mriquestions.com/spoiled-gre-parameters.html
Sflash=M0*sin(FA)*(1-E1);
Sflash=Sflash./(1-cos(FA)*E1);

%dutycycle factor: integral of T2* exponential from TE to TR*DutyCycle
dc_fac=T2star*(exp(-TE/T2star)-exp(-(TE+TR*DutyCycle)/T2star));
dc_fac=dc_fac./(DutyCycle*TR); %normalization

Sflash=Sflash.*dc_fac;
%phase evolution: not relavant for SNR
Sflash=Sflash.*exp(-1i*2*pi*freqOffset*TE);

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