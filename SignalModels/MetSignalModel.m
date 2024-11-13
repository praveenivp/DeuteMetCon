function [Msig_all,dc_fac]=MetSignalModel(Metabolites,TE,PhaseCyles,TR,dfreq,FA,Type,DutyCycle)
% [Msig_all]=SignalModel(Metabolites,TE,PhaseCyles,TR,dfreq,FA,Type)
% FA : flipangle in radians
% T1,T2,TE,TR are in seconds
% dfreq is frequency offset in Hz(meant for B0 map)
% Phasecyles in radians
% Species is a struct array with all relavent properties of your chemical compound
% Species(1)=struct('T1_s',300e-3,'T2_s',200e-3,'freq_shift_Hz',0,'name','water');
% TYPE: {'bSSFP','GRE','FISP'}
% Dutycycle : 0-1 scalar or fucntion handle as function of TR
if(~exist('DutyCycle','var'))
    DutyCycle=1;
end
dc_fac=zeros(length(Metabolites),length(TR));

% first vectorize along dfreq which can get quite slow with 3D field map
Msig_all=zeros(length(Metabolites),length(TE),length(PhaseCyles),length(TR),length(dfreq),length(FA));
for cMb=1:length(Metabolites)
    freqOffset=Metabolites(cMb).freq_shift_Hz+dfreq;
    T1=Metabolites(cMb).T1_s;
    T2=Metabolites(cMb).T2_s;
    T2star=Metabolites(cMb).T2star_s;
    
    for cTE=1:length(TE)
        for cPC=1:length(PhaseCyles)
            for cTR=1:length(TR)
                %if Dutycycle is a fucntion handle
                if(~isnumeric(DutyCycle))
                    cDutyCycle=DutyCycle(TR(cTR));  
                else
                    cDutyCycle=DutyCycle;
                end
                dc_fac(cMb,cTR)=cDutyCycle;
                for cFA=1:length(FA)

                    switch(Type)
                        case 'bSSFP'
                            phaseOffset=Freq2Phase(freqOffset,TR(cTR));
                            [Msig_all(cMb,cTE,cPC,cTR,:,cFA)]=...
                                bSSFP_profile_Ganter2(FA(cFA),T1,T2,TE(cTE),TR(cTR),PhaseCyles(cPC),phaseOffset);
                            % Caculate duty cycle factor: assuming T2* decay calc area under T2* relaxation curve from 0 to TR*DC                         
                            dcf=T2star*(exp(-min(TE)/T2star)-exp(-(min(TE)+TR(cTR)*cDutyCycle)/T2star));
                            dc_fac(cMb,cTR)= dcf./(TR(cTR)); %normalization
                        case 'bSSFP2'
                            phaseOffset=Freq2Phase(freqOffset,TR(cTR));
                            [Msig_all(cMb,cTE,cPC,cTR,:,cFA)]=...
                                bSSFP_profile_Ganter3(FA(cFA),T1,T2,TE(cTE),TR(cTR),PhaseCyles(cPC),phaseOffset);
                            % Caculate duty cycle factor: assuming T2 decay calc area under T2 relaxation curve from 0 to TR*DC                         
                            dcf=T2star*(exp(-min(TE)/T2star)-exp(-(min(TE)+TR(cTR)*cDutyCycle)/T2star));
                            dc_fac(cMb,cTR)= dcf./(TR(cTR)); %normalization
                        case {'FISP','SSFP-FID'}
                            [Msig_all(cMb,cTE,cPC,cTR,:,cFA)]=...
                                FISP(FA(cFA),T1,T2,T2star,TE(cTE),TR(cTR),freqOffset,cDutyCycle);
                            % Caculate duty cycle factor: assuming T2 decay calc area under T2 relaxation curve from 0 to TR*DC                         
                            dcf=T2star*(exp(-min(TE)/T2star)-exp(-(min(TE)+TR(cTR)*cDutyCycle)/T2star));
                            dc_fac(cMb,cTR)=dcf./(TR(cTR)); %normalization
                        case 'bSSFP-peters'
                            [Msig_all(cMb,cTE,cPC,cTR,:,cFA)]=...
                                bSSFP_peters(FA(cFA),T1,T2,TE(cTE),TR(cTR),freqOffset);
                            % Caculate duty cycle factor: assuming T2 decay calc area under T2 relaxation curve from TE to TR*DC                         
                            dcf=T2*(exp(-min(TE)/T2)-exp(-(min(TE)+TR(cTR)*cDutyCycle)/T2));
                            dc_fac(cMb,cTR)= dcf./(TR(cTR)); %normalization
                        case {'GRE','FLASH','GRE-peters'}
                            [Msig_all(cMb,cTE,cPC,cTR,:,cFA)]=...
                                GRESignal(FA(cFA),T1,T2star,TE(cTE),TR(cTR),freqOffset);
                            % % %dutycycle factor: integral of T2* exponential from TE to TR*DutyCycle
                            dcf=T2star*(exp(-min(TE)/T2star)-exp(-(min(TE)+TR(cTR)*cDutyCycle)/T2star));
                            dc_fac(cMb,cTR)=dcf./(TR(cTR)); %normalization
                        otherwise
                            error('unknown mode: {''bSSFP'',''GRE'',''bSSFP-peters'',''GRE-peters'',''FISP''}')

                    end
                end
            end
        end
    end
end

%DC=0 returns NaNs

Msig_all(isnan(Msig_all))=0;
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


% %dutycycle factor: integral of T2 exponential from TE to TR*DutyCycleDutyCycle
% DutyCycle=(TR-4e-3)/TR;
% dc_fac=T2*(exp(-TE/T2)-exp(-(TE+TR*DutyCycle)/T2));
% dc_fac=dc_fac./(DutyCycle*TR) %normalization
% bSSFP=bSSFP*dc_fac;

%% Read out signal at t=TE

bSSFP = M0*bSSFP.*exp(-TE./T2).*exp(1i*off_resonance*((TE)/TR));
end


function bSSFP=bSSFP_profile_Ganter3(flip,T1,T2,TE, TR,phi,off_resonance)
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

M0    = 1;
E1    = exp(-TR./T1);
E2    = exp(-TR./T2);

C     = E2.*(E1-1).*(1+cos(flip));
D     = (1-E1.*cos(flip))-(E1-cos(flip)).*E2.^2;
bSSFP = (1-E1).*sin(flip).*(1-E2.*exp(-complex(0,1)*(-phi+off_resonance)))./(C.*cos(-phi+off_resonance)+D);

%% Read out signal at t=TE

bSSFP = M0*bSSFP.*exp(-TE./T2).*exp(complex(0,-1)*(-1*off_resonance)*(TE/TR));
end


function bSSFP=bSSFP_peters(flip,T1,T2,TE, TR,freqOffset)
% Equation 2 from  DOI: 10.1002/mrm.28906

M0 = 1;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);


num=M0*exp(-(TR/2)/T2)*(1-E1).*sin(flip);

bSSFP=num./(1- (E1-E2)*cos(flip)-E1*E2);

%dutycycle factor: integral of T2* exponential from TE to TR*DutyCycleDutyCycle
% DutyCycle=0.8;
% dc_fac=T2*(exp(-TE/T2)-exp(-(TE+TR*DutyCycle)/T2));
% dc_fac=dc_fac./(DutyCycle*TR); %normalization
% bSSFP=bSSFP*dc_fac;

%phase evolution: not relavant for SNR
bSSFP=bSSFP.*exp(1i*2*pi*freqOffset*TE);
end


function [Sflash]=GRESignal(FA,T1,T2star,TE,TR,freqOffset)
%     sig=GREsignalmodel(flip,T1,T2star,TE,TR,freqOffset)

M0    = 1;
E1    = exp(-TR./T1);
% https://mriquestions.com/spoiled-gre-parameters.html
Sflash=M0*sin(FA)*(1-E1);
Sflash=Sflash./(1-cos(FA)*E1);
Sflash=Sflash.*exp(-TE/T2star);


%phase evolution: not relavant for SNR
Sflash=Sflash.*exp(1i*2*pi*freqOffset*TE);

end

function [Sfisp]=FISP(FA,T1,T2,T2star,TE,TR,freqOffset,DutyCycle)
%A. Oppelt, R. Graumann, H. Barfuss, H. Fischer, W. Hertl and W. Schajor. FISP: A new fast MRI sequence. Electromedica, 3: 15, 1986.

M0 = 1;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);
% DOI: 10.1002/mrm.10410
r=sqrt((1-E2.^2)./ ...
    ((1-E1*cos(FA)).^2-(E2.^2).*(E1-cos(FA)).^2));

Sfisp=M0*(sin(FA)/(1+cos(FA)))*(1-(E1-cos(FA)).*r);

% DOI: 10.1016/0022-2364(89)90083-8: infite series skipping
% https://labs.dgsom.ucla.edu/file/134807/M229_Lecture2_PulseSeqGRE.pdf
% p=1-E1*cos(flip)-(E2.^2).*(E1-cos(flip));
% q=E2.*(1-E1)*(1+cos(flip));

%phase evolution and T2 decay
Sfisp=Sfisp.*exp(1i*2*pi*freqOffset*TE).*exp(-TE./T2);

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