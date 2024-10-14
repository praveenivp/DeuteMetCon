function metabolites=getMetaboliteStruct(DataSelect,freqOffset)
% metabolites=getMetaboliteStruct(DataSelect,freqOffset=0)
% simple function to get Metabolite structure with T1/T2 and chemical
% shifts
%
%Dataselect: one of {'phantom','invivo1','invivo2','Roig7T','Peters12T'}
%freqOffset - new water frequency in Hz

if(~exist("freqOffset",'var')), freqOffset=0; end

DataSelect=validatestring(DataSelect,{'phantom','invivo1','invivo2','invivo3','Roig7T','Peters12T'});

B0=9.3879; %[T]
gammaH2=6.536; % [MHz/T]

switch(DataSelect)
    case 'phantom'
        %new phantom: /ptmp/pvalsala/deuterium/20240102_new2Hphantom
        met_name={'water','glucose','glutamic acid','lactate'};

        T1=[444.30, 69.53, 157.65, 244.38]*1e-3; %s
        T2=[272.75,52.09,114.31,245.17]*1e-3; %s
        %CSI dataset='20240813_spectral/meas_MID00857_FID14528_rpcsi_fid_Stan_res15_6_optimal.dat'; see plotT2star.m
        freq_shift=[-3.4080 -63.5005 -150.2146 -205.1276]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[77.9295   18.8752   14.0095   15.9714]*1e-3; % s
    case 'invivo1'
        % for invivo we have all values except lactate

        %invivo: /ptmp/pvalsala/deuterium/EAZW-GUMK/proc
        met_name={'water','glucose','Glx','lactate/lipid'};

        T1=[351.12 95.51  161.41 154.04]*1e-3;%s
        T2=[52.93 44.96 80.10 96.47]*1e-3;%s
        %CSI dataset: M412 % see plotT2star.m
        freq_shift=[-0.6647  -55.7530 -144.9229 -205.0421]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[24.3812   17.0160   25.2516   23.8431]*1e-3; %s
    case 'invivo2'
        %invivo: /ptmp/pvalsala/deuterium/DA77-F3UY
        met_name={'water','glucose','Glx','lactate/lipid'};
        T1=[343.58 94.21 169.32 263.46]*1e-3;%s
        T2=[80.75 36.83 109.33 117.42]*1e-3;%s
        %CSI dataset: M695 % see plotT2star.m
        freq_shift=[9.0229 -46.9185 -135.4832 -198.4703]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[20.5969 14.2529 22.2796 25.2529]*1e-3; %s

    case 'invivo3'
        %/ptmp/pvalsala/deuterium/HOSJ-D6P2/proc/T1T2
        met_name={'water','glucose','Glx','lactate/lipid'};
        T1=[362.21 64.10  155.23 208.37]*1e-3;%s
        T2=[41.65 55.59 109.13 98.31]*1e-3;%s
        %CSI dataset: M997 TR36 % see plotT2star.m
        freq_shift=[1.5423 -53.1530 -142.0912 -202.1538]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[21.0052 13.9556 22.3202 31.8583]*1e-3; %s

    case 'Roig7T'
%         DOI: 10.1002/mrm.29439
        met_name={'water1','water2','glucose','Glx','lactate/lipid'};
        freq_shift_ppm=[4.7 ,4.7,3.8,2.3 1.5];
        freq_shift=(freq_shift_ppm-freq_shift_ppm(1))*B0*gammaH2;
        T1=[496 283 66 149 159]*1e-3;%s
        T2=[412 29 44 53 159/2]*1e-3;%s

    case 'Peters12T'
%         DOI: 10.1002/mrm.28906
        met_name={'water','glucose','lactate'};
        freq_shift_ppm=[4.8,3.7,1.3];

        freq_shift=(freq_shift_ppm-freq_shift_ppm(1))*B0*gammaH2;
        T1=[320 64 297]*1e-3;%s
        T2=[12 32 61]*1e-3;%s
        T2Star=[8 13 15]*1e-3; %s 15.2 T measurement

end
%use the new water reference
if(freqOffset~=0)
    freq_shift=freq_shift-freq_shift(1)+freqOffset;
end

for i=1:length(met_name)
    metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift(i),'name',met_name{i}, ...
        'T2star_s',T2Star(i));
end


%print if no output arguments
if(nargout==0)
disp(table(1e3*[metabolites.T1_s;]',1e3*[metabolites.T2_s;]',1e3*[metabolites.T2star_s;]',freq_shift_ppm',[metabolites.freq_shift_Hz;]','VariableNames',{'T1 [ms]','T2 [ms]','T2* [ms]','chemical shift [ppm]','offset @9.4T [Hz] '},'RowNames',{metabolites.name}))
clear metabolites;
end
end