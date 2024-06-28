function metabolites=getMetaboliteStruct(DataSelect,freqOffset)
% metabolites=getMetaboliteStruct(DataSelect,freqOffset=0)
% simple function to get Metabolite structure with T1/T2 and chemical
% shifts
%
%Dataselect: one of {'phantom','invivo1','invivo2','Roig7T','Peters12T'}
%freqOffset - new water frequency in Hz

if(~exist("freqOffset",'var')), freqOffset=0; end

DataSelect=validatestring(DataSelect,{'phantom','invivo1','invivo2','Roig7T','Peters12T'});

B0=9.3879; %[T]
gammaH2=6.536; % [MHz/T]

switch(DataSelect)
    case 'phantom'
        %new phantom: /ptmp/pvalsala/deuterium/20240102_new2Hphantom
        met_name={'water','glucose','glutamic acid','lactate'};
        %         median(cell2mat(cellfun(@(x)[x.center1,x.center2,x.center3,x.center4],fitf,'UniformOutput',false)),1)
        freq_shift=[2.68  -63.83 -152.88 -213.58]; %Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T1=[444.30, 69.53, 157.65, 244.38]*1e-3; %s
        T2=[272.75,52.09,114.31,245.17]*1e-3; %s
        %         FWHM_hz=median(cell2mat(cellfun(@(x)[x.gamma1,x.gamma2,x.gamma3,x.gamma4],fitf,'UniformOutput',false)),1);
        FWHM_hz=[7.0280   23.7173   13.8651    6.1811]; %Hz
        T2Star=1./(pi*FWHM_hz);  %s
    case 'invivo1'
        % for invivo we have all values except lactate

        %invivo: /ptmp/pvalsala/deuterium/EAZW-GUMK/proc
        met_name={'water','glucose','Glx','lactate/lipid'};
        freq_shift=[2.0703  -51.2287 -138.3513 -195.4183]; % Hz from T1 data
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T1=[351.12 95.51  161.41 154.04]*1e-3;%s
        T2=[52.93 44.96 80.10 96.47]*1e-3;%s
        FWHM_hz =[15.3010   34.9108   20.1527   21.4540];% Hz from T1 data
        T2Star=1./(pi*FWHM_hz);  %s
    case 'invivo2'
        %invivo: /ptmp/pvalsala/deuterium/DA77-F3UY
        met_name={'water','glucose','Glx','lactate/lipid'};


        freq_shift=[7.5280  -44.7456 -132.6490 -196.8160]; %Hz from T1 data
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T1=[343.58 94.21 169.32 263.46]*1e-3;%s
        T2=[80.75 36.83 109.33 117.42]*1e-3;%s
        FWHM_hz =[22.0644   37.6203   18.9945    15.1671]; %Hz from T1 data
        T2Star=1./(pi*FWHM_hz);  %s

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