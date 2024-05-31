function metabolites=getMetaboliteStruct(DataSelect,freqOffset)
% metabolites=getMetaboliteStruct(DataSelect,freqOffset=0)
% simple function to get Metabolite structure with T1/T2 and chemical
% shifts
%
%Dataselect: one of {'phantom','invivo1','invivo2'}
%freqOffset - new water frequency in Hz

if(~exist("freqOffset",'var')), freqOffset=0; end

DataSelect=validatestring(DataSelect,{'phantom','invivo1','invivo2'});


switch(DataSelect)
    case 'phantom'
        %new phantom: /ptmp/pvalsala/deuterium/20240102_new2Hphantom
        met_name={'water','glucose','glutamic acid','lactate'};
%         median(cell2mat(cellfun(@(x)[x.center1,x.center2,x.center3,x.center4],fitf,'UniformOutput',false)),1)
        freq_shift=[2.68  -63.83 -152.88 -213.58]; %Hz

        T1=[444.30, 69.53, 157.65, 244.38]*1e-3; %s
        T2=[272.75,52.09,114.31,245.17]*1e-3; %s
        %         FWHM_hz=median(cell2mat(cellfun(@(x)[x.gamma1,x.gamma2,x.gamma3,x.gamma4],fitf,'UniformOutput',false)),1);
        FWHM_hz=[7.0280   23.7173   13.8651    6.1811]; %Hz
        T2Star=1000./(pi*FWHM_hz);  %ms
    case 'invivo1'
        % for invivo we have all values except lactate
           
        %invivo: /ptmp/pvalsala/deuterium/EAZW-GUMK/proc
        met_name={'water','glucose','Glx','lactate/lipid'};
        freq_shift=[2.0703  -51.2287 -138.3513 -195.4183]; % Hz from T1
        T1=[351.12 95.51  161.41 154.04]*1e-3;%s
        T2=[52.93 44.96 80.10 96.47]*1e-3;%s
        FWHM_hz =[15.3010   34.9108   20.1527   21.4540];% Hz from T1
        T2Star=1000./(pi*FWHM_hz);  %ms
    case 'invivo2'
        %invivo: /ptmp/pvalsala/deuterium/DA77-F3UY
        met_name={'water','glucose','Glx','lactate/lipid'};
        
 
        freq_shift=[7.5280  -44.7456 -132.6490 -196.8160]; %Hz from T1
        T1=[343.58 94.21 169.32 263.46]*1e-3;%s
        T2=[80.75 36.83 109.33 117.42]*1e-3;%s
        FWHM_hz =[22.0644   37.6203   18.9945    15.1671]; %Hz from T1
        T2Star=1000./(pi*FWHM_hz);  %ms
end
%use the new water reference
if(freqOffset~=0)
    freq_shift=freq_shift-freq_shift(1)+freqOffset;
end

for i=1:4 
    metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift(i),'name',met_name{i});
end


end