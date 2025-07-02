function metabolites=getMetaboliteStruct(DataSelect,freqOffset)
% metabolites=getMetaboliteStruct(DataSelect,freqOffset=0)
% simple function to get Metabolite structure with T1/T2 and chemical
% shifts
%
%Dataselect: one of {'phantom','invivo'(default),'sub-01','sub-02','sub-03','sub-04','sub-05','Roig7T','Peters12T'}
%freqOffset - new water frequency in Hz
if(~exist('DataSelect','var')), DataSelect='invivo';end % 9T fata
DataSelect=validatestring(DataSelect,{'phantom','invivo','sub-01','sub-02','sub-03','sub-04','sub-05','Roig7T','Peters12T'});


B0=9.3879; %[T]
gammaH2=6.536; % [MHz/T]

switch(DataSelect)
    case 'phantom'
        %new phantom: /ptmp/pvalsala/deuterium/20240102_new2Hphantom
        met_name={'water','glucose','glutamic acid','lactate'};
        T1=[449.8212,60.3493,162.9081,241.4547]*1e-3; %s
        T1_CI=[47.9218,23.5754,34.9505,24.2403]*1e-3; %s diff(CI95)/2
T2=[271.8672,55.5163,108.2194,264.3029]*1e-3; %s
T2_CI=[9.9736,7.5819,15.1837,31.2136]*1e-3; %s diff(CI95)/2
        ImagingFreq=61359114.0000; %Hz
        freq_shift=[2.68  -63.83 -152.88 -213.58]; %Hz has less residue

        %CSI dataset='20240813_spectral/meas_MID00857_FID14528_rpcsi_fid_Stan_res15_6_optimal.dat'; see plotT2star.m
       freq_shift=[0.0000  -64.4075 -153.6070 -212.9762]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[74.7374   22.0817   53.4715   75.7936]*1e-3; % s
        T2Star_std=[38.9049    2.8619   14.3097   28.7368]*1e-3; % s

    case 'invivo' % average all subjects: see last session to generate it
        met_name={'water','glucose','Glx','lactate/lipid'};

T1=[370.1214, 77.6572, 161.7127, 136.6657]*1e-3; %s 
T1_std=[26.6857, 8.3173, 6.7779, 21.5591]*1e-3; %s 
T2=[29.2700, 30.0429, 58.8559, 59.1913]*1e-3; %s 
T2_std=[3.5711, 7.7941, 15.2235, 58.2815]*1e-3; %s 
T2Star=[20.8950, 14.3771, 22.7435, 23.8578]*1e-3; %s 
T2Star_std=[2.1065, 1.5861, 1.4496, 5.3674]*1e-3; %s 
freq_shift=[0.0000, -54.9112, -144.2952, -205.6203]; %Hz 
freq_shift_std=[0.0000, 0.2879, 0.9122, 1.6424]; %Hz 
        water2_fac=mean([12.18,12.27,10.91,8.37,8.74]);
        water2_T2=mean([331.47, 365.04,340.79,317.66,200]);

        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;


    case 'sub-01' %subject 1
        %/ptmp/pvalsala/deuterium/HOSJ-D6P2/proc/T1T2
        met_name={'water','glucose','Glx','lactate/lipid'};
        T1=[360.5184,86.8361,158.0989,109.9415]*1e-3; %s
        T1_CI=[5.5403,9.4113,3.6082,18.5691]*1e-3; %s diff(CI95)/2
T2=[28.9490,35.9793,65.6832,59.6700]*1e-3; %s
T2_CI=[NaN,6.1671,7.4955,18.2072]*1e-3; %s diff(CI95)/2
%T2(91.05%) =28.9490 and T2(8.95%) =200.0044
        ImagingFreq=61359053.0000; %Hz
        %CSI dataset: M997 TR36 % see plotT2star.m
        freq_shift=[1.5423 -53.1530 -142.0912 -202.1538]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[21.0052 13.9556 22.3202 31.8583]*1e-3; %s

    case 'sub-02' %subject 2
        %/ptmp/pvalsala/deuterium/H4
        met_name={'water','glucose','Glx','lactate/lipid'};
        T1=[417.2924,64.3759,173.4299,142.9516]*1e-3; %s
        T1_CI=[55.7675,4.2608,15.5232,40.3232]*1e-3; %s diff(CI95)/2
T2=[31.6194,34.5575,73.5354,157.5316]*1e-3; %s
T2_CI=[NaN,4.1237,14.9912,131.8628]*1e-3; %s diff(CI95)/2
%T2(86.43%) =31.6194 and T2(13.57%) =300.1139
        ImagingFreq=61359047.0000; %Hz
        %CSI dataset: M997 TR36 % see plotT2star.m
        freq_shift=[ -1.0121 -55.5991 -145.1440 -207.7830]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[19.0796 12.7327 21.4821 17.5425]*1e-3; %s

    case 'sub-03' %subject 3
        %/ptmp/pvalsala/deuterium/H4
        met_name={'water','glucose','Glx','lactate/lipid'};
        T1=[361.3593,77.4438,157.5635,121.6826]*1e-3; %s
        T1_CI=[6.8110,4.4586,4.7049,9.6327]*1e-3; %s diff(CI95)/2
T2=[23.8321,17.3393,33.4110,7.7039]*1e-3; %s
T2_CI=[NaN,3.8072,21.2009,3.6314]*1e-3; %s diff(CI95)/2
%T2(91.41%) =23.8321 and T2(8.59%) =321.9429
        ImagingFreq=61359054.0000; %Hz
        %CSI dataset: M997 TR36 % see plotT2star.m
        freq_shift=[ 3.7134 -51.1761 -139.8969 -203.9604]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[19.4120 13.9282 22.3839 20.7921]*1e-3; %s
%         T2* std  [ms]  9.6127 11.3821 48.6879 67.1582
    case 'sub-04'
        % for invivo we have all values except lactate
T2=[33.2114,34.5932,58.6195,44.0343]*1e-3; %s
T2_CI=[NaN,4.3507,10.3322,16.6538]*1e-3; %s diff(CI95)/2
%T2(94.82%) =33.2114 and T2(5.18%) =212.3880
        %invivo: /ptmp/pvalsala/deuterium/EAZW-GUMK/proc
        met_name={'water','glucose','Glx','lactate/lipid'};
        T1=[351.2728,81.6462,157.6897,143.1148]*1e-3; %s
        T1_CI=[19.3575,7.5787,4.6537,23.5573]*1e-3; %s diff(CI95)/2

        ImagingFreq=61359045.0000; %Hz
        %CSI dataset: M412 % see plotT2star.m
        freq_shift=[-0.6647  -55.7530 -144.9229 -205.0421]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[24.3812   17.0160   25.2516   23.8431]*1e-3; %s
    case 'sub-05' %subject 5
        %invivo: /ptmp/pvalsala/deuterium/DA77-F3UY
        met_name={'water','glucose','Glx','lactate/lipid'};
        T1=[360.1643,77.9841,161.7815,165.6381]*1e-3; %s
        T1_CI=[7.2214,1.9944,5.9712,46.3589]*1e-3; %s diff(CI95)/2
T2=[28.7383,27.7451,63.0305,27.0167]*1e-3; %s
T2_CI=[NaN,8.0573,28.5160,25.3087]*1e-3; %s diff(CI95)/2
%T2(86.92%) =28.7383 and T2(13.08%) =239.6685
        freq_shift=[3.8518,-51.4442,-141.9904,-201.7315]; %Hz median
        ImagingFreq=61359063.0000; %Hz
        %CSI dataset: M695 % see plotT2star.m
        %         freq_shift=[9.0229 -46.9185 -135.4832 -198.4703]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[20.5969 14.2529 22.2796 25.2529]*1e-3; %s

    case 'Roig7T'
        %         DOI: 10.1002/mrm.29439
        met_name={'water1','water2','glucose','Glx','lactate/lipid'};
        freq_shift_ppm=[4.7 ,4.7,3.8,2.3 1.5];
        freq_shift=(freq_shift_ppm-freq_shift_ppm(1))*B0*gammaH2;
        T1=[496 283 66 149 159]*1e-3;%s
        T2=[412 29 44 53 159/2]*1e-3;%s
        T2Star=[20.8950,20.8950, 14.3771, 22.7435, 23.8578]*1e-3; %s % 9.4T
    case 'Peters12T'
        %         DOI: 10.1002/mrm.28906
        met_name={'water','glucose','lactate'};
        freq_shift_ppm=[4.8,3.7,1.3];

        freq_shift=(freq_shift_ppm-freq_shift_ppm(1))*B0*gammaH2;
        T1=[320 64 297]*1e-3;%s
        T2=[12 32 61]*1e-3;%s
        T2Star=[10 15 15]*1e-3; %s 15.2 T measurement

    case 'Feyter4T' %humans 10.1126/sciadv.aat7314 
        met_name={'water1','water2','glucose','Glx'};
        freq_shift_ppm=[4.7 ,4.7,3.8,2.3];
        freq_shift=(freq_shift_ppm-freq_shift_ppm(1))*B0*gammaH2;
        T1=[346 346 67 139 159]*1e-3;%s
        T2=[26 351 42 53 159/2]*1e-3;%s

end
%use the new water reference
if(exist("freqOffset",'var'))
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




%% plotting mean
% all_sub={'sub-01','sub-02','sub-03','sub-04','sub-05'};
% metabolites_all=cellfun(@(x) getMetaboliteStruct(x,0),all_sub,'UniformOutput',false);
% 
% fprintf('T1=[%.4f, %.4f, %.4f, %.4f]*1e-3; %%s \n',...
%     1e3*mean(cell2mat(cellfun (@(x)[x.T1_s],metabolites_all,'UniformOutput',false)')));
% fprintf('T1_std=[%.4f, %.4f, %.4f, %.4f]*1e-3; %%s \n',...
%     1e3*std(cell2mat(cellfun (@(x)[x.T1_s],metabolites_all,'UniformOutput',false)')));
% fprintf('T2=[%.4f, %.4f, %.4f, %.4f]*1e-3; %%s \n',...
%     1e3*mean(cell2mat(cellfun (@(x)[x.T2_s],metabolites_all,'UniformOutput',false)')));
% fprintf('T2_std=[%.4f, %.4f, %.4f, %.4f]*1e-3; %%s \n',...
%     1e3*std(cell2mat(cellfun (@(x)[x.T2_s],metabolites_all,'UniformOutput',false)')));
% 
% fprintf('T2Star=[%.4f, %.4f, %.4f, %.4f]*1e-3; %%s \n',...
%     1e3*mean(cell2mat(cellfun (@(x)[x.T2star_s],metabolites_all,'UniformOutput',false)')));
% fprintf('T2Star_std=[%.4f, %.4f, %.4f, %.4f]*1e-3; %%s \n',...
%     1e3*std(cell2mat(cellfun (@(x)[x.T2star_s],metabolites_all,'UniformOutput',false)')));
% 
% fprintf('freq_shift=[%.4f, %.4f, %.4f, %.4f]; %%Hz \n',...
%     mean(cell2mat(cellfun (@(x)[x.freq_shift_Hz],metabolites_all,'UniformOutput',false)')));
% fprintf('freq_shift_std=[%.4f, %.4f, %.4f, %.4f]; %%Hz \n',...
%     std(cell2mat(cellfun (@(x)[x.freq_shift_Hz],metabolites_all,'UniformOutput',false)')));
