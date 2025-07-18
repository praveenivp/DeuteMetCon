function metabolites=getMetaboliteStruct(DataSelect,freqOffset)
% metabolites=getMetaboliteStruct(DataSelect,freqOffset=0)
% simple function to get Metabolite structure with T1/T2, T2 star and chemical
% shifts
%
%EXAMPLES:
% to print table use without output arguments
%getMetaboliteStruct('invivo')
%met_Struct=getMetaboliteStruct('phantom');
%
%Dataselect: one of {'sub-01'}    {'sub-02'}    {'sub-03'}    {'sub-04'}    {'sub-05'}    {'phantom'}    {'invivo'}    {'invivo9T'}    {'Roig7T'}    {'Peters12T'}
%freqOffset - new water frequency in Hz

DataSelect=validatestring(DataSelect,[compose('sub-%02d', 1:5),{'phantom','invivo','invivo9T','Roig7T','Peters12T'}]);

B0=9.3879; %[T]
gammaH2=6.536; % [MHz/T]

switch(DataSelect)
    case 'phantom' %@9.4T
        %see: DeuteMetCon/scripts/process_all_relaxometry.m
        met_name={'water','glucose','glutamic acid','lactate'};
        T1=[458.7212 67.5265 160.5774 319.3940]*1e-3; %s
        T1_std=[4.3498 10.3023 20.8912 41.5633]*1e-3; %s
        T2=[304.1309 56.8033 156.8406 254.8299]*1e-3; %s
        T2_std=[3.6967 8.7555 17.7546 35.3229]*1e-3; %s
        ImagingFreq=61359114.0000; %Hz
        freq_shift=[2.68  -63.83 -152.88 -213.58]; %Hz has less residue

        %CSI dataset='20240813_spectral/meas_MID00857_FID14528_rpcsi_fid_Stan_res15_6_optimal.dat'; see plotT2star.m
        freq_shift=[0.0000  -64.4075 -153.6070 -212.9762]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[74.7374   22.0817   53.4715   75.7936]*1e-3; % s
        T2Star_std=[38.9049    2.8619   14.3097   28.7368]*1e-3; % s


    case [compose('sub-%02d', 1:5),{'invivo','invivo9T'}] %@9.4T
        met_name={'water','glucose','Glx','lactate/lipid'};
        subIdx=str2double(DataSelect(end-1:end));

        % see: DeuteMetCon/scripts/process_all_relaxometry.m
        % data from process_T1_inv
        T1 = [
            372.1884   71.1968  142.7747   67.7676
            400.0000   61.7577  134.3558   44.5676
            342.7002   65.8658  137.0337   78.7666
            358.2341   66.6929  144.9041   86.7971
            365.0854   63.3548  139.2259   75.5844
            ]; % in ms

        T1_std = [
            3.3597    1.4793    2.1845    4.3329
            42.1997    5.0280   13.8865   10.8930
            49.8865    8.9086   14.2134   16.6408
            0.0059    0.0121    0.0074    0.0011
            2.8877    1.5630    2.9416   20.4450
            ] ; % in ms

        %data from process_T2.m
        % water_T2,Glu_T2,Glx_t2,lipid_T2,Water2_T2 in ms and water2_fac
        T2 = [
            28.1280   48.6978   98.8150  107.8537  308.0978    0.1693
            29.4282   39.7355   99.1927  171.5556  333.7481    0.2411
            35.0754   30.4258   91.2093  158.4283  244.6254    0.1016
            29.3883   45.9762   96.0668  143.1346  277.8036    0.1822
            28.0758   40.4198  105.1521  107.9723  303.4141    0.2154
            ];
        %standard devitaiton of above with numerical hessian matrix
        T2_std = [
            0.5342    1.8976    3.4090   11.8423   16.0746    0.0059
            0.6030    1.0390    4.0714   16.4482   12.4866    0.0056
            44.7541  201.3917   16.0118   22.8937  181.6521    0.2994
            0.5321    1.6238    3.4144   15.1081   12.6802    0.0057
            0.4265    0.9315    3.4549    9.8602   10.7013    0.0035
            ];



        % CSI datasets see plotT2star.m
        csi_datasets = {'M997', 'M997', 'M997', 'M412', 'M695'};
        pseudo={'HOSJ','H4','I3BL','DA77','EAZW'};
        T2Star = [
            21.0052   13.9556   22.3202   31.8583;
            19.0796   12.7327   21.4821   17.5425;
            19.4120   13.9282   22.3839   20.7921;
            24.3812   17.0160   25.2516   23.8431;
            20.5969   14.2529   22.2796   25.2529
            ] * 1e-3; % Convert to seconds
        ImagingFreq = [61359053; 61359047; 61359054; 61359045; 61359063]; % [Hz]
        freq_shift = [
            1.5423   -53.1530  -142.0912  -202.1538;
            -1.0121   -55.5991  -145.1440  -207.7830;
            3.7134   -51.1761  -139.8969  -203.9604 ;
            -0.6647   -55.7530  -144.9229  -205.0421;
            3.8518   -51.4442  -141.9904  -201.7315;
            ]; % [Hz]

        %slice or take mean!
        if(isnan(subIdx))
            T1=mean(T1,1)*1e-3;
            T2=mean(T2(:,1:4),1)*1e-3; % ms -> s
            T2Star=mean(T2Star,2);
            freq_shift=mean(freq_shift-freq_shift(:,1),1);
        else
            T1=T1(subIdx,:)*1e-3;
            T2=T2(subIdx,1:4)*1e-3;
            T2Star=T2Star(subIdx,:);
            freq_shift=freq_shift(subIdx,:);
        end
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;

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




