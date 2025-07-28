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
        T2=[275.0866 50.5532 100.3155 260.7906]*1e-3; %s
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
        T1_allsub = [
            366.7831   88.4618  159.8666  101.5152
            424.4544   73.4691  170.2451   66.7176
            369.4601   77.2360  160.4962  101.1598
            353.0484   81.8572  157.2119  155.0612
            360.1567   75.1441  161.8468  158.7093
            ]; % in ms

        T1_std_allsub = [
            2.8037    2.6811    5.0098    8.0059
            3.8719    1.6209    6.9902    7.6678
            2.8241    1.5624    5.0157    7.1220
            2.2168    2.1873    3.8841   12.9800
            2.8606    2.2686    6.0944   23.8939
            ] ; % in ms

        %data from process_T2.m
        % water_T2,Glu_T2,Glx_t2,lipid_T2,Water2_T2 in ms and water2_fac
        T2_allsub = [
            28.1168   48.6983   98.8198  107.8784  307.6862    0.1695
            29.4018   39.7353   99.3155  169.5491  332.9553    0.2415
            35.0685   30.4256   91.2100  158.4159  244.3296    0.1018
            29.3695   45.9765   96.0736  143.1393  277.2598    0.1825
            28.0233   40.0779  105.3374  109.4906  302.7273    0.2168
            ];
        %standard devitaiton of above with numerical hessian matrix
        T2_std_allsub = [
            0.5330    1.9541    3.6833   12.1516   16.3777    0.0060
            0.6226    1.1172    4.4191   16.4863   12.3983    0.0057
            1.7423    0.7610    3.5450   22.9434   47.5489    0.0146
            0.6150    1.5767    3.2768   15.8339   15.2681    0.0071
            0.4357    0.9374    3.5570    5.4268    9.6186    0.0043
            ];



        % CSI datasets see plotT2star.m
        csi_datasets = {'M997', 'M997', 'M997', 'M412', 'M695'};
        pseudo={'HOSJ','H4','I3BL','DA77','EAZW'};
        T2Star_allsub = [
            21.0052   13.9556   22.3202   31.8583;
            19.0796   12.7327   21.4821   17.5425;
            19.4120   13.9282   22.3839   20.7921;
            24.3812   17.0160   25.2516   23.8431;
            20.5969   14.2529   22.2796   25.2529
            ] * 1e-3; % Convert to seconds
        ImagingFreq_allsub = [61359053; 61359047; 61359054; 61359045; 61359063]; % [Hz]
        freq_shift_allsub = [
            1.5423   -53.1530  -142.0912  -202.1538;
            -1.0121   -55.5991  -145.1440  -207.7830;
            3.7134   -51.1761  -139.8969  -203.9604 ;
            -0.6647   -55.7530  -144.9229  -205.0421;
            3.8518   -51.4442  -141.9904  -201.7315;
            ]; % [Hz]

        %slice or take mean!
        if(isnan(subIdx))
            T1=mean(T1_allsub,1)*1e-3;
            T2=mean(T2_allsub(:,1:4),1)*1e-3; % ms -> s
            T2Star=mean(T2Star_allsub,2);
            freq_shift=mean(freq_shift_allsub-freq_shift_allsub(:,1),1);
        else
            T1=T1_allsub(subIdx,:)*1e-3;
            T2=T2_allsub(subIdx,1:4)*1e-3;
            T2Star=T2Star_allsub(subIdx,:);
            freq_shift=freq_shift_allsub(subIdx,:);
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




