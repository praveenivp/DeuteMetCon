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
        T1= [488.1285 64.2258 210.3209 266.4113]*1e-3; %s
        T1_std=[11.4140 8.2099 28.9405 31.9601]*1e-3; %s
        T2=[275.4814 60.4826 122.4650 257.8497]*1e-3; %s
        T2_std=[1.2859 8.5585 10.4005 13.8038]*1e-3; %s
        ImagingFreq=61359114.0000; %Hz
        freq_shift=[2.68  -63.83 -152.88 -213.58]; %Hz has less residue

        %CSI dataset='20240813_spectral/meas_MID00857_FID14528_rpcsi_fid_Stan_res15_6_optimal.dat'; see plotT2star.m
        freq_shift=[0.0000  -64.4075 -153.6070 -212.9762]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[74.7374   22.0817   53.4715   75.7936]*1e-3; % s
        T2Star_std=[38.9049    2.8619   14.3097   28.7368]*1e-3; % s


         % T2Star=[39.7374   17.0817   22.4715   17.7936]*1e-3; % s

    case [compose('sub-%02d', 1:5),{'invivo','invivo9T'}] %@9.4T
        met_name={'water','glucose','Glx','lactate/lipid'};
        subIdx=str2double(DataSelect(end-1:end));

        % see: DeuteMetCon/scripts/process_all_relaxometry.m
        % data from process_T1_inv
        T1_allsub = [

  365.5809   78.6077  154.0151   60.1361
  421.0525   68.2742  158.5058   31.2551
  365.4739   71.2702  152.3236   76.4359
  352.0298   72.6661  154.6728  129.3595
  357.7670   64.1627  158.6263   52.0203
            ]; % in ms

        T1_std_allsub = [
    2.9805    2.6137    3.9621    7.0862
    3.7796    1.7737    7.0019    7.1100
    3.0743    1.7957    5.4118   13.8731
    2.5172    1.8411    3.7039   20.8335
    2.8330    1.9057    7.1643   33.2767
            ] ; % in ms

        %data from process_T2.m
        % water_T2,Glu_T2,Glx_t2,lipid_T2,Water2_T2 in ms and water2_fac
        T2_allsub = [
   28.1167   48.6983   98.8190  107.8799  307.6817    0.1695
   29.4005   39.7347   99.2027  171.5739  332.9842    0.2415
   35.0676   30.4257   91.2105  158.4321  244.2887    0.1018
   29.3691   45.9766   96.0720  143.1650  277.2470    0.1825
   28.4935   40.6221  104.7313   84.9386  316.1225    0.2080
    
            ];
        %standard devitaiton of above with numerical hessian matrix
        T2_std_allsub = [
    0.7910    2.9509    5.3958   18.1988   23.8167    0.0085
    1.2500    2.2902    8.6840   34.1884   24.7705    0.0111
    5.0645    2.4310    9.6807   69.1844  259.7372    0.0692
    1.0222    2.4055    5.2101   25.1964   23.4493    0.0115
    1.4086    3.2306   11.4029   16.8651   35.4439    0.0148
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
            T2Star=mean(T2Star_allsub,1);
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




