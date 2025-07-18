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
        %new phantom: /ptmp/pvalsala/deuterium/20240102_new2Hphantom
        met_name={'water','glucose','glutamic acid','lactate'};
        T1=[449.8212,60.3493,162.9081,241.4547]*1e-3; %s
        T1_CI=[47.9218,23.5754,34.9505,24.2403]*1e-3; %s diff(CI95)/2
        T2=[270.8520,56.7103,112.8461,247.8881]*1e-3; %s
        T2_CI=[7.1701,7.4427,12.0632,23.1362]*1e-3; %s diff(CI95)/2
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

        % data from process_T1_inv
        T1 = [
            360.5184   86.8361  158.0989  109.9415;
            417.2924   64.3759  173.4299  142.9516;
            361.3593   77.4438  157.5635  121.6826;
            351.2728   81.6462  157.6897  143.1148;
            360.1643   77.9841  161.7815  165.6381
            ] * 1e-3; % Convert to seconds

        T1_CI = [
            5.5403    9.4113    3.6082   18.5691;
            55.7675   4.2608   15.5232   40.3232;
            6.8110    4.4586    4.7049    9.6327;
            19.3575   7.5787    4.6537   23.5573;
            7.2214    1.9944    5.9712   46.3589
            ] * 1e-3; % Convert to seconds

        %data from process_T2.m
        T2 = [
            29.9663   40.4833   93.4319   50.1270  340.7986
            34.2895   38.2420  103.6166  248.7363  317.6638
            36.5687   37.2845   81.0161   32.0723  200.0000
            31.4127   35.2399   97.6876   69.8372  331.4707
            34.5786   32.2199  111.4885   65.1398  365.0476
            ] * 1e-3; % Convert to seconds

        T2_CI = [
            6.7563    5.2629   14.5070   15.3652  0;
            25.2735   3.4229   14.7962  134.0744  0;
            7.8778    6.9314   31.6491   14.2051  0;
            8.3531    4.2345   22.9319   31.9338  0;
            11.3065   5.7484   35.9891   54.9027  0
            ] * 1e-3; % Convert to seconds

        water2_percent=[
            10.9100         0;
            8.3700         0;
            8.7400         0;
            12.1800         0;
            12.2700         0;
            ]; % Water 2 compartment and standard deviation in %

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
            T1=mean(T1,1);
            T2=mean(T2(:,1:4),1);
            T2Star=mean(T2Star,2);
            freq_shift=mean(freq_shift-freq_shift(:,1),1);
        else
            T1=T1(subIdx,:);
            T2=T2(subIdx,1:4);
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




%% plotting mean
% all_sub={'invivo1','invivo2','invivo3','invivo4','invivo5'};
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
