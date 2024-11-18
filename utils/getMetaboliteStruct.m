function metabolites=getMetaboliteStruct(DataSelect,freqOffset)
% metabolites=getMetaboliteStruct(DataSelect,freqOffset=0)
% simple function to get Metabolite structure with T1/T2 and chemical
% shifts
%
%Dataselect: one of {'phantom','invivo1','invivo2','Roig7T','Peters12T'}
%freqOffset - new water frequency in Hz

DataSelect=validatestring(DataSelect,{'phantom','invivo','invivo1','invivo2','invivo3','invivo4','Roig7T','Peters12T'});

B0=9.3879; %[T]
gammaH2=6.536; % [MHz/T]

switch(DataSelect)
    case 'phantom'
        %new phantom: /ptmp/pvalsala/deuterium/20240102_new2Hphantom
        met_name={'water','glucose','glutamic acid','lactate'};

        T1=[456.6750,60.3171,164.8149,244.6637]*1e-3; %s
        T1_CI=[50.5780,23.6860,35.0992,27.6961]*1e-3; %s diff(CI95)/2
        T2=[268.7716,57.1611,111.5127,221.9746]*1e-3; %ms
        T2_CI=[8.2778,8.8415,21.7709,30.2800]*1e-3; %ms diff(CI95)/2
        ImagingFreq=61359114.0000; %Hz

        %CSI dataset='20240813_spectral/meas_MID00857_FID14528_rpcsi_fid_Stan_res15_6_optimal.dat'; see plotT2star.m
        freq_shift=[-3.4080 -63.5005 -150.2146 -205.1276]; % Hz
        freq_shift=[+0.5697 -70.2368 -154.065 -213.966]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[77.9295   18.8752   14.0095   15.9714]*1e-3; % s

    case 'invivo' % average all invivo: see last session to generate it
        met_name={'water','glucose','Glx','lactate/lipid'};

        T1=[383.3281, 80.4467, 157.6869, 120.2125]*1e-3; %ms
        T1_std=[35.0646, 10.9546, 6.0578, 19.3594]*1e-3; %ms
        T2=[36.9348, 33.8720, 99.0049, 83.2087]*1e-3; %ms
        T2_std=[2.3520, 1.2332, 12.8110, 51.0891]*1e-3; %ms
        T2Star=[21.7471, 14.7950, 23.0429, 28.2032]*1e-3; %ms
        T2Star_std=[1.7666, 1.4873, 1.4726, 4.2597]*1e-3; %ms
        freq_shift=[0.0000, -54.9437, -144.3418, -204.3382]; %Hz
        freq_shift_std=[0.0000, 0.2991, 1.0427, 0.8900]; %Hz


        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
    case 'invivo1'
        % for invivo we have all values except lactate

        %invivo: /ptmp/pvalsala/deuterium/EAZW-GUMK/proc
        met_name={'water','glucose','Glx','lactate/lipid'};

        T1=[354.3195,81.4397,158.4985,119.8924]*1e-3; %s
        T1_CI=[17.8741,10.3948,6.4426,29.0646]*1e-3; %s diff(CI95)/2
        T2=[37.6223,33.3514,94.5104,55.7864]*1e-3; %ms
        T2_CI=[5.3364,3.7084,18.4937,16.6814]*1e-3; %ms diff(CI95)/2
        ImagingFreq=61359045.0000; %Hz
        %CSI dataset: M412 % see plotT2star.m
        freq_shift=[-0.6647  -55.7530 -144.9229 -205.0421]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[24.3812   17.0160   25.2516   23.8431]*1e-3; %s
    case 'invivo2'
        %invivo: /ptmp/pvalsala/deuterium/DA77-F3UY
        met_name={'water','glucose','Glx','lactate/lipid'};
        T1=[368.5119,85.1623,152.0021,136.2558]*1e-3; %s
        T1_CI=[14.1349,5.8384,7.8401,71.7550]*1e-3; %s diff(CI95)/2
        T2=[33.6801,35.7076,117.4483,60.3702]*1e-3; %s
        T2_CI=[5.5088,16.8677,42.7416,59.7657]*1e-3; %s diff(CI95)/2
        freq_shift=[3.8518,-51.4442,-141.9904,-201.7315]; %Hz median
        ImagingFreq=61359063.0000; %Hz
        %CSI dataset: M695 % see plotT2star.m
        %         freq_shift=[9.0229 -46.9185 -135.4832 -198.4703]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[20.5969 14.2529 22.2796 25.2529]*1e-3; %s

    case 'invivo3'
        %/ptmp/pvalsala/deuterium/HOSJ-D6P2/proc/T1T2
        met_name={'water','glucose','Glx','lactate/lipid'};
        T1=[376.3679,90.2514,154.4121,93.0732]*1e-3; %s
        T1_CI=[16.5739,8.0609,5.2190,23.5062]*1e-3; %s diff(CI95)/2
        T2=[37.1612,33.3840,87.8637,56.8920]*1e-3; %s
        T2_CI=[5.8750,4.8476,18.6411,16.0938]*1e-3; %s diff(CI95)/2
        ImagingFreq=61359053.0000; %Hz
        %CSI dataset: M997 TR36 % see plotT2star.m
        freq_shift=[1.5423 -53.1530 -142.0912 -202.1538]; % Hz
        freq_shift_ppm=(freq_shift-freq_shift(1))./(B0*gammaH2)+4.7;
        T2Star=[21.0052 13.9556 22.3202 31.8583]*1e-3; %s

    case 'invivo4'
        %/ptmp/pvalsala/deuterium/H4
        met_name={'water','glucose','Glx','lactate/lipid'};
        T1=[434.1131,64.9336,165.8349,131.6284]*1e-3; %s
        T1_CI=[59.1279,5.0090,16.6350,39.3717]*1e-3; %s diff(CI95)/2
        T2=[39.2756,33.0449,96.1974,159.7863]*1e-3; %s
        T2_CI=[6.5651,2.7752,18.2076,64.7964]*1e-3; %s diff(CI95)/2
        ImagingFreq=61359047.0000; %Hz
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

% all_sub={'invivo1','invivo2','invivo3','invivo4'};
% metabolites_all=cellfun(@(x) getMetaboliteStruct(x,0),all_sub);
%
% fprintf('T1=[%.4f, %.4f, %.4f, %.4f]*1e3; %%ms \n',...
%     1e3*mean(cell2mat(cellfun (@(x)[x.T1_s],metabolites_all,'UniformOutput',false)')));
% fprintf('T1_std=[%.4f, %.4f, %.4f, %.4f]*1e3; %%ms \n',...
%     1e3*std(cell2mat(cellfun (@(x)[x.T1_s],metabolites_all,'UniformOutput',false)')));
% fprintf('T2=[%.4f, %.4f, %.4f, %.4f]*1e3; %%ms \n',...
%     1e3*mean(cell2mat(cellfun (@(x)[x.T2_s],metabolites_all,'UniformOutput',false)')));
% fprintf('T2_std=[%.4f, %.4f, %.4f, %.4f]*1e3; %%ms \n',...
%     1e3*std(cell2mat(cellfun (@(x)[x.T2_s],metabolites_all,'UniformOutput',false)')));
% 
% fprintf('T2star=[%.4f, %.4f, %.4f, %.4f]*1e3; %%ms \n',...
%     1e3*mean(cell2mat(cellfun (@(x)[x.T2star_s],metabolites_all,'UniformOutput',false)')));
% fprintf('T2star_std=[%.4f, %.4f, %.4f, %.4f]*1e3; %%ms \n',...
%     1e3*std(cell2mat(cellfun (@(x)[x.T2star_s],metabolites_all,'UniformOutput',false)')));
% 
% fprintf('freq_shift=[%.4f, %.4f, %.4f, %.4f]; %%Hz \n',...
%     mean(cell2mat(cellfun (@(x)[x.freq_shift_Hz],metabolites_all,'UniformOutput',false)')));
% fprintf('freq_shift_std=[%.4f, %.4f, %.4f, %.4f]; %%Hz \n',...
%     std(cell2mat(cellfun (@(x)[x.freq_shift_Hz],metabolites_all,'UniformOutput',false)')));
