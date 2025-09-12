
addpath(genpath('/ptmp/pvalsala/Packages/mapVBVD'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'));
addpath(genpath('/ptmp/pvalsala/Packages/pulseq'));
%% calc all T2

clear fitResults_T2_all;
for i=1:5
    MeasPath=sprintf('/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/sub-%02d',i);
    metabolites=getMetaboliteStruct('invivo');
process_T2_fmincon;
assert(strcmp(MeasPath,sprintf('/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/sub-%02d',i)))%check Measpath is unaltered
fitResults_T2_all{i}=fitResults_T2;
end

%% format for getmetstuct
T2_allsub=cell2mat(cellfun(@(x)[x.T2_ms, x.water2_T2_ms, x.water2_fac],fitResults_T2_all,'UniformOutput',false)')
T2_std_allsub=cell2mat(cellfun(@(x)[x.T2_std, x.water2_T2_std_ms, x.water2_fac_std],fitResults_T2_all,'UniformOutput',false)')
mean(T2_allsub)
std(T2_allsub)

%% phnatom
    MeasPath=sprintf('/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/phantom');
    metabolites=getMetaboliteStruct('phantom');
process_T2_fmincon;
assert(strcmp(MeasPath,'/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/phantom'))%check Measpath is unaltered

fitResults_T2_phantom=fitResults_T2;
%% T1
clear fitResults_T1_all;
for i=1:5
    MeasPath=sprintf('/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/sub-%02d',i);
    metabolites=getMetaboliteStruct('invivo');
process_T1_inv_fmincon;
assert(strcmp(MeasPath,sprintf('/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/sub-%02d',i)))%check Measpath is unaltered
fitResults_T1_all{i}=fitResults_T1;
end

%% format for getmetstuct
T1_allsub=cell2mat(cellfun(@(x)[x.T1_ms],fitResults_T1_all,'UniformOutput',false)')
T1_std_allsub=cell2mat(cellfun(@(x)[abs(x.T1_std)],fitResults_T1_all,'UniformOutput',false)')
mean(T1_allsub)
std(T1_allsub)

%% phnatom
    MeasPath=sprintf('/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/phantom');
    % MeasPath=pwd;
    metabolites=getMetaboliteStruct('phantom');
process_T1_inv_fmincon;
assert(strcmp(MeasPath,'/ptmp/pvalsala/deuterium/dataForPublication/Relaxometry/phantom'))%check Measpath is unaltered
fitResults_T1_phantom=fitResults_T1;

%%
save('fitresults_jul30_final.mat','fitResults_T1_phantom','fitResults_T2_phantom','fitResults_T1_all','fitResults_T2_all')
%% for table
clc
fprintf('%d ± %d\n',[round([mean(T1_allsub),fitResults_T1_phantom.T1_ms]); ...
                    round([std(T1_allsub),fitResults_T1_phantom.T1_std])])

%make water2 frac first
T2_ordered=[T2_allsub(:,6)*100,T2_allsub(:,[1 5 2 3 4]),];

fprintf('\n')

fprintf('%d ± %d\n',[round([mean(T2_ordered),fitResults_T2_phantom.T2_ms]); ...
                    round([std(T2_ordered),fitResults_T2_phantom.T2_std])])
