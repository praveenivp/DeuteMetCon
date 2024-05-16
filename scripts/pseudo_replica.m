%% load data and         

MeasPath='/ptmp/pvalsala/deuterium/20240424_phamtomorient';
sn=fullfile(MeasPath,'TWIX');
pn=fullfile(MeasPath,'proc/CSI','data_denoised');

 fn='*rpcsi*.dat';
MyFolderInfo = dir(fullfile(sn,fn));

met_name={'water','glucose','Glx','lactate'};
freq_shift_WGX=[3.9 -58.7  -152 -200];
T1=[432.88 69.65 147.6 190.64]*1e-3;%s
T2=[287 65.88 124.5 180]*1e-3;%s


clear metabolites;
for i=1:length(freq_shift_WGX)
    metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift_WGX(i),'name',met_name{i});
end

mcobj.metabolites=metabolites;

CSI_filename=fullfile(sn,MyFolderInfo(i).name);
        twix=mapVBVD_CSI(CSI_filename,'rmos');
        vecSize=twix{1}.hdr.MeasYaps.sSpecPara.lVectorSize;
        mcobj=MetCon_CSI(twix,'ZeroPadSize',[1 1 1 0]*0.5,'Solver','IDEAL', ...
            'doPhaseCorr','none','doDenosing',0,'phaseoffset',[0 0],'doCoilCombine','wsvd','metabolites',metabolites);%, ...
        mcobj.mask=ones(size(mcobj.mask),'logical');

        mcobj.performMetcon();


%% pseudo_replica test


data=squeeze(mcobj.img);
metcon_baseline=mcobj.Metcon;

Nrep=128;

metcon_all=zeros([size(metcon_baseline) Nrep]);
solverobj=mcobj.SolverObj; % only works for IDEAL and phaseonly modes
parfor ii=1:Nrep

noise_white = complex(randn(size(data)),randn(size(data))); 
noisy_data=data+noise_white;

% mcobj.img=noisy_data;
% mcobj.mask=ones(size(mcobj.mask),'logical');
% mcobj.performMetcon();
% metcon_all(:,:,:,:,ii)=mcobj.Metcon;

metcon_all(:,:,:,:,ii)=solverobj'*noisy_data;
end


rep_dim = ndims(metcon_all);

g = std(abs(metcon_all + max(abs(metcon_all),[],rep_dim)),[],rep_dim); %Measure variation, but add offset to create "high snr condition"
g(g < eps) = 1;
% for ii=1:4
% g(:,:,:,ii)=imgaussfilt3(g(:,:,:,ii),1);
% end
snr = mean(metcon_all,ndims(metcon_all))./g;