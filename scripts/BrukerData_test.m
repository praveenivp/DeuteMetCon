%% plot simulated and measured bSSFP phase cycle profiles in a single ROI
pn='/ptmp/pvalsala/deuterium/20220413_111704_Deuterium_SSFPtest_1_10';

addpath('/ptmp/pvalsala/deuterium/20220413_111704_Deuterium_SSFPtest_1_10/proc/Bruker_reader_Rolf')
addpath(genpath('/ptmp/pvalsala/MATLAB'))
addpath(genpath('/ptmp/pvalsala/Packages/DeuteMetCon'));
addpath(genpath('/ptmp/pvalsala/Packages/OXSA'))

%% CSI sequence
csi_para=BrReadParams({'PVM_FovCm','PVM_SpatResol','PVM_NAverages','PVM_RepetitionTime', ...
    'NumEchoes','PVM_EchoTime','EchoSpacing','PVM_Fov','PVM_FrqWork','PVM_EncMatrix', ...
    'PVM_SliceThick','PVM_EffSWh','PVM_SPackArrReadOffset','PVM_SPackArrPhase1Offset', ...
    'PVM_SPackArrPhase2Offset','PVM_SPackArrSliceOffset','PVM_SpecDwellTime','PVM_FrqWork'},...
         fullfile(pn,'77','method'));

params.complex = 1;
params.filter=[0,0,0,0];  % filter: 1: Hanning filter, 2: Hamming filter, 3: optimized generalized Hamming filter, 4: gaussian filter, 5: exponential filter
params.phase = 0;
params.zf = [0,0,0,0];
params.ft = [0,1,1,1];

[csi_fid]=BrReadFid(fullfile(pn,'77','rawdata.job0'));

csi_fid=reshape(csi_fid,[],csi_para.PVM_NAverages,csi_para.PVM_EncMatrix(1),csi_para.PVM_EncMatrix(2));
csi_fid=permute(csi_fid,[3,4,10,1,2,5:9] );

csi_fid=padarray(csi_fid,size(csi_fid).*[1,1,0,0,0]*0.5,0,'both');

csi_spec=(specFft(myfft(sum(csi_fid,5),1:3),4));

figure,
dt=2*(csi_para.PVM_SpecDwellTime*1e-6);
faxis=linspace(-0.5/dt,0.5/dt-1/(dt*size(csi_spec,4)),size(csi_spec,4));
plot(faxis,reshape(squeeze(abs(csi_spec(:,:,1,:,1))),prod(size(csi_spec,1:3)),[])')
% plot(faxis,squeeze(abs(csi_spec(1,16,10,1,:,1))))

%%
timeSel=7:size(csi_fid,4);
offset=csi_para.PVM_FrqWork(1);
met_struct=[struct("T1_s",0.4498212,"T2_s",0.270852,"freq_shift_Hz",-82,"name","acetate","T2star_s",0);...
    struct("T1_s",0.4498212,"T2_s",0.270852,"freq_shift_Hz",71,"name","glucose","T2star_s",0);...
    struct("T1_s",0.4498212,"T2_s",0.270852,"freq_shift_Hz",183,"name","water","T2star_s",0);];


TE_s=0:dt:dt*(size(csi_fid,4)-1)+csi_para.PVM_EchoTime*1e-3;
IDEAL_obj=IDEAL(met_struct,TE_s(timeSel),'solver','IDEAL','fm',[],'maxit',5,'SmoothFM',2);
dat=myfft(sum(csi_fid,5),1:3);

met=IDEAL_obj'*dat(:,:,:,timeSel);
as(met)
%% read b0map and echo image
% load('\\mrz10\MRZ14T\AGKS\rolf\Deuterium\20220413_111704_Deuterium_SSFPtest_1_10\B0map.mat')
% B0map=ndflip(permute(B0map,[2 1 3]),[1,2]); % check next session for orientaion

% read echo images
B0Map_para=BrReadParams({'PVM_FovCm','PVM_SpatResol','PVM_RepetitionTime','NumEchoes','PVM_EchoTime','EchoSpacing','PVM_Fov','PVM_FrqWork','PVM_EncMatrix','PVM_SliceThick','PVM_EffSWh','PVM_SPackArrReadOffset','PVM_SPackArrPhase1Offset','PVM_SPackArrPhase2Offset','PVM_SPackArrSliceOffset'},...
         fullfile(pn,'106','method'));
params.complex = 1;
params.filter=[0,0,0,0];  % filter: 1: Hanning filter, 2: Hamming filter, 3: optimized generalized Hamming filter, 4: gaussian filter, 5: exponential filter
params.phase = 0;
params.zf = [0,0,0,0];
params.ft = [1,1,1,1];

[echo_im,acqp_b0map]=BrReadImage(fullfile(pn,'106','rawdata.job0'),params);

%% imaging scan 16 phase cycles
params.complex = 1;
params.filter=[0,0,0,0];  % filter: 1: Hanning filter, 2: Hamming filter, 3: optimized generalized Hamming filter, 4: gaussian filter, 5: exponential filter
params.phase = 0;
params.zf = [0,0,0,0];
params.ft = [1,1,1,1];

[im,acqp]=BrReadImage(fullfile(pn,'114','rawdata.job0'),params);
[raw_data]=BrReadFid(fullfile(pn,'114','rawdata.job0'));
ExpPara = BrReadParams({'ExcPulse1','PVM_RepetitionTime','PVM_SpatResol','NumEchoes','PVM_EchoTime','EchoSpacing','PhaseAdvance','NumPhaseAdvance','PVM_Fov','PVM_FrqWork','PVM_EncMatrix','PVM_SliceThick','PVM_EffSWh','PVM_SPackArrReadOffset','PVM_SPackArrPhase1Offset','PVM_SPackArrPhase2Offset','PVM_SPackArrSliceOffset'},...
         fullfile(pn,'114','method'));
raw_data=reshape(raw_data,[ExpPara.PVM_EncMatrix(1),ExpPara.NumEchoes,ExpPara.PVM_EncMatrix(1),ExpPara.NumPhaseAdvance]);
raw_data=permute(raw_data,[1,3,5,2,4]);

raw_data=padarray(raw_data,size(raw_data).*[1,1,0,0,0]*0.5,0,'both');

im_me=myfft(sum(raw_data,50),1:3);

%%
gammaH2=6.536 ; % [MHz/T]
TE_s=ExpPara.EchoSpacing*(-floor(ExpPara.NumEchoes/2):floor(ExpPara.NumEchoes/2))+ExpPara.PVM_EchoTime;
TR=ExpPara.PVM_RepetitionTime;
offset= 3.47641366355123*gammaH2*14*0;%csi_para.PVM_FrqWork(1)*2;
% met_struct=[struct("T1_s",0.4498212,"T2_s",0.270852,"freq_shift_Hz",1.8*14*gammaH2-offset,"name","acetate","T2star_s",0);...
%     struct("T1_s",0.4498212,"T2_s",0.270852,"freq_shift_Hz",3.2*14*gammaH2-offset,"name","glucose","T2star_s",0);...
%     struct("T1_s",0.4498212,"T2_s",0.270852,"freq_shift_Hz",4.8*14*gammaH2-offset,"name","water","T2star_s",0);];


% (chemical shifts: deuterated HDO at 4.8 ppm, glucose at 3.8 ppm, Glx at 2.4 ppm, and acetate at 1.9 ...

IDEAL_obj=IDEAL(met_struct,TE_s,'solver','IDEAL','fm',[],'maxit',10,'SmoothFM',1);

met=IDEAL_obj'*sum(im_me,5);
as(met)




%%
[im,acqp]=BrReadImage(fullfile(pn,'','rawdata.job0'),params);
[im1]=BrReadFid(fullfile(pn,'114','rawdata.job0'));


%% plot everything for consitency check

figure(6),clf,
subplot(141),
imshow(ndflip(permute(imread('C:\Users\pvalsala\Documents\MATLAB\deuterium\20220413_phantomexp\phantom_setup.png'),[2 1 3]),1))
title('Phantom setup') 

subplot(142),imagesc(abs(echo_im(:,:,24,1))),title('1H Echo image')
annotation(gcf,'rectangle',...
    [0.348916666666667 0.358264081255771 0.108895833333333 0.316712834718373],'LineWidth',2,'Color','r');
 
subplot(143),imagesc(B0map(:,:,24)),caxis([-100 100]),title('B0 map(Hz)'),colorbar
 

 %B0map voxel locations
 [xgrid,ygrid,zgrid]=meshgrid(linspace(-0.5*B0Map_para.PVM_Fov(2),0.5*B0Map_para.PVM_Fov(2)-B0Map_para.PVM_SpatResol(2),B0Map_para.PVM_EncMatrix(2)),...
     linspace(-0.5*B0Map_para.PVM_Fov(1),0.5*B0Map_para.PVM_Fov(1)-B0Map_para.PVM_SpatResol(1),B0Map_para.PVM_EncMatrix(1)),...
     linspace(-0.5*B0Map_para.PVM_Fov(3),0.5*B0Map_para.PVM_Fov(3)-B0Map_para.PVM_SpatResol(3),B0Map_para.PVM_EncMatrix(3)));
 % 2D acquistion  voxel locations
 [im_x,im_y,im_z]=meshgrid(linspace(-0.5*ExpPara.PVM_Fov(2),0.5*ExpPara.PVM_Fov(2)-ExpPara.PVM_SpatResol(2),ExpPara.PVM_EncMatrix(2)),...
     linspace(-0.5*ExpPara.PVM_Fov(1),0.5*ExpPara.PVM_Fov(1)-ExpPara.PVM_SpatResol(1),ExpPara.PVM_EncMatrix(1)),...
     ExpPara.PVM_SPackArrSliceOffset);
 % register B0 map
 B0map_crop=interp3(xgrid,ygrid,zgrid,B0map,im_x,im_y,im_z);
 
 %if registeration works like that then half of MR community will be jobless
 B0map_crop=ndCircShift(B0map_crop,[1 -2],[1 2]);
%  as(cat(3,B0map_crop,sum(abs(im),5)))
 
 subplot(244),imagesc(B0map_crop),caxis([-100 100]),title('B0 map registered(Hz)')
 subplot(248),imagesc(sum(abs(im),5)),title('2D H2 image(all PC summed)')


%% create ROIs and plot the absoulte bSSFP profile
load all_masks.mat
%  figure,
% for i=1:6 
%  mask{i}=CreateMask(sum(abs(im(:,:,:)),3),'circle');
%  end

tube_contents={'D20','ACE','GLC+D20','0.5ACE+0.5D20','D20','1/3ACE+2/3D20'};
PC=deg2rad(ExpPara.PhaseAdvance); %deg passband? 
sp_list=[1 3 4 6 7 9];
figure
for i=1:6
bSSFP_profile=(squeeze(sum((im.*mask{i}),[1 2])));
subplot(3,3,sp_list(i)),plotSSFP(rad2deg(PC),bSSFP_profile),
title(tube_contents{i})
end
subplot(3,3,2 ),imagesc(sum(abs(im(:,:,:)),3))
subplot(3,3,5 ),imagesc(sum(abs(im(:,:,:)),3).*sum(cat(3,mask{:}),3))
 subplot(3,3,8),imagesc(sum(cat(3,mask{:}),3))

 %% plot all localized spectroscopy
sp_list=[1 3 4 6 7 9 8];
start_time=cell(7,2);
spec_list=[ 95    94       91    89  93    92  33;...
113 112 109 108 111 110 33]'; %113 16:24:12.447; % 89 : 15:02:44
figure(2)

subplot(3,3,2 ),imagesc(sum(abs(im(:,:,:)),3)),title('sample H2 image')
subplot(3,3,5 ),
% imagesc(sum(abs(im(:,:,:)),3).*sum(cat(3,mask{:}),3))
imshow(ndflip(permute(imread('C:\Users\pvalsala\Documents\MATLAB\deuterium\20220413_phantomexp\phantom_setup2.png'),[1 2 3]),10))
title('Phantom setup') 
%  subplot(3,3,8),imagesc(sum(cat(3,mask{:}),3))
for i=1:7
for j=1:2
specdata=BrReadFid(fullfile(pn,num2str(spec_list(i,j)),'rawdata.job0'));

SW_h=BrReadParams('SW_h', fullfile(pn,num2str(spec_list(i,j)),'acqp'));
Dwell_s=1/SW_h;

%get date&time
q=textscan(fopen(fullfile(pn,num2str(spec_list(i,j)),'acqp'),'r'),'%s',1,'HeaderLines',7-1,'Delimiter','\n');
start_time{i,j}=datetime(q{1}{1}(4:22),'Format','HH:mm:ss');

Spectroscopy_para = BrReadParams({'ExcPulse1','PVM_RepetitionTime','PVM_FrqWork','PVM_EffSWh','PVM_NAverages'},...
    fullfile(pn,num2str(spec_list(i,j)),'method'));
 specdata=reshape(specdata,[],Spectroscopy_para.PVM_NAverages);

subplot(3,3,sp_list(i)),hold on,plot(linspace(-1*SW_h/2,SW_h/2,size(specdata,1)),abs(mean(fftshift(fft(specdata,[],1),1),2))),
xlim([-350 350])
xlabel('Freq (Hz)')
title(strcat(num2str(spec_list(i,1)),',',num2str(spec_list(i,2))))
disp(Spectroscopy_para.PVM_FrqWork(1))
box on
end

legend(sprintf('%s',start_time{i,1}) ,sprintf('%s',start_time{i,2}))
end
title('Global spectrum(33)')
%153.96 Hz
 [start_time{:,2}]-[start_time{:,1}]
%% create ROIs and plot the absoulte bSSFP profile
% figure,
clear b0_ROI_mean b0_ROI_std
  for i=1:6 
%  mask_b0map{i}=CreateMask(B0map_crop,'circle');
b0_ROI_mean(i)=mean(B0map_crop(mask{i}));
b0_ROI_std(i)=std(B0map_crop(mask{i}));
  end

 %% plot simulted profile
freq_shift_AGW=[-116.27, 48,+155.4252 ]+21/2;%Hz
% freq_shift_AGW=[1.3 3.7 4.8]*ExpPara.PVM_FrqWork(1);%Hz
T1=[627e-3, 64e-3,320e-3 ];%s
T2=[627e-3, 32e-3,12e-3 ];%s
sp_name={'Acetate','Glucose','D20'};
clear metabolites;
for i=1:3 
    metabolites(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift_AGW(i),'name',sp_name{i});
end
EchoSel=1:ExpPara.NumEchoes;
PhSel=1:ExpPara.NumPhaseAdvance;


TE=(ExpPara.PVM_EchoTime+(-1*floor(ExpPara.NumEchoes/2):floor(ExpPara.NumEchoes/2))*ExpPara.EchoSpacing)*1e-3; %s
TR=ExpPara.PVM_RepetitionTime*1e-3; %s
 PC=deg2rad(ExpPara.PhaseAdvance(PhSel));
%  PC=deg2rad(0:1:359);
FA=deg2rad(ExpPara.ExcPulse1.Flipangle); %rad
B0=-1*b0_ROI_mean*(6.536 /42.567)+1e6*(Spectroscopy_para.PVM_FrqWork(1)- ExpPara.PVM_FrqWork(1)); %Hz
  [Msig_all]=bSSFP_sim_analytical(metabolites,TE(EchoSel),-PC+deg2rad(22),TR,B0,FA);
%   [Msig_all,Mss_all]=bSSFP_sim(metabolites,TE(EchoSel),PC,TR,B0,FA);
 figure(6),clf
  for i=1:6
 subplot(3,2,i),plotSSFP(rad2deg(PC),squeeze((Msig_all(1,:,1,:,1,i))).')
 title(sprintf('B0 (mean: %.1f/ std: %.1f)',b0_ROI_mean(i),b0_ROI_std(i)))
 end
legend(sp_name)

%% try ROI Fitting
  imc=im;
%  imc=bsxfun(@times,im,exp(1i*angle(sum(im,5))));

figure(81),clf
clear metabol_con;
for i=1:6
%   A= padarray(squeeze(Msig_all(1,:,1,:,1,i)),[0 1],1,'pre').';
   A=[ squeeze(Msig_all(1,:,1,:,1,i))].';
% A=[ squeeze(Msig_all(1,:,1,:,1,i))-mean(squeeze(Msig_all(1,:,1,:,1,i)),2)].';
bSSFP_profile=(squeeze(sum((imc.*mask{i})./sum(mask{i},[1 2]),[1 2])));
bSSFP_profile=bSSFP_profile./(max(abs(bSSFP_profile))-0*min(abs(bSSFP_profile)));
%phase correction
% phaseterm=sum(im(:,:,:,:,:),5)./abs(sum(im(:,:,:,:,:),5));
% bSSFP_profile=bSSFP_profile.*conj(mean(phaseterm(mask{i})));
% disp(rad2deg(angle(mean(phaseterm(mask{i})))))


b=[ bSSFP_profile.'];
 metabol_con(i,:)=A\b(:);


subplot(6,3,1+(i-1)*3),plot(rad2deg(PC),abs(squeeze(Msig_all(1,:,1,:,1,i))).');
title(sprintf('B0 (mean: %.1f/ std: %.1f)',b0_ROI_mean(i).*(6.536 /42.567),b0_ROI_std(i).*(6.536 /42.567)))
if(i==6);legend(sp_name); end

subplot(6,3,2+(i-1)*3),plot(rad2deg(PC),abs(bSSFP_profile)),hold on
plot(rad2deg(PC),abs(sum(A(end-15:end,:).*metabol_con(i,:),2))),

if(i==6); legend('meas','fit'); xlabel('deg');end
title(tube_contents{i})
% subplot(6,3,3+(i-1)*3),plot(rad2deg(PC),abs(A(end-15:end,:).*metabol_con(i,:)))
% title('Individual components')
subplot(6,3,3+(i-1)*3),plot(rad2deg(PC),angle(bSSFP_profile)),hold on
plot(rad2deg(PC),angle(sum(A(end-15:end,:).*metabol_con(i,:),2))),
title('phase(meas/fit)'),if(i==6); legend('meas','fit'); xlabel('deg');end
end



%% try pixel level 
B0=-1*B0map_crop(:)*(6.536 /42.567)-1e6*(Spectroscopy_para.PVM_FrqWork(1)- ExpPara.PVM_FrqWork(1)); %Hz

ims=imfilter(im,ones(3));

im1=reshape(imc,[],ExpPara.NumPhaseAdvance);
     [Msig_all]=bSSFP_sim_analytical(metabolites,TE(EchoSel),-PC,TR,B0,FA);
Msig_all_norm=zeros(size(Msig_all));
%try normalizing
for i=1:3
%     for j=1:size(Msig_all,6)
% Msig_all_norm(1,i,1,:,1,:)=Msig_all(1,i,1,:,1,:)./max(col(abs(Msig_all(1,i,1,:,1,:))));
        Msig_all_norm(1,i,1,:,1,:)=Msig_all(1,i,1,:,1,:);%./norm(squeeze(Msig_all(1,i,1,:,1,:)),2);
% %     end
end

 metabol_con=zeros(size(im1,1),3);
for i=1:size(im1,1)
A= padarray(squeeze(Msig_all_norm(1,:,1,:,1,i)),[0 0],1,'pre').';
% A=A./max(abs(A));
b=[im1(i,:)]./max(abs(im(:)));
% b=b./norm(b,2);
   metabol_con(i,:)=A\b(:);
% %  metabol_con(i,:)=abs(A)\abs(b(:));
%      metabol_con(i,:)=pinv(real(A'*A))*real(A'*b(:));
end
metabol_con1=reshape(metabol_con,size(im,1),size(im,2),size(im,3),[]);

%plot
figure,
for i=1:3
    r=[2 2 10]
    subplot(1,3,i),imagesc(abs((metabol_con1(:,:,1,i)))),colorbar,title(sp_name{i}),caxis([0 r(i)])
    colormap('jet')
end


%%
for i=1:size(im1,1)
A= padarray(squeeze(Msig_all(1,:,1,:,1,i)),[0 1],1,'pre').';

proj(i,:)=A*metabol_con(i,:)';
%  metabol_con(i,:)=abs(A)\abs(b(:));
%     metabol_con(i,:)=pinv(real(A'*A))*real(A'*b(:));
end
 as(reshape(proj(:,2:16),48,48,[]))
as(squeeze(im))


%%
figure,
plotSSFP(rad2deg(PC),squeeze(Msig_all(1,1,1,:,1,i)).')
%%


function plotSSFP(PC,prof)
plot(PC,abs(prof)),
% hold on,yyaxis('right'),plot(PC,rad2deg(angle(prof)))
end

function plotComplex(x,y)
figure,subplot(211),plot(rad2deg(x),abs(y)),title('abs')
hold on,plot(rad2deg(x),abs(sum(y,2)))
subplot(212),plot(rad2deg(x),angle(y)),title('phase'),hold on,plot(rad2deg(PC),angle(sum(A,2)))
end