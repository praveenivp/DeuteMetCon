% IDEAsolver demo
%% simulation

metabolites=[struct('freq_shift_Hz',7,'con',randi([2,10],1)/10,'name','1.3 ppm','T1_s',297e-3,'T2_s',61e-3),...
             struct('freq_shift_Hz',100,'con',randi([5,10],1)/10,'name','3.7 ppm ','T1_s',64e-3,'T2_s',32e-3),...
                struct('freq_shift_Hz',200,'con',randi([8,10],1)/10,'name','4.8 ppm','T1_s',320e-3,'T2_s',12e-3)];
 TE_s=[3.2:3:14]*1e-3; %s
% TE_s=[1:5]*1e-3; %s
MatSz=[64 64]/2;
%  fm_Hz= 40*getRandomFmap([MatSz 1]);

%% simulate some signal
 IDEALobj=IDEALsolver(metabolites,TE_s);
 sig=IDEALobj.simSig([32; 32;],20,0.01);
 fm_Hz=IDEALobj.FieldMap_Hz;
 metbol_mask=IDEALobj.experimental.Phantom;
%% do IDEAL and ordinary least squares 
 
 IDEALobj_pinv=IDEALsolver(metabolites,TE_s,0*fm_Hz,'pinv',0); %OLS
 metabol_con_pinv=IDEALobj_pinv'*sig;
 
  IDEALobj_ideal=IDEALsolver(metabolites,TE_s,0*fm_Hz,'IDEAL',300);
 metabol_con_ideal=IDEALobj_ideal'*sig;
 
 IDEALobj_ideal.poltFrame(metabol_con_ideal,metabol_con_pinv,metbol_mask,fm_Hz,IDEALobj_ideal.experimental.fm_est)
 
 %% real data test
%   sig=twix.image('');
 sigm=sum(sig,[6 ]);


 zp_PRS=[1 1 1]*1.;
pad_size=[round(size(sigm,1)*zp_PRS(1)) 0 round(size(sigm,3)*zp_PRS(2))  round(size(sigm,4)*zp_PRS(3))];
sigm_zp=padarray(sigm,pad_size,0,'both');

im=squeeze(myfft(sigm_zp,[1 3 4]));
[~,CoilWeights,NormalizationMAtrix]=adaptiveCombine(permute(sum(im(:,:,:,:,1,:),6), [2 1 3 4 5 6]));
imc_zp=squeeze(sum(CoilWeights.*permute(im, [2 1 3 4 5 6]),1));
im_phasecorr_zp=bsxfun(@times,imc_zp,exp(-1i*angle(imc_zp(:,:,:,1,1))));
 
 
 %
 corr2=permute(exp(0.5i*deg2rad(phi_vec(:)-phi_vec(1))),[2 3 4 5 1]);
   cs_ppm=[0 1.02 2.56 ]*1e-6; %D20 Glu Glx Lac/fat
freq_shift_WGX=[8 -72 -152]*1;%Hz
% freq_shift_AGW=[1.3 3.7 4.8]*ExpPara.PVM_FrqWork(1);%Hz
T1=[320e-3, 64e-3,297e-3 ];%s
T2=[200e-3, 32e-3,61e-3 ];%s
sp_name={'D20','Glucose','Glutamate'};
clear metabolites;
for i=1:3 
    metabolites2(i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift_WGX(i),'name',sp_name{i});
end
 fm_meas_Hz=imresize3(fieldmap,size(im_phasecorr_zp,1)/size(fieldmap,1))./(2*pi)*(6.536 /42.567);
 im_me=sum(im_phasecorr_zp(:,:,:,:,:).*corr2,5); % sum phase cycles
  IDEALobj2=IDEALsolver(metabolites2,TE,fm_meas_Hz.*0,'IDEAL',500);
  metabol_con2=IDEALobj2'*im_me;
  
 %
metabol_con1=permute(metabol_con2(:,:,:,1:3),[2 3 1 4 ]);
figure(25),clf,
for jj=1:3
subplot(3,1,jj)
imagesc(createImMontage(reshape(abs(metabol_con1(:,:,:,jj)),size(metabol_con1,1), size(metabol_con1,2),[]),0.5*size(metabol_con1,3)))
axis image,xlabel('Slices(4-11)'),ylabel('Metabolites'),
% yticks([8:16:16*3]),
title(metabolites2(jj).name)
% title(sprintf('M%d Deterium @%s:%s (Glucose Intake: 09:30)',twix.hdr.Config.MeasUID,acq_time(1:2),acq_time(4:5)))
colorbar
end
