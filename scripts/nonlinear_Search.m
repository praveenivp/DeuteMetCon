



%% linear method
metabolites=getMetaboliteStruct('phantom');
fn='/ptmp/pvalsala/deuterium/HOSJ-D6P2/TWIX/allData#S94Tuebingen#F16606#M1000#D230924#T095442#rpcsi_ssfp_Stan25_156mm.dat';
CSI_setting={'metabolites',metabolites,'doPhaseCorr','none','parfor',true,...
    'doCoilCombine','adapt1','doZeropad',[0.5 0.5 0.5 0],'mask',[]}; 

    mcobj_csi_ssfp=MetCon_CSI(fn,CSI_setting{:},'fm','IDEAL','Solver','pinv');

%%
    mask1=mcobj_csi_ssfp.getMask(80);

    TE=mcobj_csi_ssfp.DMIPara.TE;
    PC=mcobj_csi_ssfp.DMIPara.PhaseCycles;
       met=mcobj_csi_ssfp.metabolites;
       TR=mcobj_csi_ssfp.DMIPara.TR;
       FA=mcobj_csi_ssfp.DMIPara.FlipAngle*mcobj_csi_ssfp.DMIPara.pulseCorrectionFactor;



   im1= reshape(mcobj_csi_ssfp.img,[],length(TE)*length(PC)); 

   mask=zeros(size(mcobj_csi_ssfp.img,2:4),"logical");
   mask(:,:,30)=mask1(:,:,30); 
   im1=im1(mask(:),:);
    B0=mcobj_csi_ssfp.FieldMap(mask(:))./(2*pi)*(6.536 /42.567);
    
    mc0=mc0(mask(:),:);

    x_all=zeros(size(im1,1),7);
   for i=1:size(im1,1)
%        x=[FA,B0,Phaaseffset,amp]
   x0=[FA,B0(i),0,abs(mc0(i,:))];
   func=@(x) CostFucntion(x,double(im1(i,:).'),met,TE,PC,TR);
   lb=[FA*0.75,-40,-pi,abs(mc0(i,:))*0.1];
   ub=[FA*1.25,40,pi,abs(mc0(i,:))*2];

   x_all(i,:)=fmincon(func,x0,[],[],[],[],lb,ub);

   
   end
%%
   mc1=zeros(numel(mask),4);
   mc1(mask(:),:)=x_all(:,4:end);
   mc1=reshape(mc1,51,51,51,4);
   as(mc1)

    %%


    function Cost=CostFucntion(x,data,met,TE,PC,TR)

%     x=[FA,B0,Phaaseffset,amp]


    Msig_all=reshape(MetSignalModel(met,TE,PC,TR,x(2),x(1),'bSSFP'),length(met),[]);
    amp=x(4:end);
    Cost=(amp(:).'*(Msig_all.*exp(1i*x(3))));
    Cost=norm(Cost(:)-data,2);
    end