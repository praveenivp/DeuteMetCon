function plotfitoverview(mcobj,pxl)
            PhSel=mcobj.flags.PCSel;
            EchoSel=mcobj.flags.EchoSel;
            FA=mcobj.DMIPara.FlipAngle *0.9; %rad
            TE=mcobj.DMIPara.TE; %s
%             TE=TE-TE(1);
            TR=mcobj.DMIPara.TR; %s
            PC=mcobj.DMIPara.PhaseCycles;%-pi/2;%-2*pi; %rad
           
            %calc mask to average
            mask=zeros(size(mcobj.mask),'logical');
            mask(pxl(1),pxl(2),pxl(3))=true;
            mask=imdilate(mask,strel('sphere',2));
            fprintf('Averaging over %d voxels\n',sum(mask(:)))
%             as(mask)
            B0=mcobj.FieldMap(pxl(1),pxl(2),pxl(3))*(6.536 /42.567);%+1e6*(Spectroscopy_para.PVM_FrqWork(1)- ExpPara.PVM_FrqWork(1)); %Hz
            im1=reshape(mcobj.img(:,:,:,:,mcobj.flags.EchoSel,mcobj.flags.PCSel), ...
                [],length(mcobj.flags.EchoSel)*length(mcobj.flags.PCSel));
             im1=mean(im1(mask(:),:),1,'omitnan');
%                      [Msig_all]=bSSFP_sim_analytical(metabolites,TE(EchoSel),PC(PhSel),TR,zeros(size(B0)),FA);
            [Msig_all]=bSSFP_sim_analytical(mcobj.metabolites,TE(EchoSel),PC(PhSel),TR,B0,FA);
%               Msig_all=bsxfun(@times,Msig_all,exp(-1i*angle(Msig_all(:,:,1,1,:,:,:,:,:,:))));
%             zeros(1,length(Metabolites),length(TE),length(PhaseCyles),length(TR),length(dfreq),length(FA))
             metabol_con=zeros(size(im1,1),3);
            resi=zeros(size(im1));
            condnm=zeros(size(im1,1),1);
            
            for i=1:size(im1,1)
                A= padarray(reshape(squeeze(Msig_all(1,:,:,:,1,i)),length(mcobj.metabolites),length(PhSel)*length(EchoSel)),[0 0],1,'pre').';
                b=[ im1(i,:)];
                %     b=b./max(abs(b(:)));
                metabol_con(i,:)=A\b(:);
                resi(i,:)=b(:) -A*metabol_con(i,:).';
                condnm(i)=cond(A);
            end
            
disp(abs(metabol_con))
disp(angle(metabol_con))
fprintf('B0 is %.2f Hz\n',B0)




% plot
tt=tiledlayout(2,4);
clbar=lines(10);
%abs
nexttile(1)
hold on
imdata=reshape(im1,length(mcobj.flags.EchoSel),length(mcobj.flags.PCSel)).';
plot(PC,abs(imdata)),hold on
map = get(gca,'ColorOrder');
fit=A*metabol_con(1,:).';
set(gca,'ColorOrderIndex',1);
fit=reshape(fit,length(mcobj.flags.EchoSel),length(mcobj.flags.PCSel)).';
plot(PC,abs(fit),LineStyle="--"),hold on
% set(gca, 'ColorOrder',cat(1,map,map))

nexttile(1+4)
hold on
plot(PC,angle(imdata)),hold on
set(gca,'ColorOrderIndex',1);
fit=A*metabol_con(1,:).';
fit=reshape(fit,length(mcobj.flags.EchoSel),length(mcobj.flags.PCSel)).';
plot(PC,angle(fit),LineStyle="--"),hold on


for ii=1:3
%abs
nexttile(ii+1)
imdata=reshape(Msig_all(1,ii,:),length(mcobj.flags.EchoSel),length(mcobj.flags.PCSel)).';
plot(PC,abs(imdata)),hold on
% fit=imdata*metabol_con(1,ii).';
% fit=reshape(fit,length(mcobj.flags.EchoSel),length(mcobj.flags.PCSel)).';
% plot(PC,abs(fit),LineStyle="--"),hold on
title(mcobj.metabolites(ii).name)

nexttile(ii+4+1)
hold on

plot(PC,angle(imdata)),hold on

% fit=imdata*metabol_con(1,ii).';
% fit=reshape(fit,length(mcobj.flags.EchoSel),length(mcobj.flags.PCSel)).';
% plot(PC,angle(fit),LineStyle="--"),hold on

end

end