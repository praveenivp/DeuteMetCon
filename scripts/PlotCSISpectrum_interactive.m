%% figure to plot CSI spectrum
spec=abs(squeeze(fftshift(fft(mc.img,[],5),5)));
% TODO : perform nice phase correction
faxis=calcFreqAxis(mc.DMIPara.dwell,size(spec,4));
anat=sum(spec,4);
%% plot profiles

voxel=[10,10,30];
figure(25),clf,
set(gcf,'color','w','Position',[383 197 1312 918])
tt=tiledlayout(3,4,"TileSpacing","compact","Padding","compact");
im_h=nexttile(tt,6,[2 2 ]);
%  imagesc(abs(im_pc1(:,:,voxel(3),1))),colormap(gca,'gray'),axis image,%clim([0 300])
imagesc(abs(anat(:,:,voxel(3),1))),colormap(gca,'gray'),axis image,%clim([0 300])
title('Anatomy')


tile_Nr=[9 5 1 2 3 4  8  12];
color_all=[lines(7);0 0.5 0.5];
clear voxel_all;
% load('voxel_all.mat')


for i=1:length(tile_Nr)


  pt_handle=drawpoint(im_h,"color",color_all(i,:));
   voxel=round([pt_handle.Position(2) pt_handle.Position(1) voxel(3)]);
 voxel_all{i}=voxel;
% voxel=voxel_all{i};
 voxel_mask=zeros(size(anat));
 voxel_mask(voxel(1),voxel(2),voxel(3))=1;
 voxel_mask=imdilate(voxel_mask,strel('disk',6,0));
voxel_mask=voxel_mask>0;
 


rectangle(im_h,'Position',[[voxel(2) voxel(1)]-[1 1]  2 2],'EdgeColor',color_all(i,:),'LineWidth',5)

%
spec_curr=squeeze(spec(voxel(1),voxel(2),voxel(3),:));
%averaging
% im_pc1_col=reshape(im_pc1,[],216);
% prof=mean(im_pc1_col(voxel_mask(:),:).',2);




%centereing
%  [~,minidx]=min(smooth(abs(prof)));
%  prof=circshift(prof(22:end),-1*minidx+97,1);

% 
% prof=myfft(prof,1);
% faxis=linspace(-0.5/1.68,0.5/1.68,216);
% pc_deg=(0:400/size(im_pc1,4):400-400/size(im_pc1,4));
  ax_handle=nexttile(tt,tile_Nr(i));
% ax_handle=nexttile(tt,1,[1 4]);

% plot(pc_deg   ,abs(prof),'LineWidth',1.9);
grid on
hold on
plot(faxis,spec_curr,'LineWidth',1.9);
box on
% title_str={'gray matter','white matter','muscle','skin','CSF','thalamus','white matter 2','sinus'};

xlabel('frequency [Hz]'),ylabel(' amplitude [a.u]')
title('spectrum')
%  set(gca,'ColorOrder',linspecer(size(H0field,2)));
xlim([-200 200])
 set(ax_handle,'XColor',color_all(i,:),'YColor',color_all(i,:),'LineWidth',1.5)


end


function ndMat=ndflip(ndMat,dim)
for i=dim
    ndMat=flip(ndMat,i);
end

end

function [im]=myfft(kdata1,dim,dofftshift,sz)
% [im]=myfft(kdata1,dim,shift,sz)
% dim is array 
% dofftshift is boolean(fftshifgt or not)
% sz final size(tail zeropadding)
%praveenivp
switch (nargin)
    case 1
        dim=1;
        dofftshift=true;
        sz=size(kdata1);
    case 2
        dofftshift=true;
        sz=size(kdata1);
    case 3      
        sz=size(kdata1);
end

im=kdata1;
if(dofftshift)
    for i=dim
        im=fftshift(fft(ifftshift(im,i),sz(i),i),i)./sqrt(sz(i));
    end
else
    for i=dim
        im=fft(im,sz(i),i);
    end
end

end