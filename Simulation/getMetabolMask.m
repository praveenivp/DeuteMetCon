function [met_mask,mask_labels,id_labels]=getMetabolMask(MatSz,nMet) 
% met_mask=getMetabolMask(MatSz) 
% function to get overlapping 2D binary mask like venn diagrams
% 
%met_mask=getMetabolMask(128,3); 
% [met_mask,mask_labels,id_labels]=getMetabolMask(64,4) ;
%
%
% praveeen.ivp


[xx,yy]=meshgrid(ceil(MatSz/-2):ceil(MatSz/2)-1);
met_mask=zeros(MatSz*MatSz,nMet);
mask_labels=zeros(MatSz,MatSz);
id_labels={'A','B','C','A+B','A+C','A+D','B+C','B+D','C+D','A+B+C','A+B+D','A+C+D','B+C+D','A+B+C+D'};
if(MatSz==1),met_mask=ones(MatSz*MatSz,nMet);
    return;
end

if(nMet==4)

MyEllipse=@(xx,yy,cx,cy,theta) 4*((xx-cx)*cos(theta)+(yy-cy)*sin(theta)).^2+...
    1*((xx-cx)*sin(theta)-(yy-cy)*cos(theta)).^2 ...
    < round(MatSz/2.5).^2;

%center parameter
para2=round(MatSz/3);

%four ellipses
met_mask(MyEllipse(xx(:),yy(:),-0.4*para2,-0.4*para2,pi/4),1)=1;

met_mask(MyEllipse(xx(:),yy(:),0,0,pi/4),2)=1;

met_mask(MyEllipse(xx(:),yy(:),0,0,-pi/4),3)=1;
met_mask(MyEllipse(xx(:),yy(:),0.4*para2,-0.4*para2,-pi/4),4)=1;


else

    MyCircle=@(xx,yy,cx,cy) ((xx-cx)).^2+...
    ((yy-cy)).^2 ...
    < round(MatSz/3).^2;

%center parameter
para2=0.4*round(MatSz/3);

%four ellipses
met_mask(MyCircle(xx(:),yy(:),-para2,para2),1)=1;

met_mask(MyCircle(xx(:),yy(:),para2,para2),2)=1;

met_mask(MyCircle(xx(:),yy(:),0,-para2),3)=1;


end
met_mask=reshape(met_mask,MatSz,MatSz,nMet);
%  as(sum(met_mask,3))


% genrate mask label image {A,B,C,....} to {1,2,3,....}

idx_text=regexprep(id_labels,strsplit('A B C D +',' '),strsplit('1 2 3 4 ,',' '));

for cLab=1:length(id_labels)
    evalc(sprintf('cMask=sum(met_mask(:,:,[%s]),3)',idx_text{cLab}));
    mask_labels(cMask>count(id_labels{cLab},'+'))=cLab;
end

end

%
