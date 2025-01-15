function [kdata]=myifft(im,dim,dofftshift,sz)
% [im]=myfft(kdata1,dim,shift,sz)
% dim is array 
% dofftshift is boolean(fftshifgt or not)
% sz final size(tail zeropadding)
%praveenivp
switch (nargin)
    case 1
        dim=1;
        dofftshift=true;
        sz=size(im);
    case 2
        dofftshift=true;
        sz=size(im);
    case 3    
        sz=size(im);  
end

kdata=im;
if(dofftshift)
    for i=dim
        kdata=fftshift(ifft(ifftshift(kdata,i),sz(i),i),i).*sqrt(sz(i));
    end
else
    for i=dim
        kdata=ifft(kdata,sz(i),i);
    end
end

end