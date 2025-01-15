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