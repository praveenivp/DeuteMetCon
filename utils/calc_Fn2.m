function Fn = calc_Fn2(img_series,phi_vec,Np)
%function Fn = calc_Fn(img_series,phi_vec,Np)
%
%Calculation of higher-order SSFP signals from a series of phase-cycled bSSFP measurements 
%following Zur et al., Motion-Insensitive, Steady-State Free Precession Imaging, 1990
%
%INPUT
%
%image data (img_series)
%RF phase cycle vector in rad (phi_vec)
%highest mode number to be calculated (Np)
%
%OUTPUT
%
%Fn (higher order modes n)

shape = size(img_series); % number of measurements is assumend to be the last dimension
nPC=length(phi_vec);
assert(shape(end)==nPC,'phase cycle dimension should be the last dim')
img_series=reshape(img_series,[],nPC);


%!!!BE AWARE OF THE + or - MINUS SIGN!!!
% vectorized version of rheule's
p = -Np:1:Np;

Fn_basis=exp(complex(0,1)*phi_vec(:)*(p(:).'));
Fn=img_series*Fn_basis;
Fn=reshape(Fn,[shape(1:end-1) length(p)])/nPC;
end