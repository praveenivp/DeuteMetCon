function omega = my_hann(n)
    n(abs(n)>pi)=pi;
    omega=0.5.*(cos(n)+1);
    %+0.5warning('Hannfilter changed+0.5 ')
end