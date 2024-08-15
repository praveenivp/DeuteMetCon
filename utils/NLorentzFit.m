function [fitf,gof,foptions]=NLorentzFit(xdata,ydata,freq,foptions)
%  [fitf,gof,foptions]=NLorentzFit(xdata,ydata,freq)
% xdata: spectral axis (in Hz)
% freq : same unit as xdata
% foptions: fitoptions output
%ADapted from:
% fit n Lorentzians to xydata with automated guessing of peak locations
% Erik Bauch, May 2017
% github.com/ebauch/matlab
npeaks=length(freq);
if(exist("foptions",'var'))
    foptions=fitoptions(foptions);
else
    for ii=1:npeaks
        [~,idx(ii)]=min(abs(xdata-freq(ii)));
    end

    freqamp=(ydata(idx));
    freqampUL=1.5*abs(ydata(idx));
    freqampLL=-1.5*abs(ydata(idx));
    offset=0;
    gammas=(ones(npeaks,1).*18)'; %linewidth


    startpoints = [freqamp(:)' freq(:)' gammas(:)' offset ];
    UL=[freqampUL(:)' freq(:)'+10 4*gammas(:)' offset+max(ydata) ];
    LL=[freqampLL(:)' freq(:)'-10 -4*gammas(:)' offset-max(ydata)];
    % intial parameters [amplitudes centers linewidths offset]


    foptions = fitoptions('Method','NonlinearLeastSquares',...
        'Algorithm','Trust-Region',...                %'Levenberg-Marquardt', 'Gauss-Newton', or 'Trust-Region'. The default is 'Trust-Region'              'MaxIter',1000,...
        'Robust','LAR',... % LAR, Bisquare or off
        'Startpoint', startpoints,...
        'MaxFunEvals', 100 ,...
        'MaxIter', 500,...
        'TolFun', 10^-6,...
        'Upper',UL,'Lower',LL);
end


% build n-peak Lorentzian model function
model = 'A1/(1 + ((x - center1)/(gamma1/2))^2) '; % first Lorentzian
for i = 2 : npeaks
    model = [model sprintf(' + A%d/(1 + ((x - center%d)/(gamma%d/2))^2)', i, i, i)];
end
model = [model '+ offset'];

ftype = fittype(model, 'options', foptions);

[fitf, gof] = fit(xdata(:),ydata(:),ftype);
end