function [fitf,gof,foptions]=NLorentzFit(xdata,ydata,freq)
%  [fitf,gof,foptions]=NLorentzFit(xdata,ydata,freq)
% xdata: spectral axis (in Hz)
% freq : same unit as xdata
%ADapted from:
% fit n Lorentzians to xydata with automated guessing of peak locations
% Erik Bauch, May 2017
% github.com/ebauch/matlab

npeaks=length(freq);
for ii=1:npeaks
[~,idx(ii)]=min(abs(xdata-freq(ii)));
end

freqamp=(ydata(idx));
freqampUL=1.5*abs(ydata(idx));
freqampLL=-1.5*abs(ydata(idx));
offset=0;
gammas=(ones(npeaks,1).*18)'; %linewidth
% fitting

% build n-peak Lorentzian model function
model = 'A1/(1 + ((x - center1)/(gamma1/2))^2) '; % first Lorentzian
%             startpoints = [freqamp(1) freq(1) gammas(1) ];
%             UL=[freqampUL(1) freq(1)+20 1e3*gammas(1) ];
%             LL=[freqampLL(1) freq(1)-20 1e-1*gammas(1)];
for i = 2 : npeaks
   
    model = [model sprintf(' + A%d/(1 + ((x - center%d)/(gamma%d/2))^2)', i, i, i)];

%             startpoints = [startpoints freqamp(i) freq(i) gammas(i) ];
%             UL=[UL freqampUL(i) freq(i)+20 1e3*gammas(i)] ;
%             LL=[LL freqampLL(i) freq(i)-20 1e-1*gammas(i)];
    
end

 model = [model '+ offset'];
%                 startpoints = [startpoints offset];
%             UL=[UL 10];
%             LL=[LL -10];
%            

            startpoints = [freqamp(:)' freq(:)' gammas(:)' offset ];
            UL=[freqampUL(:)' freq(:)'+10 4*gammas(:)' offset+max(ydata) ];
            LL=[freqampLL(:)' freq(:)'-10 -4*gammas(:)' offset-max(ydata)];
            % intial parameters [amplitudes centers linewidths offset]                        
            

            
            foptions = fitoptions('Method','NonlinearLeastSquares',...
                  'Algorithm','Trust-Region',...                %'Levenberg-Marquardt', 'Gauss-Newton', or 'Trust-Region'. The default is 'Trust-Region'              'MaxIter',1000,...
                  'Robust','off',... % LAR, Bisquare or off
                  'Startpoint', startpoints,...
                  'MaxFunEvals', 5000 ,...
                  'MaxIter', 5000 ,...
                  'TolFun', 10^-8,...
                  'Upper',UL,'Lower',LL);
                         
            ftype = fittype(model, 'options', foptions);
%             coeffnames(ftype) % gives the order of the parameters 
            
            [fitf, gof] = fit(xdata(:),ydata(:),ftype);
%             figure,plot(xdata,ydata),hold on,plot(xdata,fitf(xdata))
end       