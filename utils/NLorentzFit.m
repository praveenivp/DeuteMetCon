function [fitf,gof,foptions]=NLorentzFit(xdata,ydata,freq)
%  [fitf,gof,foptions]=NLorentzFit(xdata,ydata,freq)
% xdata: spectral axis (in Hz)
% freq : same unit as xdata
%ADapted from:
% fit n Lorentzians to xydata with automated guessing of peak locations
% Erik Bauch, May 2017
% github.com/ebauch/matlab


npeaks=length(freq);

[~,idx(1)]=min(abs(xdata-freq(1)));
[~,idx(2)]=min(abs(xdata-freq(2)));
[~,idx(3)]=min(abs(xdata-freq(3)));
freqamp=abs(ydata(idx));
freqampUL=[1; 1;1]*max(ydata);
freqampLL=[0; 0; 0];
offset=0;
gammas=([4 4 4].*8)'; %linewidth
% fitting

% build n-peak Lorentzian model function
model = 'A1/(1 + ((x - center1)/(gamma1/2))^2) + offset'; % first Lorentzian
for i = 2 : npeaks
   
    model = [model sprintf(' + A%d/(1 + ((x - center%d)/(gamma%d/2))^2)', i, i, i)];
    
end
           
            % intial parameters [amplitudes centers linewidths offset]                        
            
            startpoints = [freqamp' freq' gammas' offset];
            UL=[freqampUL' freq'+20 1e3*gammas' offset+10];
            LL=[freqampLL' freq'-20 1e-1*gammas' offset-10];
            
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