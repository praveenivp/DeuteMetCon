function [fid_extrap] = fidExtrp(fids,ext_size)
%--------------------------------------------------------------------------
% Linear predicion of the missing points of FID signal measured with
% CSI-FID sequence.
%--------------------------------------------------------------------------
% Usage:
%       [fid_extrap] = fidExtrp(fid,ext_size)
% where:
%              fid - original FID signal;
%                  [Time x other dims]    
%         ext_size - number of points to be predicted. This depends from
%                    acquisition duration, vector size and delay time of
%                    the CSI-FID sequence;
%       fid_extrap - final fid signal with begin points extrapolated
%--------------------------------------------------------------------------
% Function calculates the missing points in the FID signal measured with
% CSI-FID sequence. Function should be used with non-complex data.
% Prediction of the missing points can be done separately for the real
% and imaginary part of the FID.
%
% Adopted from the code on:
%                                     https://gist.github.com/tobin/2843661
%
% Thread about linear prediction is available under:
%
%                                     http://dsp.stackexchange.com/a/110/64
%
% Adopted by:                                                     GCH, 2013
% Last changes:
% 15.04.2015 - GCH - adopted to ToolboxCSI v2.
% support multidim
%--------------------------------------------------------------------------


fidSz=size(fids);
fids=reshape(fids,fidSz(1),[]);

fid_extrap=zeros([size(fids) 2]);
for i=1:size(fids,2)
    for j=1:2

        if(j==1)
            fid=real(fids(:,i));
        else
            fid=imag(fids(:,i));
        end

        LPC_start = length(fid);          % Starting point for LPC
        LPC_order = length(fid)-1;        % Order of LPC autoregression
        LPC_datax = length(fid)+ext_size; % Size of expadned data set

        % FID signal is rotated so that the begin goes to the end and the LPC
        % algorithm can predict the end
        fid_rot = rot90(fid,2);

        % Autoregression with Burg method
        % from: Signal Processing Toolbox / Parametric Modelling
        LPC_arb = arburg(fid_rot,LPC_order);
        %LPC_arb = lpc(fid_rot,LPC_order);

        % Below is the expanded dataset defined:
        fid_new = zeros(LPC_datax,1);
        fid_new(1:LPC_start) = fid_rot;

        % Initial FID is runned through the filter to get the LPC coeficients
        [~, zf] = filter(-[0 LPC_arb(2:end)], 1, fid_rot(1:LPC_start));

        % Filter is used as an IIR (infinite impulse response) to extrapolate the
        % missing data points
        fid_new((LPC_start+1):LPC_datax) = filter([0 0], -LPC_arb, zeros(LPC_datax-LPC_start,1), zf);

        % Final result must be rotated back in order to match the imput time series
        fid_0 = rot90(fid_new,2);
        fid_extrap(:,i,j) = fid_0(1:LPC_start);

    end
end
fid_extrap=fid_extrap(:,:,1)+1i*fid_extrap(:,:,2);
fid_extrap=reshape(fid_extrap,fidSz);

end