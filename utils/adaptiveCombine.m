function [recon,wfull,norm]=adaptiveCombine(im,bs,modeSVD,modeSliBySli,st,rn,donorm)
%% Combines Multi-Coil images
% 
% Examples:
% combinedImage=adaptiveCombine2(CoilImages);
% [combinedImage,CoilWeights,NormalizationMAtrix]=adaptiveCombine2(CoilImages,BlockSize);
%
% INPUTS:
% im   :  Uncombined coil image/volume
%         3D or 4D complex double matrix, coil should be the first dimension
% bs    : Block Size of kernel in all physical dimension
%         Optional parameter: default is min(9,MatrixSize in that dimension)
%         Integer Array has length ndims(im)-1 
% modeSVD: 
%         flag to use SVD for coil weighting(much better for phase correction!)
%         Optional Parameter: default is true
%         Boolean: true or false
% modesliBySli: 
%         flag to treat third dimension as slice
%         Optional Parameter: default is false(treated as volume)
%         Boolean: true or false
%st:      Intepolation step size
%         Optional parmeter: deafult is ones(1,3)
%         Integer Array has length ndims(im)-1
%rn:      
%         Noise correlation matrix
%         Optional Parameter: default is 1
%         1 or 2D matrix of size [Ncoil x Ncoil]
% donorm: Output image normalization flag
%         Optional parameter: default is false
%         Boolean: true or false(default)
%
%OUTPUTS:
% recon: Combined coil image/volume
% wi   : Coil weights used for coil combination
% norm : Normalization matrix
%
%   See also adaptiveCombine2
%
%
%   Adaptive recon based on Walsh et al.
%   Walsh DO, Gmitro AF, Marcellin MW.
%   Adaptive reconstruction of phased array MR imagery. 
%   Magn Reson Med. 2000 May;43(5):682-90.
%
%    and
%
%   Mark Griswold, David Walsh, Robin Heidemann, Axel Haase, Peter Jakob. 
%   The Use of an Adaptive Reconstruction for Array Coil Sensitivity Mapping and Intensity Normalization, 
%   Proceedings of the Tenth  Scientific Meeting of the International Society for Magnetic Resonance in Medicine pg 2410 (2002)

if ~exist('modeSVD','var')
    modeSVD = true; %much better for phase correction!
end

if ~exist('modeSliBySli','var')
    modeSliBySli =  false;
end
if ~exist('donorm','var')
    donorm=0;
end
if ~exist('rn','var')
    rn = 1;
end

sz = size(im); 
nc = sz(1);
n  = ones(1,3);
n(1:numel(sz)-1) = sz(2:end);

if ~exist('st','var')
    st = min(1,n);
end
nsmall = round(n./st);

% in 3D können wir blocksize reduzieren, da generell mehr pixel fuer 
% statistik zur verfuegung stehen (b_y*b_x*bs_z statt bs_y*bs_x)
if ~exist('bs','var') || isempty(bs)
    bs = min(7,n);
    if n(3) > 1
        bs(3) = min(3,n(3));
    end
end
if modeSliBySli
    bs(3) = 1;
end

if modeSVD
    nc_svd = min(min(12,max(9,floor(nc/2))),nc);
else
    nc_svd = nc;
end
    
if false %exist('adaptiveCombine_mex','file')
    if ~modeSVD
        nc_svd = 0;
    end
    [wsmall,normsmall] = adaptiveCombine_mex(im,rn,bs,st,nc_svd,modeSliBySli);
else
    inv_rn = rn\eye(size(rn,1));
    cnt    = 0;
    wsmall = zeros(nc,prod(nsmall),class(im));
    
    maxcoil = 1;
    if ~modeSliBySli
        if modeSVD
            [~,~,V]=svd(im(:,:)*im(:,:)');
            V = V(:,1:nc_svd);
        else
            V = 1;
            [~,maxcoil] = max(sum(abs(im(:,:)),2));
        end
    end
    
    for z=1:nsmall(3)
        if modeSliBySli
            if modeSVD
                tmp = reshape(im(:,:,:,z),nc,[]);
                [~,~,V]=svd(tmp*tmp');
                V = V(:,1:nc_svd);
            else
                V = 1;
                [~,maxcoil] = max(sum(sum(abs(im(:,:,:,z)),2),3));
            end
        end
        for y=1:nsmall(2)
            for x=1:nsmall(1)
                cnt = cnt+1;
                ix = [x,y,z];
                imin = max(st.*ix-floor(bs/2)  , 1);
                imax = min(st.*ix+ceil (bs/2)-1, n);
                m1 = V'*reshape(im(:,imin(1):imax(1),imin(2):imax(2),imin(3):imax(3)),nc,[]);
                m  = m1*m1';              %Calculate signal covariance
                [v,d]   = eig(inv_rn*m);  %Eigenvector with max eigenval gives
                [~,ind] = max(diag(d));   %the correct combination coeffs.
                tmp = v(:,ind);
                tmp = V*tmp; %transform back to original coil space
                %Correct phase based on coil with max intensity
                tmp = tmp*exp(-1j*angle(tmp(maxcoil)));
                wsmall(:,cnt) = conj(tmp)/(tmp'*inv_rn*tmp);
            end
        end
    end
    wsmall = reshape(wsmall,[nc,nsmall]);
    if (donorm || nargout > 2)
        %This is the normalization proposed in the abstract 
        normsmall = reshape(sum(abs(wsmall)).^2,nsmall);
    end
end

if (prod(st) ~= 1)
    %Now have to interpolate these weights up to the full resolution. This is done separately for
    %magnitude and phase in order to avoid 0 magnitude pixels between +1 and -1 pixels.
    [x,y,z] = ndgrid(1:(nsmall(1)-1)/(n(1)-1):nsmall(1),1:(nsmall(2)-1)/(n(2)-1):nsmall(2),1:(nsmall(3)-1)/(n(3)-1):nsmall(3));
    wfull   = zeros(n,'like',im);
    % permute wsmall for faster access
    wsmall  = permute(wsmall,[2 3 4 1]);
    for c=1:nc
        if n(3) > 1
            wfull(c,:,:,:) = interpn(abs(wsmall(:,:,:,c)),x,y,z,'linear') .* exp(1j.*interpn(angle(wsmall(:,:,:,c)),x,y,z,'nearest'));
        else % for 2D scans
            wfull(c,:,:,:) = interpn(abs(wsmall(:,:,:,c)),x,y,'linear')   .* exp(1j.*interpn(angle(wsmall(:,:,:,c)),x,y,'nearest'));
        end
    end
    if (donorm || nargout > 2)
        if nz > 1
            norm = interpn(normsmall,x,y,z,'linear');
        else % for 2D scans
            norm = interpn(normsmall,x,y,'linear');
        end
    end
    clear x y z;
else
    wfull  = wsmall;
    if (donorm || nargout > 2)
        norm  = normsmall;
    end
end
clear wsmall normsmall;

recon = squeeze(sum(wfull.*im));   %Combine coil signals. 

if donorm
    recon = recon.*norm;   
end

% You should carefully read the following terms and conditions before installing or using the 
% software. Unless you have entered into a separate written license agreement with 
% Universit�t W�rzburg providing otherwise, installation or use of the software indicates your 
% agreement to be bound by these terms and conditions. 
% 
% Use of the software provided with this agreement constitutes your acceptance of these terms. 
% If you do NOT agree to the terms of this agreement, promptly remove the software together 
% with all copies from your computer. User's use of this software is conditioned upon compliance 
% by user with the terms of this agreement. 
% 
% Upon ordering, downloading, copying, installing or unencrypting any version of the software, you
% are reaffirming that you agree to be bound by the terms of this agreement. 
% 
% License to use 
% 
% Universit�t W�rzburg grants to you a limited, non-exclusive, non-transferable and non-assignable 
% license to install and use this software for research purposes. Use of this software for any 
% diagnostic imaging procedure is strictly forbidden.
% 
% License to distribute 
% 
% Please feel free to offer the non-commercial version of this software on any website, CD, or 
% bulletin board, demonstrate the non-commercial version of the software and its capabilities, or 
% give copies of the non-commercial version of the software to other potential users, so that others 
% may have the opportunity to obtain a copy for use in accordance with the license terms contained
% here. 
% 
% You agree you will only copy the non-commercial version of the software in whole with this 
% license and all delivered files, but not in part. 
% 
% Termination 
% 
% This license is effective until terminated. You may terminate it at any point by destroying 
% the software together with all copies of the software. 
% 
% If you have acquired a non-commercial version, the license granted herein shall automatically 
% terminate if you fail to comply with any term or condition of this Agreement. 
% 
% Also, Universit�t W�rzburg has the option to terminate any license granted herein if you fail 
% to comply with any term or condition of this Agreement. 
% 
% You agree upon such termination to destroy the software together with all copies of the software.
% 
% 
% Copyright 
% 
% The software is protected by copyright law. You acknowledge that no title to the intellectual 
% property in the software is transferred to you. You further acknowledge that title and full 
% ownership rights to the software will remain the exclusive property of Universit�t W�rzburg, 
% and you will not acquire any rights to the software except as expressly set forth in this 
% license. You agree that any copies of the software will contain the same proprietary notices 
% which appear on and in the software. 
% 
% Rent, lease, loan 
% 
% You may NOT rent, lease or loan the software without first negotiating a specific license
% for that purpose with Universit�t W�rzburg. 
%     
% No warranties 
% 
% Universit�t W�rzburg does NOT warrant that the software is error free. Universit�t W�rzburg 
% disclaims all warranties with respect to the software, either express or implied, including 
% but not limited to implied warranties of merchantability, fitness for a particular purpose and 
% noninfringement of third party rights. The software is provided "AS IS." 
% 
% No liability for consequential damages 
% 
% In no event will Universit�t W�rzburg be liable for any loss of profits, business, use, or data 
% or for any consequential, special, incidental or indirect damages of any kind arising out of 
% the delivery or performance or as a result of using or modifying the software, even if 
% Universit�t W�rzburg has been advised of the possibility of such damages. In no event will 
% Universit�t W�rzburg's liability for any claim, whether in contract, negligence, tort or any 
% other theory of liability, exceed the license fee paid by you, if any. 
% The licensed software is not designed for use in high-risk activities requiring fail-safe 
% performance. Universit�t W�rzburg disclaims any express or implied warranty of fitness for 
% high-risk activities. 
% 
% Severability 
% 
% In the event of invalidity of any provision of this license, the parties agree that such 
% invalidity shall not affect the validity of the remaining portions of this license.
% 
