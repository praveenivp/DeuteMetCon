function [D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix)
% [D_noise,D_image,noise_info]=CalcNoiseDecorrMat(twix)
% function to calc noise correlation function
% D_noise is noise decorr matrix from trwi.noise
% D_image is noise decorr matrix from image data (only for 0 flipangle scans)
% noise_info : struct with dwell time and oversampling removal flag
if(iscell(twix))
    twix=twix{end};
end
NormNoiseData=false;
% CoilSel
          if isfield(twix,'noise')
                noise                = permute(twix.noise(:,:,:),[2,1,3]);
                noise                = noise(:,:)';
                R                    = cov(noise);
                R(eye(size(R,1))==1) = abs(diag(R));
                if(NormNoiseData)
                    R= R./mean(abs(diag(R)));
                    D_noise               = sqrtm(inv(R));
                else
                    scale_factor=1; %dwell time are the same
                    Rinv = inv(chol(R,'lower'));
                    D_noise = Rinv*sqrt(2)*sqrt(scale_factor);

                end
            else
                D_noise = 1;
          end

%if we have 0 flipangle data
                    if isfield(twix,'image')
                noise                = permute(twix.image.unsorted,[2,1,3:12]);
                noise                = noise(:,:)';
                R                    = cov(noise);
                R(eye(size(R,1))==1) = abs(diag(R));
                if(NormNoiseData)
                    R= R./mean(abs(diag(R)));
                    D_image               = sqrtm(inv(R));
                else
                    scale_factor=1; %dwell time are the same
                    Rinv = inv(chol(R,'lower'));
                    D_image = Rinv*sqrt(2)*sqrt(scale_factor);

                end
            else
                D_image = 1;
                    end
% only for imaging sequence
noise_info.os_flag=    (twix.image.NCol~=twix.hdr.Phoenix.sKSpace.lBaseResolution);
noise_info.dwell_s=twix.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9; %s

                    
end