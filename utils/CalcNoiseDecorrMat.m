function [D_noise,D_image]=CalcNoiseDecorrMat(twix)

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

end