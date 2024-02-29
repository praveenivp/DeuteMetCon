function sig_amp=getBandingRegion(TR,fmap,metabolites,FA)
%sig_amp=getBandingRegion(TR,fmap,metabolties1,FA)
% function to get banding region depending on fieldmap,TR and metabolite
% sig_amp=getBandingRegion(10e-3,fmap_rad_s,metabolties1,40);
%
%INPUTS:
% TR : Repetition time in s
% fmap : 1H fieldmap rad_s
% metabolites : struct
% FA : flip angle in radians




if( isobject(TR))
    mcobj=TR;
    TR=mcobj.DMIPara.TR;
    fmap=mcobj.FieldMap;
    metabolites=mcobj.metabolites;
    FA=mcobj.DMIPara.FlipAngle;

end
if(FA>2)
    warning('converting FA into radians')
    FA=deg2rad(FA);
end
if(exist('metabolties','var'))
    cs_ppm=[0 1.0256 ]*1e-6; %D20 Glu Glx Lac/fat
    freq_shift_WGX=[3.9 -58.7  -152 -215];
    T1=[432.88 69.65 147.6 190.64]*1e-3;%s
    T2=[287 65.88 124.5 180]*1e-3;%s

    sp_name={'D20','Glucose','Glutamate','Lactate'};
    for i=1:4
        metabolites(5-i)=struct('T1_s',T1(i),'T2_s',T2(i),'freq_shift_Hz',freq_shift_WGX(i),'name',sp_name{i});
    end
end


% fmap in 2H hz
fmap=(fmap/(2*pi))*(6.536 /42.567); 
mask=~(fmap==0);
[Msig_all]=bSSFP_sim_analytical(metabolites,TR/2,0,TR,fmap(mask),FA);


Msig_all=squeeze(Msig_all);
sig_amp=zeros([numel(fmap) length(metabolites)]);

for i=1:length(metabolites)
sig_amp(mask(:),i)=Msig_all(i,:);
end

sig_amp=reshape(sig_amp,[size(fmap) length(metabolites)]);


end