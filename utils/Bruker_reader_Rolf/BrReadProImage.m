function [res,params] = BrReadProImage(path)
%Reads Bruker processed data with one to four dimensions and reconstructs the
%images.
%The parameter path can be either a directory, which serves as starting
%point for the search, or a 2dseq file. If path is completely missing, a
%file search dialogue is used to find the raw data file.

if nargin == 0
    path = '';
end
% if strcmp(path, '')==1
% path = BrPickDataset(path,0,'image');
% end


    
    

% Loop through the files. All files have to have equal dimensions.
NFiles = size(path);
%if numel(NFiles) > 1
    NFiles = NFiles(1);
%end
rnames = {'RECO_size','RECO_fov', 'RECO_Transposition'};
pnames = {'ACQ_size', 'NI', 'NR', 'ACQ_phase_factor', 'ACQ_obj_order', 'SW_h','NSLICES','ACQ_echo_time','IRTime', 'ACQ_slice_offset','ACQ_GradientMatrix','ACQ_grad_matrix'};
for CntFile = 1:NFiles
    if NFiles > 1
        fid = BrRead2dSeq(path{CntFile});
        % Read the parameters
        pathstr = fileparts(path{CntFile});
    else
        fid = BrRead2dSeq(path);
        % Read the parameters
        pathstr = fileparts(path);
    end        
    recofile = [pathstr,filesep,'reco'];
    disp(['Reading parameter file ',recofile]);
    r = BrReadParams(rnames,recofile);
    acqpfile = [pathstr,filesep,'..',filesep,'..',filesep,'acqp'];
    disp(['Reading parameter file ',acqpfile]);
    p = BrReadParams(pnames,acqpfile);

    % Correct for different PV versions
    if numel(p.ACQ_GradientMatrix) == 0 && numel(p.ACQ_grad_matrix) >= 3
        p.ACQ_GradientMatrix = p.ACQ_grad_matrix;
    end
    p = rmfield(p,'ACQ_grad_matrix');
    
    if numel(p.ACQ_size)<3 
        p.ACQ_size(3)=1;
    end
    if numel(p.ACQ_size)<4 
        p.ACQ_size(4)=1;
    end
    if numel(r.RECO_size)<2 
        r.RECO_size(2)=1;
    end
     if numel(r.RECO_size)<3 
        r.RECO_size(3)=1;
     end
     if r.RECO_Transposition(1) == 1
         rr = r.RECO_size(2);
         r.RECO_size(2) = r.RECO_size(1);
         r.RECO_size(1) = rr;
     end
   
    size(fid)
    %fid =rot90(reshape(fid,[r.RECO_size(1),r.RECO_size(2),r.RECO_size(3),p.NI/p.NSLICES,p.NSLICES,p.NR]));
    fid =reshape(fid,[r.RECO_size(1),r.RECO_size(2),r.RECO_size(3),p.NI/p.NSLICES,p.NSLICES,p.NR]);  % Do we need the rot90 here? Why? Can we replace that somehow?

    if CntFile == 1
        res = fid;
        s1 = size(squeeze(fid));
    else
        s = size(squeeze(fid));
        if numel(s1) == 2 & numel(s)==2 & s1 == s
            res(:,:,CntFile) = squeeze(fid);
        else
            res{CntFile}=fid;
        end
    end
        
    
end

params = p;
params.recosize = r.RECO_size;
params.fov = r.RECO_fov;
params.transposition = r.RECO_Transposition;

return;
end
