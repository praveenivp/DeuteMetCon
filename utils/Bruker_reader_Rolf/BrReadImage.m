function [res, p] = BrReadImage(path, params)
%Reads Bruker raw data with one to four dimensions and reconstructs the
%images.
%The parameter path can be either a directory, which serves as starting
%point for the search, or a raw data file. If path is completely missing, a
%file search dialogue is used to find the raw data file.

if nargin == 0
    path = '';
end
if nargin < 2
    params.complex = 0;
    params.filter=[0,0,0,0];  % filter: 1: Hanning filter, 2: Hamming filter, 3: optimized generalized Hamming filter, 4: gaussian filter, 5: exponential filter
    params.phase = 0;
    params.zf = [0,0,0,0];
    params.ft = [1,1,1,1];
end
res = [];
p = [];
if strcmp(path, '')==1
    path = BrPickDataset(path);
    if numel(path) == 0
        return
    end
elseif exist(path,'dir')==7
 path = BrPickDataset(path);
    if numel(path) == 0
        return
    end
end    

% Loop through the files. All files have to have equal dimensions.
NFiles = size(path);
if numel(NFiles) > 1
    NFiles = NFiles(1);
end
pnames = {'ACQ_size', 'ACQ_jobs', 'NI', 'NR', 'ACQ_phase_factor', 'ACQ_rare_factor','ACQ_obj_order', 'SW_h','NSLICES','ACQ_echo_time','ACQ_spatial_size_1','ACQ_spatial_size_2','ACQ_spatial_phase_1','ACQ_spatial_phase_2','ACQ_ReceiverSelect'};
for CntFile = 1:NFiles
    if NFiles > 1
        fid = BrReadFid(path{CntFile});
        % Read the parameters
        pathstr = fileparts(path{CntFile});
    else
        if iscell(path)
            path = path{1};
        end
        fid = BrReadFid(path);
        % Read the parameters
        pathstr = fileparts(path);
    end        
    acqpfile = [pathstr,filesep,'acqp'];
    pvmfile = [pathstr,filesep,'method'];
    disp(['Reading parameter file ',acqpfile]);
    p = BrReadParams(pnames,acqpfile);
    pnames = {'PVM_EncNReceivers','PVM_EncMatrix'};
    pvm = BrReadParams(pnames, pvmfile);
    if numel(pvm.PVM_EncNReceivers) >0
        NReceivers = pvm.PVM_EncNReceivers;
    else
        NReceivers = 1;
    end
    if p.version < 360
        if p.ACQ_phase_factor == 1 && p.ACQ_rare_factor > 1
            p.ACQ_phase_factor =  p.ACQ_rare_factor;
        end
        if numel(p.ACQ_size)<2 
            p.ACQ_size(2)=1;
        end
        if numel(p.ACQ_size)<3 
            p.ACQ_size(3)=1;
        end
        if numel(p.ACQ_size)<4 
            p.ACQ_size(4)=1;
        end
        datasize = p.ACQ_size(1)/2;
        ImageSize = p.ACQ_size;
        ImageSize(1) = ImageSize(1)/2;
        numbers = [1,128,192,256,320,384,448,512,576,640,704,768,832,896,960,1024, 1536,1920, 2048, 2304, 4096];
        %numbers = [1,128,256,320,384,448,512,576,640,704,768,832,896,960,1024, 1536,1920, 2048, 2304, 4096];
        while numel(fid) ~= datasize*p.ACQ_phase_factor*p.NI*p.ACQ_size(2)/p.ACQ_phase_factor*p.ACQ_size(3)*p.ACQ_size(4)*p.NR*NReceivers
        
        %qq = find(numbers == p.ACQ_size(1)/2);
        %if numel(qq) < 1
        %    disp(['Theres something wrong with the array sizes']);
        %    disp(p.ACQ_size(1)/2);
        %    disp('');
        %else
            qq = find(numbers > datasize);
            qq = qq(1);
            datasize = numbers(qq);
            %disp(['Changing ACQ_size from ',num2str(p.ACQ_size(1)/2),' to ',num2str(numbers(qq))]);
            %p.ACQ_size(1) = numbers(qq)*2;
        %end
        end
    else             % ParaVision 360
        datasize = p.ACQ_jobs(1).ScanSize/2;
        if numel(p.ACQ_spatial_size_1) == 0
            p.ACQ_spatial_size_1 = 1;
        end
        if numel(p.ACQ_spatial_size_2) == 0
            p.ACQ_spatial_size_2 = 1;
        end        
        if datasize<128
            q = round(numel(fid)/(NReceivers*p.ACQ_phase_factor*p.NI*p.ACQ_spatial_size_1/p.ACQ_phase_factor*p.ACQ_spatial_size_2*p.NR));
            if q > datasize & (q == 64 || q == 128)
                datasize = q;
            end
        end
        ImageSize = [datasize,p.ACQ_spatial_size_1,p.ACQ_spatial_size_2,1];
    end
            
    fid =reshape(fid,[datasize,NReceivers,p.ACQ_phase_factor,p.NI,ImageSize(2)/p.ACQ_phase_factor,ImageSize(3),ImageSize(4),p.NR]);
    if p.version < 360 && datasize>p.ACQ_size(1)/2
        fid = fid(1:p.ACQ_size(1)/2,:,:,:,:,:,:,:);
    end
    if numel(p.ACQ_spatial_phase_1) > 1
    fid = permute(fid,[1,3,5,6,7,4,8,2]);
    else
        fid = permute(fid,[1,5,3,6,7,4,8,2]);
    end
    fid(:,:,:,:,:,p.ACQ_obj_order+1,:,:) = fid;
    fid = reshape(fid,[ImageSize,p.NI, p.NR,NReceivers]);
    %Mit ACQ_spatial_phase
    if numel(p.ACQ_spatial_phase_1) > 1
        [junk,ind1] = sort(p.ACQ_spatial_phase_1);
        fid = fid(:,ind1,:,:,:,:,:);
    end
    if numel(p.ACQ_spatial_phase_2) > 1 && p.ACQ_size(3) > 1
        [junk,ind2] = sort(p.ACQ_spatial_phase_2);
        fid = fid(:,:,ind2,:,:,:,:);
    end
    %Filter
    if max(params.filter)>0
        siz = ones(7,1);
        s = size(fid);
        filter = ones(s);
        ftype = find(params.filter == 3);    %optimized generalized Hamming filter
        if numel(ftype) > 0
            h = MakeHanning(s(ftype),1,1,0.58,1,0);
            r = ones(1,numel(s));
            r(ftype) = s(ftype);
            h = reshape(h,r);
            q = s;
            q(ftype) = 1;
            h = repmat(h,q);
            filter = filter.*h;
        end
        ftype = find(params.filter == 1);    % Hanning filter
        if numel(ftype) > 0
            h = MakeHanning(s(ftype),1,1,0,0,0);
            r = ones(1,numel(s));
            r(ftype) = s(ftype);
            h = reshape(h,r);
            q = s;
            q(ftype) = 1;
            h = repmat(h,q);
            filter = filter.*h;
        end
        ftype = find(params.filter == 2);    % Hamming filter
        if numel(ftype) > 0
            h = MakeHanning(s(ftype),1,1,0,1,0);
            r = ones(1,numel(s));
            r(ftype) = s(ftype);
            h = reshape(h,r);
            q = s;
            q(ftype) = 1;
            h = repmat(h,q);
            filter = filter.*h;
        end

        fid = fid.*filter;
    end
    % zerofilling
    if numel(params.zf)<4
        zf = params.zf;
        params.zf = ImageSize;
        params.zf(1:numel(zf)) = zf;
    else
        mm = params.zf<ImageSize;
        if max(mm)>0
            a = ImageSize;
            params.zf(mm) = a(mm);
        end           
    end
    if max(params.zf>ImageSize) 
        store = fid;
        fid = zeros(params.zf(1),params.zf(2),params.zf(3),params.zf(4),p.NI, p.NR,NReceivers);
        fid(1:ImageSize(1),1:ImageSize(2),1:ImageSize(3),1:ImageSize(4),:,:,:) = store;
        store = 0;
    end
    
    %Fourier transform
    ftdir = zeros(4);
    ftshift = zeros(4);
    if ImageSize(1)>1 && params.ft(1) == 1
        ftdir(1)=1;
        ftshift(1) = -floor(ImageSize(1)/2);
        ndims = 1;
    end
    if ImageSize(2)>1 && params.ft(2) == 1
        ftdir(2)=1;
        ftshift(2) = -floor(ImageSize(2)/2);
        ndims = 2;
    end
    if ImageSize(3)>1 && params.ft(3) == 1
        ftdir(3)=1;
        ftshift(3) = -floor(ImageSize(3)/2);
        ndims = 3;
    end
    if ImageSize(4)>1 && params.ft(4) == 1
        ftdir(4)=1;
        ftshift(4) =-floor(ImageSize(4)/2);
        ndims = 4;
    end

    if params.complex == 1 || params.phase == 1
        fid = circshift(fid,[ftshift(1), ftshift(2), ftshift(3), ftshift(4), 0,0,0]);
    end
    for cntrec = 1:NReceivers
        for cntnr = 1:p.NR
            for cntni = 1:p.NI
                if ftdir(4) == 1
                    if params.filter(4) == 1
                    end
                    fid(:,:,:,:,cntni,cntnr,cntrec) = circshift(fft(fid(:,:,:,:,cntni,cntnr,cntrec),[],4),[0,0,0,-round(params.zf(4)/2)]);
                end
                if ftdir(3) == 1
                    if params.filter(3) == 1
                    end
                    fid(:,:,:,:,cntni,cntnr,cntrec) = circshift(fft(fid(:,:,:,:,cntni,cntnr,cntrec),[],3),[0,0,-round(params.zf(3)/2)]);
                end
                if ftdir(2) == 1
                    fid(:,:,:,:,cntni,cntnr,cntrec) = circshift(fft(fid(:,:,:,:,cntni,cntnr,cntrec),[],2),[0,-round(params.zf(2)/2)]);
                end
                if ftdir(1) == 1
                    fid(:,:,:,:,cntni,cntnr,cntrec) = circshift(fft(fid(:,:,:,:,cntni,cntnr,cntrec),[],1),-round(params.zf(1)/2));
                end
            end
        end
    end
    %fid = reshape(fid,[params.zf(1),params.zf(2),params.zf(3),params.zf(4),p.NSLICES,p.NI/p.NSLICES, p.NR, NReceivers]);
    fid = reshape(fid,[params.zf(1),params.zf(2),params.zf(3),params.zf(4),p.NI/p.NSLICES,p.NSLICES, p.NR, NReceivers]);
    dim3 = params.zf(3);
    dim4 = params.zf(4);
    if p.NSLICES > 1
        if params.zf(3)==1
            fid = permute(fid,[1,2,6,4,5,3,7,8]);
            dim3 = p.NSLICES;
        elseif params.zf(4) == 1
            fid = permute(fid,[1,2,3,6,5,4,7,8]);
            dim4 = p.NSLICES;
        end
    end
    fid = reshape(fid,[params.zf(1),params.zf(2),dim3,dim4,p.NI/p.NSLICES* p.NR, NReceivers]);
    
    if CntFile == 1
        if params.complex == 0
            if params.phase == 0
                res = abs(fid);
                if NReceivers > 1
                    res = sqrt(sum(res.^2,6));
                end
            else
                res = angle(fid)*180/pi;
            end
        else
            res = fid;
        end
        s1 = size(squeeze(fid));
    else
        s = size(squeeze(fid));
        if numel(s1) == 2 && numel(s)==2 && s1 == s
            if params.complex == 0
                if params.phase == 0
                    fid = abs(fid);
                    if NReceivers > 1
                        fid = sqrt(sum(fid.^2,6));
                    end
                else
                    fid = angle(fid)*180/pi;
                end
            end
            res(:,:,CntFile) = squeeze(fid);
        else
            res{CntFile}=fid;
        end
    end
        
    
end



return;
end


function res = hannig(N)
res = 0.5*(1+cos(linspace(0,N+1,N+2)/(N+1)*2*pi-pi));
res = res(2:N+1);
return
end

function res = gaussian(N,sigma)
res = exp(-linspace(-N/2,N/2-1,N).^2/2/sigma);
%res = res(2:N+1);
return
end
