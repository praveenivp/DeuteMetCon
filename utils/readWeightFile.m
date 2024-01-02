function [W,W1,dat]=readWeightFile(WeightFile,MatSz)
%function [W,W1,dat]=readWeightFile(WeightFile,MatSz)
% WeightFile : from CSI sequence
%optional
% MatSize=[csiTwxi.image.NLin, csiTwix.image.NPar, csiTwix.image.NSeg];
% USAGE:
% [W,W1,dat]=readWeightFile('X:\mrdata\echtdata\studies\48\experiments\KSRI-QYZ6\EXPDATA\weighting.txt',[13 13 13])

%%
% WeightFile='X:\mrdata\echtdata\studies\48\experiments\KSRI-QYZ6\EXPDATA\weighting.txt';
dat=textscan(fopen(WeightFile,'r'),'%d: %d, %d, %d: %d, %f');

AcqIdx=dat{1}+1;


if(~exist('MatSz','var'))
MatSz=double([max(dat{2}),max(dat{3}),max(dat{4})])-...
      double([min(dat{2}),min(dat{3}),min(dat{4})])+...
      [1,  1, 1];
end

LinIdx=dat{2}+ceil(MatSz(1)/2); % first phase encoding
ParIdx=dat{3}+ceil(MatSz(2)/2);% second phase encoding
SegIdx=dat{4}+ceil(MatSz(2)/2); % Third phase encoding

W=zeros(MatSz);
W1=zeros(MatSz);
for i=1:length(LinIdx)
W(LinIdx(i),ParIdx(i),SegIdx(i))=dat{5}(i);
W1(LinIdx(i),ParIdx(i),SegIdx(i))=dat{6}(i);
end
% as(cat(4,W,W1,1./W))
end

%% this  should give same W from the weighting file
% csiTwix=mapVBVD('X:\mrdata\echtdata\studies\48\experiments\KSRI-QYZ6\TWIX\allData#S94Tuebingen#F37802#M376#D170123#T103746#rpcsi_fid.dat');
% MatSz=[csiTwix.image.NLin, csiTwix.image.NPar, csiTwix.image.NSeg];
% LinIdx=csiTwix.image.Lin-1;
% ParIdx=csiTwix.image.Par-1;
% SegIdx=csiTwix.image.Seg-1;
% 
% W_twix=zeros(MatSz);
% for i=1:length(LinIdx)
% W_twix(LinIdx(i),ParIdx(i),SegIdx(i))=W_twix(LinIdx(i),ParIdx(i),SegIdx(i))+1;
% end
% as(W_twix)
