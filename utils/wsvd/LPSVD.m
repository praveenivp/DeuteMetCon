function [Ampl,b,w,Basis] = LPSVD(y, BasisSize)
%LPSVD Compute the linear prediction of a signal based on singular value
%decomposition
%   The used algorithm is described here:
%    Improved algorithm for noniterative time-domain model fitting to 
%    exponentially damped magnetic resonance signals
%    Volume 73, Issue 3, July 1987, Pages 553–557
%   A sum of damped sinusoids is fitted to the NMR signal y
%   BasisSize = number of basis vectors


%% Transpose y if it is not a column vector
if(size(y,1)<size(y,2))
    y = y.';
end

%% Number of samples
N = length(y);

%% Hankel matrix
M = ceil(0.75*N);
X = zeros(N-M+1,M);
for i=1:N-M+1
    X(i,:) = y(i:i+M-1);
end

%% Singular value decomposition of the Hankel matrix
% X = U*A*V'.
%[U,A,Vd] = svd(X);
%V = Vd';
[U,~,~] = svd(X);

%% Number of signal components
if(nargin<2)
    K = 10;
    warning('Using 10 basis vectors as default')
else
    K = BasisSize;
end

%% Truncate matrices
%Ak = A(1:K,1:K);
%Vk = V(:,1:K);
Uk = U(:,1:K);

%% Discard bottom and top row of U, respectively
Ub = Uk(1:end-1,:);
Ut = Uk(2:end,:);

%% Calculate Z prime (Eq. 12 in paper mentionned above)
Zp = Ub\Ut;

%% Diagonalise matrix
[~,D] = eig(Zp);

%% Poles zk = exp(-kb+1i*wk) are given by the eigenvalues of Z prime
zk = diag(D);

%% The damping factors b and frequencies w are given by
b = -real(log(zk))';
w = imag(log(zk))';

%% Construct basis
t = [0:length(y)-1]';
Basis = exp(bsxfun(@times,-b+1i*w,t));

%% Determine amplitude by least-square fit
Ampl = (Basis\y).';


