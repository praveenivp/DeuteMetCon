function A=ndCircShift(A,K,dim)
% Circshift for mulidimensional matrix
% A=ndCircShift(A,K,dim)
% A=ndCircShift(A,[1 -1],[1 3])
for i=1:length(dim)
A=circshift(A,K(i),dim(i));
end
end