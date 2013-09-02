function delta = isi2delta(isi)
% function delta = isi2delta(isi)
% 
% utility to convert spike isi function to delta
%
y = cumsum(isi);
delta([1 y+1])=1;

return
