function isi_shuf = shuffle_isi(isi)
% function to shuffle isi's
%
r=randperm(length(isi));
isi_shuf(r)=isi;

return
