function [ind,indy,h]=find_train(a,k);
% [IND,FIRST_IND,H]=FIND_TRAIN(A,K) finds the indices IND of the vector of ones and
% zeros A that are part of a K long train of ones. H is 1 if such trains
% exist in A and 0 otherwise. FIRST_IND is the index of the first index of the trains.

b=a;
for j=2:k
    b=[b(:,1:end-1);a(j:end)];
end
b = single(b);
indy=find(prod(b,1));
h=~isempty(indy);
if k>1
ind0=repmat(indy,k,1)+[0:k-1]'*ones(1,length(indy));
else
ind0=indy;
end
ind=unique(ind0(:));
indy(find(diff(indy)==1)+1)=[];
