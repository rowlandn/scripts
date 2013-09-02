function issignificative = evRanksum(alignmat, pop1, pop2, alpha)
%issignificative = evRanksum(alignmat, pop1, pop2, alpha)



M = size(alignmat,1);
pvalues = zeros(M,1);
issignificative = zeros(M,1);
for i = 1 : M
    a = alignmat(i,pop1);
    a = a(find(~isnan(a)));
    b = alignmat(i,pop2);
    b = b(find(~isnan(b)));
    if length(a) & length(b)
        [p,h] = ranksum(a,b,alpha);
        pvalues(i) = p;
        issignificative(i) = h;
        %[pvalues(i),issignificative(i),stats_ranksum(i)] = ranksum(a,b,.05);
    else
        pvalues(i) = NaN;
        issignificative(i) = NaN;
    end
    
end