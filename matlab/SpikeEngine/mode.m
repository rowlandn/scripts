function y=mode(x)
% MODE finds the mode of a sample.  The mode is the
% observation with the greatest frequency.
%
% i.e. in a sample x=[0, 1, 0, 0, 0, 3, 0, 1, 3, 1, 2, 2, 0, 1]
% the mode is the most frequent item, 0.


% IT'S NOT FANCY BUT IT WORKS
% Michael Robbins
% robbins@bloomberg.net
% michael.robbins@us.cibc.com

[b,i,j] = unique(x);
[m,k]=max(hist(j,length(b)));
y=b(k);
