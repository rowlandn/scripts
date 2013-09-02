function S=Poisson_surprise(ISIs,rate);

% S=POISSON_SURPRISE(ISI,RATE);
% This function calculates the Poisson surprise value of a sequence of inter-
% spike intervals (ISI), with a given RATE. Make sure that the units of 
% the rate are the reciprocal of the units of the intervals.

% JAG, 5/2/03 

N=length(ISIs);
T=sum(ISIs);
argum=rate*T;
facts=[];
for k=0:N-1
   facts=[facts factorial(k)];
end
S=-log(1-exp(-argum)*sum((argum).^(0:N-1)./facts));
