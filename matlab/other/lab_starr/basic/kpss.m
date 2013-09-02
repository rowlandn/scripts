% KPSS stationarity test ; asymptotic critical values 
% significance
% level          opt=0      opt=1  
%
%   .10:         .347        .119
%   .05:         .463        .146
%   .01          .739        .216
%
function test = kpss(x,w,opt)
T = size(x,1);
if nargin == 2,
   opt = 0;
end
if opt == 0,
   e = center(x);
   else
   [par, e] = ols(x,[ones(T,1) (1:T)']);
end
prod = zeros(w,1);
for j = 1:w
   prod(j) = e(j+1:T)'*e(1:T-j);
end
s2 = e'*e + 2*(1-(1:w)/(w+1))*prod;
S = cumsum(e);
test = T^(-1)*(S'*S)/s2;

