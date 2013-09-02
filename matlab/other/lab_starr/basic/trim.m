function y = trim(x, A)
% y = trim(x, A)
%
%	Trim input data x according to A*(median absolute deviation)
%	Assign NaN to outliers.
%	
%	Inputs:
%		x - data vector
%		A - multiplier 
%
%	MAD is scaled to approximate Stdev, so A*MAD robustly 
%	approximates an A*Stdev threshold for outlier
%	
% 

y = x;
md = median( x(~isnan(x)) );

thresh = A * 1.483 * median( abs( x(~isnan(x)) - md ) );

lowcut = md - thresh;
hicut = md + thresh;

outlier = find( x < lowcut | x > hicut );

y(outlier) = NaN;
