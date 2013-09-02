function [y,t] = delta2autocorr(A,lag)
% Function [y t] = delta2autocorr(A,lag_time) computes the autocorrelation
% of input delta function (using lag samples).
% The output t is the lag time.
%
% Sample call:
% 		[y, t] = delta2autocorr(A,1000);
% 
% This will calculate the autocorrelation function based on 
% delta function input stream A up
% to a maximal lag of 1000
%
% Default values:  lag = 500
%
% Written 2005/07/30 by RST

if nargin < 2,lag = 500;end

[y,t] = xcorr(A,lag,'none');		% calculation of autocorrelation
y(lag+1) = 0;						% remove mid-point coeff=1
y = 1000*y/length(A);			% calibrate to spike/sec
