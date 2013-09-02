function [y,t] = Spike_autocorr(A,lag)
% Function [y t] = Spike_autocorr(A,lag_time) computes the autocorrelation
% of input ISI data stream A (using lag samples).
% The output t is the lag time.
%
% Sample call:
% 		[B t] = autocorr(A,1000);
% 
% This will calculate the autocorrelation function based on spike time input stream A up
% to a maximal lag of 1000
%
% Default values:  lag = 500
%
% Written 12/21/99 by Thomas Wichmann
% Modified 10/13/03 by RST

if nargin < 4,sample_length = 1;end
if nargin < 3,ave_width = 1;end
if nargin < 2,lag = 500;end

B = spk_t2delta(A);					% conversion of A into delta function
[y,t] = xcorr(B,lag,'none');		% calculation of autocorrelation
y(lag+1) = 0;						% remove mid-point coeff=1
y = 1000 * y ./ length(A);			% calibrate to spike/sec
%y = y(((length(y)-1)/2):length(y));	% removal of negative half of autocorrelation
%y = mov_ave(y',ave_width);				% smoothing
