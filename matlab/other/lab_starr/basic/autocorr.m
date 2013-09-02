function [y,t] = autocorr(A,lag,ave_width,sample_length)
% Function [y t] = autocorr(A,lag_time,ave_width,sample_length) computes the autocorrelation
% of input ISI data stream A (using lag samples), and smoothes the positive half
% of it with the moveing average technique, using mov_ave as the averaging width.
% The output t is the length of time per sample (for autocorrelation plots).
%
% Sample call:
% 		[B t] = autocorr(A,1000,5,1);
% 
% This will calculate the autocorrelation function based on ISI input stream A up
% to a maximal lag of 1000 ms, smoothes the data with a 5 point moving average 
% technique and bases its t output on the assumption that the units of A were 1 ms 
% long.
%
% Default values:  lag = 500, ave_width = 1, sample length = 1
%
% Written 12/21/99 by Thomas Wichmann
% Modified 10/13/03 by RST

if nargin < 4,sample_length = 1;end
if nargin < 3,ave_width = 1;end
if nargin < 2,lag = 500;end

B = spk_t2delta(A);									% conversion of A into delta function
y = xcorr(B,lag);									% calculation of autocorrelation
%y = y(((length(y)-1)/2):length(y));					% removal of negative half of autocorrelation
y = mov_ave(y',ave_width);							% smoothing
t = 0:sample_length:(sample_length*length(y)-1);	% construction of t