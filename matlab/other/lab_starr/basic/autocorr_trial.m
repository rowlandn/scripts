function [spk_n,ind_n] = autocorr_trial(A,lag)
% function [spk_n,ind_n] = autocorr_trial(A,lag) computes the autocorrelation
% of input spike time data stream A (using lag samples)
%
%	Used to compile spike count autocorrs over multiple separate trials
%	ind_n contains count of index spikes used for each lag bin
%	
%	Assumes spike times in A are in sec
%	spk_n & ind_n are in msec
%
% Written 11/28/05 by RST
%

n_sp = length(A);

if nargin < 2,lag = 500;end
spk_n = zeros(lag,1);
ind_n = zeros(lag,1);

tmp_sp = 1+round( 1000 .* A );        % put spk time in msec starting @ 1 msec

for i=1:n_sp-1
	for j=(i+1):n_sp;
		d = tmp_sp(j)-tmp_sp(i);
		if d <= lag
			spk_n(d) = spk_n(d)+1;
		end
	end
	l = min( (tmp_sp(end)-tmp_sp(i)), lag);	% interval available for this index spike
	ind_n(1:l) = ind_n(1:l)+1;	% increment index spike count for this interval
end
