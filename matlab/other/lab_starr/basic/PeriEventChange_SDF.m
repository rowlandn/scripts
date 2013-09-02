function [s,cntl_mean,sig_thr] = PeriEventChange_SDF(data,cntl_inds,test_start,alpha,contig)
%function s = PeriEventChange_SDF(data,cntl_inds,test_start,alpha,contig)
%
%	Finds significant changes in peri-event SDF using 
%	statistics.
%
% Inputs:
%	data - SDF vector
%	cntl_inds - indices into data for control period
%	test_start - index at which to start testing
%	alpha - p-value for significance
%	contig - number of contiguous significant points required
%
% Outputs:
%	s.on_ind - index into data for onset of significant response (empty if
%				no significant change found)
%	s.off_ind - index into data for offset of significant response
%	s.sgn		- sign of the significant change
%	s.mean_change - mean change from control across duration of change (till
%		offset or end of data)
%	cntl_mean - mean control level
%	sig_thr - threshold to determine significance
%
%	RST 2005-08-10
%	
cntl_data = data(cntl_inds);
cntl_mean = mean( cntl_data );
data = data - cntl_mean;

% make testing data 
test = data( test_start:end );

p = alpha/(length(cntl_data)/contig);	% Compensate for multcomp
sig_thr = norminv(1-p) * std( cntl_data );

% Find threshold crossing indices
thr_inds = find(test > sig_thr | test < -1*sig_thr );

% Compute o'th-order differences - only o-successive indices will place
% o'th order difference o indices apart
o = contig-1;
odiff = thr_inds((1+o):end)-thr_inds(1:(end-o));

% Find first index into odiff that is o indices apart
s.on_ind = thr_inds( min(find( odiff == o )) ) + test_start-1;

if ~isempty( s.on_ind )
	s.sgn = sign(data(s.on_ind));

	% search smoothed data after established onset
	test = data(s.on_ind+1:end);
	
	% find returns to baseline 
	thr_inds = find(test < sig_thr & test > -1*sig_thr);
	
	% check if it lasts o-contiguous points
	odiff = thr_inds((1+o):end)-thr_inds(1:(end-o));
	zero_ind = thr_inds( min( find(odiff==o) )) + s.on_ind;
	
	% find any change in sign > sig_thr
	switch_ind = min(find( sign(test) ~= s.sgn & abs(test)>sig_thr )) + s.on_ind;
	
	s.off_ind = min( [zero_ind switch_ind] );
	if isempty(s.off_ind)
		s.mean_change = mean( data(s.on_ind:end));
		s.int_change = sum( data(s.on_ind:end))/1000;
	else
		s.mean_change = mean( data(s.on_ind:s.off_ind));
		s.int_change = sum( data(s.on_ind:s.off_ind))/1000;
	end

else
	s.sgn = [];
	s.off_ind = [];
	s.mean_change = [];
	s.int_change = [];
end



