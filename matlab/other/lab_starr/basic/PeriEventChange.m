function [s,cntl_mean,sig_thr] = PeriEventChange(data,cntl_inds,test_start,alpha,contig)
%function s = PeriEventChange(data,cntl_inds,test_start,alpha,contig)
%
%	Finds significant changes in peri-event data using 
%	statistics.
%
% Inputs:
%	data - vector of independent scalar values
%	cntl_inds - indices into data for control period
%	test_start - index at which to start testing
%	alpha - p-value for significance
%	contig - number of contiguous significant points
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
b = ones(1,contig) ./ contig;	% moving avg filter

cntl_data = trim(data(cntl_inds),6);
cntl_mean = nanmean( cntl_data );
data = data - cntl_mean;
data_sm = filter(b,1,data);
cntl_data_sm = trim(data_sm(cntl_inds),6);

% make testing data and matching smooted data
test = data( test_start:end );
test_sm = data_sm(test_start+contig-1:end);

p = alpha/length(cntl_data);	% Compensate for multcomp
sig_thr = norminv(1-p) * nanstd( cntl_data );

p = alpha/(length(cntl_data)/contig);	% Compensate for multcomp
sig_thr_sm = norminv(1-p) * nanstd( cntl_data_sm );

% Find threshold crossing indices
on_inds = find(test > sig_thr | test < -1*sig_thr );
on_inds_sm = find(test_sm > sig_thr_sm | test_sm < -1*sig_thr_sm );

% onset is intersection of the two
s.on_ind = min( intersect(on_inds,on_inds_sm) ) + test_start-1;

if ~isempty( s.on_ind )
	s.sgn = sign(data(s.on_ind));

	% search smoothed data after established onset
	test = data(s.on_ind+1:end);
	test_sm = filter(b,1,data);
	test_sm = test_sm(s.on_ind+contig:end);
	
	% find returns to baseline
	zero_inds = find(test < sig_thr & test > -1*sig_thr);
	zero_inds_sm = find(test_sm < sig_thr_sm & test_sm > -1*sig_thr_sm);
	zero_ind = min( intersect(zero_inds,zero_inds_sm) ) + s.on_ind;
	
	% find change in sign
	switch_ind = min(find( sign(test) ~= s.sgn )) + s.on_ind;
	
	s.off_ind = min( [zero_ind switch_ind] );
	if ~isempty(s.off_ind)
		s.mean_change = mean( data(s.on_ind:s.off_ind));
	else
		s.mean_change = mean( data(s.on_ind:end));
	end

else
	s.sgn = [];
	s.off_ind = [];
	s.mean_change = [];
end



