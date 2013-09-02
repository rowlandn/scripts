function [sig_inds, sig_thresh] = find_freq_peaks( pow, freq, ...
		sig, npts, cntl_frq, srch_frq );
% Function to find sig peaks in a spectrum
% Inputs:
%		pow - spectrum vector
%		freq - matching frequency vector
%		sig - base alpha level for significance
%		npts - test significance of local maxima w/in npts points
%		cntl_frq - range of frequencies to use as control data
%		srch_frq - range of frequencies to search for sig peaks
%
% Outputs:
%		sig_inds - indices into pow of significant peaks
%		sig_thresh - threshold computed for significance
%
%	Written 11/28/05 by RST
%
	
	cntl_ind = find(freq>=cntl_frq(1) & freq<cntl_frq(2) );
	srch_inds = find( freq>=srch_frq(1) & freq<srch_frq(2) );
	p = sig/length(srch_inds);	% Bonferroni correction for mult comparisons

	% For spectrum - compute significance 
	sig_thresh = norminv(1-p) * std( pow(cntl_ind) ) + 1;
	
	% & find significant peaks w/in search range
	sig_inds = find_sig_peaks( pow,	sig_thresh, npts, srch_inds );
