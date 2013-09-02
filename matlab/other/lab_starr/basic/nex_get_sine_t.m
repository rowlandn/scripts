function 	[Cycle, MDCycle, Angl_mat, Torq_mat] = nex_get_sine_t( fname, start_t, end_t );
% function 	[Cycle, MDCycle, Angl_mat, Torq_mat] = nex_get_sine_t( fname, start_t, end_t );
%
% Function to return times at which sines crossed zero
%
% INPUT:
%   filename - NEX File to query
%	start_t - start time of valid data (trim earlier cycles)
%
% OUTPUT:
%   Cycle - times (in secs) at which angle variable does positive zero
%   crossing
%	MDCycle - duration of cycles (in sec)
%	Angl_mat - Angle x cycle matrix for all valid cycles
%	Torq_mat - Torq x cycle matrix matching Angl_mat
%
%	RST 2005-08-30
%
verbose = false;

if exist(fname) ~= 2
	error( ['Unable to find the file named:  ' fname ] );
end

[adfs, n, ts, fn, Angl] = nex_cont(fname,'Angl001',verbose);

if n==0
	Cycle = [];
	warning( ['Unable to find Angl001 variable in:  ' fname ] );
	return
end

[adfs, n, ts, fn, Torq] = nex_cont(fname,'Torq001',verbose);

if n==0
	Cycle = [];
	warning( ['Unable to find Torq001 variable in:  ' fname ] );
	return
end

Angl = Angl - median(Angl);

Angl_ts = (0:(length(Angl)-1))' ./ adfs;

% Find positive zero crossings in angle = starts of each cycle
Cycl_ind = find( Angl(1:end-1)<0 & Angl(2:end)>=0 ) + 1;
Cycle = Angl_ts(Cycl_ind-1);

% Remove cycles w/ a length different from the mean cycle length
Cycl_len = diff(Cycle);
[Cycl_tlen,drp] = trimsd(Cycl_len, 1);	% outliers replaced w/ NaN
MDCycle = median(Cycl_len);

Cycle(drp) = [];

% Remove cycles starting outside valid period (i.e., before DBS test or close to end of trial)
drp = find(Cycle < start_t | Cycle > end_t);
Cycle(drp) = [];

ncycl = length(Cycle);
halfcycl = round(adfs*MDCycle/2);
cycl_inds = -1*halfcycl:halfcycl;
Angl_mat = zeros( ncycl,length(cycl_inds) );
Torq_mat = zeros( ncycl,length(cycl_inds) );
for i = 1:ncycl
	Angl_mat(i,:) = Angl( cycl_inds + round(adfs*Cycle(i)) );
	Torq_mat(i,:) = Torq( cycl_inds + round(adfs*Cycle(i)) );
end
