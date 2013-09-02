function 	[DBS_on, DBS_off] = nex_get_dbs( fname );
% function 	[DBS_on, DBS_off] = nex_get_dbs( fname )
%
% Function to return times at which DBS is turned on and off
%
% INPUT:
%   filename - NEX File to query
%
% OUTPUT:
%   DBS_on - times (in secs) at which DBS turned on
%   DBS_off - times (in secs) at which DBS turned off
%
verbose = false;

if exist(fname) ~= 2
	error( ['Unable to find the file named:  ' fname ] );
end

[n, DBSt] = nex_ts(fname,'DBSt',verbose);

if DBSt==0
	DBS_on = [];
	DBS_off = [];
	return
end

% Find indices for long gaps between shocks
DBS_inds = find( diff(DBSt) > 1);

% Make arrays of time DBS went 'on' and 'off'
DBS_on = [ min(DBSt) DBSt(DBS_inds+1) ];
DBS_off = [ DBSt(DBS_inds) max(DBSt) ];

% Check for outliers in DBS duration (i.e., overly short)
DBS_dur =  DBS_off - DBS_on ;
drp = find( DBS_dur < median(DBS_dur)-2*mad(DBS_dur));
DBS_on(drp) = [];
DBS_off(drp) = [];
OFF_dur =  DBS_on(2:end) - DBS_off(1:end-1) ;
drp = find( OFF_dur < median(OFF_dur)-2*mad(OFF_dur));
DBS_on(drp) = [];
DBS_off(drp) = [];
