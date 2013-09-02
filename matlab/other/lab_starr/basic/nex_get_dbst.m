function 	DBSt = nex_get_dbst( fname )
% function 	DBSt = nex_get_dbst( fname )
%
% Function to return times at which each DBS stim is delivered 
%
% INPUT:
%   filename - NEX File to query
%
% OUTPUT:
%   DBSt - times (in secs) at which each stim is delivered
%
verbose = false;

if exist(fname) ~= 2
	error( ['Unable to find the file named:  ' fname ] );
end

[n, DBSt] = nex_ts(fname,'DBSt',verbose);

if DBSt==0
	DBSt = [];
	return
end

% Check for indices for short gaps between shocks
short_inds = find( diff(DBSt) < 0.001);

if ~isempty(short_inds)
	warning('Found short intervals in DBSt - correcting');
	DBSt(short_inds+1) = [];
end

