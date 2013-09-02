function 	events = nex_get_events( fname, evtname_pattern );
% function 	events = nex_get_events( fname, evtname_pattern )
%
% Function to return times at which a named event occurs
%
% INPUT:
%   filename - NEX File to query
%	evtname_pattern - pattern to match (using matlab regular experssions)
%
% OUTPUT:
%   events - structure of times (in secs) at which event occured, one per
%   event found
%
%

verbose = false;
EVENT_TYPE = 1;	% Works only on event-type timestamp

if exist(fname) ~= 2
	error( ['Unable to find the file named:  ' fname ] );
end

[nvar, varname, types] = nex_info(fname,0);

var_inds = find(types == EVENT_TYPE);
n_vars = length(var_inds);

for i=1:n_vars		% Find events that match desired name pattern
	tmpname = deblank( varname(var_inds(i),:)) ;
	if isempty( regexpi( tmpname, evtname_pattern) )
		var_inds(i) = NaN;	% If names don't match, mark as bad
	else
		[events(i).n, events(i).ts] = nex_ts(fname,varname(var_inds(i),:),verbose);
	end
	
end
