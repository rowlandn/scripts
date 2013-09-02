function 	[var_inds, var_names] = find_nex_units( fname, spkname_pattern);
% function 	[var_inds, var_names] = find_nex_units( fname, spkname_pattern)
%
% Function to find the indices of NEX timestamp variables that match the
% indicated name_pattern
%
% INPUT:
%   filename - NEX File to query
%	spkname_pattern - string  for variable name matching
%			See also regexp on string-pattern matching in matlab
% OUTPUT:
%   var_inds - indices of the variables whose names match the name pattern
%   var_names - names of the variables whose names match the name pattern
%
	
SPIKE_TYPE = 0;	% Works only on Spike-type timestamp

if exist(fname) ~= 2
	error( ['Unable to find the file named:  ' fname ] );
end

[nvar, varname, types] = nex_info(fname,0);

var_inds = find(types == SPIKE_TYPE);
n_vars = length(var_inds);

for i=1:n_vars		% Find spikes that match desired name pattern
	tmpname = deblank( varname(var_inds(i),:)) ;
	if isempty( regexpi( tmpname, spkname_pattern) )
		var_inds(i) = NaN;	% If names don't match, mark as bad
	end
end
var_inds( isnan(var_inds) ) = [];
var_names = varname( var_inds, : );
