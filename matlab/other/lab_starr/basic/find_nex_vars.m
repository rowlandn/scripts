function 	[var_inds, var_names] = find_nex_vars( fname, var_type, name_pattern);
% function 	[var_inds, var_names] = find_nex_vars( fname, var_type, name_pattern)
%
% Function to find the indices of NEX variables of var_type that match the
% indicated name_pattern
%
% INPUT:
%   filename - NEX File to query
%	var_type - types 0 - 5 for different NEX variable types
%	name_pattern - optional string  for variable name matching
%			Default = all variable names.
%			See also regexp on string-pattern matching in matlab
% OUTPUT:
%   var_inds - indices of the variables whose names match the name pattern
%   var_names - names of the variables whose names match the name pattern
%
	
SPIKE_TYPE = 0;	% Spike-type timestamp
EVENT_TYPE = 1;	% Event-type timestamp
INTERVAL_TYPE = 2;
WAVEFORM_TYPE = 3;
POPVECT_TYPE = 4;
CONTINUOUS_TYPE = 5;

if ~exist('name_pattern','var')
	spkname_pattern = '\w*';	% Accept all variable names
end

if exist(fname) ~= 2
	error( ['Unable to find the file named:  ' fname ] );
end

[nvar, varname, types] = nex_info(fname,0);

var_inds = find(types == var_type);
n_vars = length(var_inds);

for i=1:n_vars		% Find spikes that match desired name pattern
	tmpname = deblank( varname(var_inds(i),:)) ;
	if isempty( regexpi( tmpname, name_pattern) )
		var_inds(i) = NaN;	% If names don't match, mark as bad
	end
end
var_inds( isnan(var_inds) ) = [];
var_names = varname( var_inds, : );
