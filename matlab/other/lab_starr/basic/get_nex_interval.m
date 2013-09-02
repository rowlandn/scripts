function 	interval = get_nex_interval( fname, intname);
% function 	interval = get_nex_interval( fname, intname)
%
% Function to return interval data from a NEX file from variables that match the
% indicated name_pattern
%
% INPUT:
%   filename - NEX File to query
%	intname - string name of variable to return
%
% OUTPUT:
%   interval - start & stop times of the matching interval if found, empty
%				otherwise
%
	
INTERVAL_TYPE = 2;	% Works only on interval-type variable

if exist(fname) ~= 2
	error( ['Unable to find the file named:  ' fname ] );
end

[nvar, varname, types] = nex_info(fname,0);

var_ind = strmatch(intname,varname);

if types(var_ind) == INTERVAL_TYPE
	[n, interval(1), interval(2)] = nex_int(fname, intname,0);
else
	interval = [];
end


