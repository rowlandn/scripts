function [name_matrix, par_matrix, data_matrix] = sim_dataset(fnamepattern, channel, npars)
%SIM_DATASET reads genesis files from a multiparameter simulation 
% 	and assembles them into sorted matrices based on the file names. 
%
% Usage: [name_matrix, par_matrix, data_matrix] = sim_dataset(fnamepattern, channel, npars);
%
% name_matrix contains the sorted filenames.
% param_matrix contains the sorted parameter values.
% data_matrix contains the sorted traces.
%
% Example:
%	[sim1_names, sim1_pars, sim1_traces] = sim_dataset('data/*.bin', 1, 3);
%
% In this example, all files with the .bin extension within the 
%	data/ folder will be loaded and the trace in channel 1 put into 
%	the matrix sim1_traces.
%	The filenames will be placed in sim1_names, and the values for the
%	three parameters into sim1_pars.
%	Each matrix will be sorted according to the param values 
%	in increasing order by par1/par2/par3

% get list of file names using the function forAllFiles, written by
%	Cengiz Gunay. Make sure it's in your path!
[name_matrix, temp] = forAllFiles(fnamepattern, 'sprintf', 1, '%s');
clear temp;	% not useful in this case.

% Parse names into parameter values. 
% ***RULES FOR FILE NAMES***
%	1. Param values must be included in the file name.
%	2. The file name must start with the first parameter value.
%	3. All param values must be followed by '_' e.g. '100_'
%	4. All but the first param value must come after '_' e.g. '_100_'
%	5. The params will be read and assigned in the order in which they
%		appear in the file name.
% For example: 100_HzSTN_10_HzSTR_0_HzGP_042204.bin

for n = 1:size(name_matrix, 2)
	tstr = char(name_matrix(n));
	start = strfind(tstr, '/');
	if size(start, 2) == 0
		start = 1; %using current directory, so params start right away.
	else 
		start = start(size(start, 2)) + 1;
		% first param starts after last '/' in file name.
	end
	sep = strfind(tstr, '_');

	% If the path name contains any underscores, skip them
	for k = 1:size(sep, 2)
		if sep(1) <= start
			sep(1) = [];
		end
	end
	
	% Get parameter values
	k = 1;
	m = 1;
	while k <= size(sep,2)
		last = sep(k) - 1;
		parval = str2num(tstr(start:last));
		if size(parval,2) > 0
			% parval = numeric
			par_matrix(n,m) = parval;
			m = m + 1;
		end
		if m > npars
			break
		end
		k = k + 1;
		start = last + 2;
	end
end

% make sorting priority list
for k = 1:npars
	tmat(1,k) = k;
end

% sort parameter matrix and return sorting index vector
[par_matrix, index] = sortrows(par_matrix, tmat);

% sort name array to match params matrix
for n = 1:size(name_matrix, 2)
	tempvec(n) = name_matrix(index(n));
end
name_matrix = tempvec;

% load traces in sorted order
for n = 1:size(name_matrix, 2)
	data_matrix(:,n) = readgenesis(char(name_matrix(n)), channel);
end

% Change scale from volts to millivolts
data_matrix = data_matrix * 1000;

clear n;
clear k;
clear m;
clear start;
clear tstr;
clear sep;
clear last;
clear parval;
clear tmat;
clear index;
clear tempvec;

