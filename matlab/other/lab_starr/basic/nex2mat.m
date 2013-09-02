function nex2mat(show)
% nex2mat(show)
% This program converts nex files into mat files 1per spike channel

global VERBOSE
VERBOSE=false;

if ~exist('show','var')
	show = true;
end

% Edit this pattern to select time-stamp variables w/ specific names
%spkname_pattern = 'Snip\w*[abcd]';	% Kevin's pattern
spkname_pattern = '\w*';	% Accept all units

cd(uigetdir);

FileLst = dir('*.nex');
if(isempty(FileLst))
	str = pwd;
	error(['Found no NEX files in current directory - ' str ]);
end

n = 0;	% Count of units processed
% For each file found in directory...
for i=1:length(FileLst)

	fname = FileLst(i).name;
	% Find variables of interest in file
	[var_inds,var_names] = find_nex_units( fname, spkname_pattern);
	if isempty(var_inds)
		continue
	end
	
	% For each unit in a file...
	for j=1:length(var_inds)
		n = n+1;

		% get spike times
		[spk_n, spk_t] = nex_ts( fname, var_names(j,:),VERBOSE);
		spk_t = spk_t';
		
		% save name of file & unit
		unitname = ...
			strrep( fname,'.nex',['_' deblank(var_names(j,:))] );
		display(['Processing...' unitname]);
		
		% write spike times to mat & txt files
		save( [unitname '.mat'],'spk_n','spk_t');
		save( [unitname '.txt'],'spk_t','-ASCII','-TABS','-DOUBLE');
	end
end
