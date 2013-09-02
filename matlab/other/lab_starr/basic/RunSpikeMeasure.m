function RunSpikeMeasure(show)
% RunSpikeMeasure()
% This program computes measures of action potential shape from mean
% snippets extracted from NEX files
% Created by RST, 2005-08-02
%
%	Run within a directory and this program will process all spikes in all 
%	NEX files in the directory
%
%	Inputs:
%			show - if 'false', then graphical output is not generated
%			(default = true)
%	
%	Averages across the whole file

global VERBOSE
SPIKE_TYPE = 0;	% Spike-type timestamp
EVENT_TYPE = 1;	% Event-type timestamp
INTERVAL_TYPE = 2;
WAVEFORM_TYPE = 3;
POPVECT_TYPE = 4;
CONTINUOUS_TYPE = 5;

VERBOSE=false;

if ~exist('show','var')
	show = false;
end

% Edit this pattern to select time-stamp variables w/ specific names
%spkname_pattern = '\w*';	% Accept all units
spkname_pattern = 'Snip\w*[abcd]';	% Kevin's pattern

cd(uigetdir);

FileLst = dir('*.nex');
if(isempty(FileLst))
	str = pwd;
	error(['Found no NEX files in current directory - ' str ]);
end

outfid = write_text_header();	% Subfunction below

n = 0;	% Count of units processed
% For each file found in directory...
for i=1:length(FileLst)

	fname = FileLst(i).name;
	% Find variables of interest in file
	[var_inds,var_names] = find_nex_vars( fname, WAVEFORM_TYPE, spkname_pattern);
	if isempty(var_inds)
		continue
	end
	
	% For each unit in a file...
	for j=1:length(var_inds)
		n = n+1;

		% get spike waveforms
		[fs, spk_n, spk_ts, spk_wfn, spk_wf] = nex_wf( fname, var_names(j,:),VERBOSE);
		
		% save name of file & unit
		unitname = deblank( strrep(var_names(j,:),'_wf','') );
		display(['Processing...' fname ' --: ' unitname]);
		
		% get waveform stats
		spk{n} = SpikeMeasure( spk_wf, fs, fname, unitname, show);
		
		% write stats to file
		write_text(outfid, spk{n});
	end
end
fclose(outfid);
if exist('spk','var')
	save SpikeMeasure spk
end


%------------------------------------------------------
% Subfunction to write a line of data to output file
function write_text(outfid, spk)

	fprintf(outfid,'%s\t%s\t%d\t%s\t', ...
		spk.fname, spk.unit, spk.n, spk.report);
	
	fprintf(outfid,'%.3f\t%.3f\t', ...
		spk.max, spk.min );

	fprintf(outfid,'%.3f\t%.3f\t', ...
		spk.start2max, spk.start2stop );

	fprintf(outfid,'%.3f\t%.3f\n', ...
		spk.min2max, spk.area);

return

%------------------------------------------------------
% Subfunction to open output file and print header line
function outfid = write_text_header()
	fname = 'SpikeMeasure.txt';
	outfid = fopen(fname,'w');
	if(outfid == -1)
       error(['Unable to open...' fname ]);
	end
	fprintf(outfid,'filename\tunitname\tNum_Spikes\tReport\t');
	fprintf(outfid,'max\tmin\t');
	fprintf(outfid,'start2max\tstart2stop\t');
	fprintf(outfid,'min2max\tarea\n');

return

