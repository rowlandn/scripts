function RunSpikeMeasureOR(show,meta)
% RunSpikeMeasureOR(show,meta)
% This program computes measures of action potential shape from
% DDT files around times of spikes in NEX files.
%
% NEX & DDT files must be in the same directory.
% Created by RST, 2005-08-02
% Modified to read DDT files 2007-10-25
%
%	Run within a directory and this program will process all spikes in all 
%	NEX files in the directory
%
%	Inputs:
%			show - if 'false', then graphical output is not generated
%			(default = false)
%			meta - if true, then graphical output is saved as a meta-file
%	
%	Averages across the whole file

global VERBOSE
global PLOT_DUR
global PLOT_PRE

PLOT_DUR = 2.5;	% in msec
PLOT_PRE = 1.0;	% in msec
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

if ~exist('meta','var')
	meta = false;
end

% Edit this pattern to select time-stamp variables w/ specific names
spkname_pattern = '\w*';	% Accept all units

cd(uigetdir);

FileLst = dir('*.nex');
if(isempty(FileLst))
	str = pwd;
	error(['Found no NEX files in current directory - ' str ]);
end

outfid = write_text_header();	% Subfunction below

n = 0;	% Count of units processed
wf_mat = [];
% For each file found in directory...
for i=1:length(FileLst)

	fname = FileLst(i).name;
	% Find variables of interest in file
	[var_inds,var_names] = find_nex_vars( fname, SPIKE_TYPE, spkname_pattern);
	if isempty(var_inds)
		continue
	end
	
	% For each unit in a file...
	for j=1:length(var_inds)
		n = n+1;

		% get spike time-stamps
		[spk_n, spk_ts] = nex_ts( fname, var_names(j,:),VERBOSE);

		% get spike wave-forms
		[fs, spk_wf] = ddt_wf( strrep(fname,'.nex','.ddt'),var_names(j,:),spk_ts);

		% save name of file & unit
		unitname = deblank( strrep(var_names(j,:),'_wf','') );
		display(['Processing...' fname ' --: ' unitname]);
		
		% get waveform stats
		spk{n} = SpikeMeasure( spk_wf, fs, fname, unitname, show, meta);
		
		% write stats to file
		write_text(outfid, spk{n});
		
		% save mean wf to matrix
		wf_mat = [wf_mat ; spk{n}.mean'];
	end
end
fclose(outfid);
if exist('spk','var')
	save SpikeMeasure spk
	fname = 'SpikeWF.txt';
	dlmwrite(fname,wf_mat,'\t');
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

