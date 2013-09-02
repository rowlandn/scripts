function RunBurstMeasureOR(show,meta)
% RunBurstMeasureOR(show,meta)
% This program measures burst morphology from NEX files.
%	Uses Legendy suprise method to detect bursts
%	Adjust suprise level & min burst number below.
%
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
global PLOT_PRE
global PLOT_PST

global Leg_suprise
global Leg_minspk
global FS
global MIN_SPK_N

Leg_suprise = 5;
Leg_minspk = 4;
FS = 1000;
MIN_SPK_N = 500;
MIN_BURST_N = 5;
PLOT_PRE = -0.2;	% in sec
PLOT_PST = 0.2;		% in sec

SPIKE_TYPE = 0;	% Spike-type timestamp
EVENT_TYPE = 1;	% Event-type timestamp
INTERVAL_TYPE = 2;
WAVEFORM_TYPE = 3;
POPVECT_TYPE = 4;
CONTINUOUS_TYPE = 5;
fac = 1.5;
local_bst_len = 0;

VERBOSE=false;

if ~exist('show','var')
	show = true;
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
n_files = length(FileLst);

outfid = write_text_header();	% Subfunction below

n = 0;	% Count of units processed

avg_len = length(round(PLOT_PRE*1000):round(PLOT_PST*1000));
avg = nan(n_files,avg_len);
ci =  nan(n_files,avg_len);
% For each file found in directory...
for i=1:n_files

	fname = FileLst(i).name;
	% Find variables of interest in file
	[var_inds,var_names] = find_nex_vars( fname, SPIKE_TYPE, spkname_pattern);
	if isempty(var_inds)
		continue
	end
	
	% For each unit in a file...
	for j=1:length(var_inds)

		info = nex_info_rst(fname,VERBOSE);
		data_dur = info.dur;
		
		% get spike time-stamps
		[spk_n, spk_ts] = nex_ts( fname, var_names(j,:),VERBOSE);
		if spk_n<MIN_SPK_N
			display(['Skipping...' fname ' --: ' unitname '..: Fewer than ' num2str(MIN_SPK_N) ' spikes']);
			continue
		end			
		n = n+1;

		% save name of file & unit
		unitname = deblank( strrep(var_names(j,:),'_wf','') );
		display(['Processing...' fname ' --: ' unitname]);
		
		% get burst stats
		spk{n} = legendy_new3(round(1000*diff(spk_ts)), fac, FS, Leg_minspk, local_bst_len, Leg_suprise);
		
		spk{n}.fname = fname;
		spk{n}.unitname = unitname;
		% make burst-triggered averages
		if spk{n}.num_bursts>MIN_BURST_N
			[avg(n,:),ci(n,:)] = MakePeriBurstAvg(spk{n},spk_ts,data_dur,show,meta);
		else
			display(['..: Fewer than ' num2str(MIN_BURST_N) ' bursts. Skipping MorphMeas.']);
		end
		spk{n}.avg = avg(n,:);
		spk{n}.ci = ci(n,:);
		
		% write stats to file
		write_text(outfid, spk{n});
		
	end
end
fclose(outfid);
if exist('spk','var')
	save BurstMeasure spk
	fname = 'BurstAvg.txt';
	dlmwrite(fname,avg,'\t');
	fname = 'BurstCI.txt';
	dlmwrite(fname,ci,'\t');
end


%------------------------------------------------------
% Subfunction to write a line of data to output file
function write_text(outfid, spk)

	fprintf(outfid,'%s\t%s\t%d\t', ...
		spk.fname, spk.unitname, spk.num_bursts);
	
	fprintf(outfid,'%.3f\t%.3f\t', ...
		spk.mean_spikes_per_burst, spk.mean_intra_burst_frequency );

	fprintf(outfid,'%.3f\t%.3f\n', ...
		spk.proportion_spikes_in_bursts, spk.proportion_time_in_bursts);

return

%------------------------------------------------------
% Subfunction to open output file and print header line
function outfid = write_text_header()
	fname = 'BurstMeasure.txt';
	outfid = fopen(fname,'w');
	if(outfid == -1)
       error(['Unable to open...' fname ]);
	end
	fprintf(outfid,'filename\tunitname\tNum_bursts\t');
	fprintf(outfid,'SpikesPerBurst\tFreqInBurst\t');
	fprintf(outfid,'FracSpikeInBurst\tFracTimeInBurst\n');

return

