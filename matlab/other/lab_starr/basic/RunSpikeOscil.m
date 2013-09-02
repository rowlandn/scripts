function RunSpikeOscil(show)
% RunSpikeOscil()
% This program analyses spontaneous spike trains in NEX for oscillations
% Created by RST, 2005-07-30
%
%	Input:
%		show - controls whether graphical output is produced for each cell
%		  default = true.
%
%	Run within a directory and this program will process all spikes in all 
%	nex files in the directory
%	
%	Assumes spontaneous data 

% Global defines used by the other functions called by this one
global SEG_PWR
global FS
global N_SHUF
global SIG_OSC
global CNTL_FRQ
global MAX_N

global SRCH_LO
global SRCH_HI
global AUTOCORR_LAG
global VERBOSE
i=1;
j=1;
SEG_PWR = 11;	% Segment length - specified as power of 2
FS = 1000;
SIG_OSC = 0.01;	% or OSC_IND > SIG_OSC
N_SHUF = 100;	% number of ISI shuffles for control spectra
CNTL_FRQ = [300 500];	% Frq range used to compute spectral SD
% CNTL_FRQ = [270 300];
SRCH_LO = 0.25;	% low freq for search for significant oscillations
SRCH_HI = 200;	% top freq for search for significant oscillations
MAX_N = 7;		% find max's in spectra over MAX_N symmetric points (must be odd) 9 => 3.5 Hz separation
AUTOCORR_LAG = 3000;	% For constructing autocorrelogram
VERBOSE=false;

N_OscPks = 3;	% max # of significant autocorr pks reported in txt output

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

outfid = write_text_header(N_OscPks);	% Subfunction below

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
		
		% save name of file & unit
		unitname = ...
			deblank( strrep( fname,'.nex',['_' var_names(j,:)] ) );
		display(['Processing...' unitname]);
		
		% get oscillation stats
		spk{n} = SpikeOscil( spk_t, unitname, show );
		
		% write stats to file
		write_text(outfid, spk{n}, N_OscPks);
	end
end
fclose(outfid);
save SpikeOscil spk


%------------------------------------------------------
% Subfunction to write a line of data to output file
function write_text(outfid, spk, N_OscPks)

	fprintf(outfid,'%s\t%.3f\t%.3f\t', ...
		spk.name, max(spk.t)/1000, spk.sig_thresh );

	% Report significant acor peaks
	for i = 1:N_OscPks
		if i <= length(spk.sig_inds)

			fprintf(outfid,'%.3f\t%.3f\t',...
				spk.freq( spk.sig_inds(i) ),...
				spk.pow_comp( spk.sig_inds(i) ) );
		else
			fprintf(outfid,'\t\t');
		end
	end
	if length(spk.sig_inds) > N_OscPks
		fprintf(outfid,'yes\n');
	else
		fprintf(outfid,'no\n');
	end	
	return

%------------------------------------------------------
% Subfunction to open output file and print header line
function outfid = write_text_header(N_OscPks)
	fname = 'SpikeOscil.txt';
	outfid = fopen(fname,'w');
	if(outfid == -1)
       error(['Unable to open...' fname ]);
	end
	fprintf(outfid,'name\trecord_dur\tsig_thresh\t');
	for i = 1:N_OscPks
		% For max reported signif acorr pks, freq & normalized power
		fprintf(outfid,'AC%d_freq\tAC%d_pow\t',i,i);
	end
	fprintf(outfid,'AC_more?');	% Header for missed significant ACorr peaks
	fprintf(outfid,'\n');		% EOL

	return

