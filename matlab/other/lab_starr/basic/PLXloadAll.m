function PLXdata = PLXLoadAll(filename,verbose)
% function PLXdata = PLXLoadAll(filename,verbose)
% 
%	Function to load all variable in a PLX file into PLXdata structure
%	Fields in PLXdata depend on contents of PLX file
%
%	inputs:	filename - file to be opened	
%			verbose - controls feedback on progress
% 

% Only look for these possible Plexon event types.  Ignore keyboard events
evnames = { 'Unknown' 'Event001' 'Event002' 'Event003' 'Event004' 'Event005' 'Event006' ...
	'Event007' 'Event008' 'Event009' 'Event010' 'Event011' 'Event012' ...
	'Event013' 'Event014' 'Event015' 'Event016' 'Strobed' 'Start' 'Stop' };
evchan2name = [ 2:17 ones(1,240) 18:20 ];

if ~exist('verbose','var') 
    verbose=false;
end

if ~exist(filename,'file')
	error(['Unable to open the file: ' filename]);
end

[p.fileName, p.version, p.snip_fs, Cmnt, Trd, NPW, PreTresh, SpikePeakV,...
	SpikeADResBits, SlowPeakV, SlowADResBits, p.Duration, p.DateTime] = ...
	plx_information(filename);
PLXdata = p;

% Find ts & wf's
[tscounts, wfcounts, evcounts] = plx_info(filename, 1);
evinds = find(evcounts);
evinds = evinds( find(evinds<260) );

tscounts_unsorted = tscounts(1,2:end);	% unsorted timestamps
wfcounts_unsorted = wfcounts(1,2:end);	% unsorted waveforms
tsinds_unsorted = find(tscounts_unsorted);
wfinds_unsorted = find(wfcounts_unsorted);

tscounts_sorted = tscounts(2:end,2:end);	% drop 1st column (blank) & row (unsorted spikes)
wfcounts_sorted = wfcounts(2:end,2:end);	% drop 1st column (blank) & row (unsorted spikes)
tsinds = find(tscounts_sorted);
wfinds = find(wfcounts_sorted);

% Get unsorted waveforms if present
for i=1:length(wfinds_unsorted)
	chan = wfinds_unsorted(i);
	unit = 0;  % unsorted waves

	varname = ['sig0' num2str(chan,'%02d') 'U'];
	[n, npw, ts, wave] = plx_waves_v(filename, chan, unit);
	eval(['PLXdata.' varname '.ts = ts;']);
	eval(['PLXdata.' varname '.snips = wave;']);
end

% Get sorted units (& waveforms if present)
for i=1:length(tsinds)
	[unit,chan] = ind2sub(size(tscounts_sorted),tsinds(i));
	[n, ts] = plx_ts(filename, chan, unit);
	
	varname = ['sig0' num2str(chan,'%02d') char(96+unit) ];
	eval(['PLXdata.' varname '.ts = ts;']);
	
	if ismember(tsinds(i),wfinds)	% get matching waveforms if present
		[n, npw, ts, wave] = plx_waves_v(filename, chan, unit);
		eval(['PLXdata.' varname '.snips = wave;']);
	end	
end

% Get events
for i=1:length(evinds)
	chan = evinds(i);
	varname = evnames{ evchan2name( chan ) };
	[n, ts, sv] = plx_event_ts(filename, chan);

	eval(['PLXdata.' varname '.ts = ts;']);
	if sv>0
		eval(['PLXdata.' varname '.sv = sv;']);
	end
	if strcmp(varname,'Strobed')
		PLXdata.trials = dysEvent2Trial(ts, sv);
	end
end

% NOT USED AT PRESENT
% Get continuous channels
% [n,names] = plx_adchan_names(filename);
% [n_cont,freqs] = plx_adchan_freqs(filename);
% ad_fs = unique(freqs);
% if length(ad_fs) ==1
% 	PLXdata.cont_fs = ad_fs;
% else
% 	PLXdata.cont_fs = freqs;
% end
% for i=1:n_cont
% 	[adfreq, n, ts, fn, ad] = plx_ad_v(filename, i);
% 	if n>0
% 		eval(['PLXdata.' deblank(names(i,:)) ' = ad;']);
% 	end
% end

if verbose		disp('PLX FILE SUCCESSFULLY LOADED'); end




