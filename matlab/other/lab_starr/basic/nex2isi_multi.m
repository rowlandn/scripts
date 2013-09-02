function nex2isi_multi()
% nex2isi_multi()
% This program reads a NEX file created by Plexon spike sorter
% and writes ISI text files for every sorted spike
%
% Created by RST, 2002-04-05
% modified 2002-10-15 to overcome apparent error in spike sorter
% modified 2005-06-29 to process multiple spikes
% modified 2007-02-27 to create output txt in same directory as nex files


global RTPATH	% starting directory containing NEX file
global fs       % Sampling rate
fs = 1000;      % sampling rate for EMG & SDF
p = path;
SPIKE_TYPE = 0;
n = 1;
isi_batching = 0;	% used for something in Nexload, added here merely for compatibility

[nexname, nexpath] = uigetfile('*.nex', 'Select spike time file (NEX)');

Nexdir(n).name = nexname;	% More compatibility w/ Nexload
if (nexname ~= 0)
    cd(nexpath);
	Nexload		% Run Rory's script for reading a nex file
else
    error(['I can''t find the NEX file:  ' nexname ' in ' nexpath]);
end


spkvar_inds = find(Nexvar.Type == SPIKE_TYPE);
n_spkvars = length(spkvar_inds);
for i = 1:n_spkvars		% Run for each spike timestamp
	
	fname = deblank( strrep(nexname,'.nex',['_' Nexvar.Name(spkvar_inds(i),:)] ) );
	ts = Nexdata.TS{ spkvar_inds(i) } ;
	
	
	isi = 1000 .* diff(ts);
	bad = find(isi < 1);
	if length(bad)
		warning(['Found ' int2str(length(bad)) ' ISI''< 1 msec!!  Correcting...']);
		ts( bad+1 ) = [];
	end

	% Output a txt file for Labview or Matlab ISI analyses
	disp( ['Writing ISI file: ' fname ]);
	write_isi_txt(ts,fname);
	% RORY: present burstiness data
%	[L, isi_mean, range, bins, histo] = roryburst(ts)
end

