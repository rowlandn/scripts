function nex2isi()
% nex2isi()
% This program reads a NEX file created by Plexon spike sorter
% and writes an ISI text file similar to the one
% created by AlphOmega MSD
%
% Created by RST, 2002-04-05
% modified 2002-10-15 to overcome apparent error in spike sorter


global RTPATH	% starting directory containing NEX file
global fs       % Sampling rate for processed EMG & spike SDF
fudge = 5;      % Conversion to make sdf = sp/sec
fs = 1000;      % sampling rate for EMG & SDF
p = path;

[nexname, RTPATH] = uigetfile('*.nex', 'Select spike time file (NEX)');

if (nexname ~= 0)
    cd(RTPATH);
    [nvar, varname, types] = nex_info(nexname);
    if nvar > 1
        error([num2str(nvar) ' spikes in ' nexname '!!!  I don''t know how to process > 1 spike yet']);
    end
    [spk.n, spk.t] = nex_ts(nexname,varname);
else
    error(['I can''t find the NEX file:  ' nexname ' in ' RTPATH]);
end

fname = strrep(nexname,'.nex','');
isi = 1000 .* diff(spk.t);
bad = find(isi < 1);
if length(bad)
    warning(['Found ' int2str(length(bad)) ' ISI''< 1 msec!!  Correcting...']);
    spk.t( bad+1 ) = [];
    spk.n = length(spk.t);
end

% Output a txt file for Labview or Matlab ISI analyses
write_isi_txt(spk.t,fname);
% RORY: present burstiness data
[L, isi_mean, range, bins, histo] = roryburst(spk.t)
