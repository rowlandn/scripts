function [adfreq, w] = ddt_wf(filename, varname, spk_ts)
% Read waveforms from a .ddt file at spk_ts timepoints
%
% [adfreq, w] = ddt_wf(filename, varname, spk_ts)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   varname - variable name
%	spk_ts - spike timestamps
%           
% OUTPUT:
%   adfreq - sampling rate for ddt file
%   w - matrix of waveform a/d values [n nf] (in millivolts)

global PLOT_DUR
global PLOT_PRE

n = 0;
adfreq = 0;
ts = spk_ts;

if(nargin < 3)
   disp('3 input arguments are required')
   return
end

if(ischar(filename) == 0)
   disp('filename argument should be character arrays')
   return
end

if(ischar(varname) == 0)
   disp('varname argument should be character arrays')
   return
end

[nch,data_len,adfreq,data] = readddt(filename);
if nch>1
	error('ddt_wf doesn''t know how to process ddt files with >1 channel of data');
end

dt = 1/adfreq;
n_pts = round(PLOT_DUR/1000/dt);
n_pre = round(PLOT_PRE/1000/dt);
peri_inds = (-n_pre:(n_pts-n_pre-1))';
spk_inds = round(spk_ts .* adfreq);
drp = find(spk_inds+min(peri_inds)<1 | spk_inds+max(peri_inds)>data_len);
spk_inds(drp) = [];
n_spk = length(spk_inds);

ind_mat = repmat(peri_inds,1,n_spk) + repmat(spk_inds,n_pts,1);

w = data(ind_mat);
