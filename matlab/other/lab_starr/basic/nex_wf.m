function [adfreq, nf, ts, n, w] = nex_wf(filename, varname, verbose)
% nex_wf(filename, varname): Read waveform variable from a .nex file
%
% [adfreq, nf, ts, n, w] = nex_wf(filename, varname, verbose)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   varname - variable name
%
%           
% OUTPUT:
%   n - number of data points in each waveform
%   ts - array of waveform timestamps (in seconds)
%   nf - number of waveforms
%   w - matrix of waveform a/d values [n nf] (in millivolts)

n = 0;
adfreq = 0;
ts = 0;
nf = 0;
w = 0;

if(nargin < 2 | nargin >3)
   disp('2 or 3 input arguments are required')
   return
end

if(ischar(filename) == 0)
   disp('input arguments should be character arrays')
   return
end

if(ischar(varname) == 0)
   disp('input arguments should be character arrays')
   return
end

if(length(filename) == 0)
   [fname, pathname] = uigetfile('*.nex', 'Select a Nex file');
	filename = strcat(pathname, fname);
end
if ~exist('verbose','var') 
	verbose = 1;
end

fid = fopen(filename, 'r');
if(fid == 0)
   return
end

magic = fread(fid, 1, 'int32');
version = fread(fid, 1, 'int32');
comment = fread(fid, 256, 'char');
freq = fread(fid, 1, 'double');
tbeg = fread(fid, 1, 'int32');
tend = fread(fid, 1, 'int32');
nvar = fread(fid, 1, 'int32');
fseek(fid, 260, 'cof');
name = zeros(1, 64);
found = 0;
for i=1:nvar
	type = fread(fid, 1, 'int32');
	var_version = fread(fid, 1, 'int32');
	name = fread(fid, [1 64], 'char');
	offset = fread(fid, 1, 'int32');
	nf = fread(fid, 1, 'int32');
	dummy = fread(fid, 32, 'char');
	adfreq = fread(fid, 1, 'double');
	adtomv = fread(fid, 1, 'double');
	n = fread(fid, 1, 'int32');
	name = setstr(name);
	name = deblank(name);
	k = strcmp(name, deblank(varname));
	if(k == 1)
		if type ~= 3
			disp(sprintf('%s is not a waveform variable', deblank(varname)));
			return;
		end
		found = 1;
		fseek(fid, offset, 'bof');
		ts = fread(fid, [1 nf], 'int32');
		w = fread(fid, [n nf], 'int16');
		break
	end
	dummy = fread(fid, 76, 'char');
end

fclose(fid);

if found == 0
	disp('did not find variable in the file');
else
	ts = ts/freq;
	w = w*adtomv;
	if verbose
		disp(strcat('file = ', filename));
		disp(strcat('number of waveforms = ', num2str(nf)));
	end
end
