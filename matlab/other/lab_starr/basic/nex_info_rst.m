function info = nex_info_rst(filename, verbose)
% info = nex_info_rst(filename, verbose) -- read and display .nex file info
%   Enhanced by RST for extra info
%
% info = nex_info+(filename, verbose)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%	verbose - optional flag to control display of results (default=1)
% OUTPUT:
%   info.nvar - number of variables in the file
%   info.names - [nvar 64] array of variable names
%   info.types - [1 nvar] array of variable types
%           Interpretation of type values: 0-neuron, 1-event, 2-interval, 3-waveform, 
%                        4-population vector, 5-continuous variable, 6 - marker
%   info.ver - version
%   info.fs - sampling fequency
%   info.tbeg - data begin sample
%   info.tend - data end sample
%   info.dur - data duration (in sec)

if ~exist('verbose','var') 
	verbose = 1;
end

if(nargin<1 | nargin>2 )
   disp('1 or 2 inputs arguments are required')
   return
end

if(length(filename) == 0)
   [fname, pathname] = uigetfile('*.nex', 'Select a Nex file');
	filename = strcat(pathname, fname);
end

fid = fopen(filename, 'r');
if(fid == -1)
	disp('cannot open file');
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
if verbose
	disp(strcat('file = ', filename));
	disp(strcat('version = ', num2str(version)));
	disp(strcat('frequency = ', num2str(freq)));
	disp(strcat('duration (sec) = ', num2str((tend - tbeg)/freq)));
	disp(strcat('number of variables = ', num2str(nvar)));
end
names = zeros(1, 64);
for i=1:nvar
	types(i) = fread(fid, 1, 'int32');
	var_version = fread(fid, 1, 'int32');
	names(i, :) = fread(fid, [1 64], 'char');
	dummy = fread(fid, 128+8, 'char');
end
names = setstr(names);
fclose(fid);

   info.nvar = nvar;    % number of variables in the file
   info.names = names;  %- [nvar 64] array of variable names
   info.types = types;  %- [1 nvar] array of variable types
   info.ver = version;  %- version
   info.fs = freq;      %- sampling fequency
   info.sbeg = tbeg;    %- begin sample
   info.send = tend;    %- end sample
   info.dur = (tend - tbeg)/freq;   % data duration

