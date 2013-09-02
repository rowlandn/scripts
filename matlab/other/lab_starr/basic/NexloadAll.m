function Nexdata = NexLoadAll(filename,verbose)
% function Nexdata = NexLoadAll(nexfilename,verbose)
% 
%	Function to load all variable in a nex file into Nexdata structure
%	Fields in Nexdata depend on contents of nex file
%
%	inputs:	nexfilename - file to be opened	
%			verbose - controls feedback on progress
% 

% Variable types
NEURON = 0;
EVENT = 1;
INTERVAL = 2;
WAVEFORM = 3;
POPVECT = 4;
CONTIN = 5;
MARKER = 6;


if ~exist('verbose','var') 
    verbose=false;
end

if ~exist(filename,'file')
	error(['Unable to open the file: ' filename]);
end

% USE NEX AUTHOR'S MATLAB FUNCTION TO GRAB: #vars, varnames, vartypes:
[nvar varnames vartypes] = nex_info(filename,verbose);
Nexdata.nvar = nvar;
Nexdata.varnames = varnames;
Nexdata.vartypes = vartypes;

for v=1:nvar
    
	varname = deblank( varnames(v,:) );
    switch vartypes(v)	% types in file.

		case {NEURON, EVENT}
			[n, ts] = nex_ts(filename, varname, verbose);
			eval(['Nexdata.' varname '.n = n;']);
			eval(['Nexdata.' varname '.ts = ts;']);
            
        case INTERVAL
			[n, ts_left, ts_right] = nex_int(filename, varname, verbose);
			eval(['Nexdata.' varname '.n = n;']);
			eval(['Nexdata.' varname '.ts = ts_left;']);
			eval(['Nexdata.' varname '.ts = ts_right;']);
            
        case WAVEFORM
			[adfreq, nf, ts, n, w] = nex_wf(filename, varname, verbose);
			eval(['Nexdata.snip_fs = adfreq;']);
			eval(['Nexdata.' varname '.ts = ts;']);
			eval(['Nexdata.' varname '.snips = w;']);
            
        case POPVECT
            % SKIP THIS FOR NOW
            
        case CONTIN 
			[adfreq, n, ts, fn, d] = nex_cont(filename, varname, verbose);
			eval(['Nexdata.cont_fs = adfreq;']);
			eval(['Nexdata.' varname '.len = n;']);
			eval(['Nexdata.' varname '.data = d;']);
            
        case MARKER
			[n, nm, nl, ts, names, m] = nex_marker(filename, varname, verbose);
			eval(['Nexdata.' varname '.n = n;']);
			eval(['Nexdata.' varname '.ts = ts;']);
			eval(['Nexdata.' varname '.names = names;']);
			eval(['Nexdata.' varname '.marker = m;']);
			
    end
end
if verbose		disp('NEX FILE SUCCESSFULLY LOADED'); end	
