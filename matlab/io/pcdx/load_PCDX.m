function [traces, ntraces] = ...
    load_PCDX_SE(file, tracelist, channel)

% load_PCDX_SE This function loads data acquired
% in the PCDX format.  Multiple traces from a single
% channel may be loaded. 
%
% traces = load_PCDX_SE(file, tracelist, channel)
% 
% Example 1: A = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',1);
% 
% NOTE: You will need the readpcdx plug-in and the 
%       tracelist2 script (@Confucius:/Lab/matlab/alfonso)
%       
% <adelgado@biology.emory.edu>


	
	[tracenums, ntraces] = gettracelist2(tracelist);
    
	progbar = waitbar(0, 'Loading PCDX Traces...');
    
	for i = 1: ntraces
        trace = readpcdx(file, tracenums(i), channel);
        traces(1: length(trace), i) = trace;
        zoom on      
        waitbar(i/ntraces, progbar)
	end
    
    
    close(progbar)
    