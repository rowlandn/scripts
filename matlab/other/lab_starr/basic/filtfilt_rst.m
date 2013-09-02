function y = filtfilt_rst(b, a, x)
% y = filtfilt_rst(b, a, x)
% This wrapper for filtfilt pads the x array 
% before passing it on to filtfilt
% Created by RST, 2004-02-03
%	

	dl = length(x);		% Data length
	fl = length(b);		% Filter length
	
	start_ind = fl:-1:1;
	stop_ind = dl:-1:(dl-fl+1);
	
	% Pad by reflecting
	pad_ra = [ x(start_ind) x x(stop_ind) ];
	y = filtfilt(b, a, pad_ra ); 

	% Remove padding
	y( (end-fl+1):end )=[];
	y( 1:fl )=[];
