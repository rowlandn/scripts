function in = inrange(chan, range)
% function in = inrange(chan, range)
%
%	Return logic if chan is within range
%
%	RST 2006-04-05

in = chan > range(1) & chan < range(2);
