function out = outrange(chan, range)
% function out = outrange(chan, range)
%
%	Return logic if chan is outside of range
%
%	RST 2006-04-05

out = chan < range(1) | chan > range(2);
