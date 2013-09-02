function inds = getinds(chan, range)
% function inds = getinds(chan, range)
%
%	Helper to find indices of a chan withing range
%
% RST 2006-04-05

inds = find(chan>range(1) & chan<range(2));

return