function resize_all_x(handle, xmin, xmax, xincr)
%RESIZE_ALL Resize x-axes of all panes in graph window to user-specified values.
% Input arguments: figure handle, xmin, xmax, xincr.
%
%resize_all_x(handle, xmin, xmax, xincr)

if nargin < 4
	error('Not enough arguments.')
elseif nargin > 4
	error('Too many arguments.')
end

tvec = get(handle, 'Children');
num = size(tvec);	% number of child axes in figure.

rmd = mod(xmin,xincr);
if rmd == 0
	xticks = [xmin:xincr:xmax];
else
	xticks = [xmin+rmd:xincr:xmax];	% make ticks at multiples of increment.
end

for n = 1:num
	set(tvec(n),'XLim', [xmin,xmax]);
	set(tvec(n),'XTick', xticks);
end

clear tvec;
clear num;
clear rmd;
clear xticks;

return

	

