function resize_all(handle, xmin, xmax, xincr, ymin, ymax, yincr)
%RESIZE_ALL Resize axes of all panes in graph window to user-specified values.
% Input arguments: figure handle, xmin, xmax, xincr, ymin, ymax, yincr.
%

if nargin < 5
	error('Not enough arguments.')
elseif nargin > 7
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

rmd = mod(ymin,yincr);
if rmd == 0
	yticks = [ymin:yincr:ymax];
else
	yticks = [ymin+rmd:yincr:ymax];
end

for n = 1:num
	set(tvec(n),'XLim', [xmin,xmax]);
	set(tvec(n),'YLim', [ymin,ymax]);
	set(tvec(n),'XTick', xticks);
	set(tvec(n),'YTick', yticks);
end

clear tvec;
clear num;
clear rmd;
clear xticks;
clear yticks;

return

	

