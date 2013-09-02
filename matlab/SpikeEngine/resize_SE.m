function resize(handle, pane,  xmin, xmax, xincr, ymin, ymax, yincr)
%RESIZE Resize axes of specified pane in graph window to user-specified values.
% Input arguments: figure handle, pane, xmin, xmax, xincr, ymin, ymax, yincr.
%
%resize(handle, pane,  xmin, xmax, xincr, ymin, ymax, yincr)

if nargin < 6
	error('Not enough arguments.')
elseif nargin == 7
	error('Inappropriate number of arguments.')
elseif nargin > 8
	error('Too many arguments.')
end

if (nargin == 6)
	ymax = ymin
	ymin = xincr
	xincr = (xmax-xmin)/2
	yincr = (ymax-ymin)/5
end

tvec = get(handle, 'Children');
num = size(tvec, 1);	% number of child axes in figure.

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

n = size(tvec, 1)
set(tvec(1+num-pane),'XLim', [xmin,xmax]);
set(tvec(1+num-pane),'YLim', [ymin,ymax]);
set(tvec(1+num-pane),'XTick', xticks);
set(tvec(1+num-pane),'YTick', yticks);

clear tvec;
clear num;
clear rmd;
clear xticks;
clear yticks;

return

	

