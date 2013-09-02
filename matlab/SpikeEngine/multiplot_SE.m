function Hfig = multiplot(fignum, nrows, ncols, xdset, xrange, ydset, varargin)
%MULTIPLOT plots multiple data columns on a single graph, with each column occupying
%	a separate axis pane.
%
% Usage:	Hfig = multiplot(fignum, ncols, nrows, xdset, xrange, ydset, varargin);
%
% Return value is the handle for the figure.
%
% If ncols or nrows == 0, value will be automatically determined.
%
% varargin = column numbers within the matrix dset.
%
% User can specify individual columns or ranges using the : operator.
%
% Example:
%	fig1 = multiplot(1, 0, 0, time, 1:10000, datatraces, 1, 3, 7:12, 25:30);
%
%	In this example, a figure with the index number 1 will be created, its
%		handle will be stored in the variable fig1, 16 subplots will be
%		put in the figure with one trace in each subplot. X-axis data are
%		from the vector dataset 'time', while y-axis data are from the matrix
%		dataset 'datatraces.' The particular columns of 'datatraces' included
%		in the plot are 1, 3, 7-12, and 25-30, and the first 10,000 rows of
%		each column will be plotted.
%
% Created by: Jeremy R. Edgerton, November 2003.
%
% SEE ALSO: resize, resize_all

nplots = 0;		% Total number of traces to be plotted.
narg = nargin - 6;	% Number of trace index args. Keep in mind
			%	that each arg could be either a single
			%	trace or a whole range of traces.
counter = 0;		% Keep track of how many traces have been plotted.
tmpcol = ncols;
tmprow = nrows;
tmp1 = 0;
tmp2 = 0;

for n = 1:narg
	[tmp1,tmp2] = size(varargin{n});
	nplots = nplots + (tmp1*tmp2);
end

if nplots < 1
	error('No traces entered, or insufficient args.')
end

if tmprow == 0	
	% If automatically setting up panes, keep grid as square as possible,
	%	adding first a row, then a column as necessary.
	if tmpcol == 0
		tmpcol = round(sqrt(nplots));	%number of columns on graph
		if tmpcol >= sqrt(nplots)
			tmprow = tmpcol;
		else
			tmprow = tmpcol + 1;
		end
	else
		if round(nplots/tmpcol) > nplots/tmpcol
			tmprow = round(nplots/tmpcol)
		else
			tmprow = round(nplots/tmpcol) + 1
		end
	end

elseif tmpcol == 0
	if round(nplots/tmprow) > nplots/tmprow
		tmpcol = round(nplots/tmprow)
	else
		tmpcol = round(nplots/tmprow) + 1
	end
	
else
	%user has supplied both nrows and ncols, so no calculation needed.
end


Hfig = figure(fignum);
for n = 1:narg
	[tmp1, tmp2] = size(varargin{n});
	for k = 1:tmp1*tmp2
		counter = counter + 1;
		subplot(tmprow, tmpcol, counter);
		plot(xdset(xrange), ydset(xrange, varargin{n}(k)));
	end
end

% echo on
% narg
% nrow
% ncol
% nplots
% counter
% echo off

clear nplots;
clear narg;
clear counter;
clear tmpcol
clear tmprow;

return

