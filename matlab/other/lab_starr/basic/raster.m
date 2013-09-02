function y = raster(A,color,width,tickheight,start,length_of_raster,num_spaces)
% Function raster(A,width,tickheight,start,length_of_raster,num_spaces) plots a raster diagram of A, using width
% as the x-axis and tickheight as tick height.  A has to be an ISI stream (such
% as is contained in the *.BIS files), and width is measured in the time unit of
% the input stream of A.  Usually this will be in ms.  The color parameter specifies the 
% color of the raster plot.
%
% The parameters start and length_of raster define the beginning and the length
% of the plot.  To plot the entire data stream, set start = 0, and length_of raster = -1
%
% Num_spaces sets the spacing between the lines of the raster plot.
%
% Not all input parameters have to be given.  Defaults are color = black, width = 1000, tickheight = 1,
% start = 0, length_of_raster = total length of data stream.
%
% Sample call
% 		raster(B,1000,1,1000,5500)
%
% This will draw a rasterplot of B with tickheight 1.  Every 1000 ms, a new line
% is started. The plot will cover 5500 ms, starting at 1000 ms.
%
% Written 12/18/99 by Thomas Wichmann
% Modified 12/21/99, 1/10/2000

A = cumsum(A);

% Taking care of command line length
if nargin < 2,color = 'k';end
if nargin < 3,width = 1000;end
if nargin < 4,tickheight = 1;end
if nargin < 5,start = 0;end
if nargin < 6,length_of_raster = A(length(A));end
if nargin < 7,num_spaces = 1-tickheight;end

% Adjusting length of raster, if input argument was -1
if length_of_raster == -1,length_of_raster = A(length(A));end

% Calculation of start of displayed data strem
begin = threshold(A,start,'-+');
if begin == -1,begin = 1;end;

% Calculation of end of displayed data stream
finish = threshold(A,start+length_of_raster,'-+');
if finish == -1,finish = length(A)+1;end;

% Display of data stream
D = A(begin:(finish-1));								% Selection of data segment to be displayed
D = D-start;													% Adjustment so that the data segment starts with 0
B = ceil(D/width)*(num_spaces+tickheight);	% Calculation of Y-axis of plot, adjusting for line spacing
C = rem(D,width);											% Calculation of X-axis of plot
plotbar_color(C,B,tickheight,color);								% Plotting the raster plot - this is a modified version of Scott's function
axis ij
axis tight														% Arranging the axes to be tight so that plot fills the figure
