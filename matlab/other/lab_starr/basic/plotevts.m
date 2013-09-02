function h = plotevts(x, y, ticktop,color,txt)
%PLOTBAR Plot two vectors simulating a '|' line
%     h = = plotevts(x, y, ticktop,color,txt) 
%	Returns a vector of handles to the lines created.
%	Arguements:
%		x = row vector of spike times
%		y = row vector of y position of events
%		ticktop = y position of top of bars
%		color (optional specification of bar color)
%		txt = optional text to write at top of each bar
%
%	RST 2006-04-05

if ~exist('color')
	color = 'b';
end

if isempty(x)
	display('No events to plot');
	return
end

if size(x,1) > 1 %if a column vector
   x = x';
end

if size(y,1) > 1
   y = y';
end

ln_x = length(x);
newx = zeros(2,ln_x);
newy = newx;

newx(1,:) = x;
newx(2,:) = x;
newy(1,:) = y;
newy(2,:) = ticktop;

h = plot(newx,newy, 'LineWidth', 0.01,'color', color);

if exist('txt')
	hold on
	txty = newy(2,1);
	for i=1:ln_x
		text( x(1,i),txty,txt,'Color',color,...
			'HorizontalAlignment','center','VerticalAlignment','top',...
			'BackgroundColor','w','FontSize',8)
	end
end