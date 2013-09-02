function h = plotbar_color(x, y, tickheight,tick_width,color)
%PLOTBAR Plot two vectors simulating a '|' linestyle
%     h = plotbar(x,y, tickheight,tick_width,color) 
% Returns a handle to the line created.
% This function was developed by Scott in Mike Crutcher's lab, and
% was only minimally altered here by TW & RST

x = x(:);
y = y(:);

newx = zeros(1,length(x)*3);
newx(1:3:end) = x;
newx(2:3:end) = x;
newx(3:3:end) = nan;
newy = zeros(1,length(y)*3);
newy(1:3:end) = y;
newy(2:3:end) = y + tickheight;
newy(3:3:end) = nan;
h = plot(newx,newy, color,'LineWidth', tick_width);