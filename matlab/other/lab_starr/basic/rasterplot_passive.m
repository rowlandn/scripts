function rasterplot_passive(s, RSTROWLEN, tickheight,txt,xlm,color, marker)
% Create raster plot
%     
%	Arguements:
%		tickheight = height of '|' markers as fraction of 1
%		txt = text to print as title of plot
%		xlm - optional x-limits, used if defined
%		color - optional color (default = black)

spike_mat = s.raster;

if ~exist('xlm','var')
	xlm = [min(min(spike_mat)) max(max(spike_mat))];
end
if ~exist('color','var')
    f = get(gcf);
    if sum(f.Color) == 0
        color = 'w';
    else
        color = 'k';
    end
end

if ~exist('marker','var')
    marker.dir(1) = 'o';
    marker.dir(2) = '+';
end


missing = find( isnan(spike_mat(:,1)) );
spike_mat(missing,:) = [];

[nrows,maxspks] = size(spike_mat);

trial_n = zeros(nrows,maxspks);

for i = 1:nrows
	trial_n(i,:) = i;
end

% transform spike times & trial numbers into row vectors
x = reshape(spike_mat',1,[]);
y = reshape(trial_n',1,[]);

newx = zeros(2,length(x));
newx(1,:) = x;
newx(2,:) = x;
newy = zeros(2,length(y));
newy(1,:) = y;
newy(2,:) = y + tickheight;

h = plot(newx,newy, 'LineWidth', 0.01,'color', color);

for j = 1:length(s.dir)
    hold on
    mvt = s.dir(j).ts;
    accel_y = floor(mvt/RSTROWLEN)+1;
    accel_x = mvt - RSTROWLEN * floor(mvt/RSTROWLEN);
    accel_newx = zeros(2,length(mvt));
    accel_newx(1,:) = accel_x;
    accel_newx(2,:) = accel_x;
    accel_newy = zeros(2,length(mvt));
    accel_newy(1,:) = accel_y;
    accel_newy(2,:) = accel_y + tickheight;
    h = plot(accel_newx,accel_newy, 'LineWidth', 0.01,'color','k','marker', marker.dir(j));

end

hold off

ylim([tickheight nrows+tickheight]); 
x = get(gca);
ind = find( rem(x.YTick,1)==0 );
ytk = x.YTick(ind);
set(gca,'YTick',ytk);

xlim( xlm );
title(txt);
