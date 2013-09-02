function rasterplot(spike_mat,tickheight,txt,xlm,color)
% Create raster plot
%     rasterplot(spike_mat,tickheight,txt)
%	Arguements:
%		spike_mat = 2-dimensional array of spike times
%			1 row per trial
%			NaN pads variable length rows
%		tickheight = height of '|' markers as fraction of 1
%		txt = text to print as title of plot
%		xlm - optional x-limits, used if defined
%		color - optional color (default = black)

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

missing = find( isnan(spike_mat(:,1)) );
spike_mat(missing,:) = [];

[ntrials,maxspks] = size(spike_mat);

trial_n = zeros(ntrials,maxspks);

for i = 1:ntrials
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
ylim([tickheight ntrials+tickheight]); 
x = get(gca);
ind = find( rem(x.YTick,1)==0 );
ytk = x.YTick(ind);
set(gca,'YTick',ytk);

xlim( xlm );
title(txt);

