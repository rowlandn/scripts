function psth = psth_rast_slow_SE(spike_times,trial_length,binwidth,trigger,sp1,sp2,sp3)

% This function plots a raster plot and corresponding peri-stimulus time
% histogram of given spike times.  Spike times should be a cell consisting
% of spike time arrays with the row of the cell array representing the
% corresponding trace from which the spikes were discriminated (e.g.,
% output of findspikes_win_SE function). Enter trial length in ms.  
% The histogram is binned with the user-defined binwidth.  A vertical line 
% is placed at the trigger (ms) to indicate the time of the stimulus. This
% algorith uses a simple for loop to plot the raster plot and can sometimes
% take a while to finish.  Use this script over psth_rast_fast_SE when
% you have a slow neuron and there are likely trials with no spikes in
% them.
% 
% psth_rast_SE(spike_times,trial_length,binwidth,trigger,sp1,sp2,sp3)

SPIKES = cat(2,spike_times{:});
%assignin('base','SPIKES',SPIKES)
% assignin('base','spike_times',spike_times)

if nargin == 7
    subplot(sp1,sp2,sp3)
end
no_bins = round(trial_length/binwidth);
no_trials = size(spike_times,1);

[hist_SPIKES_y,hist_SPIKES_x] = hist(SPIKES,no_bins);
% assignin('base','hist_SPIKES_y',hist_SPIKES_y)
% assignin('base','hist_SPIKES_x',hist_SPIKES_x)
hist_SPIKES_normalized_freq = hist_SPIKES_y/(no_trials * (binwidth/1000));
bar(hist_SPIKES_x,hist_SPIKES_normalized_freq,'k'); hold on
psth.x = hist_SPIKES_x;
psth.y = hist_SPIKES_normalized_freq;


max_hist_SPIKES_normalized_freq = max(hist_SPIKES_normalized_freq);
% assignin('base','max_hist_SPIKES_normalized_freq',max_hist_SPIKES_normalized_freq)
axis([0 hist_SPIKES_x(end) 0 max_hist_SPIKES_normalized_freq+.1*max_hist_SPIKES_normalized_freq+no_trials]);

for i = 1:size(spike_times,1)
    if size(spike_times{i,1},2) == 0
    else
        plot(spike_times{i,1},max_hist_SPIKES_normalized_freq+.1*max_hist_SPIKES_normalized_freq+i,'k.','MarkerSize',1); hold on
    end
end
% plot([trigger trigger],[0 max_hist_SPIKES_normalized_freq+.1*max_hist_SPIKES_normalized_freq+no_trials],'k')

% SPIKES = SPIKES';
% diff_SPIKES = diff(SPIKES(:,1));
% trial_breaks = find(diff_SPIKES < 0);
% 
% %SPIKES = SPIKES(:,1)/binwidth;
% 
% for i = 1:no_trials
%     if i == 1
%         trial_1_x = SPIKES((1):(trial_breaks(1)));
%         trial_1_y = linspace(max_hist_SPIKES_normalized_freq+no_trials*1.9-(1*1.75-1),...
%             max_hist_SPIKES_normalized_freq+no_trials*1.9-(1*1.75-1),length(trial_1_x));
%         plot(trial_1_x,trial_1_y,'ko','MarkerSize',1); hold on
%     elseif i == no_trials
%         trial_no_x_str = ['trial_',num2str(i),'_x = SPIKES((trial_breaks(',...
%                 num2str(i),'-1) + 1):end);'];
%         eval(trial_no_x_str)
%         trial_no_y_str = ['trial_',num2str(i),'_y = linspace(max_hist_SPIKES_normalized_freq+no_trials*1.9-(',...
%                 num2str(i),'*1.75-1),max_hist_SPIKES_normalized_freq+no_trials*1.9-(',num2str(i),'*1.75-1),length(trial_',...
%                 num2str(i),'_x));'];
%         eval(trial_no_y_str)
%         plot_trial_no_x_str = ['plot(trial_',num2str(i),'_x,trial_',...
%                 num2str(i),'_y,''ko'',''MarkerSize'',1);'];
%         eval(plot_trial_no_x_str)
%     else
%         assignin('base','i',i)
%         assignin('base','SPIKES',SPIKES)
%         assignin('base','trial_breaks',trial_breaks)
%         assignin('base','spike_times',spike_times)
%         trial_no_x_str = ['trial_',num2str(i),'_x = SPIKES((trial_breaks(',...
%                 num2str(i),'-1) + 1):(trial_breaks(',num2str(i),')));'];
%         eval(trial_no_x_str)
%         trial_no_y_str = ['trial_',num2str(i),'_y = linspace(max_hist_SPIKES_normalized_freq+no_trials*1.9-(',...
%                 num2str(i),'*1.75-1),max_hist_SPIKES_normalized_freq+no_trials*1.9-(',num2str(i),'*1.75-1),length(trial_',...
%                 num2str(i),'_x));'];
%         eval(trial_no_y_str)
%         plot_trial_no_x_str = ['plot(trial_',num2str(i),'_x,trial_',...
%                 num2str(i),'_y,''ko'',''MarkerSize'',1);'];
%         eval(plot_trial_no_x_str)
%     end
% end
% 
% xlabel('Time (ms)')
% ylabel('spikes/s                                               Trial No.')
% y_lim = get(gca,'ylim');
% % xtick = get(gca,'XTick');
% % xtick_relabel = xtick * binwidth;
% % set(gca,'XTickLabel',xtick_relabel)
% % 
% % ytick = get(gca,'YTick');
% % closest_ytick = min(find(ytick > max_hist_SPIKES_normalized_freq));
% % set(gca,'YTick',ytick(1:closest_ytick))
% 
% %trigger_x = linspace(200,200,max_hist_SPIKES_normalized_freq+(no_trials*2));
% % trigger_x = linspace(trigger/binwidth,trigger/binwidth,max_hist_SPIKES_normalized_freq+(no_trials*2));
% % assignin('base','trigger_x',trigger_x)
% plot([trigger trigger],[0 y_lim(2)],'k')
% psth.trigger = trigger; 