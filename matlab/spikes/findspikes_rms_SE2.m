function spikes_SE = findspikes_rms_SE2(traces,spike_times,rmst,art_plottype)

%
% findspikes_rms_SE 
% Syntax: A = findspikes_rms_SE(traces,spike_times,pre,post,rms_thresh)
% Accepted spikes are plotted in green and rejected spikes in red 

pre = 10; %samples
post = 10; %samples

collect_spikes = collect_spikes_SE(traces,spike_times.times.win,pre,post,0);
%assignin('base','collect_spikes',collect_spikes)

for i = 1:size(traces,2)
clear rejected_spike_count rejected_spike_index rms_all
rejected_spike_count = 0;
rejected_spike_index = [];
    for j = 1:size(collect_spikes.spk_traces{i,1},2)
        err1 = collect_spikes.spk_avg{i,1}(1:pre)-collect_spikes.spk_traces{i,1}(1:pre,j);
        %err2 = collect_spikes.spk_avg{i,1}(25:40)-collect_spikes.spk_traces{i,1}(25:40,j);
        err2 = collect_spikes.spk_avg{i,1}(pre+1:pre+post)-collect_spikes.spk_traces{i,1}(pre+1:pre+post,j);
        err = [err1',err2']';
        sqerr = err .* err;
        rms = sqrt(mean(sqerr));
        rms_all(j,1) = rms;
        %assignin('base','rms_all',rms_all)
        %# The preceeding 5 lines calculate the RMS difference between
        %# the trace under test and the stored average trace, calculated
        %# over two regions:  First, from the beginning (baseline) to the
        %# peak (at 0.2 ms; 20 points) and second, during the long-term
        %# return to baseline (from 0.45 to 0.64 ms; points 45-64)
        %# The following "if" loop will dump traces whose RMS deviation
        %# from the average trace is above some set threshold value.
    %     Curve_max = abs(max(diff(diff(collect_spikes.spk_traces(1:(pre+1),j)))));
    %         %# This calculates the maximum of the second derivative
    %         %# (curvature) of the rising part of the trace.  The "tight"
    %         %# peaks, which rise abruptly, will show up with large values of
    %         %# this curvature at the "kink"
        
        if (rms > rmst) %| Curve_max > 1800)   %# 50 for less noisy traces
        %#################################################################
        %# End of "if" test for bad traces.  Following is block to "dump"
        %# rejected traces.  Indices of bad traces are saved for later purging
        %####################################################################
            rejected_spike_count = rejected_spike_count + 1;
            rejected_spike_index = [rejected_spike_index j];
        end
    end
    spikes_SE.rejected_spike_counts(i,1) = rejected_spike_count;
    spikes_SE.rejected_spike_index{i,1} = rejected_spike_index;
    spikes_SE.rms_values{i,1} = rms_all;
    %assignin('base','a',spikes_SE)
end

spike_count_00 = cellfun(@length,spike_times.times.win);
spike_count_01 = sum(spike_count_00);

rev_spike_times = spike_times.times.win;
rev_spike_peaks = spike_times.peaks.win;
for i = 1:size(spikes_SE.rejected_spike_index,1)
    rev_spike_times{i,1}(spikes_SE.rejected_spike_index{i,1}) = [];
    rev_spike_peaks{i,1}(spikes_SE.rejected_spike_index{i,1}) = [];
end
 
spikes_SE.spike_count_01 = spike_count_01;
spikes_SE.times_win_rms = rev_spike_times;
spikes_SE.peaks_win_rms = rev_spike_peaks;
spikes_SE.collect_spikes = collect_spikes;
% spikes_SE{:,1}(rejected_spike_index) = [];



    

%##############################################################
%# At this point, all the candidate traces have been examined, and
%# indices of the "bad" traces have been stored in "dump_handles"
%# Now prepare to re-run the discrimination analysis based on the 
%# first pass through.
%###############################################################
% 
% Goodspikes1 = X.spikes;
% 
% if ~isempty(dump_handles)
% Goodspikes1(:,dump_handles)= [];
% end
% 
% if ~isempty(Goodspikes1)
% Good_spksavg1 = (mean(Goodspikes1'))';
% end
% 
% %setup_figure;
% 
% subplot(3,1,1)
% hold on;
% grid on;
% plot(peaks_up(1,:),peaks_up(2,:),'ro')
% for j = 1:length(peaks_up)
% plot([(peaks_up(1,j)-pre*dt):dt:(peaks_up(1,j)+post*dt)], ...
%     X.spikes(1:41,j),'color',[0.8 0.8 0.8])
% end
% % ylim([-300 300]);
% % xlim([xmin xmax]);
% ylabel('Potential, [mV]')
% xlabel('Time, [sec]')
% xl = get(gca, 'XLim');
% x2_3 = xl(1)+0.67*(xl(2)-xl(1));
% %# The above two lines find the width of the total data graph, and
% %# calculate the time that is 2/3 of the way across, for display of the 
% %# "Total Time" data.
% % % text(x2_3, 280,['Ending Time: ',int2str(xmax),' sec.'],...
% % %     'FontSize',10,'Color','blue')
% % % 
% % title([filename ' Traces ' int2str(starting_trace_number) '-' ...
% %         int2str(ending_trace_number) ' Detected and Rejected Spike Times'])
% 
% subplot(3,1,2)
% hold on;
% grid on;
% 
% plot(spike_x,X.spikes,'k')
% plot(spike_x,X.spksavg,'r','LineWidth', 2)
% if ~isempty(Goodspikes1)
%     plot(spike_x,Good_spksavg1,'b','LineWidth',2)      
% end
% %ylim([-300 300]);
% ylabel('Potential, [mV]')
% xlabel('Time, [msec]')
% 
% % title([filename ' Traces ' int2str(starting_trace_number) '-'...
% %     int2str(ending_trace_number) ' Average and Detected Spike ' ...
% %     'Traces (3 ^r^d pass fit to average disc.)'])
% 
% 
% 
% 
% %# Now re-set the arrays used in the discrimination routine:
% peaks_up2 = [];
% up_dumped = 0;
% dump_handles = [];
% 
% for j = 1:length(peaks_up)
%     k = j-up_dumped;
%     subplot(3,1,2)
%     plot(spike_x,X.spikes(:,j),'g')
%     
%     %#err1 = Good_spksavg1(1:20)-X.spikes(1:20,j);
%     %#err2 = Good_spksavg1(45:64)-X.spikes(45:64,j);
%     %#err = [err1',err2']';
%     if findstr(fitwidth,'narrow')
%         err = Good_spksavg1(1:25)-X.spikes(1:25,j);
%     else
%         %# err = Good_spksavg1(1:65)-X.spikes(1:65,j);
%         err = Good_spksavg1(1:40)-X.spikes(1:40,j);
%     end
%     sqerr = err .* err;
%     rms = sqrt(mean(sqerr));
%     %# The preceeding 3 lines calculate the RMS difference between
%     %# the trace under test and the average trace calculated above,
%     %# calculated over a broad region that should include the whole 
%     %# Action Potential: (from 0 to 0.65 ms; points 0-65)
%     %# The following "if" loop will dump traces whose RMS deviation 
%     %# from the average trace is above some set threshold value.
%     
%     Curve_max = max(diff(diff(X.spikes(1:(pre+1),j))));
%         %# This calculates the maximum of the second derivative 
%         %# (curvature) of the rising part of the trace.  The "tight"
%         %# peaks, which rise abruptly, will show up with large values of
%         %# this curvature at the "kink"
%     %#if rms > 80    %# 90 for moderately noisy traces  
%     if (rms > rms2 | Curve_max > 1800) 
%         %# Here the RMS value can be much larger than that in the original
%         %# pass through, since we can be more confident about the shape of
%         %# the average trace.
%      %#################################################################
%      %# End of "if" test for bad traces.  Following is block to "dump"
%      %# rejected traces.  Indices of bad traces are saved for later purging
%      %####################################################################
%         up_dumped = up_dumped +1;
%         dump_handles = [dump_handles j];
%         
%         subplot(3,1,1)
%         plot(peaks_up(1,j),peaks_up(2,j),'ko')
%         plot([(peaks_up(1,j)-pre*dt):dt:(peaks_up(1,j)+post*dt)], ...
%                 X.spikes(1:(pre+post+1),j),'m')
%         
%         subplot(3,1,2)
%         plot(spike_x,X.spikes(:,j),'k')
%                 
%     else
%         %#If we have passed ALL the tests:
%         peaks_up2(:,k)= peaks_up(:,j);
%         %# Store the next good value in the output array.
%         subplot(3,1,1)
%         plot(peaks_up2(1,k),peaks_up2(2,k),'go')     
%     end  %# While Loop
%     %# Now re-plot the trace being examined in black (cover the green)
%     subplot(3,1,2)
%     plot(spike_x,X.spikes(:,j),'k')
% end
% 
% disp([int2str(up_dumped) ' traces rejected (second pass)'])
% 
% Goodspikes2 = X.spikes;
% if ~isempty(dump_handles)
% Goodspikes2(:,dump_handles)= [];
% end
% 
% if ~isempty(Goodspikes2)
% Good_spksavg2 = (mean(Goodspikes2'))';
% %#Spike_Times = peaks_up2(1,:);   
% else
%     Spike_Times = [];
% end
% % %########################################################################
% % %#  Run through a third time......
% % %########################################################################
% 
% 
% %# Now re-set the arrays used in the discrimination routine:
% peaks_up3 = [];
% up_dumped = 0;
% dump_handles = [];
% subplot(3,1,1)
% cla;
% hold on;
% grid on;
% plot(peaks_up(1,:),peaks_up(2,:),'ro')
% for j = 1:length(peaks_up)
% plot([(peaks_up(1,j)-pre*dt):dt:(peaks_up(1,j)+post*dt)], ...
%     X.spikes(1:41,j),'color',[0.8 0.8 0.8])
% end
% 
% for j = 1:length(peaks_up)
%     k = j-up_dumped;
%     subplot(3,1,2)
%     plot(spike_x,X.spikes(:,j),'g')
%     
%     %#err1 = Good_spksavg1(1:20)-X.spikes(1:20,j);
%     %#err2 = Good_spksavg1(45:64)-X.spikes(45:64,j);
%     %#err = [err1',err2']';
%     if findstr(fitwidth,'narrow')
%         err = Good_spksavg2(1:25)-X.spikes(1:25,j);
%     else
%         %err = Good_spksavg2(1:65)-X.spikes(1:65,j);
%     end
%     sqerr = err .* err;
%     rms = sqrt(mean(sqerr));
%     %# The preceeding 3 lines calculate the RMS difference between
%     %# the trace under test and the average trace calculated above,
%     %# calculated over a broad region that should include the whole 
%     %# Action Potential: (from 0 to 0.65 ms; points 0-65)
%     %# The following "if" loop will dump traces whose RMS deviation 
%     %# from the average trace is above some set threshold value.
%     
%     Curve_max = max(diff(diff(X.spikes(1:21,j))));
%         %# This calculates the maximum of the second derivative 
%         %# (curvature) of the rising part of the trace.  The "tight"
%         %# peaks, which rise abruptly, will show up with large values of
%         %# this curvature at the "kink"
%     %#if rms > 80    %# 90 for moderately noisy traces  
%     if (rms > (0.75*rms2) | Curve_max > 1800) 
%         %# Here the RMS value can be much larger than that in the original
%         %# pass through, since we can be more confident about the shape of
%         %# the average trace.
%      %#################################################################
%      %# End of "if" test for bad traces.  Following is block to "dump"
%      %# rejected traces.  Indices of bad traces are saved for later purging
%      %####################################################################
%         up_dumped = up_dumped +1;
%         dump_handles = [dump_handles j];
%         
%         subplot(3,1,1)
%         plot(peaks_up(1,j),peaks_up(2,j),'ko')
%         plot([(peaks_up(1,j)-pre*dt):dt:(peaks_up(1,j)+post*dt)], ...
%                 X.spikes(1:(pre+post+1),j),'m')
%         
%         subplot(3,1,2)
%         plot(spike_x,X.spikes(:,j),'k')
%                 
%     else
%         %#If we have passed ALL the tests:
%         peaks_up3(:,k)= peaks_up(:,j);
%         %# Store the next good value in the output array.
%         subplot(3,1,1)
%         plot(peaks_up3(1,k),peaks_up3(2,k),'go')     
%     end  %# While Loop
%     %# Now re-plot the trace being examined in black (cover the green)
%     subplot(3,1,2)
%     plot(spike_x,X.spikes(:,j),'k')
% end
% 
% disp([int2str(up_dumped) ' traces rejected (third pass)'])
% 
% Goodspikes3 = X.spikes;
% if ~isempty(dump_handles)
% Goodspikes3(:,dump_handles)= [];
% end
% 
% if ~isempty(Goodspikes3)
% Good_spksavg3 = (mean(Goodspikes3'))';
% Spike_Times = peaks_up3(1,:);   
% else
%     Spike_Times = [];
% end
% 
% 
% %#########################################################################
% 
% 
% % outfilename = [dirname 'Spike_Times_' filename '_' ... 
% %         int2str(starting_trace_number) '_' ...
% %         int2str(ending_trace_number)];
% % save(outfilename, 'Spike_Times') 
% 
% len =size(Goodspikes3);
% dispstr1 = ['Number of spikes saved= ' int2str(len(2))];
% disp(dispstr1);
% 
% dispstr2 = ['Number of traces rejected= ' int2str(up_dumped)];
% disp(dispstr2);
% 
% %# Re-draw middle plot to include only detected spikes and average
% %# spike traces
% subplot(3,1,2)
% cla;
% plot(spike_x,X.spksavg,'r','LineWidth', 2)
% if ~isempty(Goodspikes3)
%     plot(spike_x,Goodspikes3,'color',[0.8 0.8 0.8])
%     plot(spike_x,Good_spksavg1,'k:','LineWidth',2)
%     plot(spike_x,Good_spksavg2,'g','LineWidth', 2)
%     plot(spike_x,Good_spksavg3,'b','LineWidth', 2)
% end
% text(9, 270,dispstr1,'FontSize',10,'Color','blue')
% text(11,220,['rms2 = ' int2str(rms2)],'color','blue')
% text(2,260,'Red Trace: Average of all traces','color','red')
% text(2,220,'Blue Trace: Average of Detected Spikes','color','blue')
% 
% %# Set up a third plot to display rejected traces, and compare with the
% %# average of detected spikes
% subplot(3,1,3)
% hold on
% grid on
% if ~isempty(dump_handles)
% plot(spike_x,X.spikes(:,dump_handles),'m')
% end
% if ~isempty(Goodspikes3)
% plot(spike_x,Good_spksavg3,'b','LineWidth', 2)
% end
% text(9, 270,dispstr2,'FontSize',10,'Color','blue')
% %ylim([-300 300]);
% ylabel('Potential, [mV]')
% xlabel('Time, [msec]')
% % title([filename ' Traces ' int2str(starting_trace_number) '-'...
% %     int2str(ending_trace_number) ' Rejected Spike ' ...
% %     'Traces (3 ^r^d pass fit to average disc.)'])
% 
% 
% 
% 
% 
% 
% 
