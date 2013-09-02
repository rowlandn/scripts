function spikes_SE = findspikes_win_SE(traces,Fs,spike_discrim_params,plotit,artifact,rmst,stim_t,art_plottype)

% a = inputname(9);
% assignin('base','a',a)
% findspikes_win_SE This function performs window-based spike 
% discrimination on single peaks arising from a baseline of zero.  Traces
% must be a matrix  with the rows representing voltages at successive
% sampling times  and the columns representing different traces.  The
% sampling frequency (Fs) in  kHz is used to scale the trace and spike
% times into ms.  Peaks may be either positive or negative. Peaks must
% occur  within a user-defined window of minimum and maximum amplitude  and
% minimum and maximum time width.  These spike discrimination parameters
% should be given as a cell in the following sequence: 
% spike_discrim_params = {peak_min_amp peak_max_amp peak_min_time
% peak_max_time}.   Enter 1 for plotit if you would like to see a plot of
% the discrimination results for  the last trace.  Otherwise, enter 0.
% Output consists of a cell of spike time arrays with the row of the cell
% array representing the corresponding trace from which the spikes were
% discriminated.
%
% spike_times = findspikes_win_SE(traces,Fs,spike_discrim_params,plotit)   
%
% Example 1: A = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',4);
%            spike_times = findspikes_win_SE(A,10,{30 100 .1 2},1); 
%                          (for positive peaks)
%                             
% Example 2: A = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',4);
%            spike_times = findspikes_win_SE(A,10,{-100 -30 .1 2},1); 
%                         (for negative peaks, notice that the lower
%                          of the two amplitudes is in the peak_min_amp
%                          position)
% Example 3: A = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',4);
%            B = {-100 -30 .1 2};
%            spike_times = findspikes_win_SE(A,10,B,1); 
%             (you can also save the parameters in a variable)


% traces = loadtraces('/F/Data/Raw/viv05/viv0504c.all','28',4);
% samp_freq_Hz = 10000;
% peak_min_amp = -200;
% peak_max_amp = -40;
% peak_min_time = .1;
% peak_max_time = 1;
% category = 1;
% channel = 4;
% flag = 1;
% plotit = 0;
% filter = 1;

[peak_min_amp peak_max_amp peak_min_time peak_max_time] = deal(spike_discrim_params{:}); 

% Set up waitbar
progbar = waitbar(0, 'Finding Spikes...');

% Begin for loop to analyze each trace
for j = 1:size(traces,2)
%assignin('base','j',j)
    
% Each time the for loop runs, clear relevant variables
clear find_peaks_in_amp_window peaks_in_amp_window 
clear diff_peaks_in_amp_window find_diff_peaks_in_amp_window_above_10
clear no_of_runs runs_time size_of_runs_time runs_voltage size_of_runs_voltage
clear find_new_runs_time find_new_runs_time_indices
clear find_voltage_peaks_in_amp_window_max indices_find_voltage_peaks_in_amp_window_max
clear find_time_peaks_in_amp_window_max_corr_to_voltage_peaks
clear spike_times spike_peaks diff_peaks_in_amp_window_2
clear diff_peaks_in_amp_window_3 spike_nos run_size sum_run_size
clear length_spike_times trace
clear length_find_runs_time length_find_runs_time_indices
clear length_find_runs_time_indices_zeroes 
clear length_find_new_runs_time length_find_new_runs_time_indices
clear length_find_new_runs_time_indices_zeroes
clear length_find_runs_voltage length_find_runs_voltage_indices
clear length_find_runs_voltage_indices_zeroes 
clear indices_find_voltage_peaks_in_amp_window_max
clear add_run_str_* runs*



waitbar(j/size(traces,2), progbar)
    

    

% Assign a single trace (represented by a single column)
% to the variable 'trace' during the for loop
trace = traces(:,j);

%assignin('base','trace',trace)

% Assign time stamps in ms to each sample.  Time stamps are
% calculated from the length of each trace and the 
% sampling frequency during acquisition.  Add the time stamps 
% as a second column to 'trace'.
sampling_times = [1:length(trace)];
samp_scale_factor_ms = (Fs*1000)/1000;
%assignin('base','samp_scale_factor_ms',samp_scale_factor_ms)
trace(:,2) = trace;
trace(:,1) = sampling_times';
trace(:,1) = trace(:,1)/samp_scale_factor_ms;
%assignin('base','trace',trace)

% Locate values in 'trace' that fall between the 
% user-defined min and max amplitude.  Also, keep
% the original indices from trace.
if peak_min_amp < peak_max_amp
    find_peaks_in_amp_window = find(trace(:,2) >= ...
        peak_min_amp & trace(:,2) <= peak_max_amp);
    peaks_in_amp_window = trace(find_peaks_in_amp_window,:);
    peaks_in_amp_window(:,[2:3]) = peaks_in_amp_window;
    %assignin('base','peaks_in_amp_window',peaks_in_amp_window)
    peaks_in_amp_window(:,1) = find_peaks_in_amp_window;
    %assignin('base','find_peaks_in_amp_window',find_peaks_in_amp_window)
    %assignin('base','peaks_in_amp_window',peaks_in_amp_window)
elseif peak_min_amp > peak_max_amp
    find_peaks_in_amp_window = find(trace(:,2) <= ...
        peak_min_amp & trace(:,2) >= peak_max_amp);
    peaks_in_amp_window = trace(find_peaks_in_amp_window,:);
    peaks_in_amp_window(:,[2:3]) = peaks_in_amp_window;
    peaks_in_amp_window(:,1) = find_peaks_in_amp_window;
    %assignin('base','find_peaks_in_amp_window',find_peaks_in_amp_window)
    %assignin('base','peaks_in_amp_window',peaks_in_amp_window)
end


% % IF you want to see which peaks the algorithm found, uncomment
% % this section
% figure 
% plot(trace(:,1),trace(:,2)); hold on
% plot(peaks_in_amp_window(:,2),peaks_in_amp_window(:,3),'r.')
% line_peak_min_amp = [1:max(trace(:,1))];
% line_peak_max_amp = [1:max(trace(:,1))];
% plot(line_peak_min_amp,peak_min_amp,'r')
% plot(line_peak_max_amp,peak_max_amp,'r')
% %assignin('base','peaks_in_amp_window',peaks_in_amp_window)


% This next piece of code slows the script down
% considerably but is crucial for the discrimination between
% real spikes and stimulation artifacts or noise.  Peaks_in_amp_window 
% simply finds the values of trace that fall within the 
% user-defined min and max amp values.  Ideally, these 
% sections of the trace would rise into the window and 
% fall back out of the window such that the peak
% occurs inside the window.  However, if we encounter a very large
% stimulation artifact or noise, the trace may rise into the amplitude 
% window and then continue beyond it leaving a peak that lies somewhere 
% outside the window, which is not what we want. This piece 
% of code is intended to discriminate between these two cases so 
% that one does not have to worry about discriminating stimulation 
% artifacts or noise manually. It starts by taking the values 
% found in peaks_in_amp_window.  It then finds the time stamps 
% associated with these voltage values.  For the earliest time stamps, 
% it finds the time stamp and voltage value that occurs just prior 
% to that time stamp. And for the latest time stamp, it finds the time stamp
% and voltage value that occurs just after that time stamp.
% These new first and last time and voltage values must occur
% outside of the window on the same side to be identified as a spike. 
% If not, these new first and last time and voltage values will occur 
% on opposite sides of the window and are discarded.  

% if filter == 1
%     filter_spikes
% else
% end


% Here, we separate contiguous time stamps and 
% corresponding voltage values into what the script
% calls runs, named for a run of contiguous values
% separated by some time before the next run.  Each run
% represents a spike.  If the length of the run
% is shorter than the user defined minimum time window or 
% longer than the user defined maximum time window, the 
% spike is discarded.  Once we have the final set of runs,
% we simply take the maximum voltage value of the run 
% (for positive peaks) or the minimum of the voltage run
% (for negative peaks).  Then, we find the corresponding
% time stamps for those voltage values.

% if filter == 0
     diff_peaks_in_amp_window_3 = diff(peaks_in_amp_window(:,2));
%    assignin('base','diff_peaks_in_amp_window_3',diff_peaks_in_amp_window_3);
% else    
%    diff_peaks_in_amp_window_3 = diff(peaks_in_amp_window(:,1));
%    assignin('base','diff_peaks_in_amp_window_3',diff_peaks_in_amp_window_3);
% end

% Successive time stamps will have a difference of .1 between them.
% These time stamps indicate a spike.  Therefore, counting the number of 
% successive time stamps that have a difference greater than .1 
% between them will give you the number of spikes (or runs as they are
% initially referred to here.)
find_diff_peaks_in_amp_window_above_10 = find(diff_peaks_in_amp_window_3 > .2);
% assignin('base','find_diff_peaks_in_amp_window_above_10',...
%     find_diff_peaks_in_amp_window_above_10)
% if diff_peaks_in_amp_window_3(end) < .11
%     find_diff_peaks_in_amp_window_above_10 = ...
%         [0;find_diff_peaks_in_amp_window_above_10;length(diff_peaks_in_amp_window_3)+1];
%     %assignin('base','find_diff_peaks_in_amp_window_above_10',...
%     %find_diff_peaks_in_amp_window_above_10)
% else
%     find_diff_peaks_in_amp_window_above_10 = ...
%         [0;find_diff_peaks_in_amp_window_above_10;length(diff_peaks_in_amp_window_3)+1];
%     %assignin('base','find_diff_peaks_in_amp_window_above_10',...
%     %find_diff_peaks_in_amp_window_above_10)
% end

find_diff_peaks_in_amp_window_above_10 = ...
        [0;find_diff_peaks_in_amp_window_above_10;length(diff_peaks_in_amp_window_3)+1];
%assignin('base','find_diff_peaks_in_amp_window_above_10',...
        %find_diff_peaks_in_amp_window_above_10)
no_of_runs = length(find_diff_peaks_in_amp_window_above_10);
%assignin('base','no_of_runs',no_of_runs)

runs_time = zeros(no_of_runs,peak_max_time*10);
runs_voltage = zeros(no_of_runs,peak_max_time*10);
%assignin('base','runs_voltage',runs_voltage)
%assignin('base','runs_time',runs_time)
for i = 1:no_of_runs - 1
    %assignin('base','i',i)
    %if find_diff_peaks_in_amp_window_above_10(end-1) == find_diff_peaks_in_amp_window(end)
     single_run_time = peaks_in_amp_window((find_diff_peaks_in_amp_window_above_10(i))+1:find_diff_peaks_in_amp_window_above_10(i+1));
     %assignin('base','single_run_time',single_run_time)
%      if filter == 0
        single_run_voltage = ...
            peaks_in_amp_window((find_diff_peaks_in_amp_window_above_10(i))+1:find_diff_peaks_in_amp_window_above_10(i+1),3);
%      else
%       single_run_voltage = ...
%            peaks_in_amp_window((find_diff_peaks_in_amp_window_above_10(i))+1:find_diff_peaks_in_amp_window_above_10(i+1),2);
%      end
     %assignin('base','single_run_voltage',single_run_voltage)
     runs_time(i,1:length(single_run_time)) = single_run_time;
     %assignin('base','runs_time',runs_time)
     runs_voltage(i,1:length(single_run_voltage)) = single_run_voltage';
     %assignin('base','runs_voltage',runs_voltage)
end

runs_time = runs_time';
%assignin('base','runs_time',runs_time)
% if filter == 0
     runs_time = runs_time/10;
% end
%assignin('base','runs_time',runs_time)


% Erase any empty columns of runs_time and runs_voltage
for i = 1:size(runs_time,2)
    length_find_runs_time = length(find(runs_time(:,i)));
        if  length_find_runs_time == 0
            length_find_runs_time_indices(1,i) = 0;
        else
            length_find_runs_time_indices(1,i) = 1;
            size_of_runs_time = size(runs_time);
        end
end
%assignin('base','length_find_runs_time',length_find_runs_time)
%assignin('base','length_find_runs_time_indices',length_find_runs_time_indices)
length_find_runs_time_indices_zeroes = find(length_find_runs_time_indices == 0);

for m = length_find_runs_time_indices_zeroes(end):length_find_runs_time_indices_zeroes(1)
    runs_time(:,m) = [];
end
%assignin('base','runs_time',runs_time)

runs_voltage = runs_voltage';
%assignin('base','runs_voltage',runs_voltage)
%%assignin('base','runs_time',runs_time)
for i = 1:size(runs_voltage,2)
    length_find_runs_voltage = length(find(runs_voltage(:,i)));
        if length_find_runs_voltage == 0
            length_find_runs_voltage_indices(1,i) = 0;
        else
            length_find_runs_voltage_indices(1,i) = 1;
        end
end

length_find_runs_voltage_indices_zeroes = find(length_find_runs_voltage_indices == 0);
%assignin('base','length_find_runs_voltage',length_find_runs_voltage)
%assignin('base','length_find_runs_voltage_indices',length_find_runs_voltage_indices)
%assignin('base','length_find_runs_voltage_indices_zeroes',length_find_runs_voltage_indices_zeroes)
for o = length_find_runs_voltage_indices_zeroes(end):length_find_runs_voltage_indices_zeroes(1)
    runs_voltage(:,o) = [];
end
%%assignin('base','runs_time',runs_time)
%assignin('base','runs_voltage',runs_voltage)

% Erase any runs of runs_time that are 
% either shorter or longer than the user
% defined time window.  Use the same
% indices to eliminate columns in
% runs_voltage
for i = 1:size(runs_time,2)
    length_find_new_runs_time = length(find(runs_time(:,i)));
        if length_find_new_runs_time/samp_scale_factor_ms < peak_min_time | length_find_new_runs_time/samp_scale_factor_ms > peak_max_time
            length_find_new_runs_time_indices(1,i) = 0;
        else
            length_find_new_runs_time_indices(1,i) = 1;
        end
end
%assignin('base','length_find_new_runs_time_indices',...
    %length_find_new_runs_time_indices)

%assignin('base','length_find_new_runs_time',length_find_new_runs_time)


if size(find(length_find_new_runs_time_indices(1,:))) == [1 0]
   warndlg('No spikes found.  Please choose different time window parameters.')
end

if isempty(find(length_find_new_runs_time_indices == 0)) == 1
else
    length_find_new_runs_time_indices_zeroes = find(length_find_new_runs_time_indices == 0);
        for n = length_find_new_runs_time_indices_zeroes(end):length_find_new_runs_time_indices_zeroes(1)
            runs_time(:,n) = [];
            %assignin('base','length_find_new_runs_time_indices_zeroes',...
    %length_find_new_runs_time_indices_zeroes)
        end
end

if isempty(find(length_find_new_runs_time_indices == 0)) == 1
else
    for p = length_find_new_runs_time_indices_zeroes(end):length_find_new_runs_time_indices_zeroes(1)
            runs_voltage(:,p) = [];
            %assignin('base','length_find_new_runs_time_indices_zeroes',...
        %length_find_new_runs_time_indices_zeroes)
    end
end

%%assignin('base','runs_time',runs_time)
size_of_runs_time = size(runs_time);



if peak_min_amp < 0 & peak_max_amp < 0
    [find_voltage_peaks_in_amp_window_max,indices_find_voltage_peaks_in_amp_window_max] = min(runs_voltage);
else    
    [find_voltage_peaks_in_amp_window_max,indices_find_voltage_peaks_in_amp_window_max] = max(runs_voltage);
end

%assignin('base','indices_find_voltage_peaks_in_amp_window_max',indices_find_voltage_peaks_in_amp_window_max)
%assignin('base','find_voltage_peaks_in_amp_window_max',find_voltage_peaks_in_amp_window_max)
spike_peaks_str = ['spike_peaks_trace_',num2str(j),' = find_voltage_peaks_in_amp_window_max;'];
eval(spike_peaks_str)

% if exist(['spike_peaks_trace_',num2str(j)]) == 1
% assignin('base',['spike_peaks_trace_',num2str(j)],find_voltage_peaks_in_amp_window_max)
% end
spike_peaks = find_voltage_peaks_in_amp_window_max;
size_spike_peaks = size(spike_peaks); 
for i = 1:size(runs_time,2)
    find_time_peaks_in_amp_window_max_corr_to_voltage_peaks(1,i) = ...
        runs_time(indices_find_voltage_peaks_in_amp_window_max(i),i);
end
%assignin('base',...
%     'find_time_peaks_in_amp_window_max_corr_to_voltage_peaks', ...
%     find_time_peaks_in_amp_window_max_corr_to_voltage_peaks)
spike_times_str = ['spike_times_trace_',num2str(j),' = find_time_peaks_in_amp_window_max_corr_to_voltage_peaks;'];
%assignin('base','spike_times_str',spike_times_str)
eval(spike_times_str)
% if exist(['spike_times_trace_',num2str(j)]) == 1
% assignin('base',['spike_times_trace_',num2str(j)],find_time_peaks_in_amp_window_max_corr_to_voltage_peaks)
% end

spike_times = find_time_peaks_in_amp_window_max_corr_to_voltage_peaks;
length_spike_times = length(spike_times);
%assignin('base','spike_peaks',spike_peaks)
%assignin('base','spike_times',spike_times)
%assignin('base','length_spike_times',length_spike_times)

peaks{j,1} = spike_peaks;
spikes{j,1} = spike_times;
% %assignin('base','j',j)
% %assignin('base','spikes_SE',spikes_SE)

end

%% Close Progess Bar
close(progbar)

%% Build structure
spikes_SE.traces = traces;
spikes_SE.times.win = spikes;
spikes_SE.peaks.win = peaks;
spikes_SE.discrim_params = spike_discrim_params;


%assignin('base','trace',trace)

%Plot spike times, spike peaks and thresholds
%for last trace
if plotit == 0 & nargin == 4
elseif  plotit == 1 & nargin == 4
    figure
    if peak_min_amp > 0 & peak_max_amp > 0
        if size_spike_peaks(2) == 0
            spike_times = [0];
            plot(trace(:,2)); hold on
            warndlg('No spikes found.  Please choose different window parameters')
        else
            plot(trace(:,1),trace(:,2)); hold on
            plot(spike_times,spike_peaks,'ro')
            x_lim = get(gca,'xlim');
            plot([x_lim(1) x_lim(2)],[peak_min_amp peak_min_amp],'r')
            plot([x_lim(1) x_lim(2)],[peak_max_amp peak_max_amp],'r')
            y_lim = get(gca,'ylim');
            xlabel('ms')
            ylabel('mV')
            axis([0 max(trace(:,1)) y_lim(1) y_lim(2)+.2*y_lim(2)])
        end
    else
            plot(trace(:,1),trace(:,2)); hold on
            plot(spike_times,spike_peaks,'ro')
            %plot(spike_times,1,'ro')
            x_lim = get(gca,'xlim');
            plot([x_lim(1) x_lim(2)],[peak_min_amp peak_min_amp],'r')
            plot([x_lim(1) x_lim(2)],[peak_max_amp peak_max_amp],'r')
            y_lim = get(gca,'ylim');
            xlabel('ms')
            ylabel('mV')
            axis([0 max(trace(:,1)) y_lim(1)+.2*y_lim(1) y_lim(2)])
    end

elseif plotit == 1 & nargin == 9 & strmatch('single',art_plottype) == 1
    spikes_SE_2 = findspikes_rms_SE1(traces,spikes_SE,rmst,stim_t,art_plottype);
    spikes_SE.times.win_rms = spikes_SE_2.times_win_rms;
    spikes_SE.peaks.win_rms = spikes_SE_2.peaks_win_rms;
    spikes_SE.rms_params.rejected_spike_counts = spikes_SE_2.rejected_spike_counts;
    spikes_SE.rms_params.rejected_spike_indices = spikes_SE_2.rejected_spike_index;
    spikes_SE.rms_params.rms_values = spikes_SE_2.rms_values;
    spikes_SE.rms_params.rmst = spikes_SE_2.rmst;
    
    figure
    subplot(3,1,1)
    for i = 1:size(spikes_SE_2.collect_spikes.spk_traces,1)
        plot(spikes_SE_2.collect_spikes.spk_traces{i,1}); hold on
    end
    title('All Detected Spikes')

    subplot(3,1,2)
    for i = 1:size(spikes_SE.rms_params.rejected_spike_indices,1)
        plot(spikes_SE_2.collect_spikes.spk_traces{i,1}(:,spikes_SE.rms_params.rejected_spike_indices{i,1}),'r')
        hold on
    end
    title([num2str(sum(spikes_SE.rms_params.rejected_spike_counts)),' rejected spikes out of ',num2str(spikes_SE_2.spike_count_01)])
    % assignin('base','spike_times',spike_times)
    
    subplot(3,1,3)
    plot(traces(:,end),'k'); hold on
    plot(spikes_SE.times.win{end}*10,spikes_SE.peaks.win{end},'go')
    plot(spikes_SE.times.win{end}(spikes_SE.rms_params.rejected_spike_indices{end})*10,spikes_SE.peaks.win{end}(spikes_SE.rms_params.rejected_spike_indices{end}),'ro')
elseif strmatch('all',art_plottype) == 1
    spikes_SE_2 = findspikes_rms_SE1(traces,spikes_SE,rmst,stim_t,art_plottype);
    spikes_SE.times.win_rms = spikes_SE_2.times_win_rms;
    spikes_SE.peaks.win_rms = spikes_SE_2.peaks_win_rms;
    spikes_SE.rms_params.rejected_spike_counts = spikes_SE_2.rejected_spike_counts;
    spikes_SE.rms_params.rejected_spike_indices = spikes_SE_2.rejected_spike_index;
    spikes_SE.rms_params.rms_values = spikes_SE_2.rms_values;    
    spikes_SE.rms_params.rmst = spikes_SE_2.rmst;
    
    figure
    plot(spikes_SE.traces); hold on
    for i = 1:size(spikes_SE.traces,2)
        plot(spikes_SE.times.win_rms{i,1}*10, spikes_SE.peaks.win_rms{i,1},'go')
        plot(spikes_SE.times.win{i,1}(spikes_SE.rms_params.rejected_spike_indices{i,1})*10,spikes_SE.peaks.win{i,1}(spikes_SE.rms_params.rejected_spike_indices{i,1}),'ro')
    end
end



clear d* f* i* r* n* *a   sampling_times single_run_time single_run_voltage size_of_runs_time
clear size_of_runs_voltage size_spike_peaks
clear add_run_str_2 ans cat_str diff_peaks_in_amp_window
clear diff_peaks_in_amp_window_2 diff_peaks_in_amp_window_3
clear find_diff_peaks_in_amp_window_2 
clear find_diff_peaks_in_amp_window_above_10 find_peaks_in_amp_window
clear indices_find_voltage_peaks_in_amp_window
clear j k length_find_new_runs_time length_find_runs_time
clear length_find_runs_time_indices length_spike_times
clear m peaks_in_amp_window





    