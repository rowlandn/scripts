function gauss_rd = glr_rd_SE(spike_times, Fs, max_time, flag, reference_period, response_period, gauss_width, min_dur, stdev)

% glr_rd_SE This function is used in conjunction with the gaussian local
% rate coding algorithm. Once a set of spike trains has been transformed 
% into an average gaussian trace, this function is used to search for 
% significant deviations of the trace from the mean of a specified period of 
% the trace. For example, in a typical stimulation paradigm, a cell is allowed 
% to fire spontaneously from time zero to a specified time, such as one second.  
% After 1 second of spontaneous activity, some stimulus occurs.  The response 
% of the cell is then measured for some time after the stimulus, for example, 
% for an additional 2.5 seconds. To detect signficant stimulus-evoked responses 
% in the average spike train, the spike train would again first be convolved with 
% gaussian curves of the appropriate width.  Significant deviations of the 
% post-stimulus period, from 1 to 2.5 seconds, from the mean of the spontaneous
% period, 0 to 1 s, would then be calculated and plotted.  
%
% gauss_tracerd = glr_rd_SE(spike_times, gauss_width, Fs, max_time, flag, ...
%                             reference_period, response_period, min_dur, stdev)
% 
%      spike_times = same as for glr_SE
%      gauss_width = same as for glr_SE
%               Fs = same as for glr_SE
%         max_time = same as for glr_SE
%             flag = integer identifier for spike train
% reference_period = significant deviations of the post-stimulus period will
%                    be calculated from the mean of this period in ms
%  response_period = significant deviations of this period will be calculated
%                    from the mean of the reference period
%          min_dur = minimum duration in ms that the trace must exceed the signficance 
%                    threshold to be counted as a response
%            stdev = significance threshold that the trace must exceed to be counted as 
%                    a significant response
% 
% Example 1: A = load_PCDX_SE('/Raw/viv05/viv0518d.all','2-51',4);
%            spike_times = findspikes_win_SE(A,10,{-200,-60,.1,2},1);
%            sig_responses = glr_rd_SE(spike_times, 5, 10, 2500, 1, [0 1000], [1000 2500], 20, 3);

% Change reference_period(1) == 0 to reference_period(1) == 1
% User is trying to start at the beginning of the trace (0 ms), but
% the script thinks it is a zero index, which is not possible


if reference_period(1) == 0
    reference_period(1) = .1;
end

% Scale parameters
reference_period = reference_period*Fs;
response_period = response_period*Fs;

% Set up waitbar
progbar = waitbar(0, 'Convolving...');

for i = 1:size(spike_times,1)
    
    waitbar(i/length(spike_times), progbar)
    
    % Check desired length of signal
    if max(spike_times{i,1}) > max_time
    warndlg(['Desired maximum signal length too short.  Please',...
        ' enter a desired maximum signal length that exceeds the',...
        ' maximum spike time.']);
    break
    end
    
    % Scale parameters into data points
    spike_index = round(spike_times{i,1}*Fs);
    max_index = max_time*Fs;
    spike_train = zeros(max_index,1);
    spike_train(spike_index) = 1;


    % Calculate and scale gaussian curve
    sigma = gauss_width*Fs*10;
    time_range = [-100*gauss_width:100*gauss_width]*10;
    gauss_curve = (1./(sqrt(2*pi)*sigma))*(exp(-time_range.^2/(2*sigma^2)));
    gauss_curve = gauss_curve/trapz(gauss_curve);
    
    glr = conv(spike_train, gauss_curve)*10000;
    m_cutoff = fix((length(glr)-length(spike_train))/2);
    glr = glr(1+m_cutoff:max_index+m_cutoff);
    
    gauss_mean(:,i) = glr;
    
end

% Calculate average gaussian trace and scale time axis
mean_gauss = mean(gauss_mean,2);
mean_gauss(:,2) = [1:size(mean_gauss,1)]'/Fs;
gauss_rd = mean_gauss;
%assignin('base','gauss_rd',gauss_rd)

% Close Progess Bar
close(progbar)

% Calculate mean and std of reference period
mean_reference_period = mean(gauss_rd(reference_period(1):reference_period(2),1));
mean_plus_std = mean_reference_period + stdev*std(gauss_rd(reference_period(1):reference_period(2),1));
mean_minus_std = mean_reference_period - stdev*std(gauss_rd(reference_period(1):reference_period(2),1));

% Find significant responses
find_sig_responses = find(gauss_rd(response_period(1):response_period(2),1) > mean_plus_std | gauss_rd(response_period(1):response_period(2),1) < mean_minus_std);
assignin('base','find_sig_responses',find_sig_responses)
diff_find_sig_responses = diff(find_sig_responses);
assignin('base','diff_find_sig_responses',diff_find_sig_responses)
%no_responses = [0 find(diff_find_sig_responses ~= 1)'];
if isempty(find(diff_find_sig_responses ~= 1)) == 1
    no_responses = [0 find_sig_responses(1)];
else
    no_responses = [0 find(diff_find_sig_responses ~= 1)'];
end
assignin('base','no_responses',no_responses)

if no_responses == 0
elseif size(no_responses,2) == 2
    sig_responses = find_sig_responses;
else
    sig_responses = zeros(max(diff(no_responses)),size(no_responses,2)-1);
    for i = 1:size(no_responses,2)-1
        if i == 1
            assignin('base','i',i)
            sig_responses(1:size(find_sig_responses(1:no_responses(i+1)),1),i) = find_sig_responses(1:no_responses(i+1));
            %assignin('base','sig_responses',sig_responses)
        elseif i == size(no_responses,2)-1
            sig_responses(1:size(find_sig_responses(no_responses(i)+1:end),1),i) = find_sig_responses(no_responses(i)+1:end);
            %assignin('base','sig_responses',sig_responses)
        else
            sig_responses(1:size(find_sig_responses(no_responses(i)+1:no_responses(i+1)),1),i) = find_sig_responses(no_responses(i)+1:no_responses(i+1));
            %assignin('base','sig_responses',sig_responses)
        end
    end
end

sig_responses = sig_responses/Fs + response_period(1)/Fs;
%assignin('base','sig_responses',sig_responses)
sig_responses_start = sig_responses(1,:)';
[max_i,max_j] = max(sig_responses-response_period(1)/Fs);
for i = 1:size(max_j,2)
    sig_responses_end(i,1) = sig_responses(max_j(i),i);
end

%% If the end of a response is over the maximum time of the trial,
%% cut it off
find_responses_over = find(sig_responses_end > max_time);
if size(find_responses_over,1) == 0
else
sig_responses_end(find_responses_over) = max_time;
end
%assignin('base','sig_responses_end',sig_responses_end)
%%% Assign statistics to sig_responses_stats_all
%% Flag value
sig_responses_stats_all(:,1) = linspace(flag,flag,size(no_responses,2)-1)';

%% Gaussian width
sig_responses_stats_all(:,2) = linspace(gauss_width,gauss_width,size(no_responses,2)-1)';
%% Standard deviation criterion
sig_responses_stats_all(:,3) = linspace(stdev,stdev,size(no_responses,2)-1)';
%% Minimum duration criterion
sig_responses_stats_all(:,4) = linspace(min_dur,min_dur,size(no_responses,2)-1)';
%% Start of reference period
sig_responses_stats_all(:,5) = linspace(reference_period(1)/Fs,reference_period(1)/Fs,size(no_responses,2)-1)';
%% End of reference period
sig_responses_stats_all(:,6) = linspace(reference_period(2)/Fs,reference_period(2)/Fs,size(no_responses,2)-1)';
%% Mean of reference period
sig_responses_stats_all(:,7) = linspace(mean_reference_period,mean_reference_period,size(no_responses,2)-1)';
%% Start of analysis period
sig_responses_stats_all(:,8) = linspace(response_period(1)/Fs,response_period(1)/Fs,size(no_responses,2)-1)';
%% End of analysis period
sig_responses_stats_all(:,9) = linspace(response_period(2)/Fs,response_period(2)/Fs,size(no_responses,2)-1)';
%% Start of response
sig_responses_stats_all(:,10) = sig_responses_start;
%% End of response
sig_responses_stats_all(:,11) = sig_responses_end;
%% Mean rate of response
for i = 1:size(no_responses,2)-1
assignin('base','i',i)
    sig_responses_stats_all(i,12) = mean(gauss_rd(sig_responses_stats_all(i,10)*10:sig_responses_stats_all(i,11)*10,1));
end
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)
%% Duration of response 
sig_responses_stats_all(:,13) = sig_responses_stats_all(:,11)-sig_responses_stats_all(:,10);
%% Amplitude of response
sig_responses_stats_all(:,14) = (sig_responses_stats_all(:,12)-sig_responses_stats_all(:,7))./sig_responses_stats_all(:,7);
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)


%%% Filter out responses which don't meet the mininum duration criteria
find_short_responses = find(sig_responses_stats_all(:,13) < min_dur);
sig_responses_stats_all(find_short_responses,:) = [];
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)
end
 
%%% Plot
%figure;
max_gaussrd = max(gauss_rd(:,1));
plot(gauss_rd(:,2),gauss_rd(:,1)); hold on
plot([gauss_rd(1,2) gauss_rd(end,2)],[mean_plus_std mean_plus_std],'r')
plot([gauss_rd(1,2) gauss_rd(end,2)],[mean_minus_std mean_minus_std],'r')
if no_responses == 0
else
    for i = 1:size(sig_responses_stats_all,1)
        plot([sig_responses_stats_all(i,10) sig_responses_stats_all(i,11)], [max_gaussrd+.05*max_gaussrd max_gaussrd+.05*max_gaussrd],'k')
    end
end
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
max_gaussrd = max(gauss_rd(:,1));
axis([0 max_time 0 max_gaussrd + .5*max_gaussrd])


%%% Assign output
gauss_rd.gauss_rd = mean_gauss;
gauss_rd.mean_plus_std = mean_plus_std;
gauss_rd.mean_minus_std = mean_minus_std;
if no_responses == 0
    gauss_rd.stats = linspace(0,0,14);
else
    gauss_rd.stats = sig_responses_stats_all;
end

