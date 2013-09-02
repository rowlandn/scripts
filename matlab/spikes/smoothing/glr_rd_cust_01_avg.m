function cat_sig_responses = glr_rd_cust_01_avg(spike_times,Fs,max_time,reference_period,response_period,stds,sp1,sp2,sp3)

% glr_rd_cust_01_avg This is a customized version of the glr_rd_SE
% function.
% This function looks for three different kinds of responses in the gaussian
% trace: an early excitation, a late excitation and an intervening inhibition.
%
% The early excitation is searched for using the following parameters:
% gaussian width = 1 ms
% must occur between 0 and 50 ms after stimulus
% minimum duration = 1 ms
% standard deviation criterion = 3 STD
%  
% The inhibition is searched for using the following parameters:
% gaussian width = 5 ms
% start of response must occur between 0 and 100ms after stimulus
% minimum duration = 20 ms
% standard deviation criterion = 3 STD
%  
% The late excitation is searched for using the following parameters:
% gaussian width = 20 ms
% start of response must occur between 50 and 1000ms after stimulus
% minimum duration = 40 ms
% standard deviation criterion = 3 STD
%
% One combined trace of gaussian width 5ms is printed out at 
% sp1,sp2,sp3, w/ the sle, inh and lle indicated by magenta, 
% green and yellow markers, respectively
%
% cat_sig_responses = glr_rd_cust_01_avg(spike_times,Fs,max_time,reference_period,response_period,sp1,sp2,sp3)
%
% Example: cat_sig_responses = glr_rd_cust_01_avg(DCN_spike_times,10,2500,[1 1000],[1000 2500],2,2,1);
% 
% This returns a structure with separate arrays for the sle, inh and lle responses.
% The columns of each array signifies the following:
%
% 1 - flag value: 1 for sle, 2 for inh, 3 for lle
% 2 - user specified gaussian width
% 3 - user specified standard deviation criterion
% 4 - user specified minimum duration
% 5 - user specified start of reference period
% 6 - user specified end of reference period
% 7 - mean rate during reference period
% 8 - user specified start of response period
% 9 - user specified end of response period
% 10 - start of significant response
% 11 - end of significant response
% 12 - mean rate during response period
% 13 - duration of significant response
% 14 - amplitude of significant response
%
% If no significant responses are found for a particular response component,
% a 1x14 array of zeros is returned.

clear gauss_width min_dur stdev sig_responses_stats_all gauss_rd
clear sig_responses_start sig_responses_end
clear no_responses sig_responses
%%%%%%%%%%%%%%%%%%%% SLE %%%%%%%%%%%%%%%%%%%%
flag = 1;
gauss_width = 1;
min_dur = 1;
stdev = stds(1);

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
gauss_rd_1ms = mean_gauss;
%assignin('base','gauss_rd_1ms',gauss_rd_1ms)

% Close Progess Bar
close(progbar)

% Calculate mean and std of reference period
mean_reference_period = mean(gauss_rd_1ms(reference_period(1):reference_period(2),1));
mean_plus_std_sle = mean_reference_period + stdev*std(gauss_rd_1ms(reference_period(1):reference_period(2),1));
mean_minus_std_sle = mean_reference_period - stdev*std(gauss_rd_1ms(reference_period(1):reference_period(2),1));

% Find significant responses
find_sig_responses = find(gauss_rd_1ms(response_period(1):response_period(2),1) > mean_plus_std_sle | gauss_rd_1ms(response_period(1):response_period(2),1) < mean_minus_std_sle);
%assignin('base','find_sig_responses',find_sig_responses)
diff_find_sig_responses = diff(find_sig_responses);
%assignin('base','diff_find_sig_responses',diff_find_sig_responses)

if isempty(find_sig_responses) == 1
    no_responses = 0;
elseif isempty(find_sig_responses) == 0 & isempty(find(diff_find_sig_responses ~= 1)) == 1
    no_responses = [find_sig_responses(1)];
else
    no_responses = [0 find(diff_find_sig_responses ~= 1)'];
end
%assignin('base','no_responses',no_responses)

if no_responses == 0
    sig_responses = 0;
elseif size(no_responses,2) == 1
    sig_responses = find_sig_responses(1):find_sig_responses(end);
    sig_responses = sig_responses';
else
    sig_responses = zeros(max(diff(no_responses)),size(no_responses,2));
    % assignin('base','sig_responses',sig_responses)
    for i = 1:size(no_responses,2)
        if i == 1
            %assignin('base','i',i)
            sig_responses(1:size(find_sig_responses(1:no_responses(i+1)),1),i) = find_sig_responses(1:no_responses(i+1));
            %assignin('base','sig_responses',sig_responses)
        elseif i == size(no_responses,2)
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
sig_responses_stats_all(1:size(no_responses,2),1) = linspace(1,1,size(no_responses,2))';

%% Flag value
sig_responses_stats_all(:,1) = linspace(flag,flag,size(no_responses,2))';
%% Gaussian width
sig_responses_stats_all(:,2) = linspace(gauss_width,gauss_width,size(no_responses,2))';
%% Standard deviation criterion
sig_responses_stats_all(:,3) = linspace(stdev,stdev,size(no_responses,2))';
%% Minimum duration criterion
sig_responses_stats_all(:,4) = linspace(min_dur,min_dur,size(no_responses,2))';
%% Start of reference period
sig_responses_stats_all(:,5) = linspace(reference_period(1)/Fs,reference_period(1)/Fs,size(no_responses,2))';
%% End of reference period
sig_responses_stats_all(:,6) = linspace(reference_period(2)/Fs,reference_period(2)/Fs,size(no_responses,2))';
%% Mean of reference period
sig_responses_stats_all(:,7) = linspace(mean_reference_period,mean_reference_period,size(no_responses,2))';
%% Start of analysis period
sig_responses_stats_all(:,8) = linspace(response_period(1)/Fs,response_period(1)/Fs,size(no_responses,2))';
%% End of analysis period
sig_responses_stats_all(:,9) = linspace(response_period(2)/Fs,response_period(2)/Fs,size(no_responses,2))';
%% Start of response
sig_responses_stats_all(:,10) = sig_responses_start;
%% End of response
sig_responses_stats_all(:,11) = sig_responses_end;
%% Mean rate of response
for i = 1:size(no_responses,2)
%assignin('base','i',i)
    sig_responses_stats_all(i,12) = mean(gauss_rd_1ms(sig_responses_stats_all(i,10)*10:sig_responses_stats_all(i,11)*10,1));
end
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)
%% Duration of response 
sig_responses_stats_all(:,13) = sig_responses_stats_all(:,11)-sig_responses_stats_all(:,10);
%% Amplitude of response
sig_responses_stats_all(:,14) = (sig_responses_stats_all(:,12)-sig_responses_stats_all(:,7))./sig_responses_stats_all(:,7);
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)


%%% Filter out responses which don't meet the mininum duration criteria
find_short_responses = find(sig_responses_stats_all(:,13) < min_dur);
sig_responses_stats_sle = sig_responses_stats_all;
sig_responses_stats_sle(find_short_responses,:) = [];
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)

 
% Eliminate inhibitions
find_sig_responses_stats_sle_inh = find(sig_responses_stats_sle(:,14) < 0);
sig_responses_stats_sle(find_sig_responses_stats_sle_inh,:) = [];

% Keep only responses fitting boundary criteria
find_sig_responses_stats_sle_start = find(sig_responses_stats_sle(:,10) > 1000 & sig_responses_stats_sle(:,11) < 1050);
sig_responses_stats_sle = sig_responses_stats_sle(find_sig_responses_stats_sle_start,:);
size_sig_responses_stats_sle = size(sig_responses_stats_sle);
if size_sig_responses_stats_sle(1) == 0
    sig_responses_stats_sle = [1 0 0 0 0 0 0 0 0 0 0 0 0 0];
end

%Keep only first response
if size_sig_responses_stats_sle(1) >= 2
    sig_responses_stats_sle = sig_responses_stats_sle(1,:);
end
%assignin('base','sig_responses_stats_sle',sig_responses_stats_sle)

% Descale parameters
reference_period = reference_period/Fs;
response_period = response_period/Fs;

clear gauss_width min_dur stdev sig_responses_stats_all 
clear sig_responses_start sig_responses_end
clear no_responses sig_responses

%%%%%%%%%%%%%%%%%%%% INH %%%%%%%%%%%%%%%%%%%%
flag = 2;
gauss_width = 5;
min_dur = 20;
stdev = stds(2);

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
gauss_rd_5ms = mean_gauss;
%assignin('base','gauss_rd_5ms',gauss_rd_5ms)

% Close Progess Bar
close(progbar)

% Calculate mean and std of reference period
mean_reference_period = mean(gauss_rd_5ms(reference_period(1):reference_period(2),1));
%assignin('base','mean_reference_period',mean_reference_period)
mean_plus_std_inh = mean_reference_period + stdev*std(gauss_rd_5ms(reference_period(1):reference_period(2),1));
mean_minus_std_inh = mean_reference_period - stdev*std(gauss_rd_5ms(reference_period(1):reference_period(2),1));

% Find significant responses
%find_sig_responses = find(gauss_rd_5ms(response_period(1):response_period(2),1) > mean_plus_std_inh | gauss_rd_5ms(response_period(1):response_period(2),1) < mean_minus_std_inh);
find_sig_responses = find(gauss_rd_5ms(response_period(1):response_period(2),1) > mean_plus_std_inh | gauss_rd_5ms(response_period(1):response_period(2),1) < mean_minus_std_inh);
% assignin('base','gauss_rd_5ms',gauss_rd_5ms)
% assignin('base','find_sig_responses',find_sig_responses)
diff_find_sig_responses = diff(find_sig_responses);
%assignin('base','diff_find_sig_responses',diff_find_sig_responses)
%no_responses = [0 find(diff_find_sig_responses ~= 1)'];
if isempty(find_sig_responses) == 1
    no_responses = 0;
elseif isempty(find_sig_responses) == 0 & isempty(find(diff_find_sig_responses ~= 1)) == 1
    no_responses = [find_sig_responses(1)];
else
    no_responses = [0 find(diff_find_sig_responses ~= 1)'];
end
%assignin('base','no_responses',no_responses)

if no_responses == 0
    sig_responses = 0;
elseif size(no_responses,2) == 1
    sig_responses = find_sig_responses(1):find_sig_responses(end);
    sig_responses = sig_responses';
else
    sig_responses = zeros(max(diff(no_responses)),size(no_responses,2));
    for i = 1:size(no_responses,2)
        if i == 1
            %assignin('base','i',i)
            sig_responses(1:size(find_sig_responses(1:no_responses(i+1)),1),i) = find_sig_responses(1:no_responses(i+1));
            %assignin('base','sig_responses',sig_responses)
        elseif i == size(no_responses,2)
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
sig_responses_stats_all(1:size(no_responses,2),1) = linspace(1,1,size(no_responses,2))';

%% Flag value
sig_responses_stats_all(:,1) = linspace(flag,flag,size(no_responses,2))';
%% Gaussian width
sig_responses_stats_all(:,2) = linspace(gauss_width,gauss_width,size(no_responses,2))';
%% Standard deviation criterion
sig_responses_stats_all(:,3) = linspace(stdev,stdev,size(no_responses,2))';
%% Minimum duration criterion
sig_responses_stats_all(:,4) = linspace(min_dur,min_dur,size(no_responses,2))';
%% Start of reference period
sig_responses_stats_all(:,5) = linspace(reference_period(1)/Fs,reference_period(1)/Fs,size(no_responses,2))';
%% End of reference period
sig_responses_stats_all(:,6) = linspace(reference_period(2)/Fs,reference_period(2)/Fs,size(no_responses,2))';
%% Mean of reference period
sig_responses_stats_all(:,7) = linspace(mean_reference_period,mean_reference_period,size(no_responses,2))';
%% Start of analysis period
sig_responses_stats_all(:,8) = linspace(response_period(1)/Fs,response_period(1)/Fs,size(no_responses,2))';
%% End of analysis period
sig_responses_stats_all(:,9) = linspace(response_period(2)/Fs,response_period(2)/Fs,size(no_responses,2))';
%% Start of response
sig_responses_stats_all(:,10) = sig_responses_start;
%% End of response
sig_responses_stats_all(:,11) = sig_responses_end;
%% Mean rate of response
for i = 1:size(no_responses,2)
%assignin('base','i',i)
    sig_responses_stats_all(i,12) = mean(gauss_rd_5ms(sig_responses_stats_all(i,10)*10:sig_responses_stats_all(i,11)*10,1));
end
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)
%% Duration of response 
sig_responses_stats_all(:,13) = sig_responses_stats_all(:,11)-sig_responses_stats_all(:,10);
%% Amplitude of response
sig_responses_stats_all(:,14) = (sig_responses_stats_all(:,12)-sig_responses_stats_all(:,7))./sig_responses_stats_all(:,7);
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)


%%% Filter out responses which don't meet the mininum duration criteria
find_short_responses = find(sig_responses_stats_all(:,13) < min_dur);
sig_responses_stats_inh = sig_responses_stats_all;
sig_responses_stats_inh(find_short_responses,:) = [];
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)

 
% Eliminate excitations
find_sig_responses_stats_inh_exc = find(sig_responses_stats_inh(:,14) > 0);
sig_responses_stats_inh(find_sig_responses_stats_inh_exc,:) = [];

% Keep only responses fitting boundary criteria
find_sig_responses_stats_inh_start = find(sig_responses_stats_inh(:,10) > 1000 & sig_responses_stats_inh(:,10) < 1100);
sig_responses_stats_inh = sig_responses_stats_inh(find_sig_responses_stats_inh_start,:);
size_sig_responses_stats_inh = size(sig_responses_stats_inh);
if size_sig_responses_stats_inh(1) == 0
    sig_responses_stats_inh = [2 0 0 0 0 0 0 0 0 0 0 0 0 0];
end

%Keep only first response
if size_sig_responses_stats_inh(1) >= 2
    sig_responses_stats_inh = sig_responses_stats_inh(1,:);
end
%assignin('base','sig_responses_stats_inh',sig_responses_stats_inh)

% Descale parameters
reference_period = reference_period/Fs;
response_period = response_period/Fs;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear gauss_width min_dur stdev sig_responses_stats_all 
%clear gauss_rd_5ms
clear sig_responses_start sig_responses_end
clear no_responses sig_responses

%%%%%%%%%%%%%%%%%%%% LLE %%%%%%%%%%%%%%%%%%%%
flag = 3;
gauss_width = 20;
min_dur = 40;
stdev = stds(3);

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
gauss_rd_20ms = mean_gauss;
%assignin('base','gauss_rd_20ms',gauss_rd_20ms)

% Close Progess Bar
close(progbar)

% Calculate mean and std of reference period
mean_reference_period = mean(gauss_rd_20ms(reference_period(1):reference_period(2),1));
mean_plus_std_lle = mean_reference_period + stdev*std(gauss_rd_20ms(reference_period(1):reference_period(2),1));
mean_minus_std_lle = mean_reference_period - stdev*std(gauss_rd_20ms(reference_period(1):reference_period(2),1));

% Find significant responses
find_sig_responses = find(gauss_rd_20ms(response_period(1):response_period(2),1) > mean_plus_std_lle | gauss_rd_20ms(response_period(1):response_period(2),1) < mean_minus_std_lle);
%assignin('base','find_sig_responses',find_sig_responses)
diff_find_sig_responses = diff(find_sig_responses);
%assignin('base','diff_find_sig_responses',diff_find_sig_responses)

if isempty(find_sig_responses) == 1
    no_responses = 0;
elseif isempty(find_sig_responses) == 0 & isempty(find(diff_find_sig_responses ~= 1)) == 1
    no_responses = [find_sig_responses(1)];
else
    no_responses = [0 find(diff_find_sig_responses ~= 1)'];
end
%assignin('base','no_responses',no_responses)

if no_responses == 0
    sig_responses = 0;
elseif size(no_responses,2) == 1
    sig_responses = find_sig_responses(1):find_sig_responses(end);
    sig_responses = sig_responses';
else
    sig_responses = zeros(max(diff(no_responses)),size(no_responses,2));
    for i = 1:size(no_responses,2)
        if i == 1
            %assignin('base','i',i)
            sig_responses(1:size(find_sig_responses(1:no_responses(i+1)),1),i) = find_sig_responses(1:no_responses(i+1));
            %assignin('base','sig_responses',sig_responses)
        elseif i == size(no_responses,2)
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
sig_responses_stats_all(1:size(no_responses,2),1) = linspace(1,1,size(no_responses,2))';

%% Flag value
sig_responses_stats_all(:,1) = linspace(flag,flag,size(no_responses,2))';
%% Gaussian width
sig_responses_stats_all(:,2) = linspace(gauss_width,gauss_width,size(no_responses,2))';
%% Standard deviation criterion
sig_responses_stats_all(:,3) = linspace(stdev,stdev,size(no_responses,2))';
%% Minimum duration criterion
sig_responses_stats_all(:,4) = linspace(min_dur,min_dur,size(no_responses,2))';
%% Start of reference period
sig_responses_stats_all(:,5) = linspace(reference_period(1)/Fs,reference_period(1)/Fs,size(no_responses,2))';
%% End of reference period
sig_responses_stats_all(:,6) = linspace(reference_period(2)/Fs,reference_period(2)/Fs,size(no_responses,2))';
%% Mean of reference period
sig_responses_stats_all(:,7) = linspace(mean_reference_period,mean_reference_period,size(no_responses,2))';
%% Start of analysis period
sig_responses_stats_all(:,8) = linspace(response_period(1)/Fs,response_period(1)/Fs,size(no_responses,2))';
%% End of analysis period
sig_responses_stats_all(:,9) = linspace(response_period(2)/Fs,response_period(2)/Fs,size(no_responses,2))';
%% Start of response
sig_responses_stats_all(:,10) = sig_responses_start;
%% End of response
sig_responses_stats_all(:,11) = sig_responses_end;
%% Mean rate of response
for i = 1:size(no_responses,2)
%assignin('base','i',i)
    sig_responses_stats_all(i,12) = mean(gauss_rd_20ms(sig_responses_stats_all(i,10)*10:sig_responses_stats_all(i,11)*10,1));
end
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)
%% Duration of response 
sig_responses_stats_all(:,13) = sig_responses_stats_all(:,11)-sig_responses_stats_all(:,10);
%% Amplitude of response
sig_responses_stats_all(:,14) = (sig_responses_stats_all(:,12)-sig_responses_stats_all(:,7))./sig_responses_stats_all(:,7);
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)


%%% Filter out responses which don't meet the mininum duration criteria
find_short_responses = find(sig_responses_stats_all(:,13) < min_dur);
sig_responses_stats_lle = sig_responses_stats_all;
sig_responses_stats_lle(find_short_responses,:) = [];
%assignin('base','sig_responses_stats_all',sig_responses_stats_all)

 
% Eliminate inhibitions
find_sig_responses_stats_lle_inh = find(sig_responses_stats_lle(:,14) < 0);
sig_responses_stats_lle(find_sig_responses_stats_lle_inh,:) = [];

% Keep only responses fitting boundary criteria
find_sig_responses_stats_lle_start = find(sig_responses_stats_lle(:,10) > 1050 & sig_responses_stats_lle(:,10) < 2000);
sig_responses_stats_lle = sig_responses_stats_lle(find_sig_responses_stats_lle_start,:);
size_sig_responses_stats_lle = size(sig_responses_stats_lle);
if size_sig_responses_stats_lle(1) == 0
    sig_responses_stats_lle = [3 0 0 0 0 0 0 0 0 0 0 0 0 0];
end

%Keep only first response
if size_sig_responses_stats_lle(1) >= 2
    sig_responses_stats_lle = sig_responses_stats_lle(1,:);
end
%assignin('base','sig_responses_stats_lle',sig_responses_stats_lle)

% Descale parameters
reference_period = reference_period/Fs;
response_period = response_period/Fs;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat_sig_responses.sig_responses.sle = sig_responses_stats_sle;
cat_sig_responses.sig_responses.inh = sig_responses_stats_inh;
cat_sig_responses.sig_responses.lle = sig_responses_stats_lle;
cat_sig_responses.gauss_traces.one_ms = gauss_rd_1ms;
cat_sig_responses.gauss_traces.five_ms = gauss_rd_5ms;
cat_sig_responses.gauss_traces.twenty_ms = gauss_rd_20ms;
cat_sig_responses.std.one_ms.plus = mean_plus_std_sle;
cat_sig_responses.std.one_ms.minus = mean_minus_std_sle;
cat_sig_responses.std.five_ms.plus = mean_plus_std_inh;
cat_sig_responses.std.five_ms.minus = mean_minus_std_inh;
cat_sig_responses.std.twenty_ms.plus = mean_plus_std_lle;
cat_sig_responses.std.twenty_ms.minus = mean_minus_std_lle;
cat_sig_responses.periods.reference = reference_period;
cat_sig_responses.periods.response = response_period;

 
%%% Plot
if nargin == 6
else
    figure
    subplot(sp1,sp2,sp3)
    max_gauss_rd_5ms = max(gauss_rd_5ms(:,1));
    plot(gauss_rd_5ms(:,2),gauss_rd_5ms(:,1),'.','MarkerSize',1); hold on
    plot([reference_period(1) reference_period(2)],[mean_plus_std_inh mean_plus_std_inh],'r:')
    plot([response_period(1) response_period(2)],[mean_plus_std_inh mean_plus_std_inh],'r')
    plot([reference_period(1) reference_period(2)],[mean_minus_std_inh mean_minus_std_inh],'r:')
    plot([response_period(1) response_period(2)],[mean_minus_std_inh mean_minus_std_inh],'r')

    if max_gauss_rd_5ms+.05*max_gauss_rd_5ms > mean_plus_std_inh
        y_max = max_gauss_rd_5ms+.05*max_gauss_rd_5ms;
    elseif mean_plus_std_inh > max_gauss_rd_5ms+.05*max_gauss_rd_5ms 
        y_max = mean_plus_std_inh+.1*mean_plus_std_inh;
    end

    if cat_sig_responses.sig_responses.sle(2) == 0
        text(-1,y_max+y_max*.05,'X','FontSize',6,'Color',[1 0 1])
    else
        plot([cat_sig_responses.sig_responses.sle(1,10) cat_sig_responses.sig_responses.sle(1,11)], [y_max+y_max*.05 y_max+y_max*.05],'m','LineWidth',2)
    end

    if cat_sig_responses.sig_responses.inh(2) == 0
        text(-1,y_max+y_max*.15,'X','FontSize',6,'Color',[0 1 0])
    else
        plot([cat_sig_responses.sig_responses.inh(1,10) cat_sig_responses.sig_responses.inh(1,11)], [y_max+y_max*.15 y_max+y_max*.15],'g','LineWidth',2)
    end

    if cat_sig_responses.sig_responses.lle(2) == 0
        text(-1,y_max+y_max*.25,'X','FontSize',6,'Color',[0 1 1])
    else
        plot([cat_sig_responses.sig_responses.lle(1,10) cat_sig_responses.sig_responses.lle(1,11)], [y_max+y_max*.25 y_max+y_max*.25],'c','LineWidth',2)
    end
    
    max_gauss_rd_5ms = max(gauss_rd_5ms(:,1));
    axis([0 max_time 0 max_gauss_rd_5ms + .5*max_gauss_rd_5ms])
    set(gca,'Box','off')
end




