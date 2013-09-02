% IntraOpSSEP()
% This function imports an apm file and extracts ecog data and stimulation
% trigger timestamps.  It creates peri-stim ecog voltage levels from each
% ecog channel.

% 7/22/2008 SS: IntraOpSSEP currently loads GL4k data from APM files.
% Use ReadEventData.m created by FHC once it is updated, to import
% data directly from GLR patient files without saving each snapshot to APM.  

% 8/26/2008 SS: normalized SSEP recordings that are montaged from adjact ecog contact
% pairs are processed in the following order:
%               1. find difference between adjacent contact pairs
%               2. find the mean of difference
%               3. normalize mean of difference

clear all
close all
colordef(0,'white');
plotfig1 = 1; % toggle prompt for plotting figure 1, which shows standard SSEP figure with reference to e6
% rename_variables = 0; % RLM: TOGGLE prompt for manually naming emg channels
% ecog_channel_append = '_ecog'; %  SS: default name appendage for ecog channel 


%% Define variables
% Fs = 1000;          % sampling frequency
Fs = 5000;            % 5000 Hz required to capture SSEP
PRE_STIM = -0.05;     % pre-stim period in sec
PST_STIM = 0.05;     % post-stim period in sec
STIM_TRIG_AMP = 5.5;  % amplitude of stimulus trigger

%% Import GL4k data

% AB: GL4K aux ADC values are usigned. When the AUX inputs are configured
% as bipolar, Need to substract midpoint, in order to obtain a bipolar value.
% When configures as unipolar, do not substract anything
%auxADCzero = 8192;  % For bipolar input configuration
auxADCzero = 0;   % For unipolar input configuration


while true

    [fn d]=uigetfile('*.apm','FHC APM files (*.apm)'); % choose abf file from CD or proper directory
    cd(d);
    try
        try
            apmdata;
            disp('Data is already loaded, type ''clear apmdata'' to force reloading');
        catch
            apmdata=APMReadData([d fn]);
        end;
        % Rename the variables according the way abfml.dll was originally doing it
        n_spike = 10;  % hardwire # spike channels 
        n_aux = length(apmdata.aux);
        %vn=strrep(fn,'.abf','');   % This could alter the file name assuming that the substring '.abf' occurs in the actual file name as well, not only as the extension
        vn=strrep(fn,'.apm','');   % This could alter the file name assuming that the substring '.abf' occurs in the actual file name as well, not only as the extension
        episode=1;  % This is for ABF file format compatibility
        
        % First deal with spike channels
        spike_channel = [];
        for i=1:n_spike
            if i > length(apmdata.channels) || isempty(apmdata.channels(i).continuous)
                spike_channel(i).data = [];
                spike_channel(i).time = [];
            else
                spike_channel(i).data = apmdata.channels(i).voltage_calibration(1) * apmdata.channels(i).continuous.';
                %spike_channel(i).time = apmdata.channels(i).time_calibration * 1e-6 * (apmdata.channels(i).start_continuous - apmdata.channels(i).start_trial(1) + (1:length(apmdata.channels(i).continuous)) - 1); % A time vector aligned with start_trial
                spike_channel(i).time = apmdata.channels(i).start_continuous  + (1:length(apmdata.channels(i).continuous)) - 1;
                % Keep samples between start-of-trial and end-of-trial (if any, when manually stopping recording)
                ix = find( spike_channel(i).time < apmdata.channels(i).start_trial(1) );
                if ~isempty(ix)
                    spike_channel(i).data(ix) = [];
                    spike_channel(i).time(ix) = [];
                   disp(sprintf(' Channel %d, trimmed %d leading samples before start-of-recording',i,length(ix)));
                end;
                try     % Check if recording has been stopped manually (the end_trial field is present)
                    apmdata.channels(i).end_trial(1);
                    if (apmdata.channels(i).end_trial(1) > apmdata.channels(i).start_trial(1))
                        ix = find( spike_channel(i).time >= apmdata.channels(i).end_trial(1) );
                        if ~isempty(ix)
                            spike_channel(i).data(ix) = [];
                            spike_channel(i).time(ix) = [];
                            disp(sprintf(' Channel %d, trimmed %d trailing samples after end-of-recording',i,length(ix)));
                            % SS: instead of adding a trail of zeroes to
                            % channels that did not receive the last data
                            % packet, trim samples from other channels to
                            % match data lengths
                            trim = 256-length(ix);
                            spike_channel(i).data = spike_channel(i).data(1:end-trim);
                            spike_channel(i).time = spike_channel(i).time(1:end-trim);
                            disp(sprintf(' Channel %d, trimmed %d samples from end to match data lengths with channels that did not receive its last data package',i, trim));
                            %pause(5);
%                         else
%                             % The end-of-trial code has been received, but file saving has been promptly turned off by the application, such that the last data packet is not saved any more
%                             ix = (length(spike_channel(i).data):(apmdata.channels(i).end_trial(1) - apmdata.channels(i).start_trial(1))) + 1;
%                             % Grow the data and time vectors as needed
%                             spike_channel(i).data(ix) = 0;
%                             spike_channel(i).time(ix) = spike_channel(i).time(end) + (1:length(ix));
                        end;
                    end;
                %catch
                %    disp('No end of trial found.');
                end;
                % Align time vector with the start of trial and convert to seconds.
                spike_channel(i).time = apmdata.channels(i).time_calibration * 1e-6 * (spike_channel(i).time - spike_channel(i).time(1));
                spike_channel(i).sampfreq = round(1e6/apmdata.channels(i).time_calibration);
                %[apmdata.channels(i).end_trial(1) - apmdata.channels(i).start_trial(1) + 1 length(spike_channel(i).data) length(apmdata.channels(i).continuous)]
            end;
        end;


        % Next, deal with aux channels
        j=0;
        aux_channel = [];
        for i=1:n_aux
            for k=1:length(apmdata.aux(i).input)
                if (~isempty(apmdata.aux(i).input(k).continuous))
                    j=j+1;
                    aux_channel(j).data = apmdata.aux(i).voltage_calibration(1)*(apmdata.aux(i).input(k).continuous-auxADCzero);
                    aux_channel(j).time = apmdata.aux(i).time_calibration * 1e-6 * (apmdata.aux(i).input(k).start_continuous - apmdata.aux(i).start_trial(1) + (1:length(apmdata.aux(i).input(k).continuous)) - 1);
                    aux_channel(j).sampfreq = round(1e6/apmdata.aux(i).time_calibration);
                    % Data is aligned with the gate-on/start-of-trial timestamp, may result in negative numbers as data recording may start before receiving the gate pulse
                    % Also, recording on the aux channels may continue a little bit beyond the point where the spike recording has stopped. Remove that interval as well.
                    ix = find(( aux_channel(j).time < 0 ) | (aux_channel(j).time > spike_channel(1).time(end)));
                    %For pt Diamond 4/13/10, use the line below to
                        % adjust for spike channel 2 recording
                    %ix = find(( aux_channel(j).time < 0 ) | (aux_channel(j).time > spike_channel(2).time(end)));
                    if ~isempty(ix)
                        aux_channel(j).data(ix) = [];
                        aux_channel(j).time(ix) = [];
                    end;
                end;
            end;
        end;
     
        n_aux = j;  % Update the number of channels based on what channels acutally held data
        
        %clear apmdata
        break;
    catch
        disp('Could not load data file, reason:');
        disp(lasterr);
        s=input('Retry (y/n) [y]','s');
    	if (upper(s)=='N')
            error('Canceled by user.');
        end;
    end;
end

for i =1:n_aux
    n_spike =n_spike+1;
    spike_channel(n_spike) = aux_channel(i);
end;

h = figure;
j = 0;
idx = [];
for i = 1:n_spike
    if (~isempty(spike_channel(i).data))
        j = j+1; % count number of spike channels with useful data
        idx = [idx i]; % store spike channel number
    end
end

k = 0;
for i = idx
    k = k+1;
    subplot(j,1,k)
    plot(spike_channel(i).time,spike_channel(i).data);
    title(['Channel #' num2str(i)]);
end

%% Select ecog and stim trigger channel
ecog_ch = input('Enter channel #s for good ECOG (format: [1 2 ...]): ');
trig_ch = input('Enter channel # for stimluation trigger voltage: ');
% ecog_ch = [3 4 5 6 7];
% trig_ch = 14;
%% process ecog channels

n_ecog = length(ecog_ch);
ecog = struct('contact_pair',{},'rest_time',{},'active_time',{}); %SS: ecog structure used with ecogfft script
% ecog data should have equal number of data points for all ecog channels.  Some
% channels occasionally have slightly fewer data points, as a result of not
% receiving the last data packet at the end of a recording.
% To create raw ecog data of uniform length, find the
% ecog channel with the minimum number of ecog data points, and trim
% all other channels to match its length.
len=[];
for i =1:n_ecog
    len = [len length(spike_channel(ecog_ch(i)).data)];
end
min_len = min(len);
for i = 1:n_ecog
    ecog(1).contact_pair(i).raw_ecog_signal...
        = spike_channel(ecog_ch(i)).data(1:min_len);
    ecog(1).contact_pair(i).chan_num=ecog_ch(i);
end

% disp(['Processing ' num2str(n_ecog) ' Ecog channels']);
% ecog_chan = [];
% time = [];
% j = 1;
% for i = ecog_ch
% 
%     if (isempty(time) || (length(time) == length(spike_channel(i).time)))
%         ecog(1).contact_pair(j).raw_ecog_signal = spike_channel(i).data;
%         ecog(1).contact_pair(j).chan_num = i;
% %        ecog_chan(j,:) = spike_channel(i).data;                    % For compatibility with the original script; works only when the sampling rates on all emg channels is the same
% %         ecog_names(j,:) = spike_names(i,:);                        % For compatibility with the original script
%         time = spike_channel(i).time;                               % For compatibility with the original script
%         j = j + 1;
%     else
%         disp(sprintf('WARNING: ECOG channel %d does not have the same sampling frequency like the other channels! Channel will be skipped!',i));
%     end;
% end


%% Process stim trigger voltage channel

% find stim trig timestamps using getevts
time = spike_channel(trig_ch).time;
trig_chan = spike_channel(trig_ch).data;
MARGIN = 0.35;  % 0.35 is a very wide range that handles large
% fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
% STIM_TRIG_OFF =  [0-MARGIN 0+MARGIN];
STIM_TRIG_OFF =  [4.06-MARGIN 4.06+MARGIN]; % fowler12SSEP, Danchuk01SSEP
STIM_TRIG_ON = [STIM_TRIG_AMP-MARGIN STIM_TRIG_AMP+MARGIN];
% STIM_TRIG_OFF =  [4.06-MARGIN 4.06+MARGIN]; % trig output voltage from 8-29-2008 troubleshoot
[stim_trig_time,stim_trig_val] = getevtsSSEP(time, trig_chan, STIM_TRIG_OFF,STIM_TRIG_ON);
% % stim_trig_time = [1 1.5 2 2.5 3 3.5 4];
% for adamsRSSEP
% stim_trig_time=[10.5610000000000,11.0111500000000,11.4612500000000,11.9114000000000,12.3615500000000,12.8118500000000,13.2620000000000,13.7121000000000,14.1622500000000,14.6124000000000,15.0625500000000,15.5127000000000,15.9628000000000,16.4129500000000,16.8631000000000,17.3132500000000,17.7633500000000,18.2135000000000,18.6636500000000,19.1137500000000,19.5639000000000,20.0140500000000,20.4641500000000,20.9143000000000,21.3644000000000,21.8145500000000,22.2647000000000,22.7148000000000,23.1649500000000,23.6151000000000,24.0652500000000,24.9655500000000,25.4157000000000,25.8658000000000,26.3159500000000,26.7661000000000,27.2163000000000,27.6664000000000,28.1165500000000,28.5667000000000,29.0168500000000,29.4669500000000,29.9171000000000,30.3672500000000,30.8173500000000,31.2675000000000,31.7176500000000,32.1678000000000,32.6179500000000,33.0680500000000,33.5182000000000,33.9683500000000,34.4185000000000,34.8686500000000,35.3188000000000,35.7689000000000,36.2190500000000,36.6692000000000,37.1193500000000,37.5694500000000,38.0196000000000,38.4697500000000,38.9199000000000,39.3700500000000,39.8201500000000,40.2703000000000,40.7203000000000,41.1704000000000,41.6206000000000,42.0707000000000,42.5210000000000,42.9711500000000,43.4212500000000,43.8714000000000,44.3215500000000,44.7717000000000,45.2218000000000,45.6718500000000,46.1224500000000,46.5725500000000,47.0229000000000,47.4728500000000,47.9232000000000,48.3733500000000,48.8235000000000,49.2736000000000,49.7237500000000,50.1739000000000,50.6240000000000,51.0741500000000,51.5242500000000,51.9744000000000,52.4245000000000,52.8746000000000,53.3247500000000,53.7748500000000,54.2250000000000,57.8256500000000,58.2758000000000,58.7259500000000,59.1761000000000,59.6262000000000,60.0769500000000,60.5271000000000,60.9772500000000,61.4273500000000,61.8775000000000,62.3276500000000,62.7777500000000,63.2279000000000,63.6780500000000,64.1282000000000,64.5783000000000,65.0284500000000,65.4786000000000,65.9287500000000,66.3789000000000,66.8290000000000,67.2791500000000,67.7293000000000,68.1794500000000,68.6296000000000,69.0797500000000,69.5298500000000,69.9800000000000,70.4301500000000,70.8803000000000,71.3304000000000,71.7805500000000,72.2307000000000,72.6808500000000,73.1309500000000,73.5811000000000,74.0312500000000,74.4813500000000,74.9315000000000,75.3816500000000,75.8318000000000,76.2819000000000,76.7320500000000,77.1822000000000,77.6323500000000,78.0824500000000,78.5326000000000,78.9827500000000,79.4328500000000,79.8830000000000,80.3331500000000,80.7833000000000,81.2334500000000,81.6835500000000,82.1337500000000,82.5838500000000,83.0340000000000,83.4841000000000,83.9343000000000,84.3844000000000,84.8345500000000,85.2847000000000,85.7348000000000,86.1849500000000,86.6351000000000,87.0852500000000,87.5353500000000,87.9855000000000,88.4356000000000,88.8857500000000,89.3359000000000;];
%% Peri-stim analysis

% ecog.contact_pair.raw_ecog_signal contains raw ecog data
% stim_trig_time contains stimulus trigger timestamps
% find a way to import both using ReadEventData

num_contact_pair = length(ecog.contact_pair);
stim_trig_time = stim_trig_time(stim_trig_time > abs(PRE_STIM)); % discard stim trig times at the beginning of recording that do not have adequate pre-stim period
num_stim = length(stim_trig_time);
t = 1000*(PRE_STIM:1/Fs:PST_STIM); % multiplied by 1000 to change to msec scale

% initialize
stim_epoch = zeros(num_stim,length(t),num_contact_pair);
diff_stim_epoch = zeros(num_stim,length(t),num_contact_pair);

for i = 1:num_contact_pair
    
    % parse each stim epoch from each contact pair
    for j = 1:num_stim
        tmp1 = int32((stim_trig_time(j)+PRE_STIM) * Fs);  % int32 used to keep index in integer format
        tmp2 = int32((stim_trig_time(j)+PST_STIM) * Fs);
        % since int32 rounds to the next closest integer, there may be some
        % cases where the parsing indeces are off by 1.  fix them on the fly.
        d = tmp2 - tmp1; 
        % lengt(t)-d must equal 1 for the parsing to fit stim_epoch array
        if length(t)-d == 0 
            tmp2 = tmp2-1;
        elseif length(t)-d == 2
            tmp2 = tmp2+1;
        end
        stim_epoch(j,:,i) = ecog.contact_pair(i).raw_ecog_signal(tmp1:tmp2);
    end
end

% ------------Plot 1, referenced to contact 6----------------
if plotfig1
    mean_stim_epoch = mean(stim_epoch);

    hf1 = figure;

    % for i = 1:num_contact_pair
    %     subplot(num_contact_pair,1,i);
    %     plot(t,mean_stim_epoch(1,:,i));
    % %     set(gca, 'ylim', [-20 20]);
    %     set(gca, 'ylim', [-1 1]);
    %     ylm = ylim;
    %     hold on; plot([0 0],[ylm(1) ylm(2)],'k--'); hold off;
    %     title(['Ecog contact #' num2str(ecog.contact_pair(i).chan_num-2) '- #6']);
    %     ylabel('Voltage (uV)');
    %     if i == num_contact_pair
    %         xlabel('Time (msec)');
    %     end
    % end

    hold on
    t_text = PRE_STIM*800;  % specify where along the time axis to place ecog pair text
    t_ind = find(t == t_text);  % find index at t_text
    for i = 1:num_contact_pair
        SSEPmin = min(mean_stim_epoch(1,:,i)); % max value of average
        SSEPmax = max(mean_stim_epoch(1,:,i)); % min value of average
       C = num_contact_pair - i; % constant added to stack the waves for comparison such that the first contact pair is plotted at the top
         z = (mean_stim_epoch(1,:,i)-SSEPmin)/(SSEPmax-SSEPmin) + C; % normalization step
%          z = (mean_stim_epoch(1,:,i))+C;
        plot(t,z);
        if ecog.contact_pair(i).chan_num == 8 %SSEP recorded in Ch.8 is typically Erb's Point
            text(t_text,z(t_ind)+0.25,'EP');
        else
            text(t_text,z(t_ind)+0.25,...
                ['e' num2str(ecog.contact_pair(i).chan_num - 2) '-6']);
        end
    end
    xlabel('Time (msec)');
    ylabel('normalized SSEPs');
    title(strrep(fn,'.apm',''));
    ylm = ylim;
    plot([0 0],[ylm(1) ylm(2)],'k--')
    hold off
end


%------------plot 2, montaged w/ adjacent contact pairs---------------------

for i = 1:num_contact_pair
%     % parse each stim epoch from each contact pair
%     for j = 1:num_stim
%         tmp1 = int32((stim_trig_time(j)+PRE_STIM) * Fs);  % int32 used to keep index in integer format
%         tmp2 = int32((stim_trig_time(j)+PST_STIM) * Fs);
%         stim_epoch(j,:,i) = ecog.contact_pair(i).raw_ecog_signal(tmp1:tmp2);
%     end
    % find difference between adjacent contact pairs
    if i == 1
        continue;
    else
        % pt Solis ecog data, input ecog signal and reference contacts were
% accidentally flipped.  Correct for flipped ecog signal amplitude.
        %diff_stim_epoch(:,:,i-1) = stim_epoch(:,:,i) - stim_epoch(:,:,i-1);
        diff_stim_epoch(:,:,i-1) = stim_epoch(:,:,i-1) - stim_epoch(:,:,i);        
        if i == num_contact_pair
            % pt Solis ecog data, input ecog signal and reference contacts were
% accidentally flipped.  Correct for flipped ecog signal amplitude.
            %diff_stim_epoch(:,:,i) = -stim_epoch(:,:,i);
            diff_stim_epoch(:,:,i) = stim_epoch(:,:,i);
            
        end
    end
end

mean_diff_stim_epoch = mean(diff_stim_epoch);

hf2 = figure;

% for i = 1:num_contact_pair
%     subplot(num_contact_pair,1,i);
%     plot(t,mean_diff_stim_epoch(1,:,i));
% %     set(gca, 'ylim', [-20 20]);
%     set(gca, 'ylim', [-30 30]);
%     ylm = ylim;
%     hold on; plot([0 0],[ylm(1) ylm(2)],'k--'); hold off;
%     title(['Ecog contact #' num2str(ecog.contact_pair(i).chan_num -2)...
%         ' - #' num2str(ecog.contact_pair(i).chan_num -1)]);
%     ylabel('Voltage (uV)');
%     if i == num_contact_pair
%         xlabel('Time (msec)');
%     end
%     
% end

hold on
t_text = PRE_STIM*800;  % specify where along the time axis to place ecog pair text
t_ind = find(t == t_text);  % find index at t_text

for i = 1:num_contact_pair
    SSEPmin = min(mean_diff_stim_epoch(1,:,i)); % max value of average
    SSEPmax = max(mean_diff_stim_epoch(1,:,i)); % min value of average
    
%     % normalized SSEP
    C = num_contact_pair - i; % constant added to stack the waves for comparison such that the first contact pair is plotted at the top
    z = (mean_diff_stim_epoch(1,:,i)-SSEPmin)/(SSEPmax-SSEPmin) + C; 
    % non-normalized SSEP
%     C = (SSEPmax-SSEPmin)*(num_contact_pair-i);    
%     z = mean_diff_stim_epoch(1,:,i) + C; % non-normalized SSEP

    plot(t,z);
    if ecog.contact_pair(i).chan_num == 8 %SSEP recorded in Ch.8 is typically Erb's Point
        text(t_text,z(t_ind)+0.1,'EP');
    else
        text(t_text,z(t_ind)+0.1,...
        ['e' num2str(ecog.contact_pair(i).chan_num - 2)...
        '-' num2str(ecog.contact_pair(i).chan_num-1)]);
    end
end
xlabel('Time (msec)');
ylabel('normalized SSEPs');
title(strrep(fn,'.apm',''));
ylm = ylim;
plot([0 0],[ylm(1) ylm(2)],'k--')
hold off

