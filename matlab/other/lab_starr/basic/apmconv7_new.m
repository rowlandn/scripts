% function apmconv7
% apmconv7()
% This function imports an .apm file collected by Guideline 4000 system
%   
% Input:
%       1) .apm file
% Outputs:
%       1) "fname.ddt" = a ddt file which is native to Plexon, for Offline Sorter
%       2) "fname_a.ddt" = an OPTIONAL ddt file containing accelerometer
%       for accel deflection "sorting"
%       3) "fname.mat" = EMG data
%       4) "fname_a.mat" = accel data
%       5) "fname_ecog.mat" = OPTIONAL ecog data
%
%
% RLM: created by Rob Turner, last edited 2003.10.03 by Rory Michaelis
%
% AB: modified to use APMReadData.m Matlab script for loading GL4K data in APM format
%     Note: Accelerometer data is not saved any more along with the EMG
%     data in a .mat file, since the sampling frequencies may be different
%     of the two type of channels. The accel data is saved in a separate
%     .mat file with the same suffix as the accel DDT file.
% - Andrei Barborica, FHC Inc, 2007
% 
% SS (11-18-2008): modified so that user has the option to input channel #s
% for spike, task-voltage, EMG, and accel channels if output files for such data
% are desired.  Otherwise, the input can be left blank to avoid creating
% undesired output files.

% ALC 6/19/09: modified such that all rest timestamps are preserved,
% including the last timestamp, so that the number of rest timestamps
% should always be one more than the number of active timestamps. This is
% done so that the coherence code can use the final rest timestamp in order
% to determine the end point of the final active epoch.
clear all
close all
colordef(0,'white');

% RLM: ##################################################################
rename_variables = 0; % RLM: TOGGLE prompt for manually naming emg channels
write_accel_ddt = 1;  % RLM: TOGGLE writing of LAST analog channel as accelerometer DDT
ecog_channel_append = '_ecog'; %  SS: default name appendage for ecog channel 
accel_channel_append = '_Accel'; % RLM: default name appendage for Accelerometer channel (in *.MAT)
accel_ddt_append = '_a'; % RLM: name appendage for accelerometer DDT (no extension!!!)
ipad_channel_append = '_Ipad'; % RLM: default name appendage for Accelerometer channel (in *.MAT)
ipad_ddt_append = '_i'; % RLM: name appendage for accelerometer DDT (no extension!!!)
% RLM: ##################################################################

% AB: additional paths
%addpath('../../../Guideline4000/MatlabScripts/DataImport/');   % Need to
%add APMDataImport.m script to the path

% AB: GL4K aux ADC values are usigned. When the AUX inputs are configured
% as bipolar, Need to substract midpoint, in order to obtain a bipolar value.
% When configures as unipolar, do not substract anything
%auxADCzero = 8192;  % For bipolar input configuration
auxADCzero = 0;   % For unipolar input configuration

% fn_dir = strrep(which('abfconv2'),'abfconv2.m','');   % Store function location for later use

%% import and process GL4k channels

while true

    %d='../../../Data/'; % AB: test dir
    %fn='Wav11.apm'; % AB: test file name
    %fn='TESTUCSF.apm'; % AB: test file name
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
        episode=1;  % This is for ABF file format conpatibility
        
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

% Include all aux channels below spike channels, since they may hold good EMG data, recorded directly using an electrode (not through a preamp)
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

%% select, process, and save spike channel
% spike_ch = 0; %initialize
% while isempty(spike_ch) || min(spike_ch) < 1 || max(spike_ch) > n_spike
%         spike_ch = input('Enter channel #s for spike data, otherwise leave blank (format: [1 2 ...]): ');
% %     spike_ch = 1;
% end
spike_ch = input('Enter channel #s for spike data, otherwise leave blank (format: [1 2 ...]): ');

if ~isempty(spike_ch)
    for i=spike_ch
        if (length(spike_ch)>1)
            outputname = sprintf('%s%d',vn,i); % Use var name + channel number
        else
            outputname = vn;    % Use var name
        end;
    	disp(['Writing SPIKE channel to:  ' outputname]);
    	ddtdata = (32767/(max(abs(spike_channel(i).data))))*spike_channel(i).data; % +/-2047 value is from DDT format 12bit res.    
        writeddt(ddtdata,spike_channel(i).sampfreq,outputname,16,'Creating ddt file.');
       
    end
else
    outputname=vn;
end

Fs = spike_channel(1).sampfreq;

%% select EcCOG/LFP/EMG/ACCEL channels
ecog_ch = input('Enter channel #s for ECOG and LFP (format: [3 4...8]): ');
evt_ch = input('Enter channel # for ECOG event voltage (usually 14), otherwise leave blank: ');
emg_ch = input('Enter channel #s for desired EMG, otherwise leave blank (format: [9 10 ...]): ');
ipad_ch = input('Enter channel # for ipad event(usually 11), otherwise leave blank: ');
accel_ch = input('Enter channel # for desired ACCEL, otherwise leave blank: ');
M1_ch = input('Enter channel # for M1 contact (usually 4 or 5): '); % added 6/19/08 to be carried over for later programs


%% rename aux channels if desired
% RLM: ##################################################################
% NAME CHANNELS, 64 characters per channel as per NEX channel variables.

for i = [ecog_ch emg_ch accel_ch]
    if i == ecog_ch
        varname=[num2str(i) ecog_channel_append];
    elseif i==accel_ch % make "n_emg" or 5 ?
        varname=[num2str(i) accel_channel_append];
    else
        varname=num2str(i);
    end
    if rename_variables
        % RLM: perhaps draw small channel plot to help in naming?
        disp(sprintf('\nNON-Spike Channel Number:  %d  is named:\n\t''%s''\n', i, varname));
        question='WOULD YOU LIKE TO CHANGE THIS ?\nIf so, type "y" and press ENTER. If not, just press ENTER\n # > ';
        choice=input(question, 's');
        while ~strcmpi(choice, '')
            if strcmpi(choice, 'y')
                varname = input('\nType new VARIABLE NAME now. REMEMBER, ONLY 64 CHARACTERS WILL BE WRITTEN!\n # > ', 's');
                choice = '';
            else
                choice=input(question, 's');
            end
        end
    end
    varname=deblank(varname);
    varnamesize=length(varname);
    if varnamesize < 64
        empty = 64 - varnamesize;
        empty = char(zeros(1, empty));
        varname = char([varname empty]);
    else
        varname = varname(1:64);
    end
    spike_names(i,:) = varname;
end
% RLM: ##################################################################

%% process ecog/LFP channels
n_ecog = length(ecog_ch);
ecog = struct('contact_pair',{},'rest_time',{},'active_time',{}); %SS: ecog structure used with ecogPSD script

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
minlen = min(len);
for i = 1:n_ecog
    ecog(1).contact_pair(i).raw_ecog_signal...
        = spike_channel(ecog_ch(i)).data(1:minlen); % cut to match lengths
    ecog(1).contact_pair(i).chan_num=ecog_ch(i);
end


%% parse rest/active epoch timestamps
if ~isempty(evt_ch)
    % parse ecog rest/active timestamps using getevts
    time = spike_channel(evt_ch).time;
    evt_chan = spike_channel(evt_ch).data;
    MARGIN=0.1;
%After adjusting the aux channels, a std button press voltage is as
%follows:  
REST = [6.54-MARGIN 6.54+MARGIN];
ACTIVE =  [6.17-MARGIN 6.17+MARGIN];

    [rest_time,rest_val] = getevts(time, evt_chan, REST, ACTIVE); % getevts expects row vector time inputs and column vector voltage input
    [active_time,active_val] = getevts(time,evt_chan,ACTIVE,REST);
    

    try
        % Throw out undesirably close timestamps (< 0.1 sec) within rest_time and active_time.
        % This results from large voltage fluctuations that sometimes happens when task button is pressed or unpressed
        idx = find(diff(rest_time)>= 0.1) + 1; % offset of 1 is added for indexing later
        idx = [1 idx]; % add 1st element back in
        rest_time = rest_time(idx);
        idx = find(diff(active_time)>0.1) + 1;
        idx = [1 idx];
     
        active_time = active_time(idx);

        % force ecog epochs to begin with rest
        if rest_time(1) > active_time(1)
            rest_time = [0 rest_time]; % add time '0' as the first rest epoch timestamp
        end
        % force ecog epochs to end with active
        if rest_time(end) > active_time(end)
%             rest_time = rest_time(1:end-1); % remove last rest epoch timestamp
            rest_tmp = rest_time(end); % store the last rest ts in rest_tmp for later
            rest_time = rest_time(1:end-1); % remove last rest epoch timestamp
        end
        % throw out rest/active epoch pair timestamps that are too close.  for
        % ecogPSD analysis, timestamps must be greater than 2.5 sec apart
        if length(rest_time) == length(active_time)
            idx = find((active_time - rest_time)>=2.5);
%             idx = find((active_time - rest_time)>2);
            rest_time = rest_time(idx);
            active_time = active_time(idx);
        else
            error('Rest/active epochs must be paired.  Check task-voltage channel');
        end
        
        % Add the last rest time back in because it's necessary for
        % Andrea's data analysis
        if exist('rest_tmp')
            rest_time(end+1) = rest_tmp;
        end
        
        % display warning if there are < 5 rest/active epoch pairs
        if length(rest_time) < 5
            warning('There are less than 5 rest/active epoch pairs!');
        end

        ecog(1).rest_time = rest_time;
        ecog(1).active_time = active_time;
    catch
        error('Check ecog-event voltage channel and make sure the right channel was selected');
    end
end


%% save ecog data
outputname = [outputname '.mat'];
disp(['Writing ECOG  channels to:  ' outputname(1:end-4) ecog_channel_append '.mat']);
% Include time vector if imported
save( [outputname(1:end-4) ecog_channel_append '.mat'], 'ecog','n_ecog','M1_ch');


%% process and save emg channel
if ~isempty(emg_ch)
    n_emg = length(emg_ch);
    disp(['Processing ' num2str(n_emg) ' EMG channels']);
    len = [];
    for i = emg_ch
        len = [len length(spike_channel(i).data)];
    end
    minlen = min(len);
    emg_chan = zeros(n_emg,minlen);
    j=1;
    for i = emg_ch
        tmp = spike_channel(i).data;
        tmp = tmp(1:minlen);    % cut emg data to match minimum length
        tmp = detrend(tmp);     % Remove offset & linear trend
%         tmp = abs(tmp);         % Computer absolute value
        tmp = tmp;         % Computer absolute value

        emg_chan(j,:) = tmp;
        emg_names(j,:) = spike_names(i,:);
        j = j + 1;
    end
    time = spike_channel(emg_ch(1)).time(1:minlen); %pull out time vector, cut to match minimum length
    disp(['Writing EMG  channels to:  ' outputname]);
    save( outputname, 'n_emg', 'emg_chan', 'time', 'emg_names'); % RLM: added channel names
end

%% process and save accel channel
if ~isempty(accel_ch)
    n_accel = length(accel_ch);
    disp(['Processing ' num2str(n_accel) ' Accel channels']);
    accel_chan = [];
    time = [];
    j = 1;
    for i = accel_ch
        if (isempty(time) || (length(time) == length(spike_channel(i).time)))
            accel_chan(j,:) = spike_channel(i).data;                    % For compatibility with the original script; works only when the sampling rates on all emg channels is the same
            accel_names(j,:) = spike_names(i,:);                        % For compatibility with the original script
            time = spike_channel(i).time;                               % For compatibility with the original script
            j = j + 1;
        else
            disp(sprintf('WARNING: ACCEL channel %d does not have the same sampling frequency like the other channels! Channel will be skipped!',i));
        end;
    end

    disp(['Writing ACCEL  channels to:  ' outputname(1:end-4) accel_ddt_append '.mat']);
    % Include time vector if imported
    save( [outputname(1:end-4) accel_ddt_append '.mat'], 'n_accel', 'accel_chan', 'time', 'accel_names'); % RLM: added channel names

    % RLM: ##################################################################
    % AB: As originally written, the code is going to fail if more than one accel channel is present
    if write_accel_ddt
        ddtdata = 32767/max(abs(spike_channel(accel_ch).data)) * spike_channel(accel_ch).data; % fill 16 bit integral spectrum
        writeddt(ddtdata,spike_channel(accel_ch).sampfreq,[outputname(1:end-4) accel_ddt_append],16,'ACCELEROMETER DDT created by abfconv4.m');
    end
    % RLM: ##################################################################
end


%% process and save ipad channel
if ~isempty(ipad_ch)
    n_ipad = length(ipad_ch);
    disp(['Processing ' num2str(n_ipad) ' Ipad channels']);
    ipad_chan = [];
    time = [];
    j = 1;
    for i = ipad_ch
        if (isempty(time) || (length(time) == length(spike_channel(i).time)))
            ipad_chan(j,:) = spike_channel(i).data;                    % For compatibility with the original script; works only when the sampling rates on all emg channels is the same
            ipad_names(j,:) = spike_names(i,:);                        % For compatibility with the original script
            time = spike_channel(i).time;                               % For compatibility with the original script
            j = j + 1;
        else
            disp(sprintf('WARNING: Ipad channel %d does not have the same sampling frequency like the other channels! Channel will be skipped!',i));
        end;
    end

    disp(['Writing Ipad  channels to:  ' outputname(1:end-4) ipad_ddt_append '.mat']);
    % Include time vector if imported
    save( [outputname(1:end-4) ipad_ddt_append '.mat'], 'n_ipad', 'ipad_chan', 'time', 'ipad_names'); % RLM: added channel names

end
Fs = spike_channel(3).sampfreq;
%%
save( [outputname(1:end-4) 'ecog_all.mat'], 'ecog','n_ecog','M1_ch','emg_chan','accel_chan','ipad_chan', 'time','Fs');
clear all
