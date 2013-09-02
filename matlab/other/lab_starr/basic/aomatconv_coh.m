function aomatconv_coh
% converts mat file from AlphaOmega Mapfile converter program to a format
% that can be run on existing code for ecogPSD analysis of GL4k data

% Created by S.Shimamoto (12/15/2008)

%% load data
%load .mat file 
[filename pathname]=uigetfile('*.mat','Select .mat file');
cd(pathname);
load(filename);
%% Downsampling
% for pt John, monopolar ecog channels were recorded at 6K instead of 1.5K.
% CECOG1 = downsample(CECOG1,4);
% CECOG2 = downsample(CECOG2,4);
% CECOG3 = downsample(CECOG3,4);
% CECOG4 = downsample(CECOG4,4);
% CECOG5 = downsample(CECOG5,4);

%% make sure sampling frequencies are equal for ecog/LFP/task recordings, then
% use CECOG1_KHz varible as the sampling rate
% Comment out if performing downsampling above
if isequal(CECOG1_KHz,CECOG2_KHz,CECOG3_KHz,CECOG4_KHz,CECOG5_KHz,CLFP1_KHz)
    Fs = CECOG1_KHz;
    Fs = Fs*1000;                           % multiply by 1000 for KHz->Hz conversion
else
    error('Sampling frequencies are different for ecog/LFP/task button data');
end

%% OBSOLETE - snip artificial zeros from the beginning and end of each variable
% % ecog/lfp/task recordings all have a string of artificial zeros in the
% % beginning and end of recording.  snip them out.
% inda=find(CAN_IN1);
% inde1=find(CECOG1);
% inde2=find(CECOG2);
% inde3=find(CECOG3);
% inde4=find(CECOG4);
% inde5=find(CECOG5);
% indl=find(CLFP1);
% % use indeces to pull out data segment between the first and
% % last non-zero data points
% CAN_IN1 = CAN_IN1(inda(1):inda(end));
% CECOG1  = CECOG1(inde1(1):inde1(end));
% CECOG2  = CECOG2(inde2(1):inde2(end));
% CECOG3  = CECOG3(inde3(1):inde3(end));
% CECOG4  = CECOG4(inde4(1):inde4(end));
% CECOG5  = CECOG5(inde5(1):inde5(end));
% CLFP1   = CLFP1(indl(1):indl(end));

%% force ecog/lfp data to be the same lengths for montage-ing

minlength=min([length(CECOG1) length(CECOG2) length(CECOG3)...
    length(CECOG4) length(CECOG5) length(CLFP1)]);
CECOG1  = CECOG1(1:minlength);
CECOG2  = CECOG2(1:minlength);
CECOG3  = CECOG3(1:minlength);
CECOG4  = CECOG4(1:minlength);
CECOG5  = CECOG5(1:minlength);
CLFP1   = CLFP1(1:minlength);


%% voltage calibration
C1 = 5/32768; % constant for calibration
% value 32768 is equal to 5V in the Alpha Omega system
C2 = 1e6; % conversion to microvolts
Gecog = 7000;   % ecog channel gain
Glfp = 25000;   % LFP channel gain
CAN_IN1 = CAN_IN1*C1;
CAN_IN2 = CAN_IN2*C1;
CAN_IN3 = CAN_IN3*C1;
CECOG1  = CECOG1*C1*C2/Gecog;
CECOG2  = CECOG2*C1*C2/Gecog;
CECOG3  = CECOG3*C1*C2/Gecog;
CECOG4  = CECOG4*C1*C2/Gecog;
CECOG5  = CECOG5*C1*C2/Gecog;
CLFP1   = CLFP1*C1*C2/Glfp;
% Make row vectors into column vectors
CECOG1 = CECOG1';
CECOG2 = CECOG2';
CECOG3 = CECOG3';
CECOG4 = CECOG4';
CECOG5 = CECOG5';
CLFP1 = CLFP1';
%Create matrix of raw ecog/lfp data
ecog_lfp_raw_data = [CECOG1 CECOG2 CECOG3 CECOG4 CECOG5 CLFP1];

%% Define auxiliary channels
aux1 = CAN_IN1;
aux2 = CAN_IN2;
aux3 = CAN_IN3;

%% select ECOG/LFP/EMG/ACCEL channels
% ecog_ch = input('Enter channel #s for good ECOG (format: [1 2 ...]): ');
% lfp_ch = input('Enter channel #s for good LFP (format: [1 2 ...]): ');
evt_ch = input('Enter channel # for ECOG event voltage (aux1, aux2, or aux3), otherwise leave blank: ');
% evt_ch = input('Enter channel # for ECOG event voltage: ');
% emg_ch = input('Enter channel #s for good EMG (format: [1 2 ...]): ');
% accel_ch = input('Enter channel # for ACCEL: ');

%% grab rest/active timestamps from task button
%Remember that evt_ch has already been defined as the data from the user
%chosen auxiliary channel

if ~isempty(evt_ch)
   
    nsamples = length(evt_ch);
    T = 1/(1000*CAN_IN3_KHz);  % period (sec) between samples; specifing CAN_IN3 here b/c this channel should be used for task box
% create time vector
time = linspace(0,T*nsamples,nsamples);
% time = linspace(0,T*nsamples,nsamples);
MARGIN = 0.35;  % 0.3 is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
%Old Button box values (UCSF)
% REST = [5-MARGIN 5+MARGIN];
% REST = [5.25-MARGIN 5.25+MARGIN];
%     REST = [1.69-MARGIN 1.69+MARGIN];
% ACTIVE =  [4.36-MARGIN 4.36+MARGIN];
% ACTIVE =  [4.46-MARGIN 4.46+MARGIN];
%     ACTIVE =  [5.25-MARGIN 5.25+MARGIN];
% ACTIVE =  [3.60-MARGIN 3.60+MARGIN];
 
% New button values (VA)
REST = [0.44-MARGIN 0.44+MARGIN];
ACTIVE = [5-MARGIN 5+MARGIN];

%make evt_ch and time column vectors
evt_ch = evt_ch';
time = time';

[rest_time,rest_val] = getevts(time,evt_ch,REST,ACTIVE); % getevts expects column vector inputs
[active_time,active_val] = getevts(time,evt_ch,ACTIVE,REST);

    try
        % Throw out undesirably close timestamps (< 0.1 sec) within rest_time and active_time.  
        % This results from large voltage fluctuations that sometimes happens when task button is pressed or unpressed 
        idx = find(diff(rest_time)> 0.3) + 1; % offset of 1 is added for indexing later
        idx = [1; idx]; % add 1st element back in
        rest_time = rest_time(idx);
        idx = find(diff(active_time)>0.3) + 1;
        idx = [1; idx];
        active_time = active_time(idx);
    
        % force ecog epochs to begin with rest
        if rest_time(1) > active_time(1)
            rest_time = [0; rest_time]; % add time '0' as the first rest epoch timestamp
        end
          % 9/10/08: no longer removing last rest epoch since last rest timestamps can
          % be used to calculate length of last active epoch
        % force ecog epochs to end with active
        if rest_time(end) > active_time(end)
            rest_tmp = rest_time(end); % store the last rest ts in rest_tmp for later
            rest_time = rest_time(1:end-1); % remove last rest epoch timestamp
        end

        % throw out rest/active epoch pair timestamps that are too close.  for
        % ecogfft analysis, timestamps must be greater than 3 sec apart
        if length(rest_time) == length(active_time)
          idx = find((active_time - rest_time)>=2.5);
%             idx = find((active_time - rest_time)>1.0); % 10/30/08: some of our epochs are just under 3sec. Relaxing time constraints to maximize coh data.
            rest_time = rest_time(idx);
            active_time = active_time(idx);
        else
            error('Rest/active epochs must be paired.  Check task-voltage channel');
        end

        if exist('rest_tmp')
            rest_time(end+1) = rest_tmp;
        end
    
        % display warning if there are < 5 rest/active epoch pairs
        if length(rest_time) < 5
            warning('There are less than 5 rest/active epoch pairs!');
        end

    %     ecog(1).rest_time = rest_time;  % not using the ecog structure array
    %     for apmconv7_coh
    %     ecog(1).active_time = active_time;
    catch
        error('Check ecog-event voltage channel and make sure the right channel was selected');
    end

% add rest_time and active_time to ecog_lfp_raw_data
num_row = size(ecog_lfp_raw_data,1);
rest_active_time = zeros(num_row,2); % populate 2 columns with zeros
rest_active_time(1:length(rest_time),1) = rest_time; % fill in rest time
rest_active_time(1:length(active_time),2) = active_time; % fill in active time
ecog_lfp_raw_data = [ecog_lfp_raw_data rest_active_time]; % concatenate
end

% outputname = [outputname '.mat'];
ecog_channel_append = '_ecog'; %  SS: default name appendage for ecog channel 
lfp_channel_append = '_lfp';

outputname = strrep(filename, '.mat','');
disp(['Writing ECOG/LFP  channels to:  ' outputname ecog_channel_append lfp_channel_append '.mat']);
% Include time vector if imported
% save( [outputname(1:end-4) ecog_channel_append '.mat'], 'ecog', 'n_ecog', 'ecog_chan', 'time', 'ecog_names'); % RLM: added channel names
save( [outputname ecog_channel_append lfp_channel_append '.mat'], 'ecog_lfp_raw_data');

%
% %% determine movement onset
% 
% typedata = menu('Detect movement onset using EMG/accel/task button?','EMG','accel','task button');
% 
% % determine movement onset time
% if typedata~=3
%     if typedata==1
%         % insert variable containing emg data
%         % move = ???;
%         % figure out how to calculate time vector
%         % time = ???;
%     elseif typedata==2
%         % insert variable containing accel data
%         % move = ???;
%         % figure out how to calculate time vector
%         % time = ???;
%     end
%     move_onset = DetectEMG(time,move,active_time);
%     move_onset = move_onset(~isnan(move_onset));
% else
%     move_onset = active_time;
% end

% %% Sho's aomatconv code creating a structure
% 
% % ecog = struct('contact_pair',{},'rest_time',{},'active_time',{});
% 
% ecog(1).contact_pair(1).raw_ecog_signal=CECOG1;
% ecog(1).contact_pair(2).raw_ecog_signal=CECOG2;
% ecog(1).contact_pair(3).raw_ecog_signal=CECOG3;
% ecog(1).contact_pair(4).raw_ecog_signal=CECOG4;
% ecog(1).contact_pair(5).raw_ecog_signal=CECOG5;
% ecog(1).contact_pair(6).raw_ecog_signal=CLFP1;
% 
% % ecog(1).rest_time = [0,9.94024930561930,22.6017244347558,35.7698882013466,49.1690154585078;];
% % ecog(1).active_time = [4.55196092672663,15.4029185780890,28.3325476445843,41.9810256154484,54.4593773142011;];
% 
% % use filename to create output name
% outputname=strrep(filename,'.mat','_ecog.mat');
% save(outputname, 'ecog');


