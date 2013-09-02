function aomatconv
% converts mat file from AlphaOmega Mapfile converter program to a format
% that can be run on existing code for ecogPSD analysis of GL4k data

% Created by SAS (12/15/2008)
% Edited by SAS (5/29/09) to process EMG/accel data.
% 

% SAS 5/29/09: aomatconv does NOT know how to process spike channels. For spike
% sorting single-unit MER data, use the .plx file created by Alpha Omega.

% SAS 5/29/09: aomatconv does NOT process EMG or accel data

% SAS 9/9/09: aomatconv updated to ask user-input for contact over M1
% channel a la apmconv7.m

% Input:
%       1) "fname.mat" = output file form Mapfileconverter
% Output:
%       1) "fname_ecog.mat" = ecog/LFP data


%% load data
%load .mat file 
[filename pathname]=uigetfile('*.mat','Select .mat file');
cd(pathname);
load(filename);

% % for Johnson monopolar, ecog channels were accidentally recorded with Fs=6kHz, LFP w/
% % Fs=1.5kHz.  Downsample ecog channel from 6kHZ->1.5kHz
% CECOG1=downsample(CECOG1,4);
% CECOG2=downsample(CECOG2,4);
% CECOG3=downsample(CECOG3,4);
% CECOG4=downsample(CECOG4,4);
% CECOG5=downsample(CECOG5,4); 
% CECOG1_KHz = CLFP1_KHz;
% CECOG2_KHz = CLFP1_KHz;
% CECOG3_KHz = CLFP1_KHz;
% CECOG4_KHz = CLFP1_KHz;
% CECOG5_KHz = CLFP1_KHz;

% make sure sampling frequencies are equal for ecog/LFP/task recordings, then
% use CECOG1_KHz varible as the sampling rate
if isequal(CECOG1_KHz,CECOG2_KHz,CECOG3_KHz,CECOG4_KHz,CECOG5_KHz, CLFP1_KHz)
    Fs = CECOG1_KHz;
    Fs = Fs*1000;                           % multiply by 1000 for KHz->Hz conversion
else
    error('Sampling frequencies are different for ecog/LFP data');
end

%% M1 contact
M1_ch = input('Enter channel # for M1 contact (usually 4 or 5): '); % added 9/9/09 to be carried over for later programs

%% force ecog/lfp data to be the same lengths for montage-ing

minlength=min([length(CECOG1) length(CECOG2) length(CECOG3)...
    length(CECOG4) length(CECOG5) length(CLFP1)]);
CECOG1  = CECOG1(1:minlength);
CECOG2  = CECOG2(1:minlength);
CECOG3  = CECOG3(1:minlength);
CECOG4  = CECOG4(1:minlength);
CECOG5  = CECOG5(1:minlength);
CLFP1   = CLFP1(1:minlength);

%% ecog/LFP voltage calibration
C1 = 5/32768; % constant for calibration
% note: value 32768 is equal to 5V in the Alpha Omega system
C2 = 1e6;   % constant to convert ecog/lfp channel voltage level from V->microV to match GL4k system
Gecog = 5000;   % ecog channel gain
Glfp = 25000;   % LFP channel gain

CAN_IN3 = CAN_IN3*C1;
CECOG1  = CECOG1*C1*C2/Gecog;
CECOG2  = CECOG2*C1*C2/Gecog;
CECOG3  = CECOG3*C1*C2/Gecog;
CECOG4  = CECOG4*C1*C2/Gecog;
CECOG5  = CECOG5*C1*C2/Gecog;
CLFP1   = CLFP1*C1*C2/Glfp;



%% downsample ecog/LFP 1.5kHz -> 1kHz
% sampling rate of data from AlphaOmega (1.5kHz) needs to match that of
% GL4ka data (1kHz). Use function 'resample' to resample the original data at
% 2/3 times the original sampling rate, ie. downsample 1.5kHz->1kHz
% CECOG1 = resample(CECOG1,2,3);
% CECOG2 = resample(CECOG2,2,3);
% CECOG3 = resample(CECOG3,2,3);
% CECOG4 = resample(CECOG4,2,3);
% CECOG5 = resample(CECOG5,2,3);
% CLFP1 = resample(CLFP1,2,3);

%% store processed ecog/LFP data into ecog structure array

ecog = struct('contact_pair',{},'rest_time',{},'active_time',{});

ecog(1).contact_pair(1).raw_ecog_signal=CECOG1;
ecog(1).contact_pair(2).raw_ecog_signal=CECOG2;
ecog(1).contact_pair(3).raw_ecog_signal=CECOG3;
ecog(1).contact_pair(4).raw_ecog_signal=CECOG4;
ecog(1).contact_pair(5).raw_ecog_signal=CECOG5;
ecog(1).contact_pair(6).raw_ecog_signal=CLFP1;

%% grab rest/active timestamps from task button

% ask user for aux channel with task button voltage
evt = menu('Select one','rest/active data','rest only');

if evt == 1
    
    nsamples = length(CAN_IN3); % assuming task voltage  recorded in aux ch.3
    T = 1/(1000*CAN_IN3_KHz);  % period (sec) between samples, multiplies by 1000 for Hz->kHz conversion
    % create time vector
    time = 0:T:T*(nsamples-1);
    MARGIN = 0.1; % 0.3 is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
    ACTIVE = [5-MARGIN 5+MARGIN];
     %REST = [4.35-MARGIN 4.35+MARGIN];  
% %     REST = [0.44-MARGIN 0.44+MARGIN];
%     REST = [5.25-MARGIN 5.25+MARGIN];
    %     REST = [1.69-MARGIN 1.69+MARGIN];
    % ACTIVE =  [4.36-MARGIN 4.36+MARGIN];
% %     ACTIVE =  [5-MARGIN 5+MARGIN];
%     ACTIVE =  [4.46-MARGIN 4.46+MARGIN];
%         ACTIVE =  [5.25-MARGIN 5.25+MARGIN];
    % ACTIVE =  [3.60-MARGIN 3.60+MARGIN];
    %ACTIVE =  [5-MARGIN 5+MARGIN];
    %ACTIVE = [4.35-MARGIN 4.35+MARGIN];
    %ACTIVE = [0.44-MARGIN 0.44+MARGIN];
    REST = [0.44-MARGIN 0.44+MARGIN];
    [rest_time,rest_val] = getevts(time,CAN_IN3',ACTIVE,REST); % for row vector output, input row vector time and column vector voltage
    [active_time,active_val] = getevts(time,CAN_IN3',REST,ACTIVE);
    
    
    
    try
        % Throw out undesirably close timestamps (< 0.1 sec) within rest_time and active_time.
        % This results from large voltage fluctuations that sometimes happens when task button is pressed or unpressed
        idx = find(diff(rest_time)>0.1) + 1; % offset of 1 is added for indexing later
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
            rest_tmp = rest_time(end); % save last rest timestamp for later
            rest_time = rest_time(1:end-1); % remove last rest epoch timestamp
        end
        % throw out rest/active epoch pair timestamps that are too close.  for
        % ecogfft analysis, timestamps must be greater than 2.5 sec apart
        if length(rest_time) == length(active_time)
              idx = find(abs(active_time - rest_time)>2.5);
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
        
    catch
        error('Check ecog-event voltage channel and make sure the right channel was selected');
    end
    
    % save rest/active timestamps in ecog structure array
    ecog(1).rest_time=rest_time;
    ecog(1).active_time=active_time;
    
end

%% save data
% use filename to create output name
outputname=strrep(filename,'.mat','_ecog.mat');
save(outputname, 'ecog','M1_ch','Fs');
