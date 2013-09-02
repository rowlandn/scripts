% aomatssep.m

% this m-file performs SSEP analysis on mat file generated from map format
% using Alpha Omega's mapfile converter program

%% Define variables

PRE_STIM = -0.05;     % pre-stim period in sec
PST_STIM = 0.05;     % post-stim period in sec

%% load data
%load .mat file 
[filename pathname]=uigetfile('*.mat','Select .mat file');
cd(pathname);
load(filename);

% make sure sampling frequencies are equal for all ecog channels
% as of 2/3/09, we are using sampling rate of 6kHz for both ecog and stim
% trig voltage
if isequal(CECOG_1_KHz,CECOG_2_KHz,CECOG_3_KHz,CECOG_4_KHz,CECOG_5_KHz)
    Fs = CECOG_1_KHz;
    Fs = Fs*1000;                           % multiply by 1000 for KHz->Hz conversion
else
    error('Sampling frequencies are different for ecog/LFP data');
end


%% force ecog and time vector to be the same lengths
minlength=min([length(CECOG_1) length(CECOG_2) length(CECOG_3) length(CECOG_4) length(CECOG_5)]);
CECOG_1  = CECOG_1(1:minlength);
CECOG_2  = CECOG_2(1:minlength);
CECOG_3  = CECOG_3(1:minlength);
CECOG_4  = CECOG_4(1:minlength);
CECOG_5  = CECOG_5(1:minlength);
%% voltage calibration
C1 = 5/32768; % constant for calibration
% note: value 32768 is equal to 5V in the Alpha Omega system
C2 = 1e6;   % constant to convert ecog/lfp channel voltage level from V->microV to match GL4k system
Gecog = 5000;   % ecog channel gain 
Glfp = 25000;   % LFP channel gain
% calibrate all aux channels
CAn_In__1 = CAn_In__1*C1;
CAn_In__2 = CAn_In__2*C1;
CAn_In__3 = CAn_In__3*C1;
% calibrate all ecog channels
CECOG_1  = CECOG_1*C1*C2/Gecog;
CECOG_2  = CECOG_2*C1*C2/Gecog;
CECOG_3  = CECOG_3*C1*C2/Gecog;
CECOG_4  = CECOG_4*C1*C2/Gecog;
CECOG_5  = CECOG_5*C1*C2/Gecog;

%% montage and store in ECOGALL array
% ECOGALL contains montaged ecog data for each contact in each of its columns
ECOGALL=[(CECOG_1-CECOG_2)'...
    (CECOG_2-CECOG_3)'... 
    (CECOG_3-CECOG_4)'... 
    (CECOG_4-CECOG_5)'... 
    CECOG_5'];

% ECOGALL contains raw ecog data for each contact in each of its columns
% ECOGALL=[CECOG_1'...
%     CECOG_2'... 
%     CECOG_3'... 
%     CECOG_4'... 
%     CECOG_5'];
% 
%% Process stim trigger voltage channel
% Ask user which aux channel contains stim trigger voltage recording.
% In the future, consider eliminating menu selection and hardwiring the
% aux channel number when we start using the same aux chan for every
% patient
% aux_chan = {'CAn_In__1','CAn_In__2','CAn_In__3'};
% aux_trig = aux_chan(menu('Select aux chan with stim trigger voltage recording',aux_chan));
% aux_trig = cell2mat(aux_trig);
aux_trig = cell2mat({'CAn_In__3'});

% at this point, aux_trig is a string array with name of variable
% containing stim trig voltage level
eval(['trig_chan = ' aux_trig ';']);
eval(['T = 1/(' aux_trig '_KHz*1000);']);  % period (sec) between each sample in stim trig channel
% the sampling rate in line above is multiplied by 1000 for KHz->Hz conversion
nsamples = length(trig_chan);
% using T and nsamples, create time vector of same length as trig_chan
trig_time=0:T:T*(nsamples-1);


MARGIN = 0.01;  % 0.35 is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
min_trig = mean(trig_chan)+2*std(trig_chan);
max_trig  = max(trig_chan);
STIM_TRIG_OFF =  [0-MARGIN 0+MARGIN];
% STIM_TRIG_ON = [5000-MARGIN 5000+MARGIN];
% STIM_TRIG_ON = [22260-MARGIN 22260+MARGIN];
% STIM_TRIG_ON = [3.4-MARGIN 3.4+MARGIN];
STIM_TRIG_ON = [min_trig max_trig+MARGIN];
[stim_trig_time,stim_trig_val] = getevtsSSEP(trig_time, abs(trig_chan)', STIM_TRIG_OFF,STIM_TRIG_ON);

% check that the margin is ok
diff = stim_trig_time(2:end)-stim_trig_time(1:end-1);
if nanmean(diff)> 0.6 && nanmean(diff)< 0.4 
    error('Trig SSEP different from 2 Hz');
end

% find and eliminate any stim trig times that might exceed the length of
% ecog data
tmax = (minlength/Fs)-PST_STIM;
stim_trig_time = stim_trig_time(stim_trig_time<tmax);

num_stim = length(stim_trig_time);
t = 1000*(PRE_STIM:1/Fs:PST_STIM); % multiplied by 1000 to change to msec scale

% initialize
% diff_stim_epoch contains 5 layers of 2-D array corresponding to the 5
% montaged ecog pairs. The rows of each 2-D array contain a vector of raw
% ecog data around the time of each stimulation.  
diff_stim_epoch = zeros(num_stim,length(t),5);


for i = 1:5
    
    % parse each stim epoch from each contact pair
    for j = 1:num_stim
        tmp1 = int32((stim_trig_time(j)+PRE_STIM) * Fs);  % int32 used to keep index in integer format
        tmp2 = int32((stim_trig_time(j)+PST_STIM) * Fs);
        if (tmp2-tmp1)>length(t)-1
            tmp2=tmp2-1;
        end
        diff_stim_epoch(j,:,i) = ECOGALL(tmp1:tmp2,i);
    end
end


mean_diff_stim_epoch = mean(diff_stim_epoch);

figure;
hold on
t_text = 1000*(PRE_STIM+60*T);  % specify where along the time axis to place ecog pair text
t_ind = find(t == t_text);  % find index at t_text

for i = 1:5
    SSEPmin = min(mean_diff_stim_epoch(1,:,i)); % mix value of average
    SSEPmax = max(mean_diff_stim_epoch(1,:,i)); % max value of average
    C = 5 - i; % constant added to stack the waves for comparison such that the first contact pair is plotted at the top
%     z = (mean_diff_stim_epoch(1,:,i)-SSEPmin)/(SSEPmax-SSEPmin) + C;
    z = -(mean_diff_stim_epoch(1,:,i)-SSEPmin)/(SSEPmax-SSEPmin) + C; % SSEPs recorded using AO shows N20 as an up-going potential.  invert this with a negative sign to make it down-going.
    plot(t,z); 
    text(t_text,z(t_ind)+0.1,...
        ['e' num2str(i) '-' num2str(i+1)]);
end
xlabel('Time (msec)');
ylabel('normalized SSEPs');
ylm = ylim;
plot([0 0],[ylm(1) ylm(2)],'k--')
title(strrep(filename,'.mat',''));
hold off
