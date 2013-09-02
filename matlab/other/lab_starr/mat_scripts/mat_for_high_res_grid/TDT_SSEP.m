function TDT_SSEP(name)

% this m-file performs SSEP analysis on mat file generated from map format
% using Alpha Omega's mapfile converter program

%% Define variables

PRE_STIM = -0.05;     % pre-stim period in sec
PST_STIM = 0.05;     % post-stim period in sec

%% load data
%load .mat file 
 ecog = struct('chan',{});
 aux = struct('chan',{});
 emg = struct('chan',{});
 X = struct('chan',{});
 
gain=1;
Fs=1000; % sampling freq after downsampling

curdir = cd;
dd = dir;
dd = dd(3:end);

%% load the data and store processed Aux chan data into structure array
for abc =1 : length(dd)
    if ~isempty(strfind(dd(abc).name,'htk'))
        % load the data and convert in mat.file
        file_name=dd(abc).name;
        [d,fs,dt,tc,t]=readhtk(file_name);
        
         %store the downsampled data in structure arrays
        if ~isempty(strfind(file_name,'ipad.htk'))
            aux(1).chan(1).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'accel.htk'))
            aux(1).chan(2).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'trig.htk'))
            aux(1).chan(3).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'emg.htk'))
            emg(1).chan(1).raw=resample(d,2^10,(5^5)*8);
        elseif ~isempty(strfind(file_name,'signal.htk'))
            signal(1).chan(2).raw=resample(d,2^10,(5^5)*8);
        else            
            ecog(1).contact_pair(tc).raw_ecog_signal=resample(d,2^10,5^5);
        end
    end
end

%% notch filter around 60Hz, 120Hz and 180Hz
% butterworth notch filter - model order, [low/(Fs/2) high/(Fs/2)]
    [n1_b, n1_a]=butter(3,2*[57 63]/Fs,'stop'); %60hz
    [n2_b, n2_a]=butter(3,2*[117 123]/Fs,'stop'); %120hz
    [n3_b, n3_a]=butter(3,2*[177 183]/Fs,'stop'); %180hz
    [n4_b, n4_a]=butter(3,2*[237 243]/Fs,'stop'); %60hz
    [n5_b, n5_a]=butter(3,2*[297 303]/Fs,'stop'); %120hz
    [n6_b, n6_a]=butter(3,2*[357 363]/Fs,'stop'); %180hz
    
    for k=1:length(ecog(1).contact_pair)
        ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n1_b, n1_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 60
        ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n2_b, n2_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 120
        ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n3_b, n3_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 180
        ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n4_b, n4_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 60
        ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n5_b, n5_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 120
        ecog(1).contact_pair(k).raw_ecog_signal=filtfilt(n6_b, n6_a, ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 180
    end
%% remove DC offset
 for k=1:length(ecog(1).contact_pair)
        ecog(1).contact_pair(k).raw_ecog_signal=ecog(1).contact_pair(k).raw_ecog_signal-mean(ecog(1).contact_pair(k).raw_ecog_signal); %notch out at 60
 end
for k=1:length(aux.chan)
       aux.chan(k).raw=aux.chan(k).raw-mean(aux.chan(k).raw);
end
 %% contact to remove

%  unused=input('unused contacts');
%  used = setdiff(1:length(ecog.contact_pair),unused);
%  for i = 1:28
%      j=used(i);
%      X(1).contact_pair(i).raw_ecog_signal=ecog.contact_pair(j).raw_ecog_signal;
%  end
%  ecog=X;
 
%     %% detection of bad electrodes
% figure;hold on
% for i = 1: 14
%     subplot(3,5,i)
%     plot(ecog.contact_pair(i).raw_ecog_signal)
%     title(num2str(i))
% end
% figure;hold on
% for i = 15: 28
%     subplot(3,5,i-14)
%     plot(ecog.contact_pair(i).raw_ecog_signal)
%     title(num2str(i))
% end
% bad=input('bad contacts');
% 
% close all
%% common reference
data = nan*ones(28,length(ecog.contact_pair(1).raw_ecog_signal));
for i = 1: 28
     data(i,:) = ecog.contact_pair(i).raw_ecog_signal';    
end
car=nanmean(data)/28;
for i = 1: 28
    ecog.contact_pair(i).remontaged_ecog_signal = ecog.contact_pair(i).raw_ecog_signal-car;
end

ecogall = nan*ones(28,length(ecog.contact_pair(1).raw_ecog_signal))';
for i = 1 : 28 
    ecogall(:,i) = ecog.contact_pair(i).remontaged_ecog_signal;
end

%% Process stim trigger voltage channel

% trig_chan = aux.chan(3).raw; 
T=1/Fs;% the sampling rate in line above is multiplied by 1000 for KHz->Hz conversion
% nsamples = length(trig_chan);
% % using T and nsamples, create time vector of same length as trig_chan
% trig_time=0:T:T*(nsamples-1);
% 
% 
% MARGIN = 0.01;  % 0.35 is a very wide range that handles large fluctuations in task voltage. ParseIntraOpEvents uses MARGIN = 0.05
% trig_chan=trig_chan-mean(trig_chan);
% x = trig_chan(find(trig_chan~=0));
% min_trig = mean(x(1000:end-1000))+2*std(x(1000:end-1000));
% max_trig  = max(x(1000:end-1000));
% STIM_TRIG_OFF =  [0-MARGIN 0+MARGIN];
% STIM_TRIG_ON = [min_trig max_trig+MARGIN];
% [stim_trig_time,stim_trig_val] = getevtsSSEP(trig_time, abs(trig_chan)', STIM_TRIG_OFF,STIM_TRIG_ON);
plot(aux.chan(3).raw)
STIM_TRIG_ON =input('stim threshold');
close 
inds = find(aux.chan(3).raw(15:end-15)>=STIM_TRIG_ON);
 [pos,n] = evFindGroups(inds,400,1);
stim_trig_time = inds(pos(1,:))+15;
stim_trig_time=stim_trig_time/Fs;
% check that the margin is ok
 %find active period of 1s minimum
% diff = stim_trig_time(2:end)-stim_trig_time(1:end-1);
% t = find(diff>0.4 & diff<0.6);
% stim_trig_time = stim_trig_time(t);
% if ~nanmean(diff)< 0.6 && ~nanmean(diff)> 0.4 
%     error('Trig SSEP different from 2 Hz');
% end

% find and eliminate any stim trig times that might exceed the length of
% ecog data
tmax = (length(ecog.contact_pair(1).remontaged_ecog_signal)/Fs)-PST_STIM;
stim_trig_time = stim_trig_time(stim_trig_time<tmax);
xx = find(stim_trig_time>abs(PRE_STIM));
stim_trig_time = stim_trig_time(xx);
num_stim = length(stim_trig_time);
t = 1000*(PRE_STIM:1/Fs:PST_STIM); % multiplied by 1000 to change to msec scale

% initialize
% diff_stim_epoch contains 5 layers of 2-D array corresponding to the 5
% montaged ecog pairs. The rows of each 2-D array contain a vector of raw
% ecog data around the time of each stimulation.  
diff_stim_epoch = zeros(num_stim,length(t),28);


for i = 1:28
    
    % parse each stim epoch from each contact pair
    for j = 1:num_stim
        tmp1 = int32((stim_trig_time(j)+PRE_STIM) * Fs);  % int32 used to keep index in integer format
        tmp2 = int32((stim_trig_time(j)+PST_STIM) * Fs);
        diff_stim_epoch(j,:,i) = ecogall(tmp1:tmp2,i);
    end
end


mean_diff_stim_epoch = mean(diff_stim_epoch);

figure;
hold on
t_text = 1000*(PRE_STIM+60*T);  % specify where along the time axis to place ecog pair text
t_ind = find(t == t_text);  % find index at t_text

for i = 1:14
    SSEPmin = min(mean_diff_stim_epoch(1,:,i)); % mix value of average
    SSEPmax = max(mean_diff_stim_epoch(1,:,i)); % max value of average
    C = i; % constant added to stack the waves for comparison such that the first contact pair is plotted at the top
%     z = (mean_diff_stim_epoch(1,:,i)-SSEPmin)/(SSEPmax-SSEPmin) + C;
    z = -(mean_diff_stim_epoch(1,:,i)-SSEPmin)/(SSEPmax-SSEPmin) + C; % SSEPs recorded using AO shows N20 as an up-going potential.  invert this with a negative sign to make it down-going.
    plot(t,z); 
   text(t_text,C+0.1,['e' num2str(i) ]);
end
xlabel('Time (msec)');
ylabel('normalized SSEPs');
ylm = ylim;
plot([0 0],[ylm(1) ylm(2)],'k--')
title(name);
hold off
saveas(gcf,[name 'SSEP_raw1'])

figure;
hold on
t_text = 1000*(PRE_STIM+60*T);  % specify where along the time axis to place ecog pair text
t_ind = find(t == t_text);  % find index at t_text

for i = 15:28
    SSEPmin = min(mean_diff_stim_epoch(1,:,i)); % mix value of average
    SSEPmax = max(mean_diff_stim_epoch(1,:,i)); % max value of average
    C = i-14; % constant added to stack the waves for comparison such that the first contact pair is plotted at the top
%     z = (mean_diff_stim_epoch(1,:,i)-SSEPmin)/(SSEPmax-SSEPmin) + C;
    z = -(mean_diff_stim_epoch(1,:,i)-SSEPmin)/(SSEPmax-SSEPmin) + C; % SSEPs recorded using AO shows N20 as an up-going potential.  invert this with a negative sign to make it down-going.
    plot(t,z); 
   text(t_text,C+0.1,['e' num2str(i) ]);
end
xlabel('Time (msec)');
ylabel('normalized SSEPs');
ylm = ylim;
plot([0 0],[ylm(1) ylm(2)],'k--')
title(name);
hold off
saveas(gcf,[name 'SSEP_raw2'])