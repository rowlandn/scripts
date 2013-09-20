function Determine_events

%% load file
[filename pathname]=uigetfile('*_ecog.mat','Select .mat file containing ecog/lfp raw data');
cd(pathname);
load([pathname filename]);

assignin('base','pathname',pathname)
assignin('base','filename',filename)

%% Determine channels

ipad_chan = aux.chan(1).raw;
accel_chan = aux.chan(2).raw;
button_chan = aux.chan(2).raw;
%button_chan = aux.chan(2).raw;
% emg_chan1 = emg.chan(1).raw;
% emg_chan2 = emg.chan(2).raw;
% emg_chan3 = emg.chan(3).raw;

%% Detect Ipad sound-event to sunchronize with the ecog
% % filter the signal
% Fs=2000
% lowcut = 2*999/Fs;
% highcut = 2*999.99999/Fs;
% wn=[lowcut,highcut];
% b=fir1(100, wn);
% signal_filtered = filtfilt(b,1,ipad_chan);
% 
% % Find signal above threshold
% THRESH = mean(signal_filtered)+2*std(signal_filtered);
% sup_idx=find(signal_filtered > THRESH);
% [pos,n] = evFindGroups(sup_idx,5,50);
% idx = pos(1,:);
% ipad_ON =sup_idx(idx);
% ecog.ipad_ON_time = int32(ipad_ON/2);
% 
% aux.chan(1).raw =resample(signal_filtered,1,2);
% aux.chan(2).raw=resample(aux.chan(2).raw,1,2);
% % aux.chan(3).raw=resample(aux.chan(3).raw,1,2);
% emg.chan(1).raw =resample(emg.chan(1).raw,1,2);
% emg.chan(2).raw =resample(emg.chan(2).raw,1,2);
% % emg.chan(3).raw =resample(emg.chan(3).raw,1,2);
% ecog.contact_pair(1).raw_ecog_signal=resample(ecog.contact_pair(1).raw_ecog_signal,1,2);
% ecog.contact_pair(2).raw_ecog_signal=resample(ecog.contact_pair(2).raw_ecog_signal,1,2);
% ecog.contact_pair(3).raw_ecog_signal=resample(ecog.contact_pair(3).raw_ecog_signal,1,2);
% ecog.contact_pair(4).raw_ecog_signal=resample(ecog.contact_pair(4).raw_ecog_signal,1,2);
% ecog.contact_pair(5).raw_ecog_signal=resample(ecog.contact_pair(5).raw_ecog_signal,1,2);
% ecog.contact_pair(6).raw_ecog_signal=resample(ecog.contact_pair(6).raw_ecog_signal,1,2);
% Fs=1000;
% save(name,'Fs','ecog','emg','aux','-append')
%% Take the timing of each events

%load test_ipad_import

n_steps = description{2};
n_steps = str2num(n_steps(length(n_steps)))-1;

n_trials = description{end-1};
n_trials = str2num(n_trials(length(n_trials)))+1;

n_trials = 20;
% Sound
[sound]=strfind(description,'Sound');
sound = cellfun('isempty',sound);
sound_ON_idx = find(sound==0);
sound_ON = timestamp(sound_ON_idx)';

% trial start
[trial]=strfind(description,'Trial');
trial = cellfun('isempty',trial);
trial_ON_idx = find(trial==0);
trial_ON = timestamp(trial_ON_idx)';

% Rest epoch Beg Fixation point ON red dot
[rest]=strfind(description,'Scene: 0 - Rest');
rest = cellfun('isempty',rest);
rest_ON_idx = find(rest==0);
rest_ON = timestamp(rest_ON_idx)';

% Rest epoch End Fixation point OFF red dot
[rest]=strfind(description,'Timer: 0');
rest = cellfun('isempty',rest);
rest_OFF_idx = find(rest==0);
rest_OFF = timestamp(rest_OFF_idx)';

% Rest epoch Error Fixation error by mvt
[rest]=strfind(description,'Incorrect touch: 0');
rest = cellfun('isempty',rest);
rest_error_idx = find(rest==0);
rest_error = timestamp(rest_error_idx)';

% Preparation epoch Beg ON Cue ON blue dot
[prep]=strfind(description,'Scene: 1 - Preparation');
prep = cellfun('isempty',prep);
prep_ON_idx = find(prep==0);
prep_ON = timestamp(prep_ON_idx)';

% Preparation epoch End Cue OFF blue dot
[prep]=strfind(description,'Timer: 1');
prep = cellfun('isempty',prep);
prep_OFF_idx = find(prep==0);
prep_OFF = timestamp(prep_OFF_idx)';

% Preparation epoch error by mvt
[prep]=strfind(description,'Incorrect touch: 1');
prep = cellfun('isempty',prep);
prep_error_idx = find(prep==0);
prep_error = timestamp(prep_error_idx)';

% Movements
target_ON = nan*ones(n_steps,n_trials)';
target_idx = nan*ones(n_steps,n_trials)';

touch = nan*ones(n_steps,n_trials)';
notouch = nan*ones(n_steps,n_trials)';
incorrect_touch = nan*ones(n_steps,n_trials)';
touch = nan*ones(n_steps,n_trials)';

% Find 'Scene : movement'
[target]=strfind(description,'Movement');
target = cellfun('isempty',target);
target = find(target==0);
target = target(2:end);

for i = 1:n_steps
    i
    % Target1 ON
    target_idx(:,i) = target(i:n_steps:end);
    target_ON(:,i) = timestamp(target(i:n_steps:end));
end

% Find the touch time
touch_idx = (target_idx) +1;
touch_ON = timestamp(touch_idx);

% timing


Time_events = [rest_ON prep_ON prep_OFF touch_ON];

assignin('base','Time_events',Time_events)
assignin('base','ecog',ecog)


Time_events = int32(Time_events) + ecog.ipad_ON_time(1);%%%!!!! Synchronized ipad events with GS4000 or alphaomega systeme. Always a shift between both systems!!! 
sup = find(Time_events(:,end)< length(ecog.contact_pair(1).raw_ecog_signal));%% ECog signal sometime cut
Time_events = Time_events(sup,:);


ecog.trial_beg_time = Time_events(sup,1);
ecog.prep_time = Time_events(sup,2);
ecog.go_time = Time_events(sup,3);
ecog.touch_time = Time_events(sup,4);
ecog.trial_end_time = Time_events(sup,end);

save(name,'Time_events','ecog','-append')
n_ecog = size(ecog.contact_pair,2);

%% Determine mvt onset
% if size(Time_events,1)~=length(ecog.rest_time)
    Determine_move_on_off(name)
% end

%% Plot raw data

% figure
% subplot(4,1,1)
% plot(ipad_chan)
% ylim([-0.01 0.01])
% title('ipad')
% subplot(4,1,2)
% if n_ecog>20
%     M1_ch=M1_ch1;
% end
% plot(ecog.contact_pair(M1_ch).remontaged_ecog_signal)
% ylim([-500 500])
% title('ecog')
% 
% subplot(4,1,3)
% plot(accel_chan)
% ylim([2 4])
% title('accel')
% 
% % subplot(4,1,4)
% % plot(emg_chan2)
% % ylim([-2500 2500])
% 
% title('emg')
% hold on 
% for i= 1:size(Time_events,1)
%     plot([Time_events(i,1) Time_events(i,1)], [-2500 2500],'k')
%     plot([Time_events(i,2) Time_events(i,2)], [-2500 2500],'r')
%     plot([Time_events(i,3) Time_events(i,3)], [-2500 2500],'g')
% end
% % save figure
% saveas(gcf,[name(1:end-4) '_raw'],'fig');
load(name)
trials_ok = 1:size(Time_events,1);
bad_trials = input('bad_trials: ');
bad_trials = find(ecog.bad_move_on_time~=0);
if ~isempty(bad_trials)
trials_ok = setdiff(trials_ok,bad_trials);
end
n_ecog = size(ecog.contact_pair,2);
save(name,'trials_ok','n_ecog','-append')


% % Find 'Correct touch'
% [touch]=strfind(description,'Correct touch');
% touch = cellfun('isempty',touch);
% touch_idx = find(touch==0);
% touch = timestamp(touch_idx);
%
% % No touch
% [notouch1]=strfind(description,'Timer: 2');
% notouch1 = cellfun('isempty',notouch1);
% notouch1_idx = find(notouch1==0);
% notouch1 = timestamp(notouch1_idx);
%
% % Incorrect touch
% [notouch1]=strfind(description,'Timer: 2');
% notouch1 = cellfun('isempty',notouch1);
% notouch1_idx = find(notouch1==0);
% notouch1 = timestamp(notouch1_idx);
%
% for j= 1:n_trials
%     for i = 1:n_steps
%         idx = target_idx(j,i)
%         if ~isempty(strfind(description{idx+1},'Timer'))
%             notouch = strfind(description{i},'Timer');
%         elseif
%
%         elseif
%
%
%         description{i}
%     end
%
%
%
%
%
