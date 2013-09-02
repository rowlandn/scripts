function timePSD()
% Creates a time-varying power spectral density plot and a time-varying
% coherence plot aligned to movement onset.
% Movement onset for each active movement epoch is determined using
% emg/accel/task button timing data.

% Note: .mat files containing ecog/LFP data and EMG data must be in the
% same directory, otherwise warning message will prompt user to do so then
% try again.

% Created by S.Shimamoto  11/6/2008

% Revised by SS (12/17/2008) to load 'ecog' strucutre array from APMconv7 for raw
% ecog/lfp data instead of using raw_ecog_lfp double array from
% APMconv7_coh.  Script title changed from 'time_psd' to 'timePSD' for
% consistency with 'ecogPSD'

% Revised by SS (2/24/2009) to accomodate for data recorded by AlphaOmega
% system.  timePSD works on GL4k AND AlphaOmega recordings.


%% initialize

Fs=1000;    %ecog/LFP sampling rate for GL4k
NFFT=512;
NOVERLAP=462;
WINDOW=512;
FRAME_ADVANCE=WINDOW-NOVERLAP;
PRE = 2;      % time (sec) before movement onset
POST = 2;       % time (sec) after
% POST = 6;       % SAS 5/15/09: extended for testing
ADD = 1;     % add more time to increase # windows in each snippet
BL = [2 1];  % baseline period before movement onset for calculating percent power change


%% load ecog data

%load .mat file containing ecog_lfp_raw_data
[filename pathname]=uigetfile('*_ecog.mat','Select .mat file containing ecog/lfp raw data');
cd(pathname);
load([pathname filename]);

%% ask user about system used for recording data
typesystem = menu('Select the system used for recording data','Guideline 4000','Alpha Omega');
% note: typesystem = 1 (GL4k), typesystem = 2 (AO)

%% load emg/accel data

% look for previously saved _time_psd.mat files for movement onset times
% note: emg onset times are saved with the filename ending '_time_psd_emg.mat'
% accel onset times are saved with '_time_psd_accel.mat.' 

emg_files = dir(fullfile(pathname,[filename(1:11) '*_timePSD_emg.mat']));
accel_files = dir(fullfile(pathname,[filename(1:11) '*_timePSD_accel.mat']));

if ~isempty(emg_files) || ~isempty(accel_files)
    usedata = menu('Previously saved onset times were retrieved!','Use them','Do not use them');
    % note: usedata = 1 (use previously saved onset times), usedata = 2 (do
    % NOT use them)
    
    if usedata==1
        
        % both emg and accel onset times exit
        if ~isempty(emg_files) && ~isempty(accel_files)
            typedata = menu('Both EMG and accel are available. Select one','EMG','accel');
            if typedata==1
                eval(['load ' emg_files.name ' move_onset emg_i emg_ch_names']);
                filename=strrep(filename,'_ecog','');% this takes off the _ecog ending of the filename;
                outputname = strrep(filename,'.mat',''); % outputname has no '.mat' extensions or any other misc tags; used later to save figures
            elseif typedata==2
                eval(['load ' accel_files.name ' move_onset']);
                filename=strrep(filename,'_ecog','_a');% replaces _ecog with _a
                outputname = strrep(filename,'_a.mat','');
            end
            
        % only emg onset times exist
        elseif ~isempty(emg_files) && isempty(accel_files)
            eval(['load ' emg_files.name ' move_onset emg_i emg_ch_names']);
            typedata = 1;
            filename=strrep(filename,'_ecog',''); % this takes off the _ecog ending of the filename;
            outputname = strrep(filename,'.mat','');
            
        % only accel onset times exist
        elseif isempty(emg_files) && ~isempty(accel_files)
            eval(['load ' accel_files.name ' move_onset']);
            typedata = 2;
            filename=strrep(filename,'_ecog','_a');
            outputname = strrep(filename,'_a.mat','');
        end
        
    end
    
end

%% determine movement onset

% if there are no previously saved onset times or if user wants to run
% onset detection again, run detectEMG
if (isempty(emg_files) && isempty(accel_files)) || usedata == 2

    typedata = menu('Detect movement onset using EMG/accel/task button?','EMG','accel','task button');
    % note: typedata = 1 (EMG data), typedata = 2 (accel data), typedata = 3 (task button)
    
    %load .mat file containing emg or accel data
    %note: there are separate blocks of code for loading GL4k and AO data
    
    %------------GL4k------------------
    if typesystem == 1 % GL4k recording
        try
            if typedata==1 % load emg data
                filename=strrep(filename,'_ecog',''); % this takes off the _ecog ending of the filename to load emg data
                load([pathname filename]);
                outputname = strrep(filename,'.mat','');
            elseif typedata==2 % load accel dat
                filename=strrep(filename,'_ecog','_a');% replaces _ecog with _a to load accel data
                load([pathname filename]);
                outputname = strrep(filename,'_a.mat','');
            else % task-button timestamps have already been loaded in _ecog.mat
                filename=strrep(filename,'_ecog','');
                outputname = strrep(filename,'.mat','');
            end
        catch
            disp('WARNING: matfiles must be in the same directory');
            choice = menu(['WARNING: .mat files for ecog/lfp and' sprintf('\n')...
                'emg/accel must be in the same directory.' sprintf('\n') sprintf('\n')...
                'Place files in same directory, then click button below.'], 'Try again');
            if choice == 1
                load([pathname filename]);
            end
        end
        
        if typedata==1
            % extract emg channel names
            n_emg = size(emg_chan,1);
            emg_ch_names = cell(1,n_emg);
            for i=1:n_emg
                emg_ch_names{i} = emg_names(i,:);
            end
            % select emg channel for processing movement onset
            emg_i = menu('Select emg channel',emg_ch_names);
        end
        
        % find active epoch timestamps from ecog structure
        epoch_ts = ecog.active_time;
        
        % determine movement onset time
        if typedata~=3
            if typedata==1
                move = emg_chan(emg_i,:);
            elseif typedata==2
                move = accel_chan;
            end
            move_onset = DetectEMG(time,move,epoch_ts);
            move_onset = move_onset(~isnan(move_onset));
        else
            move_onset = epoch_ts;
        end
        
        % at this point, movement onsets have been detected for GL4k data
        
    %------------AlphaOmega------------------
    elseif typesystem == 2; % AlphaOmega recording
        try
            filename=strrep(filename,'_ecog',''); % load _ecog.mat created by AOmatconv.m
            load([pathname filename]); % load .mat created by AOMapConverter program
            outputname = strrep(filename,'.mat','');
        catch
            disp('WARNING: matfiles must be in the same directory');
            choice = menu(['WARNING: .mat files for ecog/lfp and' sprintf('\n')...
                'emg/accel must be in the same directory.' sprintf('\n') sprintf('\n')...
                'Place files in same directory, then click button below.'], 'Try again');
            if choice == 1
                load([pathname filename]);
            end
        end
        
        % Note by SAS (5/29/09):
        % sampling frequency on aux channels on AO system is different from
        % GL4k (AO=1.5kHz; GL4k=1kHz).  For timePSD analysis, aux channels
        % are used to record EMG, accel, task button data.  Since such
        % non-neuronal data are used for TIMING PURPOSES (ie. to detect
        % movement onset) and NOT for rigorous signal processing (such as 
        % spectrogram or mscohere where careful consideration of Fs is necessary), 
        % AO aux channel data will NOT be downsampled and will NOT have 
        % identical sampling frequency as GL4k aux channel data.
        
        if typedata==1 % process emg data
            emg_ch_names = {'CBiceps','CDeltoid','CECR'};
            emg_i = menu('Select emg channel',emg_ch_names);
            eval(['move =' emg_ch_names{emg_i} ';']);
            eval(['T = 1/(1000*' emg_ch_names{emg_i} '_KHz);']);
        elseif typedata==2 % process accel data
            accel_ch_names= {'CAN_IN1','CAN_IN2','CAN_IN3'};
            accel_i = menu({'Select accel channel',accel_ch_names});
            eval(['move = ' accel_ch_names{accel_i} ';']);
            eval(['T = 1/(1000*' accel_ch_names{accel_i} '_KHz);']);
        end
        %create time matching time vector
        nsamples = length(move); %#ok<NODEF>
        time = 0:T:T*(nsamples-1);
        
        % find active epoch timestamps from ecog structure
        epoch_ts = ecog.active_time;
        
        % determine movement onset time
        if typedata~=3
            move_onset = DetectEMG(time,move,epoch_ts);
            move_onset = move_onset(~isnan(move_onset));
        else
            move_onset = epoch_ts;
        end
        
    end
    
    % at this point, movement onsets have been detected for AO data
end

%% montage ecog data

%extract ecog bipolar recordings from ecog_lfp_raw_data
ecog12=ecog.contact_pair(1).raw_ecog_signal-ecog.contact_pair(2).raw_ecog_signal;
ecog23=ecog.contact_pair(2).raw_ecog_signal-ecog.contact_pair(3).raw_ecog_signal;
ecog34=ecog.contact_pair(3).raw_ecog_signal-ecog.contact_pair(4).raw_ecog_signal;
ecog45=ecog.contact_pair(4).raw_ecog_signal-ecog.contact_pair(5).raw_ecog_signal;
ecog56=ecog.contact_pair(5).raw_ecog_signal;
%concatenate such that each column contains data from one channel
ecoglfpall=[ecog12' ecog23' ecog34' ecog45' ecog56'];

% if LFP data exists, include in ecoglfpall
k_LFP = menu('analyze LFP?','yes','no');
if k_LFP==1
    lfp = ecog.contact_pair(6).raw_ecog_signal;
    ecoglfpall=[ecoglfpall lfp'];
end

%% calculate normalized time-varying PSD

% calculate spectrogram for all ecog/LFP. gives the fft of each
% segment in each column, and the next column moves in time by
% window-noverlap samples

n_data_ch=size(ecoglfpall,2);
n_epochs = length(move_onset);
for i = 1:n_data_ch
    data = ecoglfpall(:,i);
    for j=1:n_epochs
        % take snippet of data around emg_onset from selected ecog/LFP channel
        % add offset to increase # windows for fft
        first = int32(Fs*(move_onset(j)-(PRE))-WINDOW/2); % WINDOW/2 offset will make spectrogram start at moveonset-PRE at appropriately centered PSD
        last = int32(Fs*(move_onset(j)+(POST+ADD))-WINDOW/2);
        snip = data(first:last);
        %calculate spectrogram of snippet
        % 3D matrix 'S' has the fft power value stored in [frequncy,time,epoch
        % number] arrangement
        S(:,:,j) = spectrogram(snip,WINDOW,NOVERLAP,NFFT,Fs); %#ok<AGROW>
        %calculate time-varying coherence if LFP exists
    end

    %find the magnitude of the fft represented in each column of S
    Smag=abs(S);

    %calculate average across all epochs
    %note: Smag contains spectrograms from each epoch in the 3rd
    %dimension. The average across all epochs are then calculated and
    %stored in the 3rd dimension of Smag_mean.  Smag_mean collects averaged
    %spectrogram for each data channel in the 3rd dimension.DO NOT GET CONFUSED!
    Smag_mean(:,:,i) = mean(Smag,3); %#ok<AGROW>
    
    % clear some variables before next loop, this is probably not necessary
    % but do it just in case
    clear data S Smag;
    
end

%setup up the frequency (faxis)and time (taxis) axes data
[nfchans,nframes] = size(Smag_mean(:,:,1));
nfchansteps = nfchans - 1;
maxfreq = Fs/2;
faxis = maxfreq*(0:nfchansteps)/nfchansteps;
t_res = FRAME_ADVANCE/Fs; % temporal resolution of spectrogram (sec)
taxis = (0:(nframes-1))* t_res;
taxis = taxis -PRE; %shift by PRE

% normalize to baseline values
if PRE<BL(1)
    error(['The baseline period for PSD plot currently begins '...
        '%d seconds before onset of movement. This value cannot be more than %d seconds '...
        'as determined by variable PRE'],BL(1),PRE);
else
    first = int32(((PRE-BL(1))/t_res)+1);
    last = int32((PRE-BL(2))/t_res);
    % to plot A with colors representing the log10 of power, uncomment this line:
    A2plot = log10(Smag_mean);
    % to plot A with colors representing raw data values, uncomment this line:
    % A2plot = Smag_mean;
    for i = 1:n_data_ch
        for j = 1:nfchans
            bl = A2plot(j,first:last,i);
            blmean = mean(bl);
            A2plot(j,:,i) = A2plot(j,:,i)/blmean; 
        end
    end
end

%% calculate time-varying transformed coherence

% calculate time-varying coherence if LFP exists
if n_data_ch==6
    lfp = ecoglfpall(:,6);
    for i=1:n_data_ch-1
        data = ecoglfpall(:,i);
        for j = 1:n_epochs
            first = int32(Fs*(move_onset(j)-(PRE))-WINDOW/2);% WINDOW/2 offset will make spectrogram start at moveonset-PRE at appropriately centered PSD
            last = int32(Fs*(move_onset(j)+(POST+ADD))-WINDOW/2);
            snip = data(first:last);
            snip_lfp = lfp(first:last);
            counter=1;
            coh_store=[];
            for k=1:nframes
                x = snip(counter:counter+WINDOW-1);
                y = snip_lfp(counter:counter+WINDOW-1);
%                 coh = mscohere(x,y,WINDOW,[],NFFT,Fs);
                coh = mscohere(x,y,[],[],NFFT,Fs);
%                 coh = mscohere(x,y,[],[],[],Fs);
                coh_store = [coh_store coh]; %#ok<AGROW>
                counter = counter+FRAME_ADVANCE;
            end
            % populate 3D matrix C with coherence for each epoch
            C(:,:,j)=coh_store; %#ok<AGROW>
        end
        % find mean across all epochs
        C_mean(:,:,i) = mean(C,3); %#ok<AGROW>
        % populate 3D matrix C_trans with transformed coherence for each
        % contact pair
        C_trans(:,:,i)=atanh(sqrt(C_mean(:,:,i))); %#ok<AGROW>
    end
end

%% plot data

% plot spectrogram for all ecog/lfp data
hf1 = figure;
val1 = min(min(min(A2plot(1:100,:,:))));
val2 = max(max(max(A2plot(1:100,:,:))));
clims1 = [val1 val2];
data_ch_names = {'e12','e23','e34','e45','e56','LFP'};

for i = 1:n_data_ch
    subplot(2,3,i);
    hold(gca,'on');
    % make the time-frequency plot
    tmp1 = A2plot(1:100,:,i); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
    set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    if i==1
        if typedata==1
%             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
            title([outputname sprintf('\n')...
                '# epochs=' num2str(n_epochs) sprintf('\n')...
                data_ch_names{i} ' aligned to EMG Ch.' emg_ch_names{emg_i}]);
        elseif typedata==2
%             filename = strrep(filename,'_a.mat',''); %delete '_a.mat' from filename
            title([outputname sprintf('\n')...
                '# epochs=' num2str(n_epochs) sprintf('\n')...
                data_ch_names{i} ' aligned to accel']);
        else
%             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
            title([outputname sprintf('\n')...
                '# epochs=' num2str(n_epochs) sprintf('\n')...
                data_ch_names{i} ' aligned to task button']);
        end
    else
        title(data_ch_names{i});
    end
%     colorbar;
end

% put a color scale indicator next to the time-frequency plot
colorbar([0.9307 0.1048 0.02354 0.8226]);


% plot coherence
if n_data_ch==6
    hf2 = figure;
    val3 = min(min(min(C_trans(1:100,:,:))));
    val4 = max(max(max(C_trans(1:100,:,:))));
    clims2 = [val3 val4];
    f=linspace(0,Fs/2,size(C_trans,1));
    for i=1:n_data_ch-1
        subplot(2,3,i);
        hold(gca,'on');
        tmp2 = C_trans(1:100,:,i);%chopping C_trans will allow the whole colobar to be represented
        f_new=f(1:100);
        imagesc(taxis,f_new,tmp2,clims2);
        %     imagesc(taxis,f,C_trans(:,:,i),clims2);
        % set the y-axis direction (YDir) to have zero at the bottom
        set(gca,'YDir','normal');
        % set xlim and ylim
        set(gca,'Xlim',[0-PRE POST]);
        set(gca,'Ylim',[0 120]);
        %plot vertical bar at movement onset
        plot([0 0],ylim,'k:');
        hold(gca,'off');
        % axis labels/title
        xlabel('time (sec)');
        ylabel('frequency (Hz)');
        if i==1
            if typedata==1
                title([outputname sprintf('\n')...
                    '# epochs=' num2str(n_epochs) sprintf('\n')...
                    data_ch_names{i} '-LFP aligned to EMG Ch.' emg_ch_names{emg_i}]);
            elseif typedata==2
                title([outputname sprintf('\n')...
                    '# epochs=' num2str(n_epochs) sprintf('\n')...
                    data_ch_names{i} '-LFP aligned to accel']);
            else
                title([outputname sprintf('\n')...
                    '# epochs=' num2str(n_epochs) sprintf('\n')...
                    data_ch_names{i} '-LFP aligned to task button']);
            end
        else
            title([data_ch_names{i} '-LFP']);
        end
    end
    % put a color scale indicator next to the time-coherence plot
    colorbar([0.9307 0.1048 0.02354 0.8226]);
end

%% save data

% save variables and figure
if typedata==1
    save([outputname '_timePSD_emg.mat'],'Smag_mean','A2plot','C_trans','faxis','taxis','move_onset','emg_i','emg_ch_names');
%     save([outputname '_timePSD_emg.mat'],'Smag_mean','faxis','taxis','move_onset','emg_i','emg_ch_names');
    saveas(hf1,[outputname '_timePSD_emg'],'fig');
    saveas(hf2,[outputname '_timeCOH_emg'],'fig');
elseif typedata==2
    save([outputname '_timePSD_accel.mat'],'Smag_mean','A2plot','C_trans','faxis','taxis','move_onset');
%     save([outputname '_timePSD_accel.mat'],'Smag_mean','faxis','taxis','move_onset');
    saveas(hf1,[outputname '_timePSD_accel'],'fig');
    saveas(hf2,[outputname '_timeCOH_accel'],'fig');
else
    save([outputname '_timePSD_button.mat'],'Smag_mean','A2plot','C_trans','faxis','taxis','move_onset');
    saveas(hf1,[outputname '_timePSD_button'],'fig');
    saveas(hf2,[outputname '_timeCOH_button'],'fig');
end

return;