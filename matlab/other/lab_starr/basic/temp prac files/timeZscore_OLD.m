%% timeZscore
%This program is designed to generate plots of z-scores of group timePSD and
%timeCOH data. 

%INPUT: timePSD_emg or timePSD_accel data from time_psd.m

%OUTPUT: none at this time (in future, will likely need excel file of
%z-score data)

%ALC 6/5/09

%Update: added a few lines for normalizing to baseline period. This is
%already done in the timePSD code, and in the future that code will output
%the normalized data, making this unnecessary here. 6/8/09ALC

%Update: now using code that has the data already normalized to baseline,
%so that normalization code has been commented out. 6/9/09ALC

    %-------------------------------------------------------
%     first = int32((0/0.05)+1); %int32(((PRE-BL(1))/t_res)+1)
%     last = int32(1/0.05); %int32((PRE-BL(2))/t_res)
    %-------------------------------------------------------
%% Set directory path for finding relevant PD and Dys files and Import data
%This gets directory info for each file within the PD folder. For each file in directory, 
%need to import transcoh struct, select resting and active coh data under M1 contact 
%and place into separate (new) rest and active structures or matrices

% for PD
pdpath = uigetdir('', 'Select directory that contains PD _timePSD_accel or _timePSD_emg files to be analyzed');
pdpath = [pdpath '\'];
cd(pdpath);

PDdir = dir('*.mat'); % selects all .mat files in directory **may want to update to make more specific to timePSD files

numPD = length(PDdir);

for i = 1:numPD
    filename = PDdir(i).name;
    load(filename);
    
    num_chan = size(A2plot, 3);
    %--------------------------------------------
%     [nfchans,nframes] = size(Smag_mean(:,:,1));
    %--------------------------------------------
    for j=1:num_chan
% PDogram will be a large structure containing z-score analysis of
%spectrogram and coherogram data further divided by contact 
% PDspect will be 3D matrix of aggregated (group) spectrogram data in the format:
%rows=freq data, col=time data, sheets = subjects
    
    % ------------------------------------------------------------------
    % temporary code to normalize each subject's data to their baseline
    % period. In the future, this will be imported along with the other
    % variables from _timePSD_data.mat
%         for k = 1:nfchans
%             bl = Smag_mean(k,first:last,i);
%             blmean = mean(bl);
%             Smag_mean(k,:,i) = Smag_mean(k,:,i)/blmean; 
%         end
    %--------------------------------------------------------------------
    PDogram.contact_pair(j).PDspect(:,:,i) = A2plot(:,:,j); 
    end
end

% BLgrand_mean = []; % initialize structure for keeping track of mean baseline values across 
                   % patients and contacts
for i = 1:num_chan
    % calculate mean and std deviation for spectrogram data (power) across
    % subjects for a given contact pair
    PDogram.contact_pair(i).meanPDspect = mean(PDogram.contact_pair(i).PDspect,3); 
    PDogram.contact_pair(i).stdPDspect = std(PDogram.contact_pair(i).PDspect,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanPD_BL = mean(PDogram.contact_pair(i).meanPDspect(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    PDogram.contact_pair(i).BLmean_value = mean(meanPD_BL);
%     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
    
    % Want to display an error message if the baseline values are somehow
    % not =1, but currently, this message is being displayed
    % inappropriately (when BLmean_value = 1)
%     if PDogram.contact_pair(i).BLmean_value ~= 1
%        disp 'The average baseline value does not equal 1!';
%     end
    
    % Calculate Z score for this contact pair
    pdZscore(:,:,i) = (PDogram.contact_pair(i).meanPDspect - ...
        PDogram.contact_pair(i).BLmean_value) ./ PDogram.contact_pair(i).stdPDspect;

    zcheck = mean(pdZscore(:,1:20,i));
    meanBLz = mean(zcheck);

end

% BLgrand_mean_value = mean(BLgrand_mean); % averages baseline values for all contact pairs together

% if BLgrand_mean_value ~= 1
%     disp 'The average baseline value does not equal 1!';
% end
%     for j = 1:numPD
%     zPDspect(:,:,i) = (PDogram.contact_pair(i).PDspect(:,:,j) - meanPDspect(:,:,i)) ./ ...
%         stdPDspect(:,:,i);
%     end


%% Plot
%faxis and taxis data are imported with time_psd data
PD2plot = pdZscore;

% plot spectrogram for all ecog/lfp data
hf1 = figure;
% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
val1 = (min(min(min(PD2plot(:,:,:))))) / 10; % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(PD2plot(:,:,:))))) / 10;
clims1 = [val1 val2];
data_ch_names = {'e12','e23','e34','e45','e56','LFP'};

for i = 1:num_chan
    subplot(2,3,i);
    hold(gca,'on');
    % make the time-frequency plot
    tmp1 = PD2plot(1:100,:,i); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
%     if i==1
% %         if typedata==1
% % %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
% %             title([outputname sprintf('\n')...
% %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% %                 data_ch_names{i} ' aligned to EMG Ch.' emg_ch_names{emg_i}]);
% %         elseif typedata==2
% % %             filename = strrep(filename,'_a.mat',''); %delete '_a.mat' from filename
% %             title([outputname sprintf('\n')...
% %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% %                 data_ch_names{i} ' aligned to accel']);
% %         else
% %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
%             title([outputname sprintf('\n')...
%                 '# epochs=' num2str(n_epochs) sprintf('\n')...
%                 data_ch_names{i} ' aligned to task button']);
% %         end
%     else
        title(data_ch_names{i});
        annotation(hf1, 'textbox','String','PD group PSD Z-scores','HorizontalAlignment','left',...
        'Linestyle','none','Position',[ 0.003469 0.97 0.1281 0.02706]);
%     end
    % put a color scale indicator next to the time-coherence plot
    colorbar([0.9307 0.1048 0.02354 0.8226]);
end
    
%% Repeat for Dystonia data
% Set directory path for finding relevant Dys files and import data
%This gets directory info for each file within the DYS folder. For each file in directory, 
%need to import Smag_mean data (eventually A2plot data) from timePSD
%output. Then compile and average across all patients in the group.

% for DYS
dyspath = uigetdir('', 'Select directory that contains DYS _timePSD_accel or _timePSD_emg files to be analyzed');
dyspath = [dyspath '\'];
cd(dyspath);

DYSdir = dir('*.mat'); % selects all .mat files in directory **may want to update to make more specific to timePSD files

numDYS = length(DYSdir);

for i = 1:numDYS
    filename = DYSdir(i).name;
    load(filename);
    
    num_chan = size(Smag_mean, 3);

    for j=1:num_chan
% PDogram will be a large structure containing z-score analysis of
%spectrogram and coherogram data further divided by contact 
% PDspect will be 3D matrix of aggregated (group) spectrogram data in the format:
%rows=freq data, col=time data, sheets = subjects
    
    % ------------------------------------------------------------------
    % temporary code to normalize each subject's data to their baseline
    % period. In the future, this will be imported along with the other
    % variables from _timePSD_data.mat
%     [nfchans,nframes] = size(Smag_mean(:,:,1));
%     first = int32((0/0.05)+1); %int32(((PRE-BL(1))/t_res)+1)
%     last = int32(1/0.05); %int32((PRE-BL(2))/t_res)
%         for k = 1:nfchans
%             bl = Smag_mean(k,first:last,i);
%             blmean = mean(bl);
%             Smag_mean(k,:,i) = Smag_mean(k,:,i)/blmean; 
%         end
%     %--------------------------------------------------------------------
    DYSogram.contact_pair(j).DYSspect(:,:,i) = A2plot(:,:,j); 
    end
end

for i = 1:num_chan
    % calculate mean and std deviation for spectrogram data (power) across
    % subjects for a given contact pair
    DYSogram.contact_pair(i).meanDYSspect = mean(DYSogram.contact_pair(i).DYSspect,3); 
    DYSogram.contact_pair(i).stdDYSspect = std(DYSogram.contact_pair(i).DYSspect,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanDYS_BL = mean(DYSogram.contact_pair(i).meanDYSspect(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    DYSogram.contact_pair(i).BLmean_value = mean(meanDYS_BL);
%     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
    
    % Want to display an error message if the baseline values are somehow
    % not =1, but currently, this message is being displayed
    % inappropriately (when BLmean_value = 1)
%     if PDogram.contact_pair(i).BLmean_value ~= 1
%        disp 'The average baseline value does not equal 1!';
%     end
    
    % Calculate Z score for this contact pair
    dysZscore(:,:,i) = (DYSogram.contact_pair(i).meanDYSspect - ...
        DYSogram.contact_pair(i).BLmean_value) ./ DYSogram.contact_pair(i).stdDYSspect;

    zcheck = mean(dysZscore(:,1:20,i));
    meanBLz = mean(zcheck);

end

%% Plot
%faxis and taxis data are imported with time_psd data
DYS2plot = dysZscore;

% plot spectrogram for all ecog/lfp data
hf2 = figure;
% val1 = (min(min(min(DYS2plot(1:100,:,:))))) / 10;
% val2 = (max(max(max(DYS2plot(1:100,:,:))))) / 10;
% clims1 = [val1 val2];
clims1 = [-3.5 3.5]; % temporarily setting colorbar limits to these values as min/max are too big
data_ch_names = {'e12','e23','e34','e45','e56','LFP'};

for i = 1:num_chan
    subplot(2,3,i);
    hold(gca,'on');
    % make the time-frequency plot
    tmp1 = DYS2plot(1:100,:,i); %chopping DYS2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
%     if i==1
% %         if typedata==1
% % %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
% %             title([outputname sprintf('\n')...
% %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% %                 data_ch_names{i} ' aligned to EMG Ch.' emg_ch_names{emg_i}]);
% %         elseif typedata==2
% % %             filename = strrep(filename,'_a.mat',''); %delete '_a.mat' from filename
% %             title([outputname sprintf('\n')...
% %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% %                 data_ch_names{i} ' aligned to accel']);
% %         else
% %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
%             title([outputname sprintf('\n')...
%                 '# epochs=' num2str(n_epochs) sprintf('\n')...
%                 data_ch_names{i} ' aligned to task button']);
% %         end
%     else
        title(data_ch_names{i});
        annotation(hf2, 'textbox','String','DYS group PSD Z-scores','HorizontalAlignment','left',...
        'Linestyle','none','Position',[ 0.003469 0.97 0.1281 0.02706]);
%     end
    % put a color scale indicator next to the time-coherence plot
    colorbar([0.9307 0.1048 0.02354 0.8226]);
end