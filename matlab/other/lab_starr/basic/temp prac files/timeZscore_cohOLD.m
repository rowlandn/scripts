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
    for j=1:num_chan-1
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
    PDogram.contact_pair(j).PDcoh(:,:,i) = C_trans(:,:,j); 
    end
end

for i = 1:num_chan-1
    % calculate mean and std deviation for spectrogram data (power) across
    % subjects for a given contact pair
    PDogram.contact_pair(i).meanPDcoh = mean(PDogram.contact_pair(i).PDcoh,3); 
    PDogram.contact_pair(i).stdPDcoh = std(PDogram.contact_pair(i).PDcoh,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanPD_BL = mean(PDogram.contact_pair(i).meanPDcoh(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    PDogram.contact_pair(i).BLmean_value = mean(meanPD_BL);
%     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
    
    % Want to display an error message if the baseline values are somehow
    % not =1, but currently, this message is being displayed
    % inappropriately (when BLmean_value = 1)
%     if PDogram.contact_pair(i).BLmean_value ~= 1
%        disp 'The average baseline value does not equal 1!';
%     end
    
    % Calculate Z score for this contact pair
    pdZscore(:,:,i) = (PDogram.contact_pair(i).meanPDcoh - ...
        PDogram.contact_pair(i).BLmean_value) ./ PDogram.contact_pair(i).stdPDcoh;

    zcheck = mean(pdZscore(:,1:20,i));
    meanBLz = mean(zcheck);

end

%% Plot
%faxis and taxis data are imported with time_psd data
PD2plot = pdZscore;
% PD2plot(PD2plot > 3) = 3; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% PD2plot(PD2plot < (-3)) = (-3);
% plot spectrogram for all ecog/lfp data
hf1 = figure;
val1 = min(min(min(PD2plot(1:100,:,:))));
val2 = max(max(max(PD2plot(1:100,:,:))));
clims1 = [val1 val2];
data_ch_names = {'e12','e23','e34','e45','e56','LFP'};

for i = 1:num_chan-1
    subplot(2,3,i);
    hold(gca,'on');
    % make the time-frequency plot
    tmp1 = PD2plot(1:100,:,i); % chopping A2plot will allow the whole colobar to be represented
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
        annotation(hf1, 'textbox','String','PD group COH Z-scores','HorizontalAlignment','left',...
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
    
    num_chan = size(C_trans, 3);

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
    DYSogram.contact_pair(j).DYScoh(:,:,i) = C_trans(:,:,j); 
    end
end

for i = 1:num_chan
    % calculate mean and std deviation for spectrogram data (power) across
    % subjects for a given contact pair
    DYSogram.contact_pair(i).meanDYScoh = mean(DYSogram.contact_pair(i).DYScoh,3); 
    DYSogram.contact_pair(i).stdDYScoh = std(DYSogram.contact_pair(i).DYScoh,0,3); 
    
    % calculate the mean of the baseline period across subjects for given
    % contact pair
    meanDYS_BL = mean(DYSogram.contact_pair(i).meanDYScoh(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
    DYSogram.contact_pair(i).BLmean_value = mean(meanDYS_BL);
%     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
    
    % Want to display an error message if the baseline values are somehow
    % not =1, but currently, this message is being displayed
    % inappropriately (when BLmean_value = 1)
%     if PDogram.contact_pair(i).BLmean_value ~= 1
%        disp 'The average baseline value does not equal 1!';
%     end
    
    % Calculate Z score for this contact pair
    dysZscore(:,:,i) = (DYSogram.contact_pair(i).meanDYScoh - ...
        DYSogram.contact_pair(i).BLmean_value) ./ DYSogram.contact_pair(i).stdDYScoh;

    zcheck = mean(dysZscore(:,1:20,i));
    meanBLz = mean(zcheck);

end


%% Plot
%faxis and taxis data are imported with time_psd data
DYS2plot = dysZscore;
DYS2plot(DYS2plot > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes
DYS2plot(DYS2plot < (-1.96)) = (-1);
DYS2plot(DYS2plot~=1 & DYS2plot~=(-1)) = 0;
% plot spectrogram for all ecog/lfp data
hf2 = figure;
val1 = min(min(min(DYS2plot(1:100,:,:))));
val2 = max(max(max(DYS2plot(1:100,:,:))));
clims1 = [val1 val2];
data_ch_names = {'e12','e23','e34','e45','e56','LFP'};

for i = 1:num_chan
    subplot(2,3,i);
    hold(gca,'on');
    % make the time-frequency plot
    tmp1 = DYS2plot(1:100,:,i); % chopping DYS2plot will allow the whole colobar to be represented
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
        annotation(hf2, 'textbox','String','DYS group COH Z-scores','HorizontalAlignment','left',...
        'Linestyle','none','Position',[ 0.003469 0.97 0.1281 0.02706]);
%     end
    % put a color scale indicator next to the time-coherence plot
    colorbar([0.9307 0.1048 0.02354 0.8226]);
end