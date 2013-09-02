function bta = bta_SE(Analog_Signal,sig_descrip,burst_times,...
                        bur_descrip,lag,Fs,plotit,sp1a,sp2a,sp3a,sp1b,sp2b,sp3b)

% bta_SE This function computes the spike-triggered average of an analog
% signal around burst onset and offset spike times, as calculated by the 
% legend_salc_SE function.  Lag (ms) defines the time of the analog signal
% before and after each spike time. Burst spike times should be supplied as
% a cell containing rows of spike time arrays which represent individual
% traces (e.g., the output of the legend_salc_SE function). Sampling
% frequency (Fs) should be given in kHz.   The burst times are also
% randomly shuffled 50 times to produce 99% confidence intervals for the
% original burst-triggered average.  If  multiple burst-triggered averages
% are to be placed on the same page, subplot (sp)  coordinates can be
% given.  If not, these can be omitted.  For plotit:
% 
%                        0 = no plot
%                        1 = onset and offset plotted side-by-side
%                        2 = onset plotted at sp1a,sp2a,sp3a and offset
%                            plotted at sp1b,sp2b,sp3b
%                        3 = onset plot only at sp1,sp2,sp3
%                        4 = offset plot only at sp1,sp2,sp3
%  
%  bta = bta_SE(Analog_Signal,sig_descrip,burst_times,bur_descrip,lag,Fs,plotit)
% 
% 
%  
%  Example 1: Analog_Signal = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',1);
%             Neuron = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',4);
%             spike_times = findspikes_win_SE(Neuron,10,-200,-40,.1,1,0);
%             burst_times = legend_salc_SE(Neuron,spike_times,10,1);
%             bta = bta_SE(Analog_Signal,'EEG',burst_times,'DCN',500,10);
%              
%  Example 2: Analog_Signal = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',1);
%             Neuron = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',4);
%             spike_times = findspikes_win_SE(Neuron,10,-200,-40,.1,1,0);
%             burst_times = legend_salc_SE(Neuron,spike_times,10,1);
%             bta = bta_SE(Analog_Signal,'EEG',burst_times,'DCN',500,10,2,2,1);
%                               (if you want to specify a certain
%                                sublot, use this notation)

% % % % Debugging code; comment out
% Analog_Signal = loadtraces_SE('/F/Data/Raw/viv05/viv0503b.all','1-10',1);
% DCN = loadtraces_SE('/F/Data/Raw/viv05/viv0503b.all','1-10',4);
% spikes_SE = findspikes_win_SE(DCN,10,-200,-40,.1,1,0);
% [burst_times,b,c] = legend_salc_SE(DCN,spikes_SE,10000,0);
% Fs = 10;
% lag = 500;



%%%%%%%%%%%%%%%%%%%%%%%%%%%% BURST ONSET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Here we find the burst onset-triggered average for the %%%%%%%%%%%%
%%%%% original analog signal and burst onset spike times.    %%%%%%%%%%%%
 

for i = 1:size(burst_times.onset,1)
    
    spike_index = round(burst_times.onset{i,1}*Fs);
    spike_train = zeros(size(Analog_Signal(:,i)));
    spike_train(spike_index) = 1;

    [orig_bta_onset(:,i) orig_m_lags_onset] = xcorr(Analog_Signal(:,i),spike_train(:,1),lag*Fs);
    
    orig_bta_onset(:,i) = orig_bta_onset(:,i)/length(burst_times.onset{i,1});
    
end

mean_orig_bta_onset = mean(orig_bta_onset,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear spike_index spike_train


%%%%% Here, we find new, shuffled burst times for each trace %%%%% 
%%%%% that we will use to make new burst-triggered averages. %%%%% 
%%%%% The new, shuffled burst-triggered averages will be used %%%%
%%%%% to calculate confidence intervals for the original burst- %%
%%%%% triggered average.                                     %%%%%

for m = 1:50
        
        %assignin('base','m',m)

        
        for l = 1:size(burst_times.onset,1)
            isi_str = ['isi_trace_',num2str(l),' = diff(burst_times.onset{',num2str(l),',1});'];
            eval(isi_str)
            randperm_str_1 = ['randperm_size_trace_',num2str(l),' = randperm(size(burst_times.onset{',num2str(l),',1},1)-1);'];
            eval(randperm_str_1)
            shuf_str_1 = ['shuf_isis_trace_',num2str(l),' = isi_trace_',num2str(l),'(randperm_size_trace_',num2str(l),');'];
            eval(shuf_str_1);       
            shuf_str_2 = ['shuf_burst_times_',num2str(m),'{',num2str(l),',1} = [1 (cumsum(shuf_isis_trace_',num2str(l),')''+1)];'];            
            eval(shuf_str_2)
        end
        
        burst_times_shuf_str_2 = ['burst_times_shuf = shuf_burst_times_',num2str(m),';'];
        eval(burst_times_shuf_str_2)
        
end

for i = 1:size(burst_times_shuf,1)
    
    spike_index = round(burst_times_shuf{i,1}*Fs);
    spike_train = zeros(size(Analog_Signal(:,i)));
    spike_train(spike_index) = 1;

    [shuf_bta_onset(:,i) shuf_m_lags_onset] = xcorr(Analog_Signal(:,i),spike_train(:,1),lag*Fs);
    
    shuf_bta_onset(:,i) = shuf_bta_onset(:,i)/length(burst_times_shuf{i,1});
    
end

mean_shuf_bta_onset_y = mean(shuf_bta_onset,2);
mean_shuf_bta_onset_x = shuf_m_lags_onset'/Fs;
std3_plus_mean_shuf_bta_onset_y = mean_shuf_bta_onset_y + 3.1*std(mean_shuf_bta_onset_y);
std3_minus_mean_shuf_bta_onset_y = mean_shuf_bta_onset_y - 3.1*std(mean_shuf_bta_onset_y);
mean_std3_plus_mean_shuf_bta_onset_y = mean(std3_plus_mean_shuf_bta_onset_y);
mean_std3_minus_mean_shuf_bta_onset_y = mean(std3_minus_mean_shuf_bta_onset_y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Plot  
% 
% 
% subplot(2,2,1)
% m_t_onset = orig_m_lags_onset/Fs;
% CI_bar_onset = patch(lag*[-1 -1 1 1],...
%      [mean_std3_plus_mean_shuf_bta_onset_y mean_std3_minus_mean_shuf_bta_onset_y mean_std3_minus_mean_shuf_bta_onset_y mean_std3_plus_mean_shuf_bta_onset_y],...
%        'w'); hold on
% set(CI_bar_onset,'facecolor',[.85 .85 .85],'edgecolor',[1 1 1])
% plot(m_t_onset, mean_orig_bta_onset,'k'); hold on
% %plot(mean_shuf_bta_onset_x,mean_shuf_bta_onset_y)
% % plot(mean_shuf_bta_onset_x,std3_plus_mean_shuf_bta_onset_y,'r')
% % plot(mean_shuf_bta_onset_x,std3_minus_mean_shuf_bta_onset_y,'r')    
% x_lim = get(gca,'xlim');
% y_lim = get(gca,'ylim');
% plot([0 0],[y_lim(1) y_lim(2)],'k')
% ylabel('Burst Triggered Average')
% xlabel('Lag (ms)')
%     
% 
% 
% bta.onset.bta_y = mean_orig_bta_onset;
% bta.onset.bta_x = orig_m_lags_onset'/Fs;
% bta.onset.std_plus_shuf_mean = std3_plus_mean_shuf_bta_onset_y;
% bta.onset.std_minus_shuf_mean = std3_minus_mean_shuf_bta_onset_y;
% bta.onset.CI_bar_x = get(CI_bar_onset,'XData');
% bta.onset.CI_bar_y = get(CI_bar_onset,'YData');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% BURST OFFSET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Here we find the burst offset-triggered average for the %%%%%%%%%%%%
%%%%% original analog signal and burst offset spike times.    %%%%%%%%%%%%
 

for i = 1:size(burst_times.offset,1)
    
    spike_index = round(burst_times.offset{i,1}*Fs);
    spike_train = zeros(size(Analog_Signal(:,i)));
    spike_train(spike_index) = 1;

    [orig_bta_offset(:,i) orig_m_lags_offset] = xcorr(Analog_Signal(:,i),spike_train(:,1),lag*Fs);
    
    orig_bta_offset(:,i) = orig_bta_offset(:,i)/length(burst_times.offset{i,1});
    
end

mean_orig_bta_offset = mean(orig_bta_offset,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear spike_index spike_train


%%%%% Here, we find new, shuffled burst times for each trace %%%%% 
%%%%% that we will use to make new burst-triggered averages. %%%%% 
%%%%% The new, shuffled burst-triggered averages will be used %%%%
%%%%% to calculate confidence intervals for the original burst- %%
%%%%% triggered average.                                     %%%%%

for m = 1:50
        
        %assignin('base','m',m)

        
        for l = 1:size(burst_times.offset,1)
            isi_str = ['isi_trace_',num2str(l),' = diff(burst_times.offset{',num2str(l),',1});'];
            eval(isi_str)
            randperm_str_1 = ['randperm_size_trace_',num2str(l),' = randperm(size(burst_times.offset{',num2str(l),',1},1)-1);'];
            eval(randperm_str_1)
            shuf_str_1 = ['shuf_isis_trace_',num2str(l),' = isi_trace_',num2str(l),'(randperm_size_trace_',num2str(l),');'];
            eval(shuf_str_1);       
            shuf_str_2 = ['shuf_burst_times_',num2str(m),'{',num2str(l),',1} = [1 (cumsum(shuf_isis_trace_',num2str(l),')''+1)];'];            
            eval(shuf_str_2)
        end
        
        burst_times_shuf_str_2 = ['burst_times_shuf = shuf_burst_times_',num2str(m),';'];
        eval(burst_times_shuf_str_2)
        
end

for i = 1:size(burst_times_shuf,1)
    
    spike_index = round(burst_times_shuf{i,1}*Fs);
    spike_train = zeros(size(Analog_Signal(:,i)));
    spike_train(spike_index) = 1;

    [shuf_bta_offset(:,i) shuf_m_lags_offset] = xcorr(Analog_Signal(:,i),spike_train(:,1),lag*Fs);
    
    shuf_bta_offset(:,i) = shuf_bta_offset(:,i)/length(burst_times_shuf{i,1});
    
end

mean_shuf_bta_offset_y = mean(shuf_bta_offset,2);
mean_shuf_bta_offset_x = shuf_m_lags_offset'/Fs;
std3_plus_mean_shuf_bta_offset_y = mean_shuf_bta_offset_y + 3.1*std(mean_shuf_bta_offset_y);
std3_minus_mean_shuf_bta_offset_y = mean_shuf_bta_offset_y - 3.1*std(mean_shuf_bta_offset_y);
mean_std3_plus_mean_shuf_bta_offset_y = mean(std3_plus_mean_shuf_bta_offset_y);
mean_std3_minus_mean_shuf_bta_offset_y = mean(std3_minus_mean_shuf_bta_offset_y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if plotit == 0
    
elseif plotit == 1
    
    % Onset
    subplot(2,2,1)
    m_t_onset = orig_m_lags_onset/Fs;
    CI_bar_onset = patch(lag*[-1 -1 1 1],...
        [mean_std3_plus_mean_shuf_bta_onset_y mean_std3_minus_mean_shuf_bta_onset_y mean_std3_minus_mean_shuf_bta_onset_y mean_std3_plus_mean_shuf_bta_onset_y],...
        'w'); hold on
    set(CI_bar_onset,'facecolor',[.85 .85 .85],'edgecolor',[1 1 1])
    plot(m_t_onset, mean_orig_bta_onset,'k'); hold on
    %plot(mean_shuf_bta_onset_x,mean_shuf_bta_onset_y)
    % plot(mean_shuf_bta_onset_x,std3_plus_mean_shuf_bta_onset_y,'r')
    % plot(mean_shuf_bta_onset_x,std3_minus_mean_shuf_bta_onset_y,'r')    
    x_lim = get(gca,'xlim');
    y_lim = get(gca,'ylim');
    plot([0 0],[y_lim(1) y_lim(2)],'k')
    % ylabel('Burst Triggered Average')
    % xlabel('Lag (ms)')
    set(gca,'FontSize',6)

    bta.onset.bta_y = mean_orig_bta_onset;
    bta.onset.bta_x = orig_m_lags_onset'/Fs;
    bta.onset.std_plus_shuf_mean = std3_plus_mean_shuf_bta_onset_y;
    bta.onset.std_minus_shuf_mean = std3_minus_mean_shuf_bta_onset_y;
    bta.onset.CI_bar_x = get(CI_bar_onset,'XData');
    bta.onset.CI_bar_y = get(CI_bar_onset,'YData');

    [min_bta_onset_y,min_bta_onset_i] = min(bta.onset.bta_y);
    [max_bta_onset_y,max_bta_onset_i] = max(bta.onset.bta_y);

    if abs(min_bta_onset_y) > abs(max_bta_onset_y)
        bta.onset.bta_peak_lag = bta.onset.bta_x(min_bta_onset_i);
        bta.onset.bta_peak_bta = min_bta_onset_y;
    elseif abs(max_bta_onset_y) > abs(min_bta_onset_y)
        bta.onset.bta_peak_lag = bta.onset.bta_x(max_bta_onset_i);
        bta.onset.bta_peak_bta = max_bta_onset_y;
    end
        
    hold on; plot(bta.onset.bta_peak_lag,bta.onset.bta_peak_bta,'ro')


    % Offset  
    subplot(2,2,2)
    m_t_offset = orig_m_lags_offset/Fs;
    CI_bar_offset = patch(lag*[-1 -1 1 1],...
        [mean_std3_plus_mean_shuf_bta_offset_y mean_std3_minus_mean_shuf_bta_offset_y mean_std3_minus_mean_shuf_bta_offset_y mean_std3_plus_mean_shuf_bta_offset_y],...
        'w'); hold on
    set(CI_bar_offset,'facecolor',[.85 .85 .85],'edgecolor',[1 1 1])
    plot(m_t_offset, mean_orig_bta_offset,'k'); hold on
    % plot(mean_shuf_bta_offset_x,mean_shuf_bta_offset_y)
    % plot(mean_shuf_bta_offset_x,std3_plus_mean_shuf_bta_offset_y,'r')
    % plot(mean_shuf_bta_offset_x,std3_minus_mean_shuf_bta_offset_y,'r')    
    x_lim = get(gca,'xlim');
    y_lim = get(gca,'ylim');
    plot([0 0],[y_lim(1) y_lim(2)],'k')
    % ylabel('Burst Triggered Average')
    % xlabel('Lag (ms)')
    set(gca,'FontSize',6)
    
    bta.offset.bta_y = mean_orig_bta_offset;
    bta.offset.bta_x = orig_m_lags_offset'/Fs;
    bta.offset.std_plus_shuf_mean = std3_plus_mean_shuf_bta_offset_y;
    bta.offset.std_minus_shuf_mean = std3_minus_mean_shuf_bta_offset_y;
    bta.offset.CI_bar_x = get(CI_bar_offset,'XData');
    bta.offset.CI_bar_y = get(CI_bar_offset,'YData');

    [min_bta_offset_y,min_bta_offset_i] = min(bta.offset.bta_y);
    [max_bta_offset_y,max_bta_offset_i] = max(bta.offset.bta_y);

    if abs(min_bta_offset_y) > abs(max_bta_offset_y)
        bta.offset.bta_peak_lag = bta.offset.bta_x(min_bta_offset_i);
        bta.offset.bta_peak_bta = min_bta_offset_y;
    elseif abs(max_bta_offset_y) > abs(min_bta_offset_y)
        bta.offset.bta_peak_lag = bta.offset.bta_x(max_bta_offset_i);
        bta.offset.bta_peak_bta = max_bta_offset_y;
    end

    hold on; plot(bta.offset.bta_peak_lag,bta.offset.bta_peak_bta,'ro')
    
elseif plotit == 2  
    
    % Onset
    subplot(sp1a,sp2a,sp3a)
    m_t_onset = orig_m_lags_onset/Fs;
    CI_bar_onset = patch(lag*[-1 -1 1 1],...
        [mean_std3_plus_mean_shuf_bta_onset_y mean_std3_minus_mean_shuf_bta_onset_y mean_std3_minus_mean_shuf_bta_onset_y mean_std3_plus_mean_shuf_bta_onset_y],...
        'w'); hold on
    set(CI_bar_onset,'facecolor',[.85 .85 .85],'edgecolor',[1 1 1])
    plot(m_t_onset, mean_orig_bta_onset,'k'); hold on
    %plot(mean_shuf_bta_onset_x,mean_shuf_bta_onset_y)
    % plot(mean_shuf_bta_onset_x,std3_plus_mean_shuf_bta_onset_y,'r')
    % plot(mean_shuf_bta_onset_x,std3_minus_mean_shuf_bta_onset_y,'r')    
    
    % ylabel('Burst Triggered Average')
    % xlabel('Lag (ms)')
    set(gca,'FontSize',6)

    bta.onset.bta_y = mean_orig_bta_onset;
    bta.onset.bta_x = orig_m_lags_onset'/Fs;
    bta.onset.std_plus_shuf_mean = std3_plus_mean_shuf_bta_onset_y;
    bta.onset.std_minus_shuf_mean = std3_minus_mean_shuf_bta_onset_y;
    bta.onset.CI_bar_x = get(CI_bar_onset,'XData');
    bta.onset.CI_bar_y = get(CI_bar_onset,'YData');

    [min_bta_onset_y,min_bta_onset_i] = min(bta.onset.bta_y);
    [max_bta_onset_y,max_bta_onset_i] = max(bta.onset.bta_y);

    if abs(min_bta_onset_y) > abs(max_bta_onset_y)
        bta.onset.bta_peak_lag = bta.onset.bta_x(min_bta_onset_i);
        bta.onset.bta_peak_bta = min_bta_onset_y;
    elseif abs(max_bta_onset_y) > abs(min_bta_onset_y)
        bta.onset.bta_peak_lag = bta.onset.bta_x(max_bta_onset_i);
        bta.onset.bta_peak_bta = max_bta_onset_y;
    end
        
    %hold on; plot(bta.onset.bta_peak_lag,bta.onset.bta_peak_bta,'ro')


    % Offset  
    subplot(sp1b,sp2b,sp3b)
    m_t_offset = orig_m_lags_offset/Fs;
    CI_bar_offset = patch(lag*[-1 -1 1 1],...
        [mean_std3_plus_mean_shuf_bta_offset_y mean_std3_minus_mean_shuf_bta_offset_y mean_std3_minus_mean_shuf_bta_offset_y mean_std3_plus_mean_shuf_bta_offset_y],...
        'w'); hold on
    set(CI_bar_offset,'facecolor',[.85 .85 .85],'edgecolor',[1 1 1])
    plot(m_t_offset, mean_orig_bta_offset,'k'); hold on
    % plot(mean_shuf_bta_offset_x,mean_shuf_bta_offset_y)
    % plot(mean_shuf_bta_offset_x,std3_plus_mean_shuf_bta_offset_y,'r')
    % plot(mean_shuf_bta_offset_x,std3_minus_mean_shuf_bta_offset_y,'r')    
    
    % ylabel('Burst Triggered Average')
    % xlabel('Lag (ms)')
    set(gca,'FontSize',6)
    
    bta.offset.bta_y = mean_orig_bta_offset;
    bta.offset.bta_x = orig_m_lags_offset'/Fs;
    bta.offset.std_plus_shuf_mean = std3_plus_mean_shuf_bta_offset_y;
    bta.offset.std_minus_shuf_mean = std3_minus_mean_shuf_bta_offset_y;
    bta.offset.CI_bar_x = get(CI_bar_offset,'XData');
    bta.offset.CI_bar_y = get(CI_bar_offset,'YData');

    [min_bta_offset_y,min_bta_offset_i] = min(bta.offset.bta_y);
    [max_bta_offset_y,max_bta_offset_i] = max(bta.offset.bta_y);

    if abs(min_bta_offset_y) > abs(max_bta_offset_y)
        bta.offset.bta_peak_lag = bta.offset.bta_x(min_bta_offset_i);
        bta.offset.bta_peak_bta = min_bta_offset_y;
    elseif abs(max_bta_offset_y) > abs(min_bta_offset_y)
        bta.offset.bta_peak_lag = bta.offset.bta_x(max_bta_offset_i);
        bta.offset.bta_peak_bta = max_bta_offset_y;
    end

    %hold on; plot(bta.offset.bta_peak_lag,bta.offset.bta_peak_bta,'ro')

    
    %%%% Adjust axes and plot line at x = 0
    %%%% Make sure you copy this code for plotit = 1
    subplot(sp1a,sp2a,sp3a)
    y_lim_onset = get(gca,'ylim');
    %assignin('base','y_lim_onset1',y_lim_onset)
    
    subplot(sp1b,sp2b,sp3b)
    y_lim_offset = get(gca,'ylim');
    %assignin('base','y_lim_offset',y_lim_offset)
    
    if y_lim_onset(2) == y_lim_offset(2)
    elseif y_lim_onset(2) < y_lim_offset(2)
        y_lim_onset(2) = y_lim_offset(2);
        %assignin('base','y_lim_onset2',y_lim_onset)
        subplot(sp1a,sp2a,sp3a)
        set(gca,'ylim',y_lim_onset)
    elseif y_lim_onset(2) > y_lim_offset(2)
        y_lim_offset(2) = y_lim_onset(2);
        subplot(sp1b,sp2b,sp3b)
        set(gca,'ylim',y_lim_offset)
    end
    
    
    if y_lim_onset(1) == y_lim_offset(1)
    elseif y_lim_onset(1) < y_lim_offset(1)
        y_lim_offset(1) = y_lim_offset(1);
        subplot(sp1b,sp2b,sp3b)
        set(gca,'ylim',y_lim_offset)
    elseif y_lim_onset(1) > y_lim_offset(1)
        y_lim_onset(1) = y_lim_offset(1);
        subplot(sp1a,sp2a,sp3a)
        set(gca,'ylim',y_lim_onset)
    end
        
    subplot(sp1a,sp2a,sp3a)
    plot([0 0],y_lim_onset,'k')
        
    subplot(sp1b,sp2b,sp3b)
    plot([0 0],y_lim_offset,'k')
    y_lim_offset_title = round(y_lim_offset(2)+.5*y_lim_offset(2));
    text(-850,y_lim_offset_title,sig_descrip)
    
          
    
        
elseif plotit == 3  %Onset only
    subplot(sp1a,sp2a,sp3a)
    m_t_onset = orig_m_lags_onset/Fs;
    CI_bar_onset = patch(lag*[-1 -1 1 1],...
        [mean_std3_plus_mean_shuf_bta_onset_y mean_std3_minus_mean_shuf_bta_onset_y mean_std3_minus_mean_shuf_bta_onset_y mean_std3_plus_mean_shuf_bta_onset_y],...
        'w'); hold on
    set(CI_bar_onset,'facecolor',[.85 .85 .85],'edgecolor',[1 1 1])
    plot(m_t_onset, mean_orig_bta_onset,'k'); hold on
    %plot(mean_shuf_bta_onset_x,mean_shuf_bta_onset_y)
    % plot(mean_shuf_bta_onset_x,std3_plus_mean_shuf_bta_onset_y,'r')
    % plot(mean_shuf_bta_onset_x,std3_minus_mean_shuf_bta_onset_y,'r')    
    x_lim = get(gca,'xlim');
    y_lim = get(gca,'ylim');
    plot([0 0],[y_lim(1) y_lim(2)],'k')
    ylabel('Burst Triggered Average')
    xlabel('Lag (ms)')
    
    bta.onset.bta_y = mean_orig_bta_onset;
    bta.onset.bta_x = orig_m_lags_onset'/Fs;
    bta.onset.std_plus_shuf_mean = std3_plus_mean_shuf_bta_onset_y;
    bta.onset.std_minus_shuf_mean = std3_minus_mean_shuf_bta_onset_y;
    bta.onset.CI_bar_x = get(CI_bar_onset,'XData');
    bta.onset.CI_bar_y = get(CI_bar_onset,'YData');

elseif plotit == 4  %Offset only
    subplot(sp1a,sp2a,sp3a)
    m_t_offset = orig_m_lags_offset/Fs;
    CI_bar_offset = patch(lag*[-1 -1 1 1],...
        [mean_std3_plus_mean_shuf_bta_offset_y mean_std3_minus_mean_shuf_bta_offset_y mean_std3_minus_mean_shuf_bta_offset_y mean_std3_plus_mean_shuf_bta_offset_y],...
       'w'); hold on
    set(CI_bar_offset,'facecolor',[.85 .85 .85],'edgecolor',[1 1 1])
    plot(m_t_offset, mean_orig_bta_offset,'k'); hold on
    % plot(mean_shuf_bta_offset_x,mean_shuf_bta_offset_y)
    % plot(mean_shuf_bta_offset_x,std3_plus_mean_shuf_bta_offset_y,'r')
    % plot(mean_shuf_bta_offset_x,std3_minus_mean_shuf_bta_offset_y,'r')    
    x_lim = get(gca,'xlim');
    y_lim = get(gca,'ylim');
    plot([0 0],[y_lim(1) y_lim(2)],'k')
    ylabel('Burst Triggered Average')
    xlabel('Lag (ms)')

    bta.offset.bta_y = mean_orig_bta_offset;
    bta.offset.bta_x = orig_m_lags_offset'/Fs;
    bta.offset.std_plus_shuf_mean = std3_plus_mean_shuf_bta_offset_y;
    bta.offset.std_minus_shuf_mean = std3_minus_mean_shuf_bta_offset_y;
    bta.offset.CI_bar_x = get(CI_bar_offset,'XData');
    bta.offset.CI_bar_y = get(CI_bar_offset,'YData');

end