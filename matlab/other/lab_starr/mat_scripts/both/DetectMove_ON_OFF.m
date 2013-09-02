function [move_onset,move_offset,bad_move_onset,bad_move_offset] = DetectMove_ON_OFF(time, signal1,go_signal,stop_signal)

% for now I am taking out signals 2 and 3
%function [move_onset,move_offset,bad_move_onset,bad_move_offset] = DetectMove_ON_OFF(time, signal1,signal2,signal3, go_signal,stop_signal)



%   Detects onset of signal, created for use in time_psd analysis
%	Input:
%       time    - vector containing time, assuming sampling rate of 1kHz
%       signal     - vector containing signal data, assuming sampling rate of 1kHz
%       epoch   - active epoch timestamps to be used for signal alignment
%	Output:
%       signal_onset   - vector containing user-defined signal onset timestamps
%
%   Created by: S.Shimamoto 11/6/2008

% Initialize variables
rest_bl = [-6 -3]; % seconds before epoch timestamp used for baseline signal activity
rest_bl2 = [3 5];% seconds after epoch timestamp used for baseline signal activity
move_on_test = [0 8]; % seconds around epoch timestamp to test signal
move_off_test = [0 3]; % seconds around epoch timestamp to test signal

n_epoch = length(go_signal);
move_onset = nan*zeros(1,n_epoch);
move_offset = nan*zeros(1,n_epoch);
bad_move_onset = nan*zeros(1,n_epoch);
bad_move_offset = nan*zeros(1,n_epoch);

c = [];

for i = 1:n_epoch
    
    stop_time=[];
    start_time=[];
    % Baseline at rest before mvt
    rest_inds = find(time>(go_signal(i)+rest_bl(1)) & time<(go_signal(i)+rest_bl(2)));
    signal = signal1*-1;
    rest_mn = mean(signal(rest_inds));
    signal = abs(signal-rest_mn);
    rest_mn = mean(signal(rest_inds));
    rest_sd = std(signal(rest_inds));
    THRESH1 = rest_mn +2*rest_sd;
    
    % plot baseline and treshold
    hold on
    %plot(time,abs(signal3)/10000,'g')
    plot(time,signal,'k')
    %plot(time,signal2,'b')
    
    xlim([go_signal(i)-5 go_signal(i)+10]);
    xlm = xlim;
    hold on
    plot(xlm,[THRESH1 THRESH1],'r:');
    title(['epoch #' num2str(i)]);
    xlabel('time (sec)');
    ylabel('signal activity');
    ylm = ylim;
    
    % Find move on
    move_on_inds = find(time>(go_signal(i)+move_on_test(1)) & time<(go_signal(i)+move_on_test(2)));
    start_ind = find(signal(move_on_inds)>THRESH1);
    xx = evFindGroups(start_ind,1,5);
    if ~isempty(xx)
        start_time = time(move_on_inds(start_ind(xx(1))));
    end
    plot([go_signal(i) go_signal(i)],ylm,'k:');
    if ~isempty(start_time)
        plot([start_time start_time],ylm, 'r')
    end
    % change move on if not correct
    [x,y] = ginput(1);
    if ~isempty(x)
        start_time = x;
%         plot([start_time start_time],ylm, 'm')
    end
    move_onset(i) = start_time;
    
    % Find move off
    xx = evFindGroups(start_ind,1,70);
    if ~isempty(xx)
        stop_time = time(move_on_inds(start_ind(xx(end))));
%         plot([stop_signal(i) stop_signal(i)],ylm,'k:');
    end
    
    if ~isempty(stop_time)
        plot([stop_time stop_time],ylm, 'r')
    end
    % change move on if not correct
    [x,y] = ginput(1);
    if ~isempty(x)
        stop_time = x;
        plot([stop_time stop_time],ylm, 'm')
    end
    move_offset(i) = stop_time;
    
    % Find inappropriate movements
    %     THRESH1 = rest_mn +8*rest_sd;
    if 1>1
        bad_move_inds = find(time>(move_offset(i-1)+3) & time<(go_signal(i)-1));
    else
        bad_move_inds = find(time>(go_signal(i)-5) & time<(go_signal(i)-1));
    end
    
    bad_ind = find(signal(bad_move_inds)>THRESH1);
    xx = evFindGroups(bad_ind,1,40);
    if ~isempty (xx)
        bad_time = time(bad_move_inds(bad_ind(xx(1))));
        bad_move_onset(i) = bad_time;
        plot([bad_time bad_time],ylm, 'b')
        % change move on if not correct
        [x,y] = ginput(1);
        if ~isempty(x)
            bad_time = x;
            plot([bad_time bad_time],ylm, 'c')
        end
        bad_move_onset(i) = bad_time;
        
        bad_time = time(bad_move_inds(bad_ind(xx(end))));
        bad_move_offset(i) = bad_time;
        plot([bad_time bad_time],ylm, 'b')
        % change move on if not correct
        [x,y] = ginput(1);
        if ~isempty(x)
            bad_time = x;
            plot([bad_time bad_time],ylm, 'c')
        end
        bad_move_offset(i) = bad_time;
    end
    close
end

return;




