function [move_onset,move_offset,bad_move_onset,bad_move_offset] = DetectMove_ON_OFF(time, signal1,signal2,signal3, go_signal,stop_signal)
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
rest_bl = [-3 -1]; % seconds before epoch timestamp used for baseline signal activity
rest_bl2 = [3 5];% seconds after epoch timestamp used for baseline signal activity
move_on_test = [-1 10]; % seconds around epoch timestamp to test signal
move_off_test = [0 5]; % seconds around epoch timestamp to test signal

n_epoch = length(go_signal);
move_onset = nan*zeros(1,n_epoch);
move_offset = nan*zeros(1,n_epoch);
bad_move_onset = nan*zeros(1,n_epoch);
bad_move_offset = nan*zeros(1,n_epoch);

c = [];

h = figure;
for i = 1:n_epoch
    
    % Baseline at rest before mvt
    rest_inds = find(time>(go_signal(i)+rest_bl(1)) & time<(go_signal(i)+rest_bl(2)));
    signal = signal1*-1;
    rest_mn = mean(signal(rest_inds));
    signal = abs(signal-rest_mn);
    rest_mn = mean(signal(rest_inds));
    rest_sd = std(signal(rest_inds));
    THRESH1 = rest_mn +3*rest_sd;
    %     THRESH2 = rest_mn -3*rest_sd;
    
    % plot baseline and treshold
    hold on
    plot(time,abs(signal3)/5,'g')
    plot(time,signal,'k')
    plot(time,signal2,'b')
    
    xlim([go_signal(i)-5 go_signal(i)+10]);
    xlm = xlim;
    hold on
    plot(xlm,[THRESH1 THRESH1],'r:');
    %     plot(xlm,[THRESH2 THRESH2],'r:');
    title(['epoch #' num2str(i)]);
    xlabel('time (sec)');
    ylabel('signal activity');
    ylm = ylim;
    
    % Find move on
    move_on_inds = find(time>(go_signal(i)+move_on_test(1)) & time<(go_signal(i)+move_on_test(2)));
    %     start_ind = find(signal(move_on_inds)>THRESH1 | signal(move_on_inds)<THRESH2); %#ok<MXFND>
    start_ind = find(signal(move_on_inds)>THRESH1); %#ok<MXFND>
    xx = evFindGroups(start_ind,1,10);
    start_time = time(move_on_inds(start_ind(xx(1))));
    
    plot([go_signal(i) go_signal(i)],ylm,'k:');
    if ~isempty(start_time)
        plot([start_time start_time],ylm, 'r')
    end
    % change move on if not correct
    [x,y] = ginput(1);
    start_time = x;
    if ~isempty(start_time)
        plot([start_time start_time],ylm, 'm')
        move_onset(i) = start_time;
    end
    
    % Baseline at rest after move
    %     rest_inds = find(time>(stop_signal(i)+rest_bl2(1)) & time<(stop_signal(i)+rest_bl2(2)));
    %     rest_mn = mean(signal(rest_inds));
    %     rest_sd = std(signal(rest_inds));
    %     THRESH1 = rest_mn +3*rest_sd;
    %     THRESH2 = rest_mn -3*rest_sd;
    %
    %     % Find move off
        move_off_inds = find(time>(stop_signal(i)+move_off_test(1)) & time<(stop_signal(i)+move_off_test(2)));
%         stop_ind = find(signal(move_off_inds)<THRESH1 & signal(move_off_inds)>THRESH2); %#ok<MXFND>
        stop_ind = find(signal(move_off_inds)<THRESH1); %#ok<MXFND>
    %     x = evFindGroups(stop_ind,1,150);
      xx = evFindGroups(stop_ind,1,150);
    stop_time = time(move_off_inds(stop_ind(xx(1))));
%     stop_time = time(move_on_inds(stop_ind(x(end))));
%     stop_time = time(move_on_inds(start_ind(xx(end))));
    plot([stop_signal(i) stop_signal(i)],ylm,'k:');
%     plot(xlm,[THRESH1 THRESH1],'g:');
    %     plot(xlm,[THRESH2 THRESH2],'g:');
    if ~isempty(stop_time)
        plot([stop_time stop_time],ylm, 'r')
    end
    % change move on if not correct
    [x,y] = ginput(1);
    stop_time = x;
    if ~isempty(stop_time)
        plot([stop_time stop_time],ylm, 'm')
        move_offset(i) = stop_time;
    end
    
    % Determine on off of inappropriate mvt
    %     % Baseline at rest before mvt
    %     rest_inds = find(time>(go_signal(i)+rest_bl(1)) & time<(go_signal(i)+rest_bl(2)));
    %     rest_mn = mean(signal(rest_inds));
    %     rest_sd = std(signal(rest_inds));
    %     THRESH1 = rest_mn +3*rest_sd;
    %     THRESH2 = rest_mn -3*rest_sd;
    
    if 1>1
        bad_move_inds = find(time>(move_offset(i-1)+1) & time<(go_signal(i)));
    else
        bad_move_inds = find(time>(go_signal(i)-5) & time<(go_signal(i)));
    end
    %     start_ind = find(signal(bad_move_inds)>THRESH1 | signal(bad_move_inds)<THRESH2); %#ok<MXFND>
    start_ind = find(signal(bad_move_inds)>THRESH1); %#ok<MXFND>
    xx = evFindGroups(start_ind,1,15);
    if ~isempty (xx)
        start_time = time(bad_move_inds(start_ind(xx(1))));
        bad_move_onset(i) = start_time;
        plot([start_time start_time],ylm, 'b')
        % change move on if not correct
        [x,y] = ginput(1);
        start_time = x;
        if ~isempty(start_time)
            plot([start_time start_time],ylm, 'c')
            bad_move_onset(i) = start_time;
        end
        
        start_time = time(bad_move_inds(start_ind(xx(2))));
        bad_move_offset(i) = start_time;
        plot([start_time start_time],ylm, 'b')
        % change move on if not correct
        [x,y] = ginput(1);
        start_time = x;
        if ~isempty(start_time)
            plot([start_time start_time],ylm, 'c')
            bad_move_offset(i) = start_time;
        end
    end
    hold off
end

close(h)

return;




