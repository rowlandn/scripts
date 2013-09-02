function betadesynch
%graph and define % beta desynchronization across groups
%% define variables for plotting
FREQ_QPSD = [4 12;...     % delta alpha band
             13 21;...   % low beta band
             22 30;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band   

FREQBANDS = {num2str(FREQ_QPSD(1,:)), num2str(FREQ_QPSD(2,:)) num2str(FREQ_QPSD(3,:)) ...
        num2str(FREQ_QPSD(4,:)) num2str(FREQ_QPSD(5,:))};
% AREAS = {'pre-motor' 'M1' 'S1' 'STN LFP'};
BRAINAREAS = {'M1' 'S1' 'STN LFP'};
LIMBAREAS = {'HAND' 'ELBOW' 'SHOULDER' 'JAW' 'FOOT' '"ARM"' '"NON-ARM"'};
YLIM = [-100 100];

%% choose channel
k = menu('Select brain area for analysis',BRAINAREAS);
l = menu('Select limb for analysis',LIMBAREAS);

%% grab PD data
% pathnamePD = uigetdir('','Select directory with PD ecogPSD.mat files');
pdpath = uigetdir('', 'Select directory that contains PD _ecogPSD files to be analyzed');
pdpath = [pdpath '\'];
cd(pdpath);

PDdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numPD = length(PDdir);
counterPD = 0;

for i=1:numPD
    filename = PDdir(i).name;
    load(filename);
    if isnan(order(k));
        continue; % skip files without contact pair over selected area
    end 
    % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
%     [allfreq subfreq] = NORMquantPSDshort(rest,active,FREQ_QPSD,order,freq,filename);
%     [subfreq] = quantDesync(rest,active,FREQ_QPSD,order,freq,filename);
    PDdesync(:,i) = ((subfreq(:,2,k)-subfreq(:,4,k))./subfreq(:,2,k))*100;
    PDdesync_abs(:,i) = ((subfreq(:,1,k)-subfreq(:,3,k))./subfreq(:,1,k))*100;
    counterPD = counterPD+1;
    clear subfreq
end
YDesync(:,1) = mean(PDdesync,2);
EDesync(:,1) = std(PDdesync,0,2);
EDesync(:,1) = EDesync(:,1)/sqrt(numPD); %calculating SEM rather than std
YDesync_abs(:,1) = mean(PDdesync_abs,2);
EDesync_abs(:,1) = std(PDdesync_abs,0,2);
EDesync_abs(:,1) = EDesync_abs(:,1)/sqrt(numPD); %calculating SEM rather than std
%% grab DYS data
% pathnamePD = uigetdir('','Select directory with PD ecogPSD.mat files');
dyspath = uigetdir('', 'Select directory that contains DYS _ecogPSD files to be analyzed');
dyspath = [dyspath '\'];
cd(dyspath);

DYSdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numDYS = length(DYSdir);
counterDYS = 0;

for i=1:numDYS
    filename = DYSdir(i).name;
    load(filename);
    if isnan(order(k));
        continue; % skip files without contact pair over selected area
    end 
    % for evaluating freq bands other than those already represented in
    % subfreq matrix, rerun quantPSDmatrix function now, with new freq
    % bands as defined above
%     [subfreq] = quantDesync(rest,active,FREQ_QPSD,order,freq,filename);
    DYSdesync(:,i) = ((subfreq(:,2,k)-subfreq(:,4,k))./subfreq(:,2,k))*100;
    DYSdesync_abs(:,i) = ((subfreq(:,1,k)-subfreq(:,3,k))./subfreq(:,1,k))*100;
    counterDYS = counterDYS+1;
    clear subfreq
end
YDesync(:,2) = mean(DYSdesync,2);
EDesync(:,2) = std(DYSdesync,0,2);
EDesync(:,2) = EDesync(:,2)/sqrt(numDYS); %calculating SEM rather than std
YDesync_abs(:,2) = mean(DYSdesync_abs,2);
EDesync_abs(:,2) = std(DYSdesync_abs,0,2);
EDesync_abs(:,2) = EDesync_abs(:,2)/sqrt(numDYS); %calculating SEM rather than std
%% grab ET data
% if ET data needs to be analyzed, change next line to "ET=true." Otherwise
% ET=false;
ET=true;
if ET
etpath = uigetdir('', 'Select directory that contains ET _ecogPSD files to be analyzed');
etpath = [etpath '\'];
cd(etpath);

ETdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numET = length(ETdir);
counterET = 0;

    for i=1:numET
        filename = ETdir(i).name;
        load(filename);
        if isnan(order(k));
            continue; % skip files without contact pair over selected area
        end 
        % for evaluating freq bands other than those already represented in
        % subfreq matrix, rerun quantPSDmatrix function now, with new freq
        % bands as defined above
%         [subfreq] = quantDesync(rest,active,FREQ_QPSD,order,freq,filename);
        ETdesync(:,i) = ((subfreq(:,2,k)-subfreq(:,4,k))./subfreq(:,2,k))*100;
        ETdesync_abs(:,i) = ((subfreq(:,1,k)-subfreq(:,3,k))./subfreq(:,1,k))*100;
        counterET = counterET+1;
        clear subfreq
    end
end
YDesync(:,3) = mean(ETdesync,2);
EDesync(:,3) = std(ETdesync,0,2);
EDesync(:,3) = EDesync(:,3)/sqrt(numET); %calculating SEM rather than std
YDesync_abs(:,3) = mean(ETdesync_abs,2);
EDesync_abs(:,3) = std(ETdesync_abs,0,2);
EDesync_abs(:,3) = EDesync_abs(:,3)/sqrt(numET); %calculating SEM rather than std
%% Plot
hf=figure;

% subplot(2,1,1);
handlesr = barEB(-1*YDesync,EDesync);
set(handlesr.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YMinorGrid','on');
title([BRAINAREAS{k} ' Movement Related Dysynchronization']);
ylabel('% desynchronization');
hPD=handlesr.bars(1);
hDYS=handlesr.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
if ET
    legend(['PD (n=' num2str(counterPD) ')'],...
        ['Dys (n=' num2str(counterDYS) ')'],...
        ['ET (n=' num2str(counterET) ')']);
    hET=handlesr.bars(3);
    set(hET,'FaceColor','m');
else
%     legend(['PD (n=' num2str(counterPD) ')'],...
%         ['Dys (n=' num2str(counterDYS) ')']);
    legend(['PD'],...
        ['Dys']);
end

hf2=figure;

% subplot(2,1,1);
handlesr = barEB(-1*YDesync_abs,EDesync_abs);
set(handlesr.ca,'XTick',[1;2;3;4;5],...
    'XTickLabel',FREQBANDS,...
    'ylim',YLIM,...
    'YMinorGrid','on');
title([BRAINAREAS{k} ' Movement Related Dysynchronization - Absolute value']);
ylabel('% desynchronization');
hPD=handlesr.bars(1);
hDYS=handlesr.bars(2);
set(hPD,'FaceColor','b');
set(hDYS,'FaceColor','g');
if ET
    legend(['PD (n=' num2str(counterPD) ')'],...
        ['Dys (n=' num2str(counterDYS) ')'],...
        ['ET (n=' num2str(counterET) ')']);
    hET=handlesr.bars(3);
    set(hET,'FaceColor','m');
else
%     legend(['PD (n=' num2str(counterPD) ')'],...
%         ['Dys (n=' num2str(counterDYS) ')']);
    legend(['PD'],...
        ['Dys']);
end
function handles = barEB(barvals,errorvals)

% barEB plots bar graph of m-by-n matrix of bar values (by calling MATLAB
% function "bar") as well as error bars (by calling "errorbar")

% note: barvals and errorvals must be the same size!

% a list of handles are listed so that the user can change the properties
% of the plot
% handles.bars: handle to bar plot
% handles.gca: handle to axis

[numgroups numbars] = size(barvals); % # of groups, # bars in each group

handles.bars = bar(barvals);
hold on;
groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    handles.errors(i) = errorbar(x, barvals(:,i), errorvals(:,i), 'k', 'linestyle', 'none');
end

handles.ca = gca;
hold off;

% return;
%% Plot Low Beta Band Only
% hf=figure;
% 
% %Beta only
% YBetaDesync = YDesync(2,:);
% EBetaDesync = EDesync(2,:);
% YBetaDesync = YBetaDesync';
% EBetaDesync = EBetaDesync';
% 
% % subplot(2,1,1);
% handlesr = barEB(YBetaDesync,EBetaDesync);
% set(handlesr.ca,...
%     'XTickLabel',{'PD';'DYS';'ET'},...
%     'ylim',YLIM,...
%     'YMinorGrid','on');
% title([BRAINAREAS{k} ' Movement Related Dysynchronization']);
% ylabel('% desynchronization');
% hPD=handlesr.bars(1);
% hDYS=handlesr.bars(2);
% set(hPD,'FaceColor','b');
% set(hDYS,'FaceColor','g');
% if ET
%     legend(['PD (n=' num2str(counterPD) ')'],...
%         ['Dys (n=' num2str(counterDYS) ')'],...
%         ['ET (n=' num2str(counterET) ')']);
%     hET=handlesr.bars(3);
%     set(hET,'FaceColor','m');
% else
% %     legend(['PD (n=' num2str(counterPD) ')'],...
% %         ['Dys (n=' num2str(counterDYS) ')']);
%     legend(['PD'],...
%         ['Dys']);
% end

% subplot(2,1,2);
% handlesa=barEB(100*Yactive,100*Eactive);
% set(handlesa.ca,'XTick',[1;2;3;4;5],...
%     'XTickLabel',FREQBANDS,...
%     'ylim',YLIM,...
%     'YMinorGrid','on');
% title([BRAINAREAS{k} ' during ' LIMBAREAS{l} ' MOVEMENT']);
% ylabel('% power');
% hPD=handlesa.bars(1);
% hDYS=handlesa.bars(2);
% set(hPD,'FaceColor','b');
% set(hDYS,'FaceColor','g');
% if ET
%     hET=handlesa.bars(3);
%     set(hET,'FaceColor','m');
% elseif OTHER
%     hOTHER=handlesa.bars(3);
%     set(hOTHER,'FaceColor','m');
% end

function handles = barEB(barvals,errorvals)

% barEB plots bar graph of m-by-n matrix of bar values (by calling MATLAB
% function "bar") as well as error bars (by calling "errorbar")

% note: barvals and errorvals must be the same size!

% a list of handles are listed so that the user can change the properties
% of the plot
% handles.bars: handle to bar plot
% handles.gca: handle to axis

[numgroups numbars] = size(barvals); % # of groups, # bars in each group

handles.bars = bar(barvals);
hold on;
groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    handles.errors(i) = errorbar(x, barvals(:,i), errorvals(:,i), 'k', 'linestyle', 'none');
end

handles.ca = gca;
hold off;

% return;