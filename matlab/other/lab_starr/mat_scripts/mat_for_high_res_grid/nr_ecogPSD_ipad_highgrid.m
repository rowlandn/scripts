function ecogPSD_ipad_highgrid(name)
% this function calculates PSD using Welch's method to plot the frequency spectrum of files containing digitized
% ecog/LFP data.

% SS 10/2/2008: the montage method assumes that GL4k Ch.3-7 contain ecog data and Ch.8
% contains non-ecog (usually LFP data or noise)

% SS 11/24/2008: updated to perform quantitative analysis of PSD using
% sub-function 'quantPSD.' See 'quantPSD' in the sub-function section below
% for more details.

% ALC 6/9/09: updated to match the previously matrix-based code
% (ecogPSD4quant). Changes include:
%1. making the green beta band on the PSD plots 13-30 (previously 8-30Hz)
%2. asks user to define the subject's diagnosis for so output can be saved 
%   in a diagnosis-specific folder for later group analysis
%3. changes to 'quantPSD' subfunction - see below
%4. calls 'exportdata' function to write the quantitative data (allfreq and
%   subfreq matrices) to an excel file

% ALC 6/19/09: updated M1_ch variable due to change in apm7conv code

%ALC 9/7/2009: Noticed that the quantpsd subfunction is
%creating allfreq and subfreq matrices with four sheets (the third
%dimension). This is unnecessary, since we're only looking at three regions
%(M1, S1, STN). I did not change this, since it'll cause problems with the
%analysis now, but if we start analysis from scratch again, we can correct
%this to clear up any confusion in the future. 
%        
%% import ecog data
% % import filename and pathname of mat file containing ecog data created by APMconv7.m
[fn pn] = uigetfile('*ecog.mat','Select _ecog.mat containing ecog/lfp data');
cd(pn);
% load ecog data
load([pn fn]);
% % remove '_ecog.mat' ending from filename
% load(name)
%fn = name;
%fn = strrep(fn,'_ecog.mat','');

assignin('base','fn',fn)
assignin('base','pn',pn)



num_contact_pair = length(ecog.contact_pair);
num_epoch = length(ecog.rest_time);
% num_epoch = trials_ok;%length(ecog.touch_time); % this is how ecog_lfp was done


%% define variables

% EPOCH_LEN = 2;  

WINDOW = 1024;        
NOVERLAP = 512;                
NFFT = 1024;                   

EPOCH_LEN = 2*Fs;  
% 
% WINDOW = 1024*Fs/1000;        
% x = floor(log2(WINDOW*1/2));
% NOVERLAP = 2^x; %512*Fs/1000;                
% NFFT = 1024*Fs/1000;                   
XLIM_SPEC = [0 150];     
YLIM_LOG = [-3 3];      % ylim for for plotting in log scale
%---frequency ranges for graphing---
FREQ_LO = [8 30];       % beta band
%FREQ_MED = [35 57];     %low gamma band
FREQ_HI = [78 100];     % high gamma band
%---frequency ranges for quantPSD analysis---
% FREQ_QPSD below skips the frequency point at 21.5 and we will include
% that point in low beta for quantPSD analysis (SAS:11/3/09)

% SAS 4/29/10: changed frequency bandwidth so that they are contiguous (ie.
% 4-13,13-22,22-31,31-55,76-100).  Note:21.48Hz point is included in low
% beta, which gives 5 tota   l data points in low gamma at frequency
% resolution of 1.95 Hz.
FREQ_QPSD = [4 13;...     % delta alpha band
             13 30;...    %beta band
             13 22;...   % low beta band
             22 31;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band   

%% Determine GS4000 or AlphaOmega, set correction factor

% system_type = input('Recording system used (1=GS4000, 2=AlphaOmega): ');
system_type =2;
if system_type==1 % Guideline 4000 correction factor
    % correction factor based on 1.95Hz frequency resolution
    % Use correction factor below to only correct spectral amplitude for frequency points above 4Hz
    crxnfactor = [1,1,1,0.3951,0.4454,0.4926,...
        0.5888,0.6422,0.6906,0.7454,0.794,0.812,0.8234,...
        0.8318,0.8389,0.8467,0.8564,0.867,0.8776,0.8876,...
        0.8964,0.9036,0.9106,0.9174,0.9242,0.9309,0.9374,...
        0.9437,0.9498,0.9558,0.9614,0.9668,0.9719,0.9767,...
        0.9811,0.9852,0.9888,0.992,0.9948,0.997,0.9988];
    x=repmat(crxnfactor,Fs/1000);
    y = [];
    for i = size(x,1)
        y=[y;x(i,:)];
    end
    y=sort(y);
    n = length(find(y==1));
    crxnfactor = [ones(1,n) y(1:end-n)];
    
elseif system_type==2 % AlphaOmega correction factor
    % correction factor based on 1.95Hz resolution
    % NOTE: there was no attenuation between 4-100Hz in AlphaOmega system.
    % Use correction factor below to correct spectral amplitude for
    % frequency points in 1.95-80Hz range
    % Use correction factor below to only correct spectral amplitude for
    % frequency points above 4 Hz.
    crxnfactor = ones(1,41);
end
%% Define disease state
%Will output ecogPSD data into folders specific to disease state
% Dx = input('Enter patient diagnosis (1=PD, 2=Dys, 3=ET, 4=Epilepsy, 5=Other): ');
Dx =1;

%% Define period to compute psd
% initialize structures that will contain all rest/active analysis
rest = struct('contact_pair',{});
prep = struct('contact_pair',{});
active = struct('contact_pair',{});

for i = 1:num_contact_pair
    for j = 1:num_epoch
        % rest
        if j==1
            start_rest = int32(ecog.ipad_ON_time(j));
        else
            start_rest = int32(ecog.rest_time(j));
        end
        end_rest = int32(ecog.prep_time(j));
%         t = start_rest + (end_rest - start_rest)/2;% psd computed using 2s of rest in the middle of that period
%         rest(1).contact_pair(i).epoch(j).remontaged_ecog_signal = ecog.contact_pair(i).remontaged_ecog_signal(t-EPOCH_LEN:t+EPOCH_LEN);
        rest(1).contact_pair(i).epoch(j).remontaged_ecog_signal = ecog.contact_pair(i).remontaged_ecog_signal(end_rest-EPOCH_LEN:end_rest);
        rest(1).contact_pair(i).epoch(j).time = [end_rest-EPOCH_LEN end_rest];
        % prep
        start_prep = int32(ecog.prep_time(j)); 
        end_prep = int32(ecog.go_time(j));
        t = start_prep + (end_prep - start_prep)/2;% psd computed using 2s of rest in the middle of that period
        prep(1).contact_pair(i).epoch(j).remontaged_ecog_signal = ecog.contact_pair(i).remontaged_ecog_signal(t-EPOCH_LEN/2:t+EPOCH_LEN/2);
        prep(1).contact_pair(i).epoch(j).time = [t-EPOCH_LEN/2 t+EPOCH_LEN/2];
        % active
        start_active = int32(ecog.active_time(j)); 
%         end_active = int32(ecog.end_time(j));
%         t = start_active + (end_active - start_active)/2;% psd computed using 2s of rest in the middle of that period
%         active(1).contact_pair(i).epoch(j).remontaged_ecog_signal = ecog.contact_pair(i).remontaged_ecog_signal(t-EPOCH_LEN:t+EPOCH_LEN);
        active(1).contact_pair(i).epoch(j).remontaged_ecog_signal = ecog.contact_pair(i).remontaged_ecog_signal(start_active:start_active+EPOCH_LEN);
        active(1).contact_pair(i).epoch(j).time = [start_active start_active+EPOCH_LEN];
    end
end


%% calculate PSD using pwelch
% perform PSD analysis using pwelch, average all the epochs for
% rest/active

for i = 1:num_contact_pair
    for j = trials_ok  %perform fft for each epoch

        [rest(1).contact_pair(i).epoch(j).PSD, rest(1).contact_pair(i).epoch(j).freq] =...
            pwelch(rest.contact_pair(i).epoch(j).remontaged_ecog_signal, WINDOW, NOVERLAP, NFFT, Fs);
        
        [prep(1).contact_pair(i).epoch(j).PSD, prep(1).contact_pair(i).epoch(j).freq] =...
            pwelch(prep.contact_pair(i).epoch(j).remontaged_ecog_signal, WINDOW, NOVERLAP, NFFT, Fs);
        
        [active(1).contact_pair(i).epoch(j).PSD , active(1).contact_pair(i).epoch(j).freq] =...
            pwelch(active.contact_pair(i).epoch(j).remontaged_ecog_signal,WINDOW,NOVERLAP,NFFT,Fs);
    end

    % calculate mean PSD across all epochs with the subfunction
    % 'mean_signal'
    rest(1).contact_pair(i).mean_PSD = mean_signal(rest.contact_pair(i),'PSD');
    prep(1).contact_pair(i).mean_PSD = mean_signal(prep.contact_pair(i),'PSD');
    active(1).contact_pair(i).mean_PSD = mean_signal(active.contact_pair(i),'PSD');

%     % Use for-loop below to correct for frequency response
%     for k=1:length(crxnfactor)
%         rest(1).contact_pair(i).mean_PSD(k) = rest(1).contact_pair(i).mean_PSD(k)/(crxnfactor(k)^2);
%         prep(1).contact_pair(i).mean_PSD(k) = prep(1).contact_pair(i).mean_PSD(k)/(crxnfactor(k)^2);
%         active(1).contact_pair(i).mean_PSD(k) =  active(1).contact_pair(i).mean_PSD(k)/(crxnfactor(k)^2);
%     end
    
    % SAS 11/24/09: use norm_idx to normalize by max power between 8-100Hz, instead of
    % the whole frequency spectrum
    norm_idx=find(rest.contact_pair(1).epoch(trials_ok(1)).freq>4 & rest.contact_pair(1).epoch(trials_ok(1)).freq<100);
   
    % normalize mean PSD to peak rest height
    rest(1).contact_pair(i).norm_mean_PSD = rest.contact_pair(i).mean_PSD/max(rest.contact_pair(i).mean_PSD(norm_idx(1):norm_idx(end)));
    prep(1).contact_pair(i).norm_mean_PSD = prep.contact_pair(i).mean_PSD/max(rest.contact_pair(i).mean_PSD(norm_idx(1):norm_idx(end)));
    active(1).contact_pair(i).norm_mean_PSD = active.contact_pair(i).mean_PSD/max(rest.contact_pair(i).mean_PSD(norm_idx(1):norm_idx(end)));
   
    % calculate log base 10 of mean PSD
    rest(1).contact_pair(i).log_mean_PSD = log10(rest.contact_pair(i).mean_PSD);
    prep(1).contact_pair(i).log_mean_PSD = log10(prep.contact_pair(i).mean_PSD);
    active(1).contact_pair(i).log_mean_PSD = log10(active.contact_pair(i).mean_PSD);
   
    %normalize log mean PSD to peak rest height
    rest(1).contact_pair(i).norm_log_mean_PSD = rest.contact_pair(i).log_mean_PSD/max(rest.contact_pair(i).log_mean_PSD(norm_idx(1):norm_idx(end)));
    prep(1).contact_pair(i).norm_log_mean_PSD = prep.contact_pair(i).log_mean_PSD/max(rest.contact_pair(i).log_mean_PSD(norm_idx(1):norm_idx(end)));
    active(1).contact_pair(i).norm_log_mean_PSD =  active.contact_pair(i).log_mean_PSD/max(rest.contact_pair(i).log_mean_PSD(norm_idx(1):norm_idx(end)));
   
    PSD_rest = psd_all(rest.contact_pair(i),'PSD');
    PSD_prep = psd_all(prep.contact_pair(i),'PSD');
    PSD_active = psd_all(active.contact_pair(i),'PSD');
    
    % calculate difference in log mean PSD between active and rest
    num_epoch = length(trials_ok);
    active.contact_pair(i).difference = active.contact_pair(i).log_mean_PSD - rest.contact_pair(i).log_mean_PSD;
    active.contact_pair(i).difference_sig = evRanksum([PSD_rest;PSD_active]', 1:num_epoch, num_epoch + 1: num_epoch+ num_epoch, 0.05);
    
    % calculate difference in log mean PSD between active and rest
    prep.contact_pair(i).difference = active.contact_pair(i).log_mean_PSD - prep.contact_pair(i).log_mean_PSD;
    prep.contact_pair(i).difference_sig = evRanksum([PSD_prep;PSD_active]', 1:num_epoch, num_epoch + 1: num_epoch+ num_epoch, 0.05);
    
%     active.contact_pair(i).difference =...
%         active.contact_pair(i).norm_log_mean_PSD - rest.contact_pair(i).norm_log_mean_PSD;

    % calculate ratio of active to rest non-normalized mean PSD
%     active.contact_pair(i).ratio =...
%         100*(active.contact_pair(i).mean_PSD./rest.contact_pair(i).mean_PSD);
end

freq = rest(1).contact_pair(1).epoch(trials_ok(1)).freq;

%% quantitative analysis of PSD
% now that all the data crunching is done and the results are stored in
% 'rest' and 'active' structure arrays, run sub-function quantPSD for
% quantitative analysis of PSD data
% [allfreq subfreq order] = quantPSD(rest,prep,active,FREQ_QPSD,freq,fn,M1_ch);

%% Save and write quantPSD data as excel spreadsheet
% If the statistical analysis is repeated measures, use the following line
% of code: 
% [data4export] = exportdata3x31(allfreq, subfreq, fn);

% If not using repeated measures, use this code: 
% [data4export] = exportdata(allfreq, subfreq, fn);

% cd(pn); %exportdata changes the path in order to write the excel data. This brings back current path.

%% Save and write quantPSD data 
% Creates .mat output file of name (filename_ecogPSD.mat). transcoh.mat
% files can then be used to analyze coherence data across groups
% if Dx==1
%     outputdir = ['C:\Cora\raw data\ECOG data\quantPSD data\PD\'];
% elseif Dx==2
%     outputdir = ['C:\Cora\raw data\ECOG data\quantPSD data\DYS\'];
% elseif Dx==3
%     outputdir = ['C:\Cora\raw data\ECOG data\quantPSD data\ET\'];
% elseif Dx==4
%     outputdir = ['C:\Cora\raw data\ECOG data\quantPSD data\Epilepsy\'];
% else
%     outputdir = ['C:\Cora\raw data\ECOG data\quantPSD data\Other\'];
% end
% outputname = [outputdir fn,'_PSD.mat'];
% disp(['Writing ecog PSD data to:  ' outputname]);
% save(outputname,'ecogPSD');
% save(outputname,'order','rest', 'active', 'freq', 'allfreq', 'subfreq');
save(fn,'rest','active','prep', 'freq','-append');

assignin('base','rest',rest)
assignin('base','prep',prep)
assignin('base','active',active)

%% create figure
XLIM_SPEC1 = [0 150];     % range of spectral frequency for plotting (zoomed out)
XLIM_SPEC2 = [0 50];% range of spectral frequency for plotting (zoomed in)
x= max(max(PSD_rest));
necog=length(ecog.contact_pair);
YLIM_LOG = [-x x];  
if necog ==5
   length_row = 5;
elseif necog ==6
   length_row = 6;
elseif necog <=32
    length_row = 14;
    M1_ch=M1_ch1;
elseif necog >= 64;
    length_row = 16;
    M1_ch=M1_ch1;
end

% % figure plot nrom PSD 
% row=1;
% for ii = 1:length_row:necog-1    
%     hf1 = figure;
%     for i = 1:length_row
%         chan=ii+i-1;
%         % 1st row subplot
%         subplot(4,4,i)
%         ha = gca;
%         hold(ha,'on');
%         plot(freq, rest.contact_pair(chan).norm_mean_PSD,freq, prep.contact_pair(chan).norm_mean_PSD,freq, active.contact_pair(chan).norm_mean_PSD,'r','LineWidth',2);
%         if i==1;
%             title([fn sprintf('\n') 'e' num2str(chan) ]); % allows title to have file name
%             xlabel('frequency (Hz)');
%             ylabel('norm PSD');
%         elseif i == M1_ch
%             title(['e' num2str(chan) ],...
%                 'FontWeight','b');
%         else
%             title(['e' num2str(chan) ]);
%         end
%         set(ha,'YLim',[0 1.3]);
%     set(ha,'XLim',[0 150]);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
%     end
%     row = row;
%     %saveas(hf1,[fn '_normpsd' num2str(ii)],'fig');
%     saveas(hf1,[fn(1:11),'fig_psdn_ecg_',fn(20:end-4)],'fig');
%     print(hf1,[fn(1:11),'pdf_psdn_ecg_',fn(20:end-4)],'-dpdf');
%     row = row+1;
% end

% figure plot log PSD 
% clim(1)=max(rest.contact_pair(chan).log_mean_PSD);
% clim(2)=min(rest.contact_pair(chan).log_mean_PSD);
row=1;
for ii = 1:length_row:necog-1    
    hf1 = figure;
    for i = 1:length_row
        chan=ii+i-1;
        % 1st row subplot
        subplot(4,4,i)
        ha = gca;
        hold(ha,'on');
        plot(freq, rest.contact_pair(chan).log_mean_PSD,freq, prep.contact_pair(chan).log_mean_PSD,freq, active.contact_pair(chan).log_mean_PSD,'r','LineWidth',2);
        if i==1;
            title([fn sprintf('\n') 'e' num2str(chan) ]); % allows title to have file name
            xlabel('frequency (Hz)');
            ylabel('log PSD');
        elseif i == M1_ch
            title(['e' num2str(chan) ],...
                'FontWeight','b');
        else
            title(['e' num2str(chan) ]);
        end
        set(ha,'YLim',[-3 3]);
    set(ha,'XLim',[0 150]);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');
    end
    row = row;
    %saveas(hf1,[fn '_logpsd' num2str(ii)],'fig');
    saveas(hf1,[fn(1:11),'fig_psdl_ecg_',fn(20:end-4)],'fig');
    print(hf1,[fn(1:11),'pdf_psdl_ecg_',fn(20:end-4)],'-dpdf');
    row = row+1;
end

% hf1 = figure;
% for i = 1:14
%     
%     subplot(3,5,i);
%     ha = gca;       % create handle for current axes
%     hold(ha,'on');
%     plot(freq, rest.contact_pair(i).norm_mean_PSD,freq, prep.contact_pair(i).norm_mean_PSD,freq, active.contact_pair(i).norm_mean_PSD,'r','LineWidth',2);
%     if i == 1
%         title([fn sprintf('\n') 'e' num2str(i)]);
%         ylabel(['normalized average PSD' sprintf('\n') 'across all epochs'],'FontSize',8);
%         xlabel('frequency (Hz)');
%         hl = legend('Rest','Prep','Active');
%         set(hl,'FontSize',5,'Box','off');
%     elseif i==M1_ch1
%             title(['e' num2str(i)],'FontWeight','b');   % if i == num_contact_pair
%     elseif i == 29 % assume 6th contact_pair contains LFP data
%             title('LFP');
%     else
%             title(['e' num2str(i) ]);   
%     end
% 
%     set(ha,'YLim',[0 1.3]);
%     set(ha,'XLim',[0 50]);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
% 
% end
% 
% % save figure
% saveas(hf1,[fn '_PSD_raw1'],'fig');
% 
% hf1 = figure;
% for i = 15:num_contact_pair
%     j=i-14;
%     subplot(3,5,j);
%     ha = gca;       % create handle for current axes
%     hold(ha,'on');
%     plot(freq, rest.contact_pair(i).norm_mean_PSD,freq, prep.contact_pair(i).norm_mean_PSD,freq, active.contact_pair(i).norm_mean_PSD,'r','LineWidth',2);
%     if i == 1
%         title([fn sprintf('\n') 'e' num2str(i)]);
%         ylabel(['normalized average PSD' sprintf('\n') 'across all epochs'],'FontSize',8);
%         xlabel('frequency (Hz)');
%         hl = legend('Rest','Prep','Active');
%         set(hl,'FontSize',5,'Box','off');
%     elseif i==M1_ch1
%             title(['e' num2str(i)],'FontWeight','b');   % if i == num_contact_pair
%     elseif i == 29 % assume 6th contact_pair contains LFP data
%             title('LFP');
%     else
%             title(['e' num2str(i)]);   
%     end
% 
%     set(ha,'YLim',[0 1.3]);
%     set(ha,'XLim',[0 50]);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     hold(ha,'off');
% 
% end
% 
% % save figure
% saveas(hf1,[fn '_PSD_raw2'],'fig');
%% create optional figures
% %this plot visualizes data in allfreq
% hf2 = figure;
% areas = {'pre-motor' 'M1' 'M1-S1' 'LFP'};
% freqbands = {'4-12' '13-21' '22-30' '31-55' '76-100'};
% for i=1:4
%     subplot(4,1,i);
%     if isnan(order(i));
%         text(0.5,0.5,...
%             [areas{i} ' not available'],...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle');
%         continue;
%     end
%     ha = gca;
%     hold(ha,'on');
%     plot(freq, rest.contact_pair(order(i)).mean_PSD,...
%         freq, active.contact_pair(order(i)).mean_PSD,'r','LineWidth',2);
%     plot(allfreq(1,1,i),allfreq(1,2,i),'o');
%     text(allfreq(1,1,i),allfreq(1,2,i),['Max Freq: ' num2str(allfreq(1,1,i)) sprintf('\n')...
%         'Power Value: ' num2str(allfreq(1,2,i))]);
%     plot(allfreq(2,1,i),allfreq(2,2,i),'or');
%     text(allfreq(2,1,i),allfreq(2,2,i),['Max Freq: ' num2str(allfreq(2,1,i)) sprintf('\n')...
%         'Power Value: ' num2str(allfreq(2,2,i))]);
%
%     set(ha,'XLim',[0 50]);
%     FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
%     title(areas{i});
%     hold(ha,'off');
% end
%
% %this plot visualizes data in subfreq
% hf3 = figure;
% x = [1 2 3 4 5];
% for i=1:4
%     subplot(4,1,i);
%     if isnan(order(i));
%         text(0.5,0.5,...
%             [areas{i} ' not available'],...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','middle');
%         continue;
%     end
%     ha=gca;
%     hold(ha,'on');
%     Y(:,1)=subfreq(:,2,i);
%     Y(:,2)=subfreq(:,4,i);
%     bar(x,Y);
%     set(ha,'XTick',[1;2;3;4;5]);
%     set(ha,'XTickLabel',freqbands);
%     title(areas{i});
%     hold(ha,'off');
% end


%% sub-functions below
%------------------------------
function out = mean_signal(contact_pair, fieldname)
num_epoch = length(contact_pair.epoch);
A = [];
for j = 1:num_epoch
    tmp = eval(['contact_pair.epoch(' num2str(j) ').' fieldname])';
    A = [A;tmp]; %#ok<AGROW> % concatenate
end
out = mean(A);
return;

function out = psd_all(contact_pair, fieldname)
num_epoch = length(contact_pair.epoch);
A = [];
for j = 1:num_epoch
    tmp = eval(['contact_pair.epoch(' num2str(j) ').' fieldname])';
    A = [A;tmp]; %#ok<AGROW> % concatenate
end
out = A;
return;


function FFT_band_fill(ha,FREQ_HI,FREQ_LO)
% FFT_band_fill fills bands of specified frequncies on axes specified by
% the handle ha

ylm = get(ha,'YLim');
x_lo = [FREQ_LO(1) FREQ_LO(2) FREQ_LO(2) FREQ_LO(1)];
x_hi = [FREQ_HI(1) FREQ_HI(2) FREQ_HI(2) FREQ_HI(1)];
y = [ylm(1) ylm(1) ylm(2) ylm(2)];
fill(x_lo,y,'g',...
    'EdgeColor','none',...
    'FaceAlpha',0.3);
fill(x_hi,y,'y',...
    'EdgeColor','none',...
    'FaceAlpha',0.3);
return;

%------------------------------
function [allfreq subfreq order] = quantPSD(rest,prep,active,FREQ_QPSD,freq,filename,M1_ch)
% quantPSD performs quantitative analysis on rest/active structure
% arrays with the given frequency ranges.  allows user to define ecog
% contact closest to M1
%*UPDATE ALC 6/9/09: no longer allows user to define M1 contact, as this is
%occurs in main function
%
% output: one .mat file with two 3D matrices and # ecog contact that is
% closest to M1.

% select ecog contact closest to M1 and use that to reassign ecog contact
% data array
% ncontact = {'1' '2' '3' '4' '5' '6'};
% contactm1 = menu('Select ecog contact closest to M1',ncontact);
% assign contact numbers for each structure relative to m1 contact
m1 = M1_ch;
% pre = contactm1+1;
% m1s1 = contactm1-2; % updated 2/13:from M1 contact #, subtract 2 instead of 1 to find S1
s1 = m1 - 2; %updated 6/9/09 for clarity

if  m1==1
    s1 = NaN;
    menu(['There is no contact pair over S1.' sprintf('\n')...
        'Click OK to continue'],'OK');
elseif m1==2
    s1 = NaN;
    menu(['There is no contact pair over S1.' sprintf('\n')...
        'Click OK to continue'],'OK');
% elseif contactm1==5
%     pre = NaN;
%     menu(['There is no contact pair over premotor.' sprintf('\n')...
%         'Click OK to continue'],'OK');
elseif m1==6
    m1 = NaN;
%     pre = NaN;
    menu(['There is no contact pair over M1 or premotor.' sprintf('\n')...
        'Click OK to continue'],'OK');
end

% look for LFP channel
num_contact_pair = length(rest.contact_pair);
if num_contact_pair == 6
    lfp = 6; % assume that 6th contact_pair in rest/active structures always contains LFP data
elseif num_contact_pair<6
    lfp = NaN;
end

order = [m1 s1 lfp];

% initialize and populate the 2x3x4 matrix 'allfreq'
% allfreq is a 2x3x4 3-D matrix in which the four 2x3 arrays correspond to the
% 4 data recording channels (premotor ecog, M1 ecog, M1-S1 ecog, and stn
% lfp).
% The 3 rows contain:
%   1st row     -   resting state data
%   2nd row     -   active movement data
%   3rd row     -   preparation data
% The 3 columns contain:
%   1st col     -   frequency at which max power occurs
%   2nd col     -   max power value
%   3rd col     -   total power across all 5 frequency bands

allfreq = zeros(3,3,4);
for i=1:3
    if isnan(order(i))
        allfreq(:,:,i) = NaN;
        continue;
    end
    rest_data = rest.contact_pair(order(i)).mean_PSD;
    prep_data = prep.contact_pair(order(i)).mean_PSD;
    active_data = active.contact_pair(order(i)).mean_PSD;
    
% 6/9/09: Need to limit rest and active data to those in freq range of interest (typically 4Hz-100Hz)
% 6/19/09: Need to further limit to exclude the frequencies between the bins (ie: exclude 60Hz)    
    rest_data_tmp = rest_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2)); %assumes we are looking at 5 freq bands, total
    prep_data_tmp = prep_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2)); %assumes we are looking at 5 freq bands, total
    active_data_tmp = active_data(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));
    freq_tmp = freq(freq>FREQ_QPSD(1,1) & freq<FREQ_QPSD(5,2));    

    [c1 i1] = max(rest_data_tmp);%finds max value of rest data and that value's index position
    allfreq(1,1,i)=freq_tmp(i1);%puts freq corresponding to max value into allfreq matrix
    allfreq(1,2,i)=c1;%puts max value into allfreq matrix

    [c2 i2] = max(prep_data_tmp);
    allfreq(2,1,i)=freq_tmp(i2);
    allfreq(2,2,i)=c2;
    
    [c3 i3] = max(active_data_tmp);
    allfreq(3,1,i)=freq_tmp(i3);
    allfreq(3,2,i)=c3;

    array1 = [];
    array2 = [];
    array3 = [];
 
    for j = 1:size(FREQ_QPSD,1)
        tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        tmp2 = prep_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        tmp3 = active_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));

        array1 = [array1 tmp1];  %#ok<AGROW>
        array2 = [array2 tmp2]; %#ok<AGROW>
        array3 = [array3 tmp3]; %#ok<AGROW>

    end
    allfreq(1,3,i)=sum(array1);
    allfreq(2,3,i)=sum(array2);
    allfreq(3,3,i)=sum(array3);
end

% initialize and populate the 5x5x4 matrix 'subfreq'
% The four 5x5 arrays of subfreq correspond to the 4 data channels
% (pre-motor,M1,M1-S1,STN LFP).  Each of the 5 rows correspond to the 5
% frequency bands defined by variable FREQ_QPSD.
% The 5 columns contain:
%   1st col     -   total power in given frequency band at rest
%   2nd col     -   power in given frequency band, divided by total power
%                   in all 5 freq bands, at rest
%   3rd col     -   total power in given frequency band during movement
%   4th col     -   power in given frequency band, divided by total power
%                   in all 5 freq bands, during movement
%   5th col     -   power during movement divided by power at rest (col 3
%                   divided by col 1)
subfreq = zeros(5,5,4);
for i=1:3
    if isnan(order(i))
        subfreq(:,:,i)=NaN;
        continue;
    end
    rest_data = rest.contact_pair(order(i)).mean_PSD;
    prep_data = prep.contact_pair(order(i)).mean_PSD;
    active_data = active.contact_pair(order(i)).mean_PSD;
    for j=1:size(FREQ_QPSD,1)
        tmp1 = rest_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        tmp2 = prep_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        tmp3 = active_data(freq>FREQ_QPSD(j,1) & freq<FREQ_QPSD(j,2));
        subfreq(j,1,i) = sum(tmp1);
        subfreq(j,2,i) = sum(tmp1)/allfreq(1,3,i);
        subfreq(j,3,i) = sum(tmp2);
        subfreq(j,4,i) = sum(tmp2)/allfreq(2,3,i);
        subfreq(j,5,i) = sum(tmp2)/sum(tmp1);
        subfreq(j,6,i) = sum(tmp3);
        subfreq(j,7,i) = sum(tmp3)/allfreq(2,3,i);
        subfreq(j,8,i) = sum(tmp3)/sum(tmp3);
    end
end
save([filename '_PSD'],'allfreq','subfreq','order');
return;