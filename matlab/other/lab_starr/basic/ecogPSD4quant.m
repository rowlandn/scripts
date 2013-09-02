%% ecogPSD4quant
% this function calculates PSD using Welch's method to plot the frequency spectrum of files containing digitized
% ecog/LFP data.

% SS 10/2/2008: the montage method assumes that GL4k Ch.3-7 contain ecog data and Ch.8
% contains non-ecog (usually LFP data or noise)

% SS 11/24/2008: updated to perform quantitative analysis of PSD using
% sub-function 'quantPSD.' See 'quantPSD' in the sub-function section below
% for more details.

% AC 02/18/09: ecogPSD4quant a modification of Sho's ecogPSD function
% created to work on _ecog_lfp.mat files for Andrea

%AC 04/06/09: updated to call exportdata function to write data from
%allfreq and subfreq tables to excel automatically. 
%% define variables       
% Fs=1000;                % samples per unit time, in this case digitization at 1000 hz rate
% Fs = 1.502403855323792e3;   % Alpha Omega ecog/lfp are recorded using different sampling rate. 
% note: in the future, consider downsampling AO data to make it as close to
% GL4k data as possible
OFFSET = 0.5;           % epoch time offset in seconds
EPOCH_LEN = 2;          % epoch length in seconds (1.792 is the length of 6 overlapping segments)
WINDOW = 512;           % segment length and Hamming window length for welch's method
NOVERLAP = 256;         % # signal samples that are common to adjacent segments for welch's method
NFFT = 512;             % length of fft
XLIM_SPEC = [0 150];     % range of spectral frequency for plotting
YLIM_LOG = [-2 2.5];      % ylim for for plotting in log scale
%---frequency ranges for graphing---
FREQ_LO = [13 30];       % beta band --changed beta low bound from 8 to 18: ALC 5/13/09
%FREQ_MED = [35 57];     %low gamma band
FREQ_HI = [78 100];     % high gamma band
%---frequency ranges for quantPSD analysis---
FREQ_QPSD = [4 12;...     % delta alpha band
             13 21;...   % low beta band
             22 30;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band   
%---variables for plotting bar graph---
FREQBANDS = {'4-12' '13-21' '22-30' '31-55' '76-100'};
% FREQBANDS = {num2str(FREQ_QPSD(1,:)), num2str(FREQ_QPSD(2,:)) num2str(FREQ_QPSD(3,:)) ...
%         num2str(FREQ_QPSD(4,:)) num2str(FREQ_QPSD(5,:))};
% AREAS = {'pre-motor' 'M1' 'S1' 'STN LFP'};
BRAINAREAS = {'M1' 'S1' 'STN LFP'};
LIMBAREAS = {'HAND' 'ELBOW' 'SHOULDER' 'JAW' 'FOOT' '"ARM"' '"NON-ARM"'};
YLIM = [0 70];
%% import ecog data
% import filename and pathname of mat file containing ecog data created by APMconv7_coh.m
[fn pn] = uigetfile('*_ecog_lfp.mat','Select .mat containing ecog/lfp data');
cd(pn);
% load ecog data
load([pn fn]);
% remove '_ecog.mat' ending from filename
fn = strrep(fn,'_ecog_lfp.mat','');
filename = fn; % exportdata function (called later) requires a "filename" variable
%% Define OR site
%Will choose appropriate Fs based on site/recording system
OR = input('Enter OR site: 1=UCSF, 2=VA : ');
if OR==1
    Fs=1000;
elseif OR==2
    Fs = 1.502403855323792e3;
end

%% Define disease state
%Will output coh data into folders specific to disease state
Dx = input('Enter patient diagnosis: 1=PD, 2=Dys, 3=ET, 4=Epilepsy, 5=Other : ');

%% Remontage Data
%such that contact 1 is referencing contact 2, 2 references 3, etc

ecog1v2 = ecog_lfp_raw_data(:,1) - ecog_lfp_raw_data(:,2);
ecog2v3 = ecog_lfp_raw_data(:,2) - ecog_lfp_raw_data(:,3);
ecog3v4 = ecog_lfp_raw_data(:,3) - ecog_lfp_raw_data(:,4);
ecog4v5 = ecog_lfp_raw_data(:,4) - ecog_lfp_raw_data(:,5);
ecog5v6 = ecog_lfp_raw_data(:,5);
LFP = ecog_lfp_raw_data(:,end-2);
% %% Downsampling for Miller case (accidentally recorded at 5K)
% ecog1v2 = downsample(ecog1v2,5);
% ecog2v3 = downsample(ecog2v3,5);
% ecog3v4 = downsample(ecog3v4,5);
% ecog4v5 = downsample(ecog4v5,5);
% ecog5v6 = downsample(ecog5v6,5);
% LFP = downsample(LFP,5);

recog = [ecog1v2 ecog2v3 ecog3v4 ecog4v5 ecog5v6 LFP];
%% Round timestamps to 3 decimal places
%AlphaOmega system retains 4 decimal places, which causes rounding errors
%when parsing rest/active epochs

time_stamps = ecog_lfp_raw_data(:,end-1:end);
time_stamps = time_stamps*1000;
time_stamps = round(time_stamps);
time_stamps = time_stamps/1000;
ecog_lfp_raw_data(:,end-1:end) = time_stamps; 

%% Parse data
% num_contact_pair = length(ecog.contact_pair);
% num_epoch = length(ecog.rest_time);

[num_row num_col]=size(ecog_lfp_raw_data); %Need to keep raw data matrix here b/c has rest and active onsets
num_contact_pair = num_col-2; %last two columns are rest and active onsets
rest_time = ecog_lfp_raw_data(:,end-1);
active_time = ecog_lfp_raw_data(:,end);
diff = rest_time - active_time;
rest_time = rest_time(diff ~= 0);
active_time = active_time(diff~=0);
if active_time(end)== 0 %usually there will be one more rest onset than active onset, in which case we will throw out any extra zeros from active_time vector
    active_time = active_time(1:end-1);
end
num_epoch = length(active_time);

if num_epoch<5
    menu(['There are only ' num2str(num_epoch) ' epochs.' sprintf('\n')...
        'Click OK to continue'],'OK');
end
% Create empty 2D matrix where data from each epoch can be
% taken from ecog_lfp_raw_data and condensed into a single column for each
% contact pair without adding excess zeros, based on knowledge of the
% length of each epoch

% epoch_length_r = zeros(1,num_epoch); %vector of length = number of epochs

%Determine length of each epoch(i) and  replace zeros in above epoch_lenghts array with length of each epoch 
% for i = 1:num_epoch
%         epoch_length_r(i) = int32(((active_time(i) - rest_time(i))- OFFSET) * Fs); 
% end 
% tot_epoch_length_r = int32(sum(epoch_length_r)); %adds all above epoch lengths together; purpose of int32 is to ensure this variable is read as an integer (has been a problem)
% all_rest = zeros(tot_epoch_length_r,num_contact_pair); %creates zeros matrix: #rows= sum of all rest epochs, #columns = #contact pairs 

% Now the same for Active epochs

% epoch_length_a = zeros(1,num_epoch);
% 
% for i = 1:num_epoch
%         epoch_length_a(i) = int32(((rest_time(i+1) - active_time(i))- OFFSET) * Fs); 
% end
% tot_epoch_length_a = int32(sum(epoch_length_a));
% all_active = zeros([tot_epoch_length_a num_contact_pair]);
    
    %4/9/09: Now limiting the epoch length to 2s, still need to create all_rest and
    %all_active matrices
all_rest = zeros([(EPOCH_LEN*Fs)*num_epoch num_contact_pair]);
all_active = zeros([(EPOCH_LEN*Fs)*num_epoch num_contact_pair]);

% Now populate those empty matrices with "raw" time domain data that is
% segregated by activity state (rest vs movement)

for i = 1:num_contact_pair
   
    % Parse each rest/active epoch from each contact pair
    
    for j = 1:num_epoch
        
    % Determine starting point and end point of each epoch, allowing for a
    %   time offset, if needed. The start and end points must refer to a
    %   row value from the raw data matrix, thus for the first rest epoch
    %   (starts at time 0) must make alternate index in case OFFSET = 0
        start_rest(j) = int32((rest_time(j) + OFFSET) * Fs);      
        if start_rest(j) == 0
            start_rest(j) = 1; 
        end
%         end_rest(j) = int32(active_time(j) * Fs - 1); %substract 1 because end of rest period is the datum just before the start of the active period that follows
        end_rest(j) = int32((start_rest(j) + (Fs * EPOCH_LEN)) - 1);
    % Fill in zeros matrix with appropriate data
        if j==1
%             all_rest(1:epoch_length_r(j),i) = ecog_lfp_raw_data(start_rest(j):end_rest(j),i);
            all_rest(1:int32(EPOCH_LEN*Fs),i) = recog(start_rest(j):end_rest(j),i);
            epoch_starts_r(j)= 1;
            epoch_ends_r(j) = int32(EPOCH_LEN*Fs); %epoch_starts and _ends is a way to check that your indexing is working properly
        else
%             next_epoch_r = int32(sum(epoch_length_r(1:(j-1))) + 1); %to determine where to place real data from epochs beyond the first one, have to tell how far down column to start later epochs
%             next_end_r = int32((next_epoch_r - 1) + epoch_length_r(j)); %add new epoch length to the previous epoch(s) length(s)
%             all_rest(next_epoch_r:next_end_r,i) = ecog_lfp_raw_data(start_rest(j):end_rest(j),i);
            next_epoch_r = int32(((EPOCH_LEN*Fs)*(j-1))+1);
            next_end_r = int32((next_epoch_r - 1)+(EPOCH_LEN*Fs));
            all_rest(next_epoch_r:next_end_r,i) = recog(start_rest(j):end_rest(j),i);
            epoch_starts_r(j) = next_epoch_r; %fills epoch start times from all_rest matrix into the epoch_starts matrix for use in coherence section
            epoch_ends_r(j) = next_end_r;
        end
       
        start_active(j) = int32((active_time(j)+ OFFSET) * Fs);
        end_active(j) = int32(start_active(j) + (Fs * EPOCH_LEN) - 1);
        
        if j==1
%             all_active(1:epoch_length_a(j),i) = ecog_lfp_raw_data(start_active(j):end_active(j),i);
            all_active(1:int32(EPOCH_LEN*Fs),i) = recog(start_active(j):end_active(j),i);
            epoch_starts_a(j) = 1;
            epoch_ends_a(j) = int32(EPOCH_LEN*Fs);
        else
%             next_epoch_a = int32(sum(epoch_length_a(1:(j-1))) + 1); 
%             next_end_a = int32((next_epoch_a - 1) + epoch_length_a(j));
%             all_active(next_epoch_a:next_end_a,i) = ecog_lfp_raw_data(start_active(j):end_active(j),i);
            next_epoch_a = int32(((EPOCH_LEN*Fs)*(j-1))+1);
            next_end_a = int32((next_epoch_a - 1)+(EPOCH_LEN*Fs));
            all_active(next_epoch_a:next_end_a,i) = recog(start_active(j):end_active(j),i);
            epoch_starts_a(j) = next_epoch_a;
            epoch_ends_a(j) = next_end_a;
        end
       
    end
    
end

%% Calculate mean PSD
% perform PSD analysis using pwelch, average all the epochs for
% rest/active

%Create empty 3D matrix where rows = signal after window averaging, 
%columns = number of comparisons (each contact vs LFP), and 
% 3rd dimesion = number of epochs
all_epoch_meanPSD_r = zeros(NFFT/2+1,num_contact_pair,num_epoch); 
all_epoch_meanPSD_a = zeros(NFFT/2+1,num_contact_pair,num_epoch);

for i = 1:num_epoch
       
    for j = 1:num_contact_pair
        
        inp1_r = all_rest(epoch_starts_r(i):epoch_ends_r(i),j); %input containing data from 1 contact during 1 epoch
        inp1_a = all_active(epoch_starts_a(i):epoch_ends_a(i),j);
        
% [PSD,F] = pwelch(x,window,noverlap,nfft,Fs)               
            [PSD,F] = pwelch(inp1_r,WINDOW,NOVERLAP,NFFT,Fs);
            all_epoch_meanPSD_r(:,j,i) = PSD;
            
            [PSD,F] = pwelch(inp1_a,WINDOW,NOVERLAP,NFFT,Fs);
            all_epoch_meanPSD_a(:,j,i) = PSD;
    end
end

% Collapse 3D matrix into 2D matrix
%This will average all the epochs together. Resulting matrix will have 
%rows = PSD, columns = ecog contacts

grand_meanPSD_r = mean(all_epoch_meanPSD_r,3);
grand_meanPSD_a = mean(all_epoch_meanPSD_a,3);
   
for i = 1:num_contact_pair
  % normalize mean PSD to peak rest height  
    norm_meanPSD_r(:,i) = grand_meanPSD_r(:,i)/max(grand_meanPSD_r(:,i));
%     norm_meanPSD_a(:,i) = grand_meanPSD_a(:,i)/max(grand_meanPSD_a(:,i));
    norm_meanPSD_a(:,i) = grand_meanPSD_a(:,i)/max(grand_meanPSD_r(:,i));
    
    % calculate log base 10 of mean PSD
    log_meanPSD_r(:,i) = log10(grand_meanPSD_r(:,i));
    log_meanPSD_a(:,i) = log10(grand_meanPSD_a(:,i));

    %normalize log mean PSD to peak rest height
    norm_log_meanPSD_r(:,i) = log_meanPSD_r/max(log_meanPSD_r);
%     norm_log_meanPSD_a(:,i) = log_meanPSD_a/max(log_meanPSD_a);
    norm_log_meanPSD_a(:,i) = log_meanPSD_a/max(log_meanPSD_r);
   
    % calculate difference in log mean PSD between active and rest
%     active.contact_pair(i).difference =...
%         active.contact_pair(i).log_mean_PSD - rest.contact_pair(i).log_mean_PSD;
%     active.contact_pair(i).difference =...
%         active.contact_pair(i).norm_log_mean_PSD - rest.contact_pair(i).norm_log_mean_PSD;

    % calculate ratio of active to rest non-normalized mean PSD
%     active.contact_pair(i).ratio =...
%         100*(active.contact_pair(i).mean_PSD./rest.contact_pair(i).mean_PSD);
end

freq = F;

%% quantitative analysis of PSD
% now that all the data crunching is done and the results are stored in
% 'rest' and 'active' structure arrays, run sub-function quantPSD for
% quantitative analysis of PSD data

% rename mean PSD matrices to satisfy quantPSDmatrix variables
rest = grand_meanPSD_r';
active = grand_meanPSD_a';
normrest = norm_meanPSD_r';
normactive = norm_meanPSD_a';
logrest = log_meanPSD_r';
logactive = log_meanPSD_a';
normlogrest = norm_log_meanPSD_r';
normlogactive = norm_log_meanPSD_a';

[allfreq subfreq order] = quantPSDmatrix(rest,active,FREQ_QPSD,freq,filename);

% %% Create structure for exporting data
% ecogPSD.rest = rest; 
% ecogPSD.active = active; 
% ecogPSD.allfreq = allfreq;
% ecogPSD.subfreq = subfreq; 
% ecogPSD.order = order; 
% ecogPSD.freq = freq;
%% Save and write quantPSD matrix data as excel spreadsheet
[data4export] = exportdata(allfreq, subfreq, filename);

cd(pn); %exportdata changes the path in order to write the excel data. This brings back current path.
%% Save and write quantPSD data 
%Creates .mat output file of name (filename_ecogPSD.mat). transcoh.mat
%files can then be used to analyze coherence data across groups
if Dx==1
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\PD\'];
elseif Dx==2
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\DYS\'];
elseif Dx==3
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\ET\'];
elseif Dx==4
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\Epilepsy\'];
else
    outputdir = ['C:\Users\Starr\Documents\ECOG data\quantPSD data\Other\'];
end
outputname = [outputdir filename,'_ecogPSD.mat'];
disp(['Writing ecog PSD data to:  ' outputname]);
% save(outputname,'ecogPSD');
save(outputname,'order','rest', 'active', 'freq', 'logrest','logactive',...
    'normrest', 'normactive','normlogrest', 'normlogactive', 'allfreq', 'subfreq');
