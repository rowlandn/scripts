% Lburst2.m
% modification of Lburst.m by Sho Shimamoto 2006.12.19 to accept txt files output from
% nex2isi_multi
%
% INPUT:  *nex.txt files output from next2isi containing a column of ISIs
% OUTPUTS: 
%   TXT file containing 3 columns:
%   MeanRate
%   ISIMean
%   L-statistic

%
% Contents of the resulting tab-delimited file should look something like the following  -
% Input_File	MeanRate	ISIMean	L_statistic
% HustL05ISPFb_nex.txt	59.933	16.685	  4
% HustL06ISPCa_nex.txt	59.158	16.904	  6
% HustL07ISPFa_nex.txt	85.441	11.704	  5
% HustL08ISPFa_nex.txt	47.771	20.933	  6
% HustL11ISPFbkiel_nex.txt	62.241	16.067	  6
% HustR16ISPFc_nex.txt	59.009	16.947	  5
% HustR18ISPCa_nex.txt	65.691	15.223	  5

%dirname = uigetdir('C:\', 'Select directory that contains TXT files to be analyzed');
%cd(dirname);

txt_files = dir('*nex.txt'); % selects TXT files that were created by nex2isi with 'nex.txt' ending as dictated in the function write_isi_txt

numfiles = length(txt_files);

fid = fopen([num2str(numfiles) 'set_burstanalysis.txt'],'w'); % OPEN OUTPUT TEXT FILE FOR WRITING
fprintf(fid,'Input_File\tMeanRate\tISIMean\tL-statistic\r\n'); % COLUMNS

for fileindex = 1:numfiles % ITERATE  THROUGH ALL TXT FILES IN DIRECTORY
    temp = importdata(txt_files(fileindex).name); %import text file, store in temp array
    isi_data = temp.data(2:end,3); % extract just the column that contains ISIs, to match Rory's version, the first element on the ISI array will be purposely excluded
    isi_num = length(isi_data); % num spikes - 1
    isi_sum = sum(isi_data); % ELAPSED TIME OF RECORDING (from first spike)
    isi_mean = mean(isi_data);
    isi_median = median(isi_data); % unused statistic
    bin_n = isi_num; % # of bins needed (each is isi_mean size)
    bins = zeros(1,bin_n); % construct empty bins structure
    
    isi_index = 1; % start with first ISI
    bins(1) = 1; % TO CATCH THE FIRST SPIKE (at t=0 of first ISI) ???????
    isi_remainder = 0; % for storing remaining isi time BEYOND current bin
    for bin_index = 1:bin_n % MOVING THROUGH THE SPIKE COUNT STORAGE BINS
        bin_current = 1; % flag to check for moving to next bin
        bin_length = isi_mean; % SET LENGTH (ms) OF CURRENT BIN
        while bin_current == 1  % STAY INSIDE THIS BIN MOVING THROUGH ISI's
            if isi_remainder == 0
                if isi_data(isi_index) <= bin_length   % if isi "fits inside" bin
                    bins(bin_index) = bins(bin_index) + 1; % increment spk count
                    bin_length = bin_length - isi_data(isi_index); % remaining in bin
                    isi_index = isi_index + 1; % move to next isi
                else
                    isi_remainder = isi_data(isi_index) - bin_length; % isi BEYOND current bin
                    bin_length = 0; % should NOT matter, as we are changing bins!
                    bin_current = 0; % tell while() loop to change bins
                end
            else
                if isi_remainder <= bin_length % if CURRENT PARTIAL ISI "fits inside" bin
                    bins(bin_index) = bins(bin_index) + 1; % increment spk count
                    bin_length = bin_length - isi_remainder; % remaining in bin
                    isi_remainder = 0; % used all of current isi
                    isi_index = isi_index + 1; % move to next isi
                else
                    isi_remainder = isi_remainder - bin_length; %isi BEYOND current bin
                    bin_length = 0; % should NOT matter, as we are changing bins!
                    bin_current = 0; % tell while() loop to change bins
                end
            end
            if isi_index > isi_num % be certain ISI index does not exceed dimensions
                bin_current = 0;
            end
        end
    end
    if isi_remainder > 0 % Catch missed spikes at the end of some ISI files
        bins(bin_index) = bins(bin_index) + 1;
    end
    % ############################################################# 
    % CHECK RESULTS FOR ACCURACY AND SENSE (debug)
    spk_total = isi_num+1; % assume a spike at beginning and one at end
    binsum = sum(bins);
    if binsum ~= spk_total
        disp('TOTAL SPIKES in BINS NOT NUMISI + 1');
        binsum
        spk_total
    end
    binmean = sum(bins)/(bin_n+1) % increment numb of bins to account for extra spike in n-1 bins
    if binmean ~= 1
        disp('MEAN of the BINS IS NOT 1!!!');
        binmean
    end
    % #############################################################     
    % CALCULATE SOME STATISTICS:
    range = [min(bins), max(bins), max(bins)-min(bins)] % create range variable
    histo = 0; % Begin Histogram of values (# of spikes) in bins
    j = 1;
    while j <= bin_n
        if (bins(j)+1) > length(histo)
            histo(bins(j)+1) = 1; % initialize that spike counter with 1
        else
            histo(bins(j)+1) = histo(bins(j)+1) + 1; % increment spike counter
        end
        j = j+1;
    end
    histo % just to PRINT the histogram to console.
    fhisto = []; % CREATE COMMA-DELIMITED "formatted histograms" for output
    for i=1:length(histo)-1 
        fhisto = [fhisto num2str(histo(i)) ','];
    end
    fhisto = [fhisto num2str(histo(end))];
    L = length(find(histo)) % L Statistic == # of DISTINCT values in bins array
    mean_rate = 1000/isi_mean; % mean firing rate Hz
    % ############################################################# 
    % Now WRITE OUT the ISI stats to a file.
    fprintf(fid, '%s\t', txt_files(fileindex).name);
    fprintf(fid, '%4.3f\t', mean_rate);
    fprintf(fid, '%4.3f\t', isi_mean);
    fprintf(fid, '%3d\r\n', L);
end
% #############################################################
fclose(fid);
% #############################################################