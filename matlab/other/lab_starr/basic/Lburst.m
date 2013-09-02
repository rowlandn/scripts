% Lburst.m
% RORY M. 2003.05.16
%
% INPUT:  *.isi TEXT files with single-column of row-delimited ISI values
% OUTPUTS: L = burst statistic for # unique spike counts throughout all bins
%   isi_mean = mean of inter-spike-intervals in current file
%   range = 1x3 array of bin range statistics [min, max, max-min]
%   histo = histogram of values taken by the bins
%   bins = complete array of bins with spike counts used in analysis
%
% Contents of the resulting tab-delimited file should look something like the following  -
% 	Input_File     L_stat  MeanISI MinSpk/bin	MxSpk/bin Histogram
% 	beranleft01.isi	3	    15.2		0		2       10,3,1
% 	beranleft03.isi	5   	21.2		0		5       100,62,29,5,0,1

isi_files = dir('*.isi');
numfiles = length(isi_files);

fid = fopen([num2str(numfiles) 'set_burstanalysis.txt'],'w'); % OPEN OUTPUT TEXT FILE FOR WRITING
fprintf(fid,'Input_File\tL_Statistic\tMeanISI\tMinSpk/bin\tMaxSpk/bin\tHistogram\r\n'); % COLUMNS

for fileindex = 1:numfiles % ITERATE  THROUGH ALL ISI FILES IN DIRECTORY
    isi_data = load(isi_files(fileindex).name);
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
    % ############################################################# 
    % Now WRITE OUT the ISI stats to a file.
    fprintf(fid, '%s\t', isi_files(fileindex).name);
    fprintf(fid, '%3d\t', L);
    fprintf(fid, '%4.3f\t', isi_mean);
%    fprintf(fid, '%4.3f\t', isi_median); % just extra info
    fprintf(fid, '%3d\t%3d\t', range(1), range(2));
    fprintf(fid, '%s\r\n', fhisto); % ROB probably wants this commented out
end
% #############################################################
fclose(fid);
% #############################################################