function [L, isi_mean, range, bins, histo] = roryburst(ts)
% roryburst.m
% INCLUDES WRITING A TEXT FILE
% OUTPUTS: L = statistic for # unique values throughout the bins
%   isi_mean = mean of inter-spike-intervals
%   range = 1x2 array range (min,max) of values found in bins
%   bins = complete array of spikes/bin
%   histo = histogram of values taken by the bins
%
% RORY M. 2003.04.10

isi_data = diff(ts);
isi_mean = mean(isi_data);
isi_n = length(isi_data);

ts_n = length(ts);
bin_n = ceil(ts(ts_n)/isi_mean);
bin_low = 0;
bin_high = bin_low+isi_mean;
bins = zeros(1,bin_n); % Stochastic process "Pt"

i=1; j=1;
spk_count = 0;
%for i = 1:ts_n
while i <= ts_n
    if ts(i) < bin_high
        spk_count = spk_count+1;
        i=i+1;
    else
        bins(j) = spk_count;
        bin_high = bin_high + isi_mean;
        spk_count = 0;
        j=j+1;
    end
end
bins(j) = spk_count % do the FINAL bin
if mean(bins) ~= 1
    error('Incorrect binning!  Should be mean == 1 spk/bin');
end

range = [min(bins), max(bins)]

histo = 0; % Begin Histogram of values in bins (Pt)
j = 1;
while j <= bin_n
    if (bins(j)+1) > length(histo)
        histo(bins(j)+1) = 1; % index counter = 1
    else
        histo(bins(j)+1) = histo(bins(j)+1) + 1; % increment counter
    end
    j = j+1;
end
histo
L = length(find(histo)) % Number of DISTINCT values in bins (Pt)

% #############################################################
% Now WRITE OUT the ISI data to a file.

fid = fopen( [fname '_burststats.txt'],'w');
fprintf(fid,'L = %3d\r\n',L);
fprintf(fid,'Mean ISI = %4.3f\r\n',isi_mean);
fprintf(fid,'Range of Spikes/Bin: Min = %3d\t Max = %3d\r\n', range);
fprintf(fid,'------------------------\r\n');
for i=1:bin_n
     fprintf(fid,'%3d\r\n',bins(i));
end

fclose(fid);
