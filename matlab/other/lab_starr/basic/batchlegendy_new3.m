% batchlegendy_new3
% Created by Sho to perform a batch job of legendy_new3 on current directory w/ surprise cutoff
% at 5.  Outputs txt file containing data relevant to our study.


% 2 lines of code commented out by Sho (8/1/2007)
% directory will be selected in Isistats line 87
% dirname = uigetdir('', 'Select directory that contains TXT files to be analyzed');
% cd(dirname);

txt_files = dir('*nex.txt'); % selects TXT files that were created by nex2isi with 'nex.txt' ending as dictated in the function write_isi_txt

numfiles = length(txt_files);

fid = fopen([num2str(numfiles) 'set_legendyanalysis.txt'],'w'); % OPEN OUTPUT TEXT FILE FOR WRITING
fprintf(fid,'Input_File\tMean Rate\tMean freq in burst\tProp spikes in bursts\r\n'); % COLUMNS

for i = 1:numfiles
    temp = importdata(txt_files(i).name); %import text file, store in temp array
    ISI = temp.data(:,3);
    Burst = legendy_new3(ISI,[],[],[],[],5); %perform analysis with surprise_cutoff = 5
    fprintf(fid, '%s\t', txt_files(i).name);
    if isempty(Burst(1).baseline_rate)
        fprintf(fid, '%4.3f\t', NaN);
        fprintf(fid, '%4.3f\t', NaN);
        fprintf(fid, '%4.3f\r\n', NaN);
    else
    fprintf(fid, '%4.3f\t', Burst(1).baseline_rate);
    fprintf(fid, '%4.3f\t', Burst(1).mean_intra_burst_frequency);
    fprintf(fid, '%4.3f\r\n', Burst(1).proportion_spikes_in_bursts);
    end
end

fclose(fid);