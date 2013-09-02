function batch_relative_pausefind2 %(minlength)
% batch_pausefind2 outputs txt file containing data relevant to our study.

%% define function
%this function analyzes pausers, asking for its input, the output of nex2isi
%example: [IPI,pausefreq,rectime]=pausefind2(data,100)

%if nargin <1
%    minlength=500;
%end

%%
dirname = uigetdir('', 'Select directory that contains TXT files to be analyzed');
cd(dirname);

txt_files = dir('*nex.txt'); % selects TXT files that were created by nex2isi with 'nex.txt' ending as dictated in the function write_isi_txt

numfiles = length(txt_files);

%%
fid = fopen([num2str(numfiles) 'roughdata_pauselength.txt'],'w'); % OPEN OUTPUT TEXT FILE FOR WRITING
fprintf(fid,'Input_File\tpauselengths\r\n'); % COLUMNS

for i = 1:numfiles;
    temp = importdata(txt_files(i).name); %import text file, store in temp array
    isi = temp.data(:,3);
    fprintf(fid, '%s\t', txt_files(i).name);
    minlength = (mean(isi)*20);

    %% calculates pause parameters
    numspikes=length(isi)+1; %calculates total number of spikes
    rectime=sum(isi)/1000; % calculates total recording time in seconds
    timestamp=cumsum(isi); %calculates individual ISI timestamps
    Idx = find(isi>minlength);
    pauselength = (isi(Idx)); % length of each pause

    if isempty(pauselength)
        fprintf(fid, '%4.3f\r\n', NaN);
        continue;
    end

    fprintf(fid, '%4.2f\t', pauselength);
    fprintf(fid, '    \r\n');


end

fclose(fid);


%%
fid = fopen([num2str(numfiles) 'roughdata_IPI.txt'],'w'); % OPEN OUTPUT TEXT FILE FOR WRITING
fprintf(fid,'Input_File\tIPIs\r\n'); % COLUMNS

for i = 1:numfiles;
    temp = importdata(txt_files(i).name); %import text file, store in temp array
    isi = temp.data(:,3);
    fprintf(fid, '%s\t', txt_files(i).name);
    minlength = (mean(isi)*20);

    %% calculates pause parameters
    numspikes=length(isi)+1; %calculates total number of spikes
    rectime=sum(isi)/1000; % calculates total recording time in seconds
    timestamp=cumsum(isi); %calculates individual ISI timestamps
    Idx = find(isi>minlength);
    pauselength = (isi(Idx)); % length of each pause
    numpause=length(pauselength);
    pausefreq = numpause./rectime;

    if isempty(pauselength)
        fprintf(fid, '%4.3f\r\n', NaN);
        continue;
    end

    if numpause == 1;
        fprintf(fid, '%4.3f\r\n', NaN);
        continue;
    end

    eptime = timestamp(Idx); % time of the end of each pause
    IPIunc = diff(eptime);
    IPIcf = pauselength(2:(end));
    IPI = IPIunc - IPIcf;
    IPI(abs(IPI)<0.001)=0; % takes all values if IPI very close to zero and rounds to zero

    fprintf(fid, '%4.2f\t', IPI);
    fprintf(fid, '   \r\n');

end

fclose(fid);


%%
fid = fopen([num2str(numfiles) 'pausefrequency.txt'],'w'); % OPEN OUTPUT TEXT FILE FOR WRITING
fprintf(fid,'Input_File\tpausefrequency\r\n'); % COLUMNS

for i = 1:numfiles;
    temp = importdata(txt_files(i).name); %import text file, store in temp array
    isi = temp.data(:,3);
    fprintf(fid, '%s\t', txt_files(i).name);
    minlength = (mean(isi)*20);

    %% calculates pause parameters
    numspikes=length(isi)+1; %calculates total number of spikes
    rectime=sum(isi)/1000; % calculates total recording time in seconds
    timestamp=cumsum(isi); %calculates individual ISI timestamps
    Idx = find(isi>minlength);
    pauselength = (isi(Idx)); % length of each pause
    numpause=length(pauselength);
    pausefreq = numpause./rectime;

    fprintf(fid, '%4.3f\t', pausefreq);
    fprintf(fid, '  \r\n');

end

fclose(fid);

%%
fid = fopen([num2str(numfiles) 'mean_pausevalues.txt'],'w'); % OPEN OUTPUT TEXT FILE FOR WRITING
fprintf(fid,'Input_File\tmeanpauselength\tmeanIPI\tpausefrequency\r\n'); % COLUMNS

for i = 1:numfiles;
    temp = importdata(txt_files(i).name); %import text file, store in temp array
    isi = temp.data(:,3);
    fprintf(fid, '%s\t', txt_files(i).name);
    minlength = (mean(isi)*20);

 %% calculates pause parameters
    numspikes=length(isi)+1; %calculates total number of spikes
    rectime=sum(isi)/1000; % calculates total recording time in seconds
    timestamp=cumsum(isi); %calculates individual ISI timestamps
    Idx = find(isi>minlength);
    pauselength = (isi(Idx)); % length of each pause
    numpause=length(pauselength);
    pausefreq = numpause./rectime;

    if isempty(pauselength)
        fprintf(fid, '%4.3f\r\n', NaN);
        continue;
    end

    if numpause == 1;
        fprintf(fid, '%4.3f\r\n', NaN);
        continue;
    end

    eptime = timestamp(Idx); % time of the end of each pause
    IPIunc = diff(eptime);
    IPIcf = pauselength(2:(end));
    IPI = IPIunc - IPIcf;
    IPI(abs(IPI)<0.001)=0; % takes all values if IPI very close to zero and rounds to zero

    IPI3 = Idx(IPI>0);
    meanpauselength = mean(pauselength);
    meanIPI = mean(IPI3);

    fprintf(fid, '%4.2f\t', meanpauselength);
    fprintf(fid, '%4.2f\t', meanIPI);
    fprintf(fid, '%4.2f\t', pausefreq);
    fprintf(fid, '  \r\n');
end
fclose(fid);