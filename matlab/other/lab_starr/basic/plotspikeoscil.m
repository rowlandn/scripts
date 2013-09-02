function plotspikeoscil()

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes in names of excel file/sheet that contains the data 
% output from spike oscil analysis (RunSpikeOscil.m), then outputs
% a bar graph that displays the data.  It is an updated version of
% plotspikeoscil2.m.
% 
% CREATED BY: Sho Shimamoto 1/18/2007
% MODIFIED: SS 3/29/2007, 7/2/2007
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


% change directory to one containing SpikeOscil excel files
% cd('C:\Documents and Settings\Sho\My Documents\Lab Documents\Data\Spreadsheets\SpikeOscil');
% note: another option is to allow user to manually select directory by
% commenting out the directory change above and uncommenting the next line:
cd(uigetdir('','Selected folder that contains SpikeOscil sheets'));


% % Create filename options for the menu function
% file = {'SpikeOscilGPi.xls' 'SpikeOscilGPe.xls' 'SpikeOscilSTN.xls'...
%     'SpikeOscilTH.xls' 'SpikeOscilPPN.xls' 'SpikeOscilSNr.xls'};
% 
% % Create the str 'filename'
% file_k = menu('Select file', file);
% filename = file{file_k};
% 
% % Create sheetname options for the menu function
% if strcmp(filename,'SpikeOscilGPi.xls')
%     sheet = {'Dyst GPi spont' 'Dyst GPi active' 'Dyst GPi task'...
%         'PD GPi spont' 'PD GPi active'...
%         'Chorea GPi spont' 'Chorea GPi active'...
%         'Midbrain tremor GPi spont' 'Midbrain tremor GPi active'...
%         'nNHP GPi spont' 'nNHP GPi task' 'dysNHP GPi spont' 'dysNHP GPi task'};
% elseif strcmp(filename,'SpikeOscilGPe.xls')
%     sheet = {'Dyst GPe spont' 'Dyst GPe active'...
%         'PD GPe spont' 'Chorea GPe spont'...
%         'Midbrain tremor GPe spont' 'nNHP GPe spont' 'nNHP GPe task'};
% elseif strcmp(filename,'SpikeOscilSTN.xls')
%     sheet = {'Dyst STN spont' 'Dyst STN active'...
%         'PD STN spont' 'PD STN active' 'PD STN task'};
% elseif strcmp(filename, 'SpikeOscilTH.xls')
%     sheet = {'Midbrain tremor TH spont' 'Midbrain tremor TH active'};
% elseif strcmp(filename, 'SpikeOscilPPN.xls')
%     sheet = {'PPN all spont' 'PPN spont narrow AP' 'PPN spont wideAP' 'PPN active' 'PPN passive'};
% elseif strcmp(filename, 'SpikeOscilSNr.xls')
%     sheet = {'Dyst SNr spont'};
% end
FileList = dir('SpikeOscil*.xls');

if(isempty(FileList))
	str = pwd;
	error(['Found no SpikeOscil Excel files in current directory - ' str ]);
end

files = struct2cell(FileList);
files = files(1,:);

% Select file to be analyzed
file_k = menu('Select file', files);

%Create user-selected filename
filename = files{file_k};

% Get Excel file info
[typ, desc] = xlsfinfo(filename);

% Delete irrelevant sheets
expr = filename(11:12); % expr is nucleus identifier
start = regexp(desc, expr);
logic = cellfun(@isempty, start, 'UniformOutput', false);
logic = ~cell2mat(logic);
sheet = desc(logic);

% Create the str 'sheetname'
sheet_k = menu('Select sheet', sheet);
sheetname = sheet{sheet_k};

% Import oscillatory data and patient disorder type from selected
% 'filename' and 'sheetname'
[data, type] = xlsread(filename, sheetname);

numunits = size(data,1); % calculate number of units being analyzed
numcols = size(data,2); %calculate number of columns in data
numsigunits = sum(data(:,3)>0);  % calculate number of units with at least one significant oscillation
percentsig = numsigunits/numunits * 100; % calculate percent of units with significant oscillation

%% 1st plot (all data)
[allsigfreq1]=extract(data);
figure1 = figure;
% edges = [1e-3 1 2 3 5 10 15 20 30 60 100 200]; % define bin edges, do not include 0 as first edge
% bin_range = {'0-1','1-2','2-3','3-5','5-10','10-15','15-20','20-30','30-60','60-100','100-200 ',''};
edges = [1e-3 3 8 13 30 60 100 200]; % define bin edges, do not include 0 as first edge
len_edges = length(edges);
bin_range = {'0-3','3-8','8-13','13-30','30-60','60-100','100-200',''};
n = 100* histc(allsigfreq1, edges)/numunits;
if strncmp('PD',sheetname,2) || ...
        strncmp('Midbrain',sheetname,8)||...
        strncmp('PPN',sheetname,3)||...
        strncmp('nNHP', sheetname,4)||...
        strncmp('dysNHP', sheetname,6)
    axes('Parent', figure1,...
        'XTickLabel',bin_range,...
        'XTick',1:len_edges);
    box('on')
    hold('all')
    bar(n);
else
    subplot1 = subplot(2,1,1,'Parent',figure1,...
        'XTickLabel',bin_range,...
        'XTick',[1:len_edges]);
    box on;
    hold on;
    bar(n);
end
xlim([0 len_edges]);
ylim([0 20]);
xlabel('bins (Hz)'); ylabel ('Percent of units studied','fontsize',10);
title([sheetname,'  #units=', int2str(numunits),  ', %sig osc=', int2str(percentsig)]);

%% 2nd plot (subtype)
% confirm that sheetname is dyst spont sheet
if strncmp('Dyst', sheetname,4)

    % initialize identity matrices for each subtype
    idio_I = false(numunits, numcols);
    scnd_I = false(numunits, numcols);
    dytp_I = false(numunits, numcols);
    dytn_I = false(numunits, numcols);
    trdv_I = false(numunits, numcols);
    meig_I = false(numunits, numcols);
    unkn_I = false(numunits, numcols);
    % iterate through each unit in type and sort each subtype into its own identity
    for i = 1:numunits
        if strcmp('adult idio',type{i+1,11});
            idio_I(i,:)=true;

        elseif strcmp('secondary',type{i+1,11});
            scnd_I(i,:)=true;
        elseif strcmp('juv dyt 1+',type{i+1,11});
            dytp_I(i,:)=true;
        elseif strcmp('juv dyt 1-',type{i+1,11});
            dytn_I(i,:)=true;
        elseif strcmp('tardive',type{i+1,11});
            trdv_I(i,:)=true;
        elseif strcmp('meige', type{i+1,11});
            meig_I(i,:)=true;
        elseif strcmp('unknown',type{i+1,11});
            unkn_I(i,:)=true;
        elseif isempty(type{i+1,11});
            unkn_I(i,:)=true;
        end
    end
    % create data matrix for each subtype, then extract just the raw data,
    % eliminating empty entries
    data_idio = data.*idio_I;
    [idio] = extract(data_idio);
    data_scnd = data.*scnd_I;
    [scnd] = extract(data_scnd);
    data_dytp = data.*dytp_I;
    [dytp] = extract(data_dytp);
    data_dytn = data.*dytn_I;
    [dytn] = extract(data_dytn);
    data_trdv = data.*trdv_I;
    [trdv] = extract(data_trdv);
    data_meig = data.*meig_I;
    [meig] = extract(data_meig);
    data_unkn = data.*unkn_I;
    [unkn] = extract(data_unkn);
    
    % find max length
    c = max([length(idio) length(scnd) length(dytp) length(dytn) length(trdv) length(meig) length(unkn)]);

    % correct for length differnce so dimensions agree for horzcat function
    idio = vertcat(idio,NaN(c-length(idio),1));
    scnd = vertcat(scnd,NaN(c-length(scnd),1));
    dytp = vertcat(dytp,NaN(c-length(dytp),1));
    dytn = vertcat(dytn,NaN(c-length(dytn),1));
    trdv = vertcat(trdv,NaN(c-length(trdv),1));
    meig = vertcat(meig,NaN(c-length(meig),1));
    unkn = vertcat(unkn,NaN(c-length(unkn),1));

    % concatenate so each subtype is organized into columns
    allsigfreq2 = horzcat(idio,scnd,dytp,dytn,trdv,meig,unkn);

% perform the same subtype analysis for chorea pts
elseif strncmp('Chorea',sheetname,6)
    
    hd_I = false(numunits, numcols);
    juv_I = false(numunits, numcols);
    unkn_I = false(numunits, numcols);
    
    for i = 1:numunits
        if strcmp('HD',type{i+1,11});
            hd_I(i,:)=true;
        elseif strcmp('juv onset',type{i+1,11});
            juv_I(i,:)=true;
        elseif isempty(type{i+1,11});
            unkn_I(i,:)=true;
        end
    end
    
    data_hd = data.*hd_I;
    [hd] = extract(data_hd);
    data_juv = data.*juv_I;
    [juv] = extract(data_juv);
    data_unkn = data.*unkn_I;
    [unkn] = extract(data_unkn);

    c=max([length(hd) length(juv) length(unkn)]);
    
    % correct for length difference so dimesions agree for horzcat function
    hd = vertcat(hd,NaN(c-length(hd),1));
    juv = vertcat(juv,NaN(c-length(juv),1));
    unkn = vertcat(unkn,NaN(c-length(unkn),1));

    % concatenate
    allsigfreq2 = horzcat(hd,juv,unkn);
else
    return;
end

% plot histogram with bins

subplot2 = subplot(2,1,2,'Parent',figure1,...
    'XTickLabel',bin_range,...
    'XTick',[1:len_edges]);
box on;
hold all;
n = 100* histc(allsigfreq2, edges)/numunits;
bar1 = bar(n,'Parent',subplot2);

% set legend and bar graph color for each disorder type
j=1;
pointer=0;
colors = [0.8471 0.1608 0;...
    0.6784 0.9216 1;...
    1 0.6941 0.3922;...
    0.749 0.749 0;...
    0.03922 0.1412 0.4157;...
    0.9529 0.8706 0.7333;...
    1 0 1];

if strncmp('Dyst', sheetname,4)
    tmp = {'adult idio' 'secondary' 'juv dyt 1+' 'juv dyt 1-' 'tardive' 'meige' 'unknown'};
elseif  strncmp('Chorea',sheetname,6) 
    tmp = {'HD' 'juv onset' 'unknown'};
end

for i = 1:length(tmp)
    if allsigfreq2(1,i) > 0 %skip subtypes that have no sig oscil
        labels{j}= tmp{i};
        pointer(j)=i;
        h(j) = bar1(i);
        j= j+1;
    else
        continue;
    end
end

for j = 1:length(labels)
    set(h(j),'FaceColor',colors(pointer(j),:),'DisplayName',labels{j});
end

legend(h,labels,'Location','Best');
xlim([0 len_edges]);
ylim([0 20]);    %set range of histogram up to 12 bins
xlabel('bins (Hz)'); ylabel ('Percent of units studied','fontsize',10);
title([sheetname,'  #units=', int2str(numunits),  ', %sig osc=', int2str(percentsig)]);




%% subfunction extract
function [allsigfreq]=extract(data)
% This subfunction takes in data array and extracts positive non-zero
% numbers and stores them in a column vector called allsigfreq

columns = size(data,2); % calculate number of column in data

if columns == 4
    temp = data(:,3); % create a temp array that contains 1 row of significant oscillation
elseif columns == 6
    temp = [data(:,3) data(:,5)]; % create a temp array that combines 2 rows of significant oscillation
elseif columns == 8
    temp = [data(:,3) data(:,5) data(:,7)]; % create a temp array that combines 3 rows of significant oscillation
end

allsigfreq = temp(temp>0); % create a column vector that contains just the significant frequencies from all units


