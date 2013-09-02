function plotrateburst()

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function allows the user to select from a menu box the name 
% of RateBurst XLS file that requires statistical
% analysis.  It also allows user to decide if subtype analysis is needed.
% The output is a TXT file with the name 'XLSfilenamedata.txt'
% that contains the results of the analysis.
%
%
% CREATED BY: Sho Shimamoto 1/18/2007
% Modified: SS 1/24/2007 , SS 1/31/2007, SS 2/6/2007, SS 3/29/2007
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


% Change current directory to one containing SpikeOscil excel files. If the location of
% the excel files changes,  the directory must be edited accordingly.
cd('C:\Documents and Settings\HP_Administrator\My Documents\Lab Documents\data\Spreadsheets\RateBurst');
% note: another option is to allow user to manually select directory by
% commenting out the directory change above and uncommenting the next line:
%cd(uigetdir);

% Create filename options for the menu function.
% The filename assignments are written under the assumption that the XLS filenames have not been changed.  
% In the case that the filenames have been changed, the next line must be edited accordingly.
file = {'RateBurstGPi.xls' 'RateBurstGPe.xls' 'RateBurstSTN.xls' 'RateBurstTH.xls'};

% Allow user to select file to be analyzed
file_k = menu('Select file', file);
% Allow user to decide if a subtype analysis is needed.  This analysis will be performed in "subtype analysis" cell.
subtype_k = menu('Subtype analysis', 'Yes', 'No');

%Create the str 'filename' for compatibility using strcmp
filename = file{file_k};

% Create sheetname options for the menu function
% The sheetname assignments are created under the assumption that the
% sheetnames have not been changed.  In the case that the sheetnames have been
% changed, the name assingments must be edited accordingly.
if strcmp(filename,'RateBurstGPi.xls')
    sheet = {'Dyst GPi spont' 'Dyst GPi active' 'Dyst GPi task'...
        'PD GPi spont' 'PD GPi active'...
        'Chorea GPi spont' 'Chorea GPi active'...
        'Midbrain tremor GPi spont' 'Midbrain tremor GPi active'...
        'nNHP GPi spont' 'nNHP GPi task' 'dysNHP GPi spont' 'dysNHP GPi task'};
elseif strcmp(filename,'RateBurstGPe.xls')
    sheet = {'Dyst GPe spont' 'Dyst GPe active'...
        'PD GPe spont' 'Chorea GPe spont'...
        'Midbrain tremor GPe spont' 'nNHP GPe spont' 'nNHP GPe task'};
elseif strcmp(filename,'RateBurstSTN.xls')
    sheet = {'Dyst STN spont' 'Dyst STN active'...
        'PD STN spont' 'PD STN active' 'PD STN task'};
elseif strcmp(filename, 'RateBurstTH.xls')
    sheet = {'Midbrain tremor TH spont' 'Midbrain tremor TH active'};
end

fid = fopen([filename(1:length(filename)-4) 'data.txt'],'w'); % open text file for writing
fprintf(fid,'type\tcount\tstatistic\tMeanRate\tISIMean\tISIStDev\tBurstIndex\tMeanRate\tMean freq in burst\tProp spikes in burst\tMeanRate\tISIMean\tL-statistic\r\n'); % columns
fprintf(fid,'all type analysis\r\n');
%% all type analysis
% Perform analysis on all types on each sheet

% initialize
d(1,1).temp=0;  
d(1,1).extract=0;
alldata(1,1).mean(1,1) = 0;
alldata(1,1).std(1,1) = 0;
alldata(1,1).numunits = 0;

for i=1:length(sheet)

    sheetname = sheet{i};

    % Import rateburst data from selected
    % 'filename' and 'sheetname' and store 
    [data(i).store, type(i).store] = xlsread(filename, sheetname);
    
    if isempty(data(i).store) %check for sheet with no data and skip
        continue
    end
    % calculate number of units
    alldata(i).numunits = size(data(i).store,1);
    % delete 'NaN' entries so that raw data can be extracted and stored in
    % structure "d" for 10 columns of imported data
    for j = 1:10
        d(j).temp =data(i).store(:,j);
        d(j).extract =d(j).temp(d(j).temp>0);
        alldata(i).mean(j) = mean(d(j).extract);
        alldata(i).std(j) = std(d(j).extract);
    end
end

% write to TXT output
for i=1:length(sheet) 
    fprintf(fid, '%s\t', sheet{i});
    if isempty(alldata(i).numunits)
        fprintf(fid, '0\r\n');
        continue;
    end
    fprintf(fid, '%3d\t', alldata(i).numunits);
    fprintf(fid, '%s\t', 'mean');
    for j = 1:9
        fprintf(fid, '%4.4f\t', alldata(i).mean(j));
    end
    fprintf(fid, '%4.4f\r\n', alldata(i).mean(10));
    fprintf(fid, '%s\t','');
    fprintf(fid, '%s\t','');
    fprintf(fid, '%s\t', 'stand dev');
    for j = 1:9
        fprintf(fid, '%4.4f\t', alldata(i).std(j));
    end
    fprintf(fid, '%4.4f\r\n', alldata(i).std(10));
end
%% Subtype analysis

% initialize
D(1,1).temp=0;  
D(1,1).extract=0;
subdata(1,1).mean(1,1)=0;
subdata(1,1).std(1,1)=0;
subdata(1,1).numunits(1,1)=0;


if subtype_k == 1 % perform subtype analysis
    % for each filename selction, store spont and active sheet indeces into
    % index sheet_i
    if strcmp(filename,'RateBurstGPi.xls')
        sheet_i = [1 2 6 7];
    elseif strcmp(filename,'RateBurstGPe.xls')
        sheet_i = [1 2 4];
    elseif strcmp(filename,'RateBurstSTN.xls')
        sheet_i = [1 2];
    end
    %---------- main loop ----------%
    %iterates through all sheets that contain subtypes
    for j = 1:length(sheet_i)

        [numunits numcols] = size(data(sheet_i(j)).store);
        
        if j <= 2

            % initialize identity matrix for each dystonia subtype
            idio_I = false(numunits, numcols);
            unkn_I = false(numunits, numcols);
            scnd_I = false(numunits, numcols);
            juvp_I = false(numunits, numcols);
            juvn_I = false(numunits, numcols);
            trdv_I = false(numunits, numcols);
            meig_I = false(numunits, numcols);

            % --------- identity loop ---------- &
            % iterate through each unit then sort each subtype into its own identity matrix
            for m = 1:numunits
                if isempty(type(sheet_i(j)).store{m+2,13});
                    unkn_I(m,:)=true;
                elseif strcmp('adult idio',type(sheet_i(j)).store{m+2,13});
                    idio_I(m,:)=true;
                elseif strcmp('unknown',type(sheet_i(j)).store{m+2,13});
                    unkn_I(m,:)=true;
                elseif strcmp('secondary',type(sheet_i(j)).store{m+2,13});
                    scnd_I(m,:)=true;
                elseif strcmp('juv dyt 1+',type(sheet_i(j)).store{m+2,13});
                    juvp_I(m,:)=true;
                elseif strcmp('juv dyt 1-',type(sheet_i(j)).store{m+2,13});
                    juvn_I(m,:)=true;
                elseif strcmp('tardive',type(sheet_i(j)).store{m+2,13});
                    trdv_I(m,:)=true;
                elseif strcmp('meige', type(sheet_i(j)).store{m+2,13});
                    meig_I(m,:)=true;
                end
            end

            % store each dystonia subtype data into the structure 'subtype' field 'data'
            subtype(j,1).data = data(sheet_i(j)).store.*idio_I;
            subtype(j,2).data = data(sheet_i(j)).store.*scnd_I;
            subtype(j,3).data = data(sheet_i(j)).store.*juvp_I;
            subtype(j,4).data = data(sheet_i(j)).store.*juvn_I;
            subtype(j,5).data = data(sheet_i(j)).store.*trdv_I;
            subtype(j,6).data = data(sheet_i(j)).store.*meig_I; 
            subtype(j,7).data = data(sheet_i(j)).store.*unkn_I;

            % define subtype name for dystonia pts
            subtype(j).name = {'adult idio' 'secondary' 'juv dyt 1+' 'juv dyt 1-' 'tardive' 'meige' 'unknown'};

            % ---------- analysis loop --------- %
            % iterate through each subtype and perform statistical analysis within loop
            for n = 1:length(subtype(j).name)
                % delete 'NaN' entries so that 10 columns of imported data stored in subtype structure can be extracted and stored in
                % structure "D," then store mean and std of each column into
                % output structure "subdata"
                % --------- column loop -------- %
                for k = 1:10
                    D(k).temp = subtype(j,n).data(:,k);
                    D(k).extract = D(k).temp(D(k).temp>0);
                    subdata(j).mean(n,k) = mean(D(k).extract);
                    subdata(j).std(n,k) = std(D(k).extract);
                end
                subdata(j).numunits(n) = length(D(1).extract);
            end
           

        elseif j > 2

            % initialize identity matrix for each chorea subtype
            hd_I = false(numunits, numcols);
            juv_I = false(numunits, numcols);
            unk_I = false(numunits, numcols);

            % ---------- identity loop --------- %
            % iterate through each unit in TYPE2 then sort each chorea subtype into its own
            % identity matrix
            for m = 1:numunits
                if strcmp('HD', type(sheet_i(j)).store{m+2,13});
                    hd_I(m,:)=true;
                elseif strcmp('juv onset', type(sheet_i(j)).store{m+2,13});
                    juv_I(m,:)=true;
                elseif isempty(type(sheet_i(j)).store{m+2,13});
                    unk_I(m,:)=true;
                elseif strcmp('unknown',type(sheet_i(j)).store{m+2,13});
                    unk_I(m,:)=true;
                end
            end
            % store each chorea subtype data into the structure 'subtype' field 'data'
            subtype(j,1).data = data(sheet_i(j)).store.*hd_I;
            subtype(j,2).data = data(sheet_i(j)).store.*juv_I;
            subtype(j,3).data = data(sheet_i(j)).store.*unk_I;

            % define subtype name for chorea pts
            subtype(j).name = {'HD' 'juv onset' 'unknown'};
            
            % ---------- analysis loop ----------- %
            % iterate through each subtype and perform statistical analysis within loop
            for n = 1:length(subtype(j).name)
                % delete 'NaN' entries so that 10 columns of imported data stored in subtype structure can be extracted and stored in
                % structure D
                % ----------- column loop ---------%
                for k = 1:10
                    D(k).temp = subtype(j,n).data(:,k);
                    D(k).extract = D(k).temp(D(k).temp>0);
                    subdata(j).mean(n,k) = mean(D(k).extract);
                    subdata(j).std(n,k) = std(D(k).extract);
                end
                subdata(j).numunits(n) = length(D(1).extract);
            end
            
        end
    end
    % write to txt output
    fprintf(fid,'subtype analysis\r\n');
    for j = 1:length(subdata)
        fprintf(fid, '%s\r\n', sheet{sheet_i(j)});
        for n = 1:length(subtype(j).name)
            fprintf(fid, '%s\t', subtype(j).name{n});
            if isempty(subdata(j).numunits(n))
                fprintf(fid,'0\r\n');
                continue;
            end
            fprintf(fid, '%3d\t', subdata(j).numunits(n));
            fprintf(fid, '%s\t', 'mean');
            for k = 1:9
                fprintf(fid, '%4.4f\t', subdata(j).mean(n,k));
            end
            fprintf(fid, '%4.4f\r\n', subdata(j).mean(n,10));
            fprintf(fid, '%s\t','');
            fprintf(fid, '%s\t','');
            fprintf(fid, '%s\t', 'stand dev');
            for k = 1:9
                fprintf(fid, '%4.4f\t', subdata(j).std(n,k));
            end
            fprintf(fid, '%4.4f\r\n', subdata(j).std(n,10));

        end
    end

end



fclose(fid);


%% plot all type error bars
X(1)=0;
Y(1)=0;
E(1)=0;    
if strcmp(filename,'RateBurstGPi.xls')
    spontsheet_i = [1 4 6 8 10 12];
elseif strcmp(filename,'RateBurstGPe.xls')
    spontsheet_i = [1 3 4 5 6];
elseif strcmp(filename,'RateBurstSTN.xls')
    spontsheet_i = [1 3];
elseif strcmp(filename,'RateBurstTH.xls')
    spontsheet_i = 1;
end

% initialize
disordernames = cell(1,length(spontsheet_i));
for j = 1:length(spontsheet_i)
    name = sheet{spontsheet_i(j)};
    disordernames{j} = name(1:length(name)-9);
end

for j = 1:length(spontsheet_i)
    X(j) = j;
    Y(j) = alldata(spontsheet_i(j)).mean(1);
    E(j) = alldata(spontsheet_i(j)).std(1)/sqrt(alldata(spontsheet_i(j)).numunits);
end

% Create figure
figure1 = figure;
axes('Parent',figure1,'XTickLabel',disordernames,...
    'XTick', X);
box('on');
hold('on');
% Create errorbar
errorbar(X,Y,E,'MarkerSize',15,'Marker','.','LineStyle','none','Color',[0 0 1]);
title([filename(10:length(filename)-4) ' spontaneous discharge rate']);
ylabel('dicharge rate (Hz)');

%% create all type boxplot
for j = 1:length(spontsheet_i)
    x(j).data = data(spontsheet_i(j)).store(:,1);
    num(j) = length(x(j).data);
end

c=max(num);

%initialize X
X = zeros(c,length(spontsheet_i));

% correct for length difference so dimesions agree for horzcat function
for j = 1:length(spontsheet_i)
    x(j).data = vertcat(x(j).data,NaN(c-length(x(j).data),1));
    X(:,j) = x(1,j).data;
end

figure2 = figure;

box('on');
boxplot(X,'symbol','.','labels',disordernames);
title([filename(10:length(filename)-4) ' spontaneous discharge rate']);
ylabel('Discharge rate (Hz)');
%% create subtype error bar
% X(1)=0;
% Y(1)=0;
% E(1)=0;    
% if strcmp(filename,'RateBurstGPi.xls')
%     spontsheet_i = [1 2 6 7];
% elseif strcmp(filename,'RateBurstGPe.xls')
%     spontsheet_i = [1 2 4];
% elseif strcmp(filename,'RateBurstSTN.xls')
%     spontsheet_i = [1 2];
% end
% 
% A(1,1).store=0;
% for j = 1:length(spontsheet_i)
%     for n = 1:length(subtype(j).name)
%         temp = subtype(j,n).data(:,1);
%         A(j,n).store = temp(temp>0);
%     end
% end
% i = 1;
% for j = 1:length(spontsheet_i)
%     for n = 1:length(subtype(j).name)
%         num(i) = length(A(j,n).store);
%         i = i+1;
%     end
% end
% c = max(num);
% X = NaN(c,length(num));
% i = 1;
% for j = 1:row
%     for n = 1:col
%         if isempty(A(j,n).store)
%             continue
%         end
%         X(:,i) = vertcat(A(j,n).store,NaN(c-length(A(j,n).store),1));
%         i = i+1;
%     end
% end
% 
% for 
% 
% 
% % Create figure
% figure3 = figure;
% axes('Parent',figure1,'XTickLabel',disordernames,...
%     'XTick', X);
% box('on');
% hold('on');
% % Create errorbar
% errorbar(X,Y,E,'MarkerSize',15,'Marker','.','LineStyle','none','Color',[0 0 1]);
% title([filename(10:length(filename)-4) ' spontaneous discharge rate for subtypes']);
% ylabel('dicharge rate (Hz)');
% end


return;






