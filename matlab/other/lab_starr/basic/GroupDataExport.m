%Group Data Export (not 3x31 version)
%% Notes
% GroupDataExport is a program to concatenate all the data from allfreq and
% subfreq matrices into a 3x31 array that will be saved as an excel
% spreadsheet. It does this for each limb file for a single patient with 
% either monopolar or bipolar recordings (not both). The excel spreadsheet 
% can then be copy/pasted into SPSS. This way there is limited time spent on 
% data entry.

%4/3/2009 ALC
%Updated 3/14/10 ALC to accomodate rest data. This code cannot yet
%accomodate KJM epilepsy data
%% Input task condition vs rest only condition
MVMTCOND = {'Movement Task' 'Rest Only'};
k = menu('Select either Movement Task or Rest Condition',MVMTCOND);
%% Call files
path = uigetdir('', 'Select directory that contains _ecogPSD files to be analyzed');
path = [path '\'];
cd(path);

if k==2;
    Ptdir = dir('*_ecogPSDrest.mat');
    numPt = length(Ptdir);
    for i = 1:numPt
        PtNames(i) = importdata(Ptdir(i).name);
    end
    for i = 1: numPt
        fn = strrep(Ptdir(i).name,'_ecogPSDrest.mat','');
        PtNames(i).name = fn;
    end

else
    Ptdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that have allfreq and subfreq variables
    numPt = length(Ptdir);
    for i = 1:numPt
        PtNames(i) = importdata(Ptdir(i).name);
    end
    for i = 1: numPt
        fn = strrep(Ptdir(i).name,'_ecogPSD.mat','');
        PtNames(i).name = fn;
    end
end
%% linear indexing and transposition of allfreq in preparation for concatenation 
   
for i = 1:numPt
    % allfreq may be 2x3x4 rather than 2x3x3, but now only need 2x3x3, so
    % will limit allfreq to those dimensions ALC 3/14/10
    PtNames(i).allfreq = PtNames(i).allfreq(:,:,1:3);
    PtNames(i).vertallfreq = vertcat(PtNames(i).allfreq(:,:,1), ...
        PtNames(i).allfreq(:,:,2), PtNames(i).allfreq(:,:,3)); % make 2x3x3 table into 6x3 matrix
end
%% Create matrix with concatenated data from allreq and subfreq
% concatenate data into variable "data4export." This is in a format that 
% can be copy/pasted directly into the SPSS spreadsheet
if k==2
    for i = 1:numPt
        PtNames(i).subfreqRM = vertcat(PtNames(i).subfreq(:,1:2,1), ...
            PtNames(i).subfreq(:,1:2,2), PtNames(i).subfreq(:,1:2,3));
    end
else
    for i = 1:numPt
        PtNames(i).subfreqRM = vertcat(PtNames(i).subfreq(:,1:2,1), ...
            PtNames(i).subfreq(:,1:2,2), PtNames(i).subfreq(:,1:2,3),...
            PtNames(i).subfreq(:,3:4,1), PtNames(i).subfreq(:,3:4,2), ...
            PtNames(i).subfreq(:,3:4,3));
        
        PtNames(i).subfreqAVR = vertcat(PtNames(i).subfreq(:,5,1), ...
            PtNames(i).subfreq(:,5,2), PtNames(i).subfreq(:,5,3));
    end
end
%% Save data to Excel
% variable "status" tells you if it saved properly ("true")
% variable "message" will give you an error message if status=false
% save the file name in the format 'PtName_Limb_Pol' where PtName = first 4
% letters of pt's last name; Limb = hand/elbow/shoulder/jaw/foot; and Pol=
% bipolar/monopolar. If the data was not collected as Recog-Rlfp-Llimb or
% Lecog-Llfp-Rlimb, make appropriate notation in filename.

% Directory for saving the excel data 
outputpath = ['C:\Users\Starr\Documents\ECOG data analysis\Excel data\'];
cd(outputpath);

if k==2
    for i = 1:numPt
        [status, message] = xlswrite([PtNames(i).name], PtNames(i).vertallfreq, 'Sheet1', 'A1');
        [status, message] = xlswrite([PtNames(i).name], PtNames(i).subfreqRM, 'Sheet1', 'D7');
    end
else
    for i = 1:numPt
    [status, message] = xlswrite([PtNames(i).name], PtNames(i).vertallfreq, 'Sheet1', 'A1');
    [status, message] = xlswrite([PtNames(i).name], PtNames(i).subfreqRM, 'Sheet1', 'D7');
    [status, message] = xlswrite([PtNames(i).name], PtNames(i).subfreqAVR, 'Sheet1', 'F37');
    end 
end