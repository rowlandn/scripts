%% Notes
% GroupDataExport is a program to concatenate all the data from allfreq and
% subfreq matrices into a 3x31 array that will be saved as an excel
% spreadsheet. It does this for each limb file for a single patient with 
% either monopolar or bipolar recordings (not both). The excel spreadsheet 
% can then be copy/pasted into SPSS. This way there is limited time spent on 
% data entry.

%4/3/2009 ALC

%% Call files
path = uigetdir('', 'Select directory that contains _ecogPSD files to be analyzed');
path = [path '\'];
cd(path);

Ptdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that have allfreq and subfreq variables
numlimbs = length(Ptdir);

for i = 1:numlimbs
    LimbNames(i) = importdata(Ptdir(i).name);
end

for i = 1: numlimbs
    fn = strrep(Ptdir(i).name,'_ecogPSD.mat','');
    LimbNames(i).name = fn;
end
%% linear indexing and transposition of allfreq in preparation for concatenation 

for i = 1:numlimbs
    LimbNames(i).lidxallfreq = LimbNames(i).allfreq(:); % make 2x3x3 table into 18x1 column vector
    LimbNames(i).lidxallfreq = LimbNames(i).lidxallfreq'; %make column vector into row vector
end
%% Create matrix with concatenated data from allreq and subfreq
% concatenate data into variable "data4export." This is in a format that 
% can be copy/pasted directly into the SPSS spreadsheet

for i = 1:numlimbs
    LimbNames(i).data4export = [LimbNames(i).lidxallfreq(1:6)...
        LimbNames(i).subfreq(1,:,1) LimbNames(i).subfreq(2,:,1)...
        LimbNames(i).subfreq(3,:,1) LimbNames(i).subfreq(4,:,1)...
        LimbNames(i).subfreq(5,:,1); LimbNames(i).lidxallfreq(7:12)...
        LimbNames(i).subfreq(1,:,2) LimbNames(i).subfreq(2,:,2)...
        LimbNames(i).subfreq(3,:,2) LimbNames(i).subfreq(4,:,2)...
        LimbNames(i).subfreq(5,:,2); LimbNames(i).lidxallfreq(13:18)...
        LimbNames(i).subfreq(1,:,3) LimbNames(i).subfreq(2,:,3)...
        LimbNames(i).subfreq(3,:,3) LimbNames(i).subfreq(4,:,3) LimbNames(i).subfreq(5,:,3)];
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

for i = 1:numlimbs
    [status, message] = xlswrite([LimbNames(i).name], LimbNames(i).data4export);
end 