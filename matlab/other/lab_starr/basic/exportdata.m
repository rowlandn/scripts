function [status, message] = exportdata(allfreq, subfreq, filename);
% Similar to the program "DataExport.m", the function "exportdata"
% will concatenate all the data from allfreq and
% subfreq matrices into the format required by SPSS. It will save this
% formatted data into an excel spreadsheet.The excel spreadsheet can 
% then be copy/pasted into SPSS. This way there is limited time spent 
% on data entry.

% The format is 51 rows x 6 columns, with a lot of empty spaces. 
% Column headings are:
% {'MaxFQ','MaxPR','TotPR','POWER','PERCT','AcVsR'} - (SPSS variable names)

%4/5/2009 ALC
%% Make 3D allfreq matrix into 2D matrix with vertcat
%Creates 6x3 array for MaxFQ, MAXPR, TotPR, such that:
%row 1 = M1, rest
%row 2 = M1, active
%row 3 = S1, rest
%row 4 = S1, active
%row 5 = STN, rest
%row 6 = STN, active

vertallfreq = vertcat(allfreq(:,:,1), allfreq(:,:,2), allfreq(:,:,3));

%% Make subfreq matrix suitable for export
%Creates 30x2 array for POWER and PERCT variables for each brain 
%region (M1, S1, or STN) and each movement condition (rest, active) and
%each frequency band of interest (delta, lo beta, hi beta, lo gamma, hi
%gamma)
%rows 1-5 = M1, rest, each freq band
%rows 6-10 = S1, rest, each freq band
%rows 11-15 = STN, rest, each freq band
%rows 16-20 = M1, active, each freq band
%rows 21-25 = S1, active, each freq band
%rows 26-30 = STN, active, each freq band

subfreqRM = vertcat(subfreq(:,1:2,1), subfreq(:,1:2,2), subfreq(:,1:2,3),...
                    subfreq(:,3:4,1), subfreq(:,3:4,2), subfreq(:,3:4,3));
                
%% Add the rest of the subfreq data
%The last column of subfreq = power during movement/power during rest
%Create 15x1 array for this variable (coded AcVsR in SPSS)
%rows 1-5 = M1, each freq band
%rows 6-10 = S1, each freq band
%rows 11-15 = STN, each freq band
 
subfreqAVR = vertcat(subfreq(:,5,1), subfreq(:,5,2), subfreq(:,5,3));

%% Save data to Excel
% variable "status" tells you if it saved properly ("true")
% variable "message" will give you an error message if status=false
% save the file name in the format 'PtName_Limb_Pol' where PtName = first 4
% letters of pt's last name; Limb = hand/elbow/shoulder/jaw/foot; and Pol=
% bipolar/monopolar. If the data was not collected as Recog-Rlfp-Llimb or
% Lecog-Llfp-Rlimb, make appropriate notation in filename.
outputpath = ['C:\Users\Starr\Documents\ECOG data analysis\Excel data\'];
cd(outputpath);
[status, message] = xlswrite(filename, vertallfreq, 'Sheet1', 'A1');
[status, message] = xlswrite(filename, subfreqRM, 'Sheet1', 'D7');
[status, message] = xlswrite(filename, subfreqAVR, 'Sheet1', 'F37');
