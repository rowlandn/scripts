function [status, message] = exportRESTdata(allfreq, subfreq, filename)
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
%Creates 3x3 array for MaxFQ, MAXPR, TotPR, such that:
%row 1 = M1, rest
%row 2 = S1, rest
%row 3 = STN, rest

vertallfreq = vertcat(allfreq(:,:,1), allfreq(:,:,2), allfreq(:,:,3));

%% Make subfreq matrix suitable for export
%Creates 15x2 array for POWER and PERCT variables for each brain 
%region (M1, S1, or STN) and each frequency band of interest 
%(delta, lo beta, hi beta, lo gamma, hi gamma)
%rows 1-5 = M1, each freq band
%rows 6-10 = S1, each freq band
%rows 11-15 = STN, each freq band

subfreqRM = vertcat(subfreq(:,:,1), subfreq(:,:,2), subfreq(:,:,3));
                
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
[status, message] = xlswrite(filename, subfreqRM, 'Sheet1', 'D4');

