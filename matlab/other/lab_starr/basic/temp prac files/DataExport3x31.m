%% Notes
% DataExport is a program to concatenate all the data from allfreq and
% subfreq matrices into a 3x31 array that will be saved as an excel
% spreadsheet. The excel spreadsheet can then be copy/pasted into SPSS.
% This way there is limited time spent on data entry.

%4/3/2009 ALC


%% linear indexing and transposition of allfreq in preparation for concatenation 
lidxallfreq = allfreq(:);

lidxallfreq = lidxallfreq';
%% Create matrix with concatenated data from allreq and subfreq
% concatenate data into variable "data4export." This is in a format that 
% can be copy/pasted directly into the SPSS spreadsheet
data4export = [lidxallfreq(1:6) subfreq(1,:,1) subfreq(2,:,1) subfreq(3,:,1) subfreq(4,:,1) subfreq(5,:,1);...
lidxallfreq(7:12) subfreq(1,:,2) subfreq(2,:,2) subfreq(3,:,2) subfreq(4,:,2) subfreq(5,:,2);...
lidxallfreq(13:18) subfreq(1,:,3) subfreq(2,:,3) subfreq(3,:,3) subfreq(4,:,3) subfreq(5,:,3)];

%% Save data to Excel
% variable "status" tells you if it saved properly ("true")
% variable "message" will give you an error message if status=false
% save the file name in the format 'PtName_Limb_Pol' where PtName = first 4
% letters of pt's last name; Limb = hand/elbow/shoulder/jaw/foot; and Pol=
% bipolar/monopolar. If the data was not collected as Recog-Rlfp-Llimb or
% Lecog-Llfp-Rlimb, make appropriate notation in filename.

[status, message] = xlswrite('PtName_Limb_Pol', data4export);