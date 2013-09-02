function [num] = spkanaldata(worksheet)
% This M-file extracts oscillatory analysis data from the SpikeOscil.xls
% file.
%
% Ideas for the future:
%   - automatically generate bar graph

% change directory to one containing the excel file
cd('C:\Documents and Settings\HP_Administrator\My Documents\Lab Documents\data\Spreadsheets\SpikeOscil');

% grab data from excel file
num = xlsread('SpikeOscil.xls',worksheet);

