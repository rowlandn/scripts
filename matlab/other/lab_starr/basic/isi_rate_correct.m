function isi_rate_correct()
% isi_rate_correct()
% This program reads an ISI file created by Labview files (bare interspike intervals)
% and creats a new ISI file with corrected intervals 
%
% Created by RST, 2002-07-29

global RTPATH	% starting directory containing input files

while 1
    p = pwd;
	d = dir;
	str = {d.name};
	[s,v] = listdlg('Name',p,...
                    'PromptString','Select files to process:',...
                    'SelectionMode','multiple',...
                    'ListString',str,...
                    'ListSize', [300 300]);
    if v == 0
        return          % Stop if nothing selected
    end
    if d(s(1)).isdir    % change dir if directory selected
        cd(d(s(1)).name)
    else
        break
    end
end

prompt  = {'Original incorrect rate (i.e. 22.05 for 22.05 kHz)', 'Correct rate (i.e., 20 kHz)'};
title   = 'Input rates';
lines= 1;
def     = {'22.050','20.00'};
answer  = inputdlg(prompt,title,lines,def);
inrate = str2num(answer{1});
outrate = str2num(answer{2});
fix = inrate/outrate;

for i = 1:length(s)
    inname = d(s(i)).name;
    outname = strrep(d(s(i)).name,'.','.fixd.');
    disp(['Converting:  ' inname ' -> ' outname]);
    inra = load(inname);
    outra = fix .* inra;
    
    save(outname,'outra','-ASCII');
end

  