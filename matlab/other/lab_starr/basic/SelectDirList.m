function Lst = SelectDirList()
% This reads a list of directories to be processed into an array of strings.
% The dir /b /s > command must be executed to make up the file list.

global RTPATH
global DIRLIST
global DIRNUM

[filename, RTPATH] = uigetfile('*.lst', 'Select file list');
if (filename ~= 0)
   cd(RTPATH);
   fid = fopen(filename, 'rt');
   
   i = 1;
   while (feof(fid) == 0)
      name = fgetl(fid);			% reads a file name from the file list
      DIRLIST{i} = name;			% inserts name into cell array 
      i = i + 1;
   end   
   
   fclose(fid);
end 

DIRNUM = i-1;

