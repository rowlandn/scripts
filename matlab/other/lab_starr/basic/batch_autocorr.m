function y = batchautocorr(bfname,lag,ave_width,sample_length)

%	This function calculates autocorrelation functions on all files listed in bfname.
% The files listed in bfname must be listed in 'bare' format.  See the following example:
%
%	!dir *.bis /b /o:n > f.txt
%	batchautocorr('f.txt',1000,10,1);
% 
%	This command sequence will first produce a file f.txt (using the MSDOS 'dir' command)
%	which contains the data file names that will then be used to open the data files and 
%	generate the autocorrelograms with a lag time of 1000 ms, uses a moving average of 10.
%  The time output is based on the assumption that the unit is sample_length ms long.  
%  For each input file an output file is generated that contains the correlogram and the time units. 
%	These files are stored in *.mat format.  The file name is preceeded by a 'AC' (e.g., the input 
%	file H0149.BIS would generate an output file ACH0149.mat.  
%
%	Written 1/12/1999.  Modifications & comments 11/10/2000.  T. Wichmann


fidi = fopen(bfname,'r');
m = 1;
while 1
   linef = fgetl(fidi);
   if ~isstr(linef), 
      break,
   end
   fprintf(1,'Currently processing: %s\r',linef(1:length(linef)-3))
   A = load(linef);
   if isempty(find(A <= 0)) % this if loop inserted because of problems with some MSD files that generate string of zeros or negative numbers.
      [acorr,time] = autocorr(A,lag,ave_width,sample_length);
      save AC acorr time;
      copyfile('AC.mat',['ACH' linef(1:length(linef)-4) '.mat']);
   end
end
fclose('all');