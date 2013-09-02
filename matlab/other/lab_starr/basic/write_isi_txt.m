function write_isi_txt(t_spk,fname)
% Function to write an text file containing ISIs
% File format emulates that produced by the MSD-Win spike sorter
%
% Enter as arguement:
%
%       t_spk = times of spikes in msec
%       fname = character array containing root 
%                   name for file to be written.
%
% Output:
%       Writes file in current directory with
%           nex.txt appended to fname
%
%   RST March 25, 2002

[i,n_spk] = size(t_spk);
isi_spk = diff(t_spk);
outname = [fname '_nex.txt'];
fid = fopen(outname,'w');
fprintf(fid, '%s\r\n',date);
fprintf(fid, 'number of channels = 1\r\n');
fprintf(fid, 'Sampling Rate (Channel# 1) = 1000\r\n');
fprintf(fid, 'START ISI 14:37:42\r\n');
fprintf(fid, 'Channel Template  ISI\r\n');
for i=1:n_spk-1 
    fprintf(fid,'1\t1\t%.3f\r\n',1000*isi_spk(i));
end
fprintf(fid,'------------------------\r\n');
fprintf(fid,'STOP  ISI 14:38:36\r\n');
fprintf(fid,'End Of The File\r\n');
fprintf(fid, '%s\r\n',date);
fclose(fid);
