function RunPeriORTask(show)
% RunPeriORTask(show)
% This program tests for significant changes in firing around movement
% onset
% Created by RST, 2005-08-22
% Revised by RST, 2006-04-06 - Use OR task events written in MAT file
%
%       Input:
%               show - controls whether graphical output is produced for each cell
%                 default = true.
%
%       Run within a directory and this program will process all spikes in all
%       nex files in the directory
%
%   Requires matching NEX & MAT files originating from abfconv_task.m
%

% Global defines used by the other functions called by this one
global PRE_MVT
global PST_MVT
global BIN_SZ
global CNTL_PER
global MAX_N
global ALPHA
global CONTIG

global SRCH_LO
global SRCH_HI
global VERBOSE

PRE_MVT = -1.5;
PST_MVT = 1;
SMOOTH = 50;    % SD of smoothing gaussian (in msec)
CNTL_PER = [-1.5 -0.5];
ALPHA = 0.05;
CONTIG = 50;    % Mean of 'CONTIG' adjacent points must be significant
JOIN = 0.2;             % Join together responses that are < JOIN sec apart
MIN_SPK_N = 500;
EVENTS = {'home' 'cue' 'stopsignal' 'lift' 'press' };
EVENTS_ENGLISH = {'Home Position' 'Movement Cue' 'Stop Signal' ...
                  'Lifting the Hand' 'Pressing the Key'};
TRIAL_TYPES = {'Go' 'Stop (succeed)' 'Stop (fail)' 'Stop (either)'};

VERBOSE=false;

N_Pks = 2;      % max # of significant changes

if ~exist('show','var')
       show = true;
end

% Edit this pattern to select time-stamp variables w/ specific names
%spkname_pattern = 'Snip\w*[abcd]';     % Kevin's pattern
spkname_pattern = '\w*';        % Accept all units

cd(uigetdir);

FileLst = dir('*.nex');
if(isempty(FileLst))
       str = pwd;
       error(['Found no NEX files in current directory - ' str ]);
end

outfid = write_text_header(N_Pks);      % Subfunction below

n = 0;  % Count of units processed
% For each file found in directory...
for i=1:length(FileLst)

       nexname = FileLst(i).name;
   [pth,fname,ext] = fileparts(nexname);
   matfile = [fname '.mat'];
   if ~exist(matfile,'file')
       error(['Unable to find matching MAT file for ' nexname]);
   end

       % Find variables of interest in file
   info = nex_info_rst(nexname,VERBOSE);
   data_start = 0;
   data_stop = info.dur;
       [spk_inds,spk_names] = find_nex_units( nexname, spkname_pattern);
       if isempty(spk_inds)
               continue
       end

       load( matfile );
       if ~exist('events','var')
               display(['Found no ''events'' variable in file: ' matfile '. Skipping.']);
               continue
       %elseif length(events.light)==1
       %        display(['Found only 1 mvt hand in file: ' fname '. Continuing']);
       end
 
       %A - Left Go with Response
       %B - Right Go with Response
       %C - Left Stop, Successful
       %D - Right Stop, Successful
       %E - Left Stop, Failed
       %F - Right Stop, Failed
       
       L_Go = find(events.type == 'A');
       R_Go = find(events.type == 'B');
       L_StopSucc = find(events.type == 'C');
       R_StopSucc = find(events.type == 'D');
       L_StopFail = find(events.type == 'E');
       R_StopFail = find(events.type == 'F');
       
       % Simple way to select alignment event
       disp('Indicate alignment event:');
       disp('  1) HOME');
       disp('  2) CUE');
       disp('  3) STOP_SIGNAL');
       disp('  4) LIFT off of button');
       disp('  5) PRESS down on button');
       
       evt_n = input('Enter a number (1-5):');
    
       eval(['tmp_evt = events.' EVENTS{evt_n} ';']);

     for type_n = 1:4  
       %disp('Indicate trial types of interest:');
       %disp('  1) GO trials')
       %disp('  2) STOP trials (successful)');
       %disp('  3) STOP trials (failed)');
       %disp('  4) STOP trials (both successful AND failed)');
       
       %type_n = input('Enter a number (1-4):');
       
       switch(type_n)
        case 1
         L_inds = L_Go;
         R_inds = R_Go;
        case 2
         L_inds = L_StopSucc;
         R_inds = R_StopSucc;
        case 3
         L_inds = L_StopFail;
         R_inds = R_StopFail;
        case 4
         try
           L_inds = [L_StopSucc L_StopFail];
           R_inds = [R_StopSucc R_StopFail];
         catch
           L_inds = [L_StopSucc; L_StopFail];
           R_inds = [R_StopSucc; R_StopFail];
         end % try, catch
       end % switch type_n
              
       Mvt(1).ts = tmp_evt(L_inds);
       Mvt(1).n = length(Mvt(1).ts);
       Mvt(2).ts = tmp_evt(R_inds);
       Mvt(2).n = length(Mvt(2).ts);

       for k=1:length(Mvt)         % Filter out events too close to data boundaries
               drp = find(Mvt(k).ts < data_start-PRE_MVT | Mvt(k).ts > data_stop-PST_MVT);
               Mvt(k).ts(drp) = [];
               Mvt(k).n = length(Mvt(k).ts);
       end

       % For each unit in a file...
       for j=1:length(spk_inds)
               n = n+1;
               spk{n}.fname = fname;
               spk{n}.align_evt = EVENTS_ENGLISH{evt_n};
               spk{n}.trial_type = TRIAL_TYPES{type_n};

               % save name of file & unit
               spk{n}.unitname = deblank( spk_names(j,:) );
               display(['Processing...' spk{n}.fname '....:' spk{n}.unitname]);

               % get spike times
               [spk_n, spk_t] = nex_ts( nexname, spk_names(j,:),VERBOSE);
               if spk_n < MIN_SPK_N
                       display(['Found only ' num2str(spk_n) ' spikes. Skipping.']);
                       continue
               end

               % for each hand of movement
               for k=1:length(Mvt)
                       spk{n}.hand(k).n_reps = Mvt(k).n;


                       [spk{n}.hand(k).histog, spk{n}.bins] = perievent_sdf(spk_t, Mvt(k).ts, ...
                               PRE_MVT, PST_MVT, SMOOTH );

                       spk{n}.hand(k).raster = perievent_raster(spk_t, Mvt(k).ts, ...
                               PRE_MVT, PST_MVT);

                       cntl_inds = find(spk{n}.bins>=CNTL_PER(1) & spk{n}.bins<CNTL_PER(2));
                       test_start = max(cntl_inds)+1;

                       [spk{n}.hand(k).chng(1),spk{n}.hand(k).cntl_mean,spk{n}.hand(k).sig_thr] = ...
                               PeriEventChange_SDF(spk{n}.hand(k).histog,cntl_inds,test_start,ALPHA,CONTIG);

                       sgn1 = spk{n}.hand(k).chng(1).sgn;
                       sgn2 = [];
                       del2nd = 0;
                       if ~isempty( spk{n}.hand(k).chng(1).off_ind )
                               test_start = spk{n}.hand(k).chng(1).off_ind;
                               [spk{n}.hand(k).chng(2),spk{n}.hand(k).cntl_mean,q] = ...
                                       PeriEventChange_SDF(spk{n}.hand(k).histog,cntl_inds,test_start,ALPHA,CONTIG);
                               sgn2 = spk{n}.hand(k).chng(2).sgn;
                       end
                       if sgn2==sgn1
                               off1_t = spk{n}.bins( spk{n}.hand(k).chng(1).off_ind ) ;
                               on2_t = spk{n}.bins( spk{n}.hand(k).chng(2).on_ind ) ;
                               if (on2_t - off1_t)< JOIN
                                       on1_ind = spk{n}.hand(k).chng(1).on_ind ;
                                       off1_ind = spk{n}.hand(k).chng(2).off_ind ;
                                       spk{n}.hand(k).chng(1).off_ind = off1_ind;
                                       his = spk{n}.hand(k).histog - spk{n}.hand(k).cntl_mean;
                                       if ~isempty(off1_ind)
                                               spk{n}.hand(k).chng(1).mean_change = ...
                                                       mean( his(on1_ind:off1_ind) );
                                               spk{n}.hand(k).chng(1).int_change = ...
                                                       sum( his(on1_ind:off1_ind) )/1000;
                                       else
                                               spk{n}.hand(k).chng(1).mean_change = ...
                                                       mean( his(on1_ind:end) );
                                               spk{n}.hand(k).chng(1).int_change = ...
                                                       sum( his(on1_ind:end) )/1000;
                                       end
                                       del2nd = 1;
                               end
                       end
                       if isempty( spk{n}.hand(k).chng(1).off_ind ) | del2nd
                               spk{n}.hand(k).chng(2).on_ind = [];
                               spk{n}.hand(k).chng(2).sgn = [];
                               spk{n}.hand(k).chng(2).off_ind = [];
                               spk{n}.hand(k).chng(2).mean_change = [];
                               spk{n}.hand(k).chng(2).int_change = [];
                       end
               end
               if( show)
                       make_figure(spk{n},N_Pks);
%pause
               end

               % write stats to file
               write_text(outfid,spk{n},N_Pks);
       end
     end
     
end
fclose(outfid);
save PeriORTask spk

%------------------------------------------------------
% Subfunction to make figure of results
function make_figure(s, N_Pks)
       bins = s.bins;

       %%%%%%%%%%%%%% Plotting
       % Set up axes
       MARGIN = 0.06;
       TOP = 1-MARGIN;         % Top margin of page
       LEFT = MARGIN;  % Left margin of page
       WIDTH = (1-3*MARGIN)/2; % give space for 3 margin widths including middle
       HEIGHT = 0.4;   % Height of histograms
       CLR = [0.25,0.25,0.25 ; 0.75,0.75,0.75];

       figure
       set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
       % Size to make it look good
       c = get(gcf);
       c.Position(2) = 275;
       c.Position(3) = 870;
       c.Position(4) = 680;
       set(gcf,'Position',c.Position);

       % Find max across hands
       for j = 1:length(s.hand)
               x(j) = max(s.hand(j).histog);
               ymax = max(x)+5;
       end

       % Plot for 2 movements
       for j=1:length(s.hand)
               mvt = s.hand(j);
               left = MARGIN + (WIDTH+MARGIN)*(j-1);
               width = WIDTH;
               height = HEIGHT;
               bottom = TOP-height;
               subplot('position',[left bottom width height]);

               h=area(bins, mvt.histog);
               set(h,'FaceColor',[0.5,0.5,0.5],'EdgeColor','k');
               xlim([min(bins) max(bins)]);
               if isnan(ymax)
                 ymax = 50;
               end
               
               ylim([0 ymax]);
               ylm = ylim;
               xlm = xlim;
               hold on
               plot(xlim,[mvt.cntl_mean mvt.cntl_mean],'k-');
               plot(xlim,[mvt.cntl_mean+mvt.sig_thr mvt.cntl_mean+mvt.sig_thr],'k:');
               plot(xlim,[mvt.cntl_mean-mvt.sig_thr mvt.cntl_mean-mvt.sig_thr],'k:');
               plot([0,0],ylm,'k-');
               xlabel('seconds');
               ylabel('spikes/sec');

               if isempty(mvt.chng(1).on_ind)
                       text( mean([max(bins) min(bins)]), ylm(2)/2, 'No sig change found',...
                                       'HorizontalAlignment','center','Color','r');
               else
                       for i = 1:2
                               chng = mvt.chng(i);

                               if ~isempty(bins(chng.on_ind))
                                       on = bins(chng.on_ind);
                                       if ~isempty(bins(chng.off_ind))
                                               off = bins(chng.off_ind);
                                       else
                                               off = xlm(2);   %If no offset found, change lasts to end
                                       end
                                       x = [ on on off off ];
                                       y = [ ylm ylm(2:-1:1)];
                                       fill(x,y,CLR(i,:),'EdgeColor','none','FaceAlpha',0.3)
                               end
                       end
               end


               if j ==1
                       str = [ 'Datafile ' s.fname ', aligned on ' ...
                               s.align_evt ' for trial type ' s.trial_type];
                       title(str,'Interpreter','none','FontSize',12,...
                               'Position',[xlm(2),ylm(2)+5,0]);
               end

               bottom = bottom-height-MARGIN;
               subplot('position',[left bottom width height]);
               rasterplot(mvt.raster,0.9,'');
       end

return


%------------------------------------------------------
% Subfunction to write a line of data to output file
function write_text(outfid, spk,N_Pks)
       bins = spk.bins;

       fprintf(outfid,'%s\t%s\t', ...
               spk.fname, spk.unitname);

       for j = 1:length(spk.hand)       % movement hand
               mvt = spk.hand(j);
               fprintf(outfid,'%d\t%.3f\t%.3f\t',...
                        mvt.n_reps, mvt.cntl_mean, mvt.sig_thr );

               for i = 1:N_Pks
                       chng = mvt.chng(i);

                       % Report significant peaks
                       if ~isempty(chng.on_ind)
                               fprintf(outfid,'%.3f\t%.3f\t%.3f\t%.3f\t',...
                                       bins(chng.on_ind),chng.sgn,chng.mean_change,chng.int_change);
                       else
                               fprintf(outfid,'-\t-\t-\t-\t');
                       end
                       if ~isempty(chng.off_ind)
                               fprintf(outfid,'%.3f\t',bins(chng.off_ind));
                       else
                               fprintf(outfid,'-\t');
                       end
               end
       end
       fprintf(outfid,'\n');
       return

%------------------------------------------------------
% Subfunction to open output file and print header line
function outfid = write_text_header(N_Pks)
       fname = 'PeriORTask.txt';
       outfid = fopen(fname,'w');
       if(outfid == -1)
      error(['Unable to open...' fname ]);
       end
       fprintf(outfid,'fname\tunitname\t');
       for j = 1:2
               fprintf(outfid,'Nreps\tcntl_mean\tsig_thresh\t');
               for i = 1:N_Pks
                       % For max reported signif acorr pks, freq & normalized power
                       fprintf(outfid,'Onset%d\tSign%d\tMeanChange%d\tIntChange%d\tOffset%d\t',i,i,i,i,i);
               end
       end
       fprintf(outfid,'\n');           % EOL

       return



function info = nex_info_rst(filename, verbose)
% info = nex_info_rst(filename, verbose) -- read and display .nex file info
%   Enhanced by RST for extra info
%
% info = nex_info+(filename, verbose)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%       verbose - optional flag to control display of results (default=1)
% OUTPUT:
%   info.nvar - number of variables in the file
%   info.names - [nvar 64] array of variable names
%   info.types - [1 nvar] array of variable types
%           Interpretation of type values: 0-neuron, 1-event, 2-interval, 3-waveform,
%                        4-population vector, 5-continuous variable, 6 - marker
%   info.ver - version
%   info.fs - sampling fequency
%   info.tbeg - data begin sample
%   info.tend - data end sample
%   info.dur - data duration (in sec)

if ~exist('verbose','var')
       verbose = 1;
end

if(nargin<1 | nargin>2 )
  disp('1 or 2 inputs arguments are required')
  return
end

if(length(filename) == 0)
  [fname, pathname] = uigetfile('*.nex', 'Select a Nex file');
       filename = strcat(pathname, fname);
end

fid = fopen(filename, 'r');
if(fid == -1)
       disp('cannot open file');
  return
end

magic = fread(fid, 1, 'int32');
version = fread(fid, 1, 'int32');
comment = fread(fid, 256, 'char');
freq = fread(fid, 1, 'double');
tbeg = fread(fid, 1, 'int32');
tend = fread(fid, 1, 'int32');
nvar = fread(fid, 1, 'int32');
fseek(fid, 260, 'cof');
if verbose
       disp(strcat('file = ', filename));
       disp(strcat('version = ', num2str(version)));
       disp(strcat('frequency = ', num2str(freq)));
       disp(strcat('duration (sec) = ', num2str((tend - tbeg)/freq)));
       disp(strcat('number of variables = ', num2str(nvar)));
end
names = zeros(1, 64);
for i=1:nvar
       types(i) = fread(fid, 1, 'int32');
       var_version = fread(fid, 1, 'int32');
       names(i, :) = fread(fid, [1 64], 'char');
       dummy = fread(fid, 128+8, 'char');
end
names = setstr(names);
fclose(fid);

  info.nvar = nvar;    % number of variables in the file
  info.names = names;  %- [nvar 64] array of variable names
  info.types = types;  %- [1 nvar] array of variable types
  info.ver = version;  %- version
  info.fs = freq;      %- sampling fequency
  info.sbeg = tbeg;    %- begin sample
  info.send = tend;    %- end sample
  info.dur = (tend - tbeg)/freq;   % data duration


function        [var_inds, var_names] = find_nex_units( fname, spkname_pattern);
% function      [var_inds, var_names] = find_nex_units( fname, spkname_pattern)
%
% Function to find the indices of NEX timestamp variables that match the
% indicated name_pattern
%
% INPUT:
%   filename - NEX File to query
%       spkname_pattern - string  for variable name matching
%                       See also regexp on string-pattern matching in matlab
% OUTPUT:
%   var_inds - indices of the variables whose names match the name pattern
%   var_names - names of the variables whose names match the name pattern
%

SPIKE_TYPE = 0; % Works only on Spike-type timestamp

if exist(fname) ~= 2
       error( ['Unable to find the file named:  ' fname ] );
end

[nvar, varname, types] = nex_info(fname,0);

var_inds = find(types == SPIKE_TYPE);
n_vars = length(var_inds);

for i=1:n_vars          % Find spikes that match desired name pattern
       tmpname = deblank( varname(var_inds(i),:)) ;
       if isempty( regexpi( tmpname, spkname_pattern) )
               var_inds(i) = NaN;      % If names don't match, mark as bad
       end
end
var_inds( isnan(var_inds) ) = [];
var_names = varname( var_inds, : );

function [nvar, names, types] = nex_info(filename, verbose)
% nex_info(filename) -- read and display .nex file info
%
% [nvar, names, types] = nex_info(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%       verbose - optional flag to control display of results (default=1)
% OUTPUT:
%   nvar - number of variables in the file
%   names - [nvar 64] array of variable names
%   types - [1 nvar] array of variable types
%           Interpretation of type values: 0-neuron, 1-event, 2-interval, 3-waveform,
%                        4-population vector, 5-continuous variable, 6 - marker

if ~exist('verbose','var')
       verbose = 1;
end

if(nargin<1 | nargin>2 )
  disp('1 or 2 inputs arguments are required')
  return
end

if(length(filename) == 0)
  [fname, pathname] = uigetfile('*.nex', 'Select a Nex file');
       filename = strcat(pathname, fname);
end

fid = fopen(filename, 'r');
if(fid == -1)
       disp('cannot open file');
  return
end

magic = fread(fid, 1, 'int32');
version = fread(fid, 1, 'int32');
comment = fread(fid, 256, 'char');
freq = fread(fid, 1, 'double');
tbeg = fread(fid, 1, 'int32');
tend = fread(fid, 1, 'int32');
nvar = fread(fid, 1, 'int32');
fseek(fid, 260, 'cof');
if verbose
       disp(strcat('file = ', filename));
       disp(strcat('version = ', num2str(version)));
       disp(strcat('frequency = ', num2str(freq)));
       disp(strcat('duration (sec) = ', num2str((tend - tbeg)/freq)));
       disp(strcat('number of variables = ', num2str(nvar)));
end
names = zeros(1, 64);
for i=1:nvar
       types(i) = fread(fid, 1, 'int32');
       var_version = fread(fid, 1, 'int32');
       names(i, :) = fread(fid, [1 64], 'char');
       dummy = fread(fid, 128+8, 'char');
end
names = setstr(names);
fclose(fid);

function [n, ts] = nex_ts(filename, varname, verbose)
% nex_ts(filename, varname): Read timestamps from a .nex file
%
% [n, ts] = nex_ts(filename, varname)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   varname - variable name
% OUTPUT:
%   n - number of timestamps
%   ts - array of timestamps (in seconds)

n = 0;
ts = 0;

if(nargin < 2 | nargin > 3 )
  disp('2 or 3 input arguments are required')
  return
end

if(ischar(filename) == 0)
  disp('input arguments should be character arrays')
  return
end

if(ischar(varname) == 0)
  disp('input arguments should be character arrays')
  return
end

if(length(filename) == 0)
  [fname, pathname] = uigetfile('*.nex', 'Select a Nex file');
       filename = strcat(pathname, fname);
end
if ~exist('verbose','var')
       verbose = 1;
end

fid = fopen(filename, 'r');
if(fid == -1)
       disp('cannot open file');
  return
end

magic = fread(fid, 1, 'int32');
version = fread(fid, 1, 'int32');
comment = fread(fid, 256, 'char');
freq = fread(fid, 1, 'double');
tbeg = fread(fid, 1, 'int32');
tend = fread(fid, 1, 'int32');
nvar = fread(fid, 1, 'int32');
fseek(fid, 260, 'cof');
name = zeros(1, 64);
found = 0;
for i=1:nvar
       type = fread(fid, 1, 'int32');
       var_version = fread(fid, 1, 'int32');
       name = fread(fid, [1 64], 'char');
       offset = fread(fid, 1, 'int32');
       n = fread(fid, 1, 'int32');
       name = setstr(name);
       name = deblank(name);
       k = strcmp(name, deblank(varname));
       if(k == 1)
               found = 1;
               fseek(fid, offset, 'bof');
               ts = fread(fid, [1 n], 'int32');
               break
       end
       dummy = fread(fid, 128, 'char');
end

fclose(fid);

if found == 0
       disp('did not find variable in the file');
else
       ts = ts/freq;
       if verbose
               disp(strcat('file = ', filename));
               disp(strcat('number of timestamps = ', num2str(n)));
       end
end

function [sdf, bin_edges] = perievent_sdf(spk_t, event_t, pre, pst, smooth)
% Construct peri-event SDF from vector of spike times
%
% Inputs
%       spk_t - vector of spike times
%       event_t - vector of event times
%       pre - time before event at which to start histogram (negative value
%       assummed)
%       pst - time following event at which to stop histogram
%       smooth - sigma of gaussian (in msec)
%
%
% Output
%       sdf - mean rate of spike events convolved w/ gaussian kernel & expressed in spikes/sec
%       bin_edges - matching vector of independent sample intervals
%
% RST 2005-08

bin_edges = pre:0.001:pst;
data_len = length(bin_edges);

raster = perievent_raster(spk_t, event_t, pre, pst);
delta = spk_t2delta( raster, data_len);
sdf = delta2sdf(delta,smooth);


function raster = perievent_raster(spk_t, event_t, pre, pst)
% PERI-EVENT RASTER
% function raster = perievent_raster(spk_t, event_t, pre, pst)
%
%       Make peri-event raster array of spike times - 1 row per trial
%       NaN pads variable length rows

n_events = length(event_t);
raster = NaN*zeros(1,1);
rst_len = size(raster,2);
for i = 1:n_events
       pret = event_t(i)+pre;
       pstt = event_t(i)+pst;
       rstt = spk_t( find( spk_t>=pret & spk_t<pstt ) ) - event_t(i);
       rsttlen = size(rstt,2);
       % Resize raster array if necessary
       if rsttlen > rst_len
               raster(:,(rst_len+1:rsttlen)) = NaN;
               rst_len = size(raster,2);
       elseif rsttlen < rst_len
               rstt(rsttlen+1:rst_len) = NaN;
       end
       raster = [ raster; rstt ];
end
% Delete 1st dummy row when finished
raster(1,:) = [];

function spk_delta = spk_t2delta(spk_t,data_len,msec)
% Function to produce a delta function from a list of spike times
%
%       function spk_delta = spk_t2delta(spk_t,data_len,msec)
%
% Enter as arguement:
%
%       t_spk = times of spikes in seconds - may be a vector or array
%                       arrays assume rows = trials
%       data_len = length of delta function to be returned (in msec - optional)
%               msec = flag TRUE indicates that input is in msec (optional, default
%                       is seconds)
%       msec = optional flag to indicate spike times are in milliseconds
%
% Returns:
%       spk_delta = msec resolution delta function
%                   (spike = 1, no-spike = 0)
%
%   RST March 25, 2002
%       Oct 03, allow option of no data len
%       RST Nov, 2003   allow array of spike times
%       RST Aug, 2005 check for and correct negative spike times
%
seconds = true;

% Remove negative offset if present
mn = min(min(spk_t));
if mn<0
       spk_t = spk_t -mn;
end

% Discover data length (in msec) if not specified
if nargin < 2
       data_len = 1+round(1000*max(nanmax(spk_t,[],2)));
end

if nargin ==3
       if msec,        seconds = false;        end
end

[n_row,width] = size(spk_t);
spk_delta = zeros(n_row,data_len); % make delta fn as long as required

for i = 1:n_row
       % Convert spike train to delta function
       if seconds
               tmp_spk = 1+round(1000*spk_t(i,:));        % put spk time in msec starting @ 1 msec
       else
               tmp_spk = 1+spk_t(i,:);
       end
       tmp_spk(find(isnan(tmp_spk))) = [];             % delete NaNs
       spk_delta(i,squeeze(tmp_spk)) = 1;              % code times of spikes as 1's
end



function [varargout] = nanmax(varargin)
%NANMAX Maximum value, ignoring NaNs.
%   M = NANMAX(A) returns the maximum of A with NaNs treated as missing.
%   For vectors, M is the largest non-NaN element in A.  For matrices, M is
%   a row vector containing the maximum non-NaN element from each column.
%   For N-D arrays, NANMAX operates along the first non-singleton
%   dimension.
%
%   [M,NDX] = NANMAX(A) returns the indices of the maximum values in A.  If
%   the values along the first non-singleton dimension contain more than
%   one maximal element, the index of the first one is returned.
%
%   M = NANMAX(A,B) returns an array the same size as A and B with the
%   largest elements taken from A or B.  Either one can be a scalar.
%
%   [M,NDX] = NANMAX(A,[],DIM) operates along the dimension DIM.
%
%   See also MAX, NANMIN, NANMEAN, NANMEDIAN, NANMIN, NANVAR, NANSTD.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.12.2.3 $  $Date: 2004/06/25 18:52:54 $

% Call [m,ndx] = max(a,b) with as many inputs and outputs as needed
[varargout{1:nargout}] = max(varargin{:});

function sdf = delta2sdf(delta,gaus_sigma)
% sdf = delta2sdf(delta,gaus_sigma)
%
% Convert delta array to mean spike density function (sdf)
%  Argument:    delta = trial X msec array,  1 denoting spike events
%                               gaus_sigma = width of gaussian (in msec; optional, default = 20)
%
%  Returns:             sdf = mean spike density function across trials
%
%       04-01-27        Revised to take delta in either column or row vector
%       Revised 2005-08-24 to make gaussian on the fly - any sigma can be used
%                       also speeded execution by first summing delta and convolving
%                       one summed-delta stream.

[ntrial spk_len] = size(delta);

% Make the filter
if ~exist('gaus_sigma','var')
       gaus_sigma = 20;
end
flt_len = 5*gaus_sigma; % usually 5*stdev is long enough to make is smooth
x = -1*flt_len:flt_len;
flt = normpdf(x,0,gaus_sigma);
% Make sure edges are smooth & scaling correct (weights should sum to 1000)
flt = flt - min(flt);
flt = 1000.*flt ./(sum(flt)*ntrial);

flt_len = (length(flt)-1)/2;

siz = size(delta);

if min(siz) == 1
       ntrial = 1;
       spk_len = max(siz);
       if siz(2) == 1
               delta = delta';
       end
else
       ntrial = siz(1);
       spk_len = siz(2);
end

% sum the delta function across trials
sdelta = sum(delta,1);

% indices into right and left ends of original delta vector
revlist1 = flt_len+1:-1:2;
revlist2 = spk_len-1:-1:spk_len-flt_len;

% build a mirror-padded sdelta
s = [ sdelta(revlist1) sdelta sdelta(revlist2) ];

% convolve w/ filter
sdf = conv(flt, s );

% remove pads
sdf(1: 2*flt_len)=[];
sdf(spk_len+1: spk_len + 2*flt_len)=[];

if siz(2) == 1
       sdf = sdf';
end

function y = normpdf(x,mu,sigma)
%NORMPDF Normal probability density function (pdf).
%   Y = NORMPDF(X,MU,SIGMA) returns the pdf of the normal distribution with
%   mean MU and standard deviation SIGMA, evaluated at the values in X.
%   The size of Y is the common size of the input arguments.  A scalar
%   input functions as a constant matrix of the same size as the other
%   inputs.
%
%   Default values for MU and SIGMA are 0 and 1 respectively.
%
%   See also NORMCDF, NORMFIT, NORMINV, NORMLIKE, NORMRND, NORMSTAT.

%   References:
%      [1] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.10.4.2 $  $Date: 2004/08/20 20:06:04 $

if nargin<1
   error('stats:normpdf:TooFewInputs','Input argument X is undefined.');
end
if nargin < 2
   mu = 0;
end
if nargin < 3
   sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
   y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
catch
   error('stats:normpdf:InputSizeMismatch',...
         'Non-scalar arguments must match in size.');
end

function [x,xlo,xup] = norminv(p,mu,sigma,pcov,alpha)
%NORMINV Inverse of the normal cumulative distribution function (cdf).
%   X = NORMINV(P,MU,SIGMA) returns the inverse cdf for the normal
%   distribution with mean MU and standard deviation SIGMA, evaluated at
%   the values in P.  The size of X is the common size of the input
%   arguments.  A scalar input functions as a constant matrix of the same
%   size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   [X,XLO,XUP] = NORMINV(P,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
%   for X when the input parameters MU and SIGMA are estimates.  PCOV is a
%   2-by-2 matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%   bounds.  XLO and XUP are arrays of the same size as X containing the lower
%   and upper confidence bounds.
%
%   See also ERFINV, ERFCINV, NORMCDF, NORMFIT, NORMLIKE, NORMPDF,
%            NORMRND, NORMSTAT.

%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, 1046pp., sections 7.1, 26.2.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.16.4.2 $  $Date: 2004/08/20 20:06:03 $

if nargin<1
   error('stats:norminv:TooFewInputs','Input argument P is undefined.');
end
if nargin < 2
   mu = 0;
end
if nargin < 3
   sigma = 1;
end

% More checking if we need to compute confidence bounds.
if nargout>2
  if nargin<4
     error('stats:norminv:TooFewInputs',...
           'Must provide covariance matrix to compute confidence bounds.');
  end
  if ~isequal(size(pcov),[2 2])
     error('stats:norminv:BadCovariance',...
           'Covariance matrix must have 2 rows and columns.');
  end
  if nargin<5
     alpha = 0.05;
  elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
     error('stats:norminv:BadAlpha',...
           'ALPHA must be a scalar between 0 and 1.');
  end
end

% Return NaN for out of range parameters or probabilities.
sigma(sigma <= 0) = NaN;
p(p < 0 | 1 < p) = NaN;

x0 = -sqrt(2).*erfcinv(2*p);
try
   x = sigma.*x0 + mu;
catch
   error('stats:norminv:InputSizeMismatch',...
         'Non-scalar arguments must match in size.');
end

% Compute confidence bounds if requested.
if nargout>=2
  xvar = pcov(1,1) + 2*pcov(1,2)*x0 + pcov(2,2)*x0.^2;
  if any(xvar<0)
     error('stats:norminv:BadCovariance',...
           'PCOV must be a positive semi-definite matrix.');
  end
  normz = -norminv(alpha/2);
  halfwidth = normz * sqrt(xvar);
  xlo = x - halfwidth;
  xup = x + halfwidth;
end

function [s,cntl_mean,sig_thr] = PeriEventChange_SDF(data,cntl_inds,test_start,alpha,contig)
%function s = PeriEventChange_SDF(data,cntl_inds,test_start,alpha,contig)
%
%       Finds significant changes in peri-event SDF using
%       statistics.
%
% Inputs:
%       data - SDF vector
%       cntl_inds - indices into data for control period
%       test_start - index at which to start testing
%       alpha - p-value for significance
%       contig - number of contiguous significant points required
%
% Outputs:
%       s.on_ind - index into data for onset of significant response (empty if
%                               no significant change found)
%       s.off_ind - index into data for offset of significant response
%       s.sgn           - sign of the significant change
%       s.mean_change - mean change from control across duration of change (till
%               offset or end of data)
%       cntl_mean - mean control level
%       sig_thr - threshold to determine significance
%
%       RST 2005-08-10
%
cntl_data = data(cntl_inds);
cntl_mean = mean( cntl_data );
data = data - cntl_mean;

% make testing data
test = data( test_start:end );

p = alpha/(length(cntl_data)/contig);   % Compensate for multcomp
sig_thr = norminv(1-p) * std( cntl_data );

% Find threshold crossing indices
thr_inds = find(test > sig_thr | test < -1*sig_thr );

% Compute o'th-order differences - only o-successive indices will place
% o'th order difference o indices apart
o = contig-1;
odiff = thr_inds((1+o):end)-thr_inds(1:(end-o));

% Find first index into odiff that is o indices apart
s.on_ind = thr_inds( min(find( odiff == o )) ) + test_start-1;

if ~isempty( s.on_ind )
       s.sgn = sign(data(s.on_ind));

       % search smoothed data after established onset
       test = data(s.on_ind+1:end);

       % find returns to baseline
       thr_inds = find(test < sig_thr & test > -1*sig_thr);

       % check if it lasts o-contiguous points
       odiff = thr_inds((1+o):end)-thr_inds(1:(end-o));
       zero_ind = thr_inds( min( find(odiff==o) )) + s.on_ind;

       % find any change in sign > sig_thr
       switch_ind = min(find( sign(test) ~= s.sgn & abs(test)>sig_thr )) + s.on_ind;

       s.off_ind = min( [zero_ind switch_ind] );
       if isempty(s.off_ind)
               s.mean_change = mean( data(s.on_ind:end));
               s.int_change = sum( data(s.on_ind:end))/1000;
       else
               s.mean_change = mean( data(s.on_ind:s.off_ind));
               s.int_change = sum( data(s.on_ind:s.off_ind))/1000;
       end

else
       s.sgn = [];
       s.off_ind = [];
       s.mean_change = [];
       s.int_change = [];
end




function rasterplot(spike_mat,tickheight,txt,xlm,color)
% Create raster plot
%     rasterplot(spike_mat,tickheight,txt)
%       Arguements:
%               spike_mat = 2-dimensional array of spike times
%                       1 row per trial
%                       NaN pads variable length rows
%               tickheight = height of '|' markers as fraction of 1
%               txt = text to print as title of plot
%               xlm - optional x-limits, used if defined
%               color - optional color (default = black)

if ~exist('xlm','var')
       xlm = [min(min(spike_mat)) max(max(spike_mat))];
end
 if ~exist('color','var')
    f = get(gcf);
    if sum(f.Color) == 0
       color = 'w';
    else
       color = 'k';
    end
 end

missing = find( isnan(spike_mat(:,1)) );
spike_mat(missing,:) = [];

[ntrials,maxspks] = size(spike_mat);

trial_n = zeros(ntrials,maxspks);

for i = 1:ntrials
       trial_n(i,:) = i;
end

% transform spike times & trial numbers into row vectors
x = reshape(spike_mat',1,[]);
y = reshape(trial_n',1,[]);

newx = zeros(2,length(x));
newx(1,:) = x;
newx(2,:) = x;
newy = zeros(2,length(y));
newy(1,:) = y;
newy(2,:) = y + tickheight;
h = plot(newx,newy, 'LineWidth', 0.01,'color', color);
if ntrials
  ylim([tickheight ntrials+tickheight]);
else
  ylim(tickheight + [0 1]);
end
x = get(gca);
ind = find( rem(x.YTick,1)==0 );
ytk = x.YTick(ind);
set(gca,'YTick',ytk);

if ~isempty(xlm)
  xlim( xlm );
else
  xlim( [0 1] );
end

title(txt);
