% This program analyses performance of tracking
% Created by RST, 11-5-99

global RTPATH		% starting directory containing dirlist file
global DIRNUM		% number of directories in list file
global DIRLIST	% list of directory names to be processed
global XPOS		% X posn, column 2 = Target
global TARG		% Target position
global TIME
global TVel
global XVel
global ncoeffs
global tsign
global nrevs
global trev
global revpos
global revtim
global pkvel
global srchrng
global mvtamp
global onset
global offset
global plothand
global fig_time
global fig_freq
global Hd
global VELTHRESH

global perror
global terror
global mamp
global mvel
global totmt
global mvtn
global assym
global task
global t_freq
global t_pow

p = path;
newp = 'c:\rst\pet\perform\matanal';
if( isempty(findstr(newp,p)) )
   addpath(newp);
end

colordef black;
scrsz = get(0,'ScreenSize');
fig_time = figure('MenuBar','none');
fig_freq = figure('MenuBar','none');

VELTHRESH = 1.0;

SelectDirList;	% List of subdirectories to process
n = 0;
while( n < DIRNUM )
   n = n + 1;

   cd([RTPATH DIRLIST{n}]);
   Lst = dir('scan*');
   if(length(Lst) > 10)
      error(['Found > 10 scan files in ' DIRLIST{n} ]);
   end
      
   outfile = MesFile(DIRLIST{n});	% Returns handle for *.mes file
   % Write header info
   fprintf(outfile,'name\ttask\tscan\ttdist\tperr\tterr\tmamp\tmvel\ttotmt\tmvtn\n');
   m = 0;
   while( m < length(Lst))
      m = m+1;
      XPOS = zeros(1,14999);
      TARG = zeros(1,14999);
      TIME = [];
      TVel = [];
      XVel = [];
      trev = [];
      revpos = [];
      revtim = [];
      mvtamp = [];
      onset = [];
      offset = [];
      scan = [];
      perror = [];
      terror = [];
      mamp = [];
      mvel = [];
      totmt = [];
      mvtn = [];
      assym = 0;
      t_freq = 0;
      t_pow = 0;
      
      disp( ['Reading... ' Lst(m).name] );
      Hd = readfile(Lst(m).name);			% XPOS & TARG are filled and Header information is returned
      if( m==1 )
         db(n).subj = Hd.name;
      end
      if( Hd.rand>0 )
         if( Hd.curs>=1 )
            task = ['step' int2str(Hd.tdist)];
         else
            task = 'stepeye';
         end
      else
         if( Hd.curs>=1 )
            task = ['trac' int2str(Hd.tdist)];
         else
            task = 'traceye';
         end
      end
      disp( ['Processing... ' Lst(m).name '...' task] );
      
      srchrng = Hd.interv/Hd.sampl;
      
      VelFilter20;						% Get velocity weighting fn
      XVel = XPOS - XPOS(1);			% Remove offset to eliminate initial transient
      TVel = TARG - TARG(1);
      XVel = 12.55 * filter(b, .1, XVel);	% Compute velocity calibrated to cm/sec
      if(Hd.rand == 0)					% TVel only really useful for trac trials
         TVel = 12.55 * filter(b, .1, TVel);
      else
         TVel = filter(b, .1, TVel);
      end
      
      % Redimension everything to eliminate vel phase lag introduced by filter
      ncoeffs = floor(length(b)/2);
      TARG(:,1+end-ncoeffs:end) = [];
      XPOS(:,1+end-ncoeffs:end) = [];
      XVel(:,1:ncoeffs) = [];
      TVel(:,1:ncoeffs) = [];
      TIME = 0:0.004:0.004*(length(TARG)-1);
      
      % Plot basic performance stuff
      figure(fig_time);
      if(not(isempty(get(gcf,'CurrentAxes'))))
         delete(gca);
      end
      if(not(isempty(get(gcf,'CurrentAxes'))))
         delete(gca);
      end
      subplot('Position',[0 0.5 1 0.5]);
      line([0 60],[0 0],'Color',[0.2 0.2 0.2]);
      hold on
      plot(TIME,TARG,':b','LineWidth',0.2);
      plot(TIME,XPOS,'y');
      hold off
      subplot('Position',[0 0 1 0.5]);
      line([0 60],[0 0],'Color',[0.2 0.2 0.2]);
      hold on
      plot(TIME,TVel,':b','LineWidth',0.2);
      plot(TIME,XVel,'y');
      hold off
      plothand = get(gcf,'Children');	% Remember which plot is which - top is #2
      
      if(Hd.rand == 0)		% Do analysis for smooth tracking condition
         TracPerf;
      end
      if(Hd.rand > 0)		% Do analysis for step tracking condition
         StepPerf;
      end
%      pause;
      
      mvtn = 90*mvtn/Hd.collecsec;	% Calculate total mvts 
      
      % Output results of analysis to individual files for each subdir
      scan = sscanf(lower(Lst(m).name),'scan%d');
      fprintf(outfile,'%s\t%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.0f\n',...
         char(Hd.name),task,scan,Hd.tdist,perror,terror,mamp,mvel,totmt,mvtn);
      
      % Also save results in structure for whole data run
      tmpstruct = struct('scan',scan,'perror',perror,'terror',terror,'mamp',mamp,...
         'mvel',mvel,'totmt',totmt,'mvtn',mvtn);
      switch task
         case 'stepeye'
            o = 1;
            p = 1;
         case 'step3'
            o = 2;
            p = 1;
         case 'step6'
            o = 3;
            p = 1;
         case 'step9'
            o = 4;
            p = 1;
         case 'step12'
            o = 5;
            p = 1;
         case 'traceye'
            o = 1;
            p = 2;
         case 'trac3'
            o = 2;
            p = 2;
         case 'trac6'
            o = 3;
            p = 2;
         case 'trac9'
            o = 4;
            p = 2;
         case 'trac12'
            o = 5;
            p = 2;
      end
      ascan(n,o,p) = scan;
      aperror(n,o,p) = perror;
      aterror(n,o,p) = terror;
      amamp(n,o,p) = mamp;
      amvel(n,o,p) = mvel;
      atotmt(n,o,p) = totmt;
      amvtn(n,o,p) = mvtn;
      aassym(n,o,p) = assym;
      at_frq(n,o,p) = t_freq;
      at_pow(n,o,p) = t_pow;
  end
   fclose(outfile);
end
cd( RTPATH );
save all_beh
