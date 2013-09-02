function PeriEvent_Ivry(show)
% PeriEvent_Ivry(show)
% This program constructs peri-event SDF around specific events
% onset
% Created by RST, 2005-12-05
%
%	Input:
%		show - controls whether graphical output is produced for each cell
%		  default = true.
%
%	Run within a directory and this program will process all spikes in all 
%	nex files in the directory
%	


% Global defines used by the other functions called by this one
global PRE_EVT
global PST_EVT
global BIN_SZ
global CNTL_PER
global MAX_N
global ALPHA
global CONTIG

global SRCH_LO
global SRCH_HI
global VERBOSE

PRE_EVT = -1.5;	% Duration of SDF pre-trigger
PST_EVT = 1;	% Duration of SDF post-trigger
SMOOTH = 50;	% SD of smoothing gaussian for spike train (in msec)
CNTL_PER = [-1.5 -0.5];
ALPHA = 0.05;
CONTIG = 50;	% Mean of 'CONTIG' adjacent points must be significant
JOIN = 0.2;		% Join together responses that are < JOIN sec apart
i=1;
j=1;
MIN_SPK_N = 500;

VERBOSE=false;

N_Pks = 2;	% max # of significant changes

if ~exist('show','var')
	show = true;
end

% Edit this pattern to select time-stamp variables w/ specific names
%spkname_pattern = 'Snip\w*[abcd]';	% Kevin's pattern
spkname_pattern = '\w*';	% Accept all units
evtname_pattern = 'Accel\w*';	% Events must start w/ "Accel"


cd(uigetdir);

FileLst = dir('*.nex');
if(isempty(FileLst))
	str = pwd;
	error(['Found no NEX files in current directory - ' str ]);
end

%outfid = write_text_header(N_Pks);	% Subfunction below

n = 0;	% Count of units processed
% For each file found in directory...
for i=1:length(FileLst)

	matname = strrep(FileLst(i).name,'.nex','.mat');
	load(matname);
	nexname = FileLst(i).name;
	% Find spike variables in nexfile - allow for more than 1 spike in file
	[spk_inds,spk_names] = find_nex_units( nexname, spkname_pattern);
	if isempty(spk_inds)
		continue
	end

	%  Evt = Detect movement events of interest;
	if isempty(Evt) 
		display(['Found no event time stamps in file: ' fname '. Skipping.']);
		continue
	elseif length(Evt)==1
		display(['Found only 1 evt in file: ' fname '. Continuing']);
	end
	for k=1:length(Evt)		% Filter out events too close to data boundaries
		drp = find(Evt(k).ts < data_start-PRE_MVT | Evt(k).ts > data_stop-PST_MVT);
		Mvt(k).ts(drp) = [];
	end	
	
	% For each unit in a file...
	for j=1:length(spk_inds)
		n = n+1;
		spk{n}.fname = fname;

		% save name of file & unit
		spk{n}.unitname = deblank( spk_names(j,:) );
		display(['Processing...' spk{n}.fname '....:' spk{n}.unitname]);
		
		% get spike times
		[spk_n, spk_t] = nex_ts( nexname, spk_names(j,:),VERBOSE);
		if spk_n < MIN_SPK_N
			display(['Found only ' num2str(spk_n) ' spikes. Skipping.']);
			continue
		end

		% for each direction of movement
		for k=1:length(Mvt)
			spk{n}.dir(k).n_reps = Mvt(k).n;

			[spk{n}.dir(k).histog, spk{n}.bins] = perievent_sdf(spk_t, Mvt(k).ts, ...
				PRE_MVT, PST_MVT, SMOOTH );

			spk{n}.dir(k).raster = perievent_raster(spk_t, Mvt(k).ts, ...
				PRE_MVT, PST_MVT);

			cntl_inds = find(spk{n}.bins>=CNTL_PER(1) & spk{n}.bins<CNTL_PER(2));
			test_start = max(cntl_inds)+1;
			
			[spk{n}.dir(k).chng(1),spk{n}.dir(k).cntl_mean,spk{n}.dir(k).sig_thr] = ...
				PeriEventChange_SDF(spk{n}.dir(k).histog,cntl_inds,test_start,ALPHA,CONTIG);
			
			sgn1 = spk{n}.dir(k).chng(1).sgn;
			sgn2 = [];
			del2nd = 0;
			if ~isempty( spk{n}.dir(k).chng(1).off_ind )
				test_start = spk{n}.dir(k).chng(1).off_ind;
 				[spk{n}.dir(k).chng(2),spk{n}.dir(k).cntl_mean,q] = ...
 					PeriEventChange_SDF(spk{n}.dir(k).histog,cntl_inds,test_start,ALPHA,CONTIG);
				sgn2 = spk{n}.dir(k).chng(2).sgn;
			end
			if sgn2==sgn1
				off1_t = spk{n}.bins( spk{n}.dir(k).chng(1).off_ind ) ;
				on2_t = spk{n}.bins( spk{n}.dir(k).chng(2).on_ind ) ;
				if (on2_t - off1_t)< JOIN
					on1_ind = spk{n}.dir(k).chng(1).on_ind ;
					off1_ind = spk{n}.dir(k).chng(2).off_ind ;
					spk{n}.dir(k).chng(1).off_ind = off1_ind;
					his = spk{n}.dir(k).histog - spk{n}.dir(k).cntl_mean;
					if ~isempty(off1_ind)
						spk{n}.dir(k).chng(1).mean_change = ...
							mean( his(on1_ind:off1_ind) );
						spk{n}.dir(k).chng(1).int_change = ...
							sum( his(on1_ind:off1_ind) )/1000;
					else
						spk{n}.dir(k).chng(1).mean_change = ...
							mean( his(on1_ind:end) );
						spk{n}.dir(k).chng(1).int_change = ...
							sum( his(on1_ind:end) )/1000;
					end
					del2nd = 1;
				end
			end
			if isempty( spk{n}.dir(k).chng(1).off_ind ) | del2nd
				spk{n}.dir(k).chng(2).on_ind = [];
				spk{n}.dir(k).chng(2).sgn = [];
				spk{n}.dir(k).chng(2).off_ind = [];
				spk{n}.dir(k).chng(2).mean_change = [];
				spk{n}.dir(k).chng(2).int_change = [];
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
fclose(outfid);
save PeriPassive spk


%------------------------------------------------------
% Subfunction to write a line of data to output file
function write_text(outfid, spk,N_Pks)
	bins = spk.bins;
	
	fprintf(outfid,'%s\t%s\t', ...
		spk.fname, spk.unitname);

	for j = 1:length(spk.dir)	% mvt directions
		mvt = spk.dir(j);
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
	fname = 'PeriPassive.txt';
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
	fprintf(outfid,'\n');		% EOL

	return

	

%------------------------------------------------------
% Subfunction to make figure of results
function make_figure(s, N_Pks)
	bins = s.bins;
	
	%%%%%%%%%%%%%% Plotting
	% Set up axes
	MARGIN = 0.06;	
	TOP = 1-MARGIN;		% Top margin of page
	LEFT = MARGIN;	% Left margin of page
	WIDTH = (1-3*MARGIN)/2;	% give space for 3 margin widths including middle
	HEIGHT = 0.4;	% Height of histograms
	CLR = [0.25,0.25,0.25 ; 0.75,0.75,0.75];

	figure
	set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
	% Size to make it look good
	c = get(gcf);
	c.Position(2) = 275;
	c.Position(3) = 870;
	c.Position(4) = 680;
	set(gcf,'Position',c.Position);

	% Find max across directions
	for j = 1:length(s.dir)
		x(j) = max(s.dir(j).histog);
		ymax = max(x)+5;
	end

	% Plot for 2 movements
	for j=1:length(s.dir)
		mvt = s.dir(j);
		left = MARGIN + (WIDTH+MARGIN)*(j-1);
		width = WIDTH;
		height = HEIGHT;      
		bottom = TOP-height;
		subplot('position',[left bottom width height]);

		h=area(bins, mvt.histog);
		set(h,'FaceColor',[0.5,0.5,0.5],'EdgeColor','k');
		xlim([min(bins) max(bins)]);
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
						off = xlm(2);	%If no offset found, change lasts to end
					end
					x = [ on on off off ];
					y = [ ylm ylm(2:-1:1)];
					fill(x,y,CLR(i,:),'EdgeColor','none','FaceAlpha',0.3)
				end
			end
		end

		
		if j ==1
			title([ s.fname ':    ' s.unitname],'Interpreter','none','FontSize',14,...
				'Position',[xlm(2),ylm(2)+5,0]);
		end

		bottom = bottom-height-MARGIN;
		subplot('position',[left bottom width height]);
		rasterplot(mvt.raster,0.9,'');
	end

return


