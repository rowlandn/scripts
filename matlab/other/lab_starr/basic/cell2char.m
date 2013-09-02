function [STR,ii,jj] = cell2char(C,intercol,lines,nanstr) 
% CELL2CHAR  Convert cell to character matrix, without the limitations of CELL2MAT or CHAR.
%    STR = CELL2CHAR(C)         Convert cell to character matrix,
%
%    STR = CELL2CHAR(C,n)       adds inter-column space of width n (default is 0).
%
%    STR = CELL2CHAR(C,n,'-')   adds a horizontal line (of '-') under the column titles.
%    You can also use '_' and '=' for this.
%
%    STR = CELL2CHAR(C,n,'|')   adds a vertical line (of '|') after the row titles.
%    You can also use ':' for this.
%
%    STR = CELL2CHAR(C,n,'-|')  adds both horizontal and vertical lines
%
%    STR = CELL2CHAR(C,n,'',nanstr)  Replace NaNs by the string nanstr, (default is 'NaN')  
%
%    [STR,ii,jj] = CELL2CHAR(...)  returns vectors containing indexes in the char matrix corres-
%    ponding to rows (ii) and columns (jj) in the original cell.
%
% Related MATLAB functions:
% See also CHAR, CELL2MAT.

% ben 20 sep 2006
%     27 aug 2009 (Cosmetic)


if nargin < 2, intercol = 0; end
if nargin < 3, lines = '';   end
if nargin < 4, nanstr = ' '; end

f = [strfind(lines,'_') strfind(lines,'-') strfind(lines,'=')...
		strfind(lines,'–') strfind(lines,'—') strfind(lines,'­') strfind(lines,'­') strfind(lines,'¯')];
horline  = lines(f);
f = [strfind(lines,'|') strfind(lines,'¦') strfind(lines,':')];
vertline = lines(f);
	
maxlines = zeros(size(C,1),1);
maxcols  = zeros(1,size(C,2));

% 1) Convert all elements of the cell to char
for j = 1 : size(C,2)
    for i = 1 : size(C,1)
		if isnumeric(C{i,j}) || islogical(C{i,j})
			if isnan(C{i,j})
				C{i,j} = char(nanstr);   % replace NaN by nanstr (insure it's char type)
			else
				C{i,j} = num2str(C{i,j}); % numbers --> strings
			end
		elseif iscell(C{i,j})
			C{i,j} = cell2char(C{i,j},1,'',nanstr); % recursive call for each sub-cell
		end
		maxlines(i) = max([maxlines(i) size(C{i,j},1)]);
		maxcols(j)  = max([maxcols(j)  size(C{i,j},2)]);
	end
end

% 2) Uniformize their size
for i = 1 : size(C,1)
	for j = 1 : size(C,2)
		% Add lines
		n = size(C{i,j},1);
		m = size(C{i,j},2);
		C{i,j} = [C{i,j}; repmat(' ',maxlines(i)-n,m)];
		% Add cols
		n = size(C{i,j},1);
		C{i,j} = [C{i,j}  repmat(' ',n,maxcols(j)-m+intercol)];
	end
end

% 3) Add title sep. lines
if ~isempty(vertline)
	C = [C(:,1)  cell(size(C,1),1)  C(:,2:end)];
	for i = 1 : size(C,1)
		if intercol > 1
			C{i,1} = C{i,1}(:,1:end+1-intercol); % remove intercol. space for aesthetics
		end
		C{i,2} = repmat([vertline repmat(' ',1,intercol)],size(C{i,1},1),1);
    end
    if strcmp(vertline,':') && ~isempty(horline)
        C{1,2}(1) = ' ';  % (cosmetic)
    end
end
if ~isempty(horline)
	C = [C(1,:); cell(1,size(C,2)); C(2:end,:)];
	for j = 1 : size(C,2)
		C{2,j} = repmat(horline,1,size(C{1,j},2));
	end
end

% 4) Convert to char matrix
STR = cell2mat(C);

% ii: row indexes in the char matrix
tot = 0;
ii = zeros(size(C,1),1);
for i = 1 : length(ii)
	ii(i) = tot + 1;
	tot = tot + size(C{i,1},1);
end
if ~isempty(horline)
	ii(2) = [];
end

% jj: col. indexes in the char matrix
tot = 0;
jj = zeros(1,size(C,2));
for j = 1 : length(jj)
	jj(j) = tot + 1;
	tot = tot + size(C{1,j},2);
end
if ~isempty(vertline)
	jj(2) = [];
end