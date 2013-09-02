function C = cell2num(C)
% M = cell2num(C)


% str --> num
for i = find(cellfun('isclass',C,'char'))
    C{i} = str2num(C{i});
    if isempty(C{i})
        C{i} = NaN;
    end
end

% cell --> matrix
M = cell2mat(C);