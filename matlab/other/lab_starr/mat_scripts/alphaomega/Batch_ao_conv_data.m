
files = dir('*mat');

for i = 1:length(files)
    filename = files(i).name;
    if isempty(strfind(filename,'ipad.mat'))
    aomatconv_ipad(filename,4)
    end
end

%%VA
files = dir('*mat');

for i = 1:length(files)
    filename = files(i).name;
    if isempty(strfind(filename,'ipad.mat'))
    aomatconv_ipad_VA(filename,5)
    end
end