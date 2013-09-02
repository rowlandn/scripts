function axis_2(action)

if nargin == 0
    %varargin = 'init'
    %varargin = 'init';
    subplot(12,10,[111:117])
    set(gca,'XTick',[],'YTick',[],'Box','off','YColor',[.8 .8 .8],'XColor',[.8 .8 .8])
% elseif nargin == 1
%      action = 'init2'
else

switch action
    case 'Temp_filename_raw'
        a = 1
        subplot(12,10,[111:117])
        plot(Temp_filename_raw(:,end),Temp_filename_raw(:,1:end-1),'k')
        subplot(10,10,[1:7,11:17,21:27,31:37,41:47,51:57,61:67,71:77,81:87])
    case '1'           
end
end