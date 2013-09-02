figure
set(gcf,'Position',[150 40 1000 700],'MenuBar','none','Name','NeuroSage ''.data'' file parser');
Parse_NS_box = uicontrol(gcf,'Style','edit','Position',[100 100 800 500],...
    'String',parse_string,'Max',10,'Min',1,'HorizontalAlignment','left');

% ,'CallBack',...
%     ['LoadTraces_Type_userinpt = {''PCDX'',''GENESIS'',''NeuroSage'',''Text''};,'...



