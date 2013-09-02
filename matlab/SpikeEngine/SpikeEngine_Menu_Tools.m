% Set up Tools menu
Tools_menu = uimenu(gcf,'Label','Tools','Position',4);

Make_movie_menu = uimenu(Tools_menu,'Label','Make Movie');

Parse_log_file_menu = uimenu(Tools_menu,'Label','Parse log file');
Parse_PCDX_parameters = {'Enter Filename:'};
Parse_PCDX_dialogbox_title = 'Parse PCDX log file';
Parse_PCDX_log_file_menu = uimenu(Parse_log_file_menu,'Label','PCDX','CallBack',...
    ['[pcdx_log_filename,pcdx_log_pathname, pcdx_log_filterindex] = uigetfile(''*.log'', ''Select a .log file'');,'...
     'parse_pcdx_log([pcdx_log_pathname,pcdx_log_filename])']);
Parse_NeuroSage_log_file_menu = uimenu(Parse_log_file_menu,'Label','NeuroSage','CallBack',...
    ['[ns_data_filename,ns_data_pathname, ns_data_filterindex] = uigetfile(''*.data'', ''Select a .data file'');,'...
     'parse_NS_log_quick([ns_data_pathname,ns_data_filename])']);

Record_macro_menu = uimenu(Tools_menu,'Label','Record macro');


