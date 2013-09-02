% Set up Transform menu
Transform_menu = uimenu(gcf,'Label','Transform','Position',5);

Filter_menu = uimenu(Transform_menu,'Label','Filter Traces');
Design_Filter_menu = uimenu(Filter_menu,'Label','Design Filter');
Fast_spikes_menu = uimenu(Filter_menu,'Label','Fast spikes - 300 to 5000Hz');
LFP_filter_menu = uimenu(Filter_menu,'Label','Low-frequency field potential (LFP) - 0 to 500Hz');

Invert_traces_menu = uimenu(Transform_menu,'Label','Invert Traces');
Resample_menu = uimenu(Transform_menu,'Label','Resample');
Subtrace_artifact_menu = uimenu(Transform_menu,'Label','Subtract artifact');
Subtrace_mean_menu = uimenu(Transform_menu,'Label','Subtract mean');

