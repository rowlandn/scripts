function [word, ts] = StripEvents(word, ts, BIT)
% function [word, ts] = StripStrobes(task, ts, bit)
%
% removes event codes & associated time stamps for 
% events that contain bit


strip_inds = find( bitget(word,BIT) );

word( strip_inds ) = [];
ts( strip_inds ) = [];

