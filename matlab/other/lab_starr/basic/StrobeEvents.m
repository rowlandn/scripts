function [word, ts] = StrobeEvents(word, ts, BIT)
% function [word, ts] = StrobeEvents(task, ts, bit)
%
% Returns event codes & associated time stamps  
% that contain BIT
% Remove strobe BIT as well


inds = find( ~bitget(word,BIT) );

word( strip_inds ) = [];
ts( strip_inds ) = [];

