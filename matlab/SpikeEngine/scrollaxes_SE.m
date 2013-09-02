% SCROLLAXES displays five buttons on the bottom left corner of the current 
% figure. The buttons pan all the plots in the current axes up, down and left
% and right by 10% of the axis range increments. The button at the center of 
% the 'arrow heads' resets the axes limits.
% 
% You can also use the arrow keys to pan with finer increments.
% 
% Use '<' and '>' to zoom out and in; use 'x' and 'X' and 'y' and 'Y'to zoom 
% out and in only the x or y axis respectively. Use 'a' to autoscale and 'g' 
% to toggle the axis grid.
%
% Usage:
%   
%   SCROLLAXES, by itself, toggles the display of the buttons
%   
%   SCROLLAXES([STEPX STEPY]) turns the buttons on or updates the steps on 
%                             x and y to be STEPX and STEPY, respectively.
%
% <adelgado@biology.emory.edu>

function scrollaxes(deltas)
load /usr/local/matlab6p5/toolbox/signal/sptoolgui/private/siggui_icons
    scrollaxes_slider = uicontrol(gcf,'Style','slider','Position',[5 100 85 20],'CallBack',[]);
        %['scrollaxes_value = get(scrollaxes_slider,''Value''); pan_rate = 1 + ((scrollaxes_value/10)*-1);,']);
         %'assignin(''base'',''pan_rate'',pan_rate);']);
    assignin('base','scrollaxes_slider',scrollaxes_slider)
    scrollaxes_value = get(scrollaxes_slider,'Value');
    assignin('base','scrollaxes_value',scrollaxes_value)
    pan_rate = 1 + ((scrollaxes_value/10)*-1);
    assignin('base','pan_rate',pan_rate)
    
    if isempty(findobj(gcf, 'Tag', 'AlfScrollerCenter'))
        
        set(gcf, 'DoubleBuffer', 'on')
        set(gcf, 'KeyPressFcn', 'scrollaxes AlfKeyPress')
        
        if nargin == 0, 
            xLim = get(gca, 'XLim');
            yLim = get(gca, 'YLim');
            deltas = [abs(xLim(2)-xLim(1))/.09 abs(yLim(2)-yLim(1))/.09];
        end
        
        % -- Scroll to the Right
        1
        SCRLRB = uicontrol('Style', 'pushbutton', 'Position', [60 25 30 30],'cdata',bmp.mousezoom);
        assignin('base','SCRLRB',SCRLRB)
        set(SCRLRB, 'String', '>', 'FontWeight', 'bold', 'FontSize', 16)
        set(SCRLRB, 'Tag', 'AlfScrollerRight', 'UserData', deltas)
        
        %while get(SCRLRB,'Value') == 1
        cb = ['dx = get(findobj(gcf, ''Tag'', ''AlfScrollerRight''), ' ...
              '''UserData''); xlim(xlim + dx(1))'];
        set(SCRLRB, 'CallBack', ['SCRLRB_value = 1; while SCRLRB_value == 1, cb, end;'])
        %set(SCRLRB, 'CallBack', ['SCRLRB_value = 1; SCRLRB_value; cb'])
        %end
        
        % -- Scroll to the Left
        SCRLLB = uicontrol('Style', 'pushbutton', 'Position', [4 25 30 30]);
        set(SCRLLB, 'String', '<', 'FontWeight', 'bold', 'FontSize', 16)
        set(SCRLLB, 'Tag', 'AlfScrollerLeft', 'UserData', deltas)
        cb = ['dx = get(findobj(gcf, ''Tag'', ''AlfScrollerLeft''), ' ...
              '''UserData''); xlim(xlim - dx(1))'];
        set(SCRLLB, 'CallBack', cb)
        
        % -- Scroll Upwards
        SCRLUB = uicontrol('Style', 'pushbutton', 'Position', [32 54 30 30]);
        set(SCRLUB, 'String', '^', 'FontWeight', 'bold', 'FontSize', 16)
        set(SCRLUB, 'Tag', 'AlfScrollerUp', 'UserData', deltas)
        cb = ['dy = get(findobj(gcf, ''Tag'', ''AlfScrollerUp''), ' ...
              '''UserData''); ylim(ylim + dy(2))'];
        set(SCRLUB, 'CallBack', cb)
        
        % -- Scroll Downwards
        SCRLDB = uicontrol('Style', 'pushbutton', 'Position', [32 2 30 30]);
        set(SCRLDB, 'String', 'v', 'FontWeight', 'bold', 'FontSize', 16)
        set(SCRLDB, 'Tag', 'AlfScrollerDown', 'UserData', deltas)
        cb = ['dy = get(findobj(gcf, ''Tag'', ''AlfScrollerDown''), ' ...
              '''UserData''); ylim(ylim - dy(2))'];
        set(SCRLDB, 'CallBack', cb)
        
        % -- Re-center
        SCRLCB = uicontrol('Style', 'pushbutton', 'Position', [32 25 30 30]);
        set(SCRLCB, 'String', 'o', 'FontWeight', 'normal', 'FontSize', 16)
        set(SCRLCB, 'Tag', 'AlfScrollerCenter', 'CallBack', 'axis auto')
        
    elseif (nargin == 1) & strcmp('double', class(deltas))
        set(findobj(gcf, 'Tag', 'AlfScrollerRight'), 'UserData', deltas)
        set(findobj(gcf, 'Tag', 'AlfScrollerLeft'), 'UserData', deltas)
        set(findobj(gcf, 'Tag', 'AlfScrollerUp'), 'UserData', deltas)
        set(findobj(gcf, 'Tag', 'AlfScrollerDown'), 'UserData', deltas)
        
    elseif (nargin == 1) & strcmp('char', class(deltas))
        panFactor  = 0.02;
        zoomFactor = 0.90;
        
        pressedKey = get(gcbf, 'CurrentCharacter');
        
        if isempty(pressedKey), return, end
        
        switch pressedKey
            case 'a',
                axis auto
            case 'g',
                grid
            case 28, % right arrow
                xLim = get(gca, 'XLim');
                xNewLim = xLim + panFactor * diff(xLim);
                set(gca, 'XLim', xNewLim)
            case 29, % left arrow
                xLim = get(gca, 'XLim');
                xNewLim = xLim - panFactor * diff(xLim);
                set(gca, 'XLim', xNewLim)
            case 30, % up arrow
                yLim = get(gca, 'YLim');
                yLimNew = yLim - panFactor * diff(yLim);
                set(gca, 'YLim', yLimNew)
            case 31, % down arrow
                yLim = get(gca,'YLim');
                yLimNew = yLim + panFactor * diff(yLim);
                set(gca, 'YLim', yLimNew)
            case {'x', 'X'},
                if pressedKey == 'X',
                    zoomFactor = 1/zoomFactor;
                end
                xLim = get(gca, 'XLim');
                xNewLim = [0 zoomFactor*diff(xLim)] + xLim(1) + ...
                          (1 - zoomFactor) * diff(xLim)/2;
                set(gca, 'XLim', xNewLim)
            case {'y', 'Y'},
                if pressedKey == 'Y',
                    zoomFactor = 1/zoomFactor;
                end
                yLim = get(gca, 'YLim');
                yLimNew = [0 zoomFactor*diff(yLim)] + yLim(1) + ...
                          (1 - zoomFactor) * diff(yLim)/2;
                set(gca, 'YLim', yLimNew)
            case {'<', '>'},
                if pressedKey == '<',
                    zoomFactor = 1/zoomFactor;
                end
                xLim = get(gca, 'XLim');
                yLim = get(gca, 'YLim');
                    
                xNewLim = [0 zoomFactor*diff(xLim)] + xLim(1) + ...
                          (1 - zoomFactor) * diff(xLim)/2;
                yLimNew = [0 zoomFactor*diff(yLim)] + yLim(1) + ...
                          (1 - zoomFactor) * diff(yLim)/2;
                
                set(gca, 'XLim', xNewLim, 'YLim', yLimNew)
        end
    else
        delete(findobj(gcf, 'Tag', 'AlfScrollerRight'))
        delete(findobj(gcf, 'Tag', 'AlfScrollerLeft'))
        delete(findobj(gcf, 'Tag', 'AlfScrollerUp'))
        delete(findobj(gcf, 'Tag', 'AlfScrollerDown'))
        delete(findobj(gcf, 'Tag', 'AlfScrollerCenter'))
        
        set(gcf, 'KeyPressFcn', '')
      % set(gcf, 'DoubleBuffer', 'off') % won't hurt...
    end
    