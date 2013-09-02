
figure;
hold on
plot(FR_active, depth_active,...
    'MarkerFaceColor',[0 0 0],'Marker','o','LineStyle','none',...
    'Color',[0 0 0]);
% plot(norm_AP_width_spont_wide,depth_spont_wide,...
%     'MarkerFaceColor',[0 0 0],'Marker','square','LineStyle','none',...
%     'Color',[0 0 0]);
ha=gca;
% set(ha,'ylim',[-1 9]);
XLIM = get(ha,'xlim');
plot([XLIM(1) XLIM(2)],[30.5 30.5],'k--',...
    [XLIM(1) XLIM(2)],[36.5 36.5],'k--');

% hl = legend('narrow','wide');
xlabel('firing rate (Hz)');
ylabel('depth with respect to obex (mm)');
title('active firing rate vs depth');
hold off