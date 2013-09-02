function plot_glr_rd_cust_01_avg(glr_rd,gauss_width,sp1,sp2,sp3)

if gauss_width == 1
    subplot(sp1,sp2,sp3)
    max_gauss_rd = max(glr_rd.gauss_traces.one_ms(:,1));
    plot(glr_rd.gauss_traces.one_ms(:,2),glr_rd.gauss_traces.one_ms(:,1),'.','MarkerSize',5); hold on
    plot([glr_rd.periods.reference(1) glr_rd.periods.reference(2)],[glr_rd.std.one_ms.plus glr_rd.std.one_ms.plus],'r:')
    plot([glr_rd.periods.response(1) glr_rd.periods.response(2)],[glr_rd.std.one_ms.plus glr_rd.std.one_ms.plus],'r')
    plot([glr_rd.periods.reference(1) glr_rd.periods.reference(2)],[glr_rd.std.one_ms.minus glr_rd.std.one_ms.minus],'r:')
    plot([glr_rd.periods.response(1) glr_rd.periods.response(2)],[glr_rd.std.one_ms.minus glr_rd.std.one_ms.minus],'r')


    if max_gauss_rd+.05*max_gauss_rd > glr_rd.std.one_ms.plus
        y_max = max_gauss_rd+.05*max_gauss_rd;
    elseif glr_rd.std.one_ms.plus > max_gauss_rd+.05*max_gauss_rd 
        y_max = glr_rd.std.one_ms.plus+.1*glr_rd.std.one_ms.plus;
    end

    if glr_rd.sig_responses.sle(2) == 0
        text(-1,y_max+y_max*.05,'X','FontSize',6,'Color',[1 0 1])
    else
        plot([glr_rd.sig_responses.sle(1,10) glr_rd.sig_responses.sle(1,11)], [y_max+y_max*.05 y_max+y_max*.05],'m','LineWidth',2)
    end

    if glr_rd.sig_responses.inh(2) == 0
        text(-1,y_max+y_max*.15,'X','FontSize',6,'Color',[0 1 0])
    else
        plot([glr_rd.sig_responses.inh(1,10) glr_rd.sig_responses.inh(1,11)], [y_max+y_max*.15 y_max+y_max*.15],'g','LineWidth',2)
    end

    if glr_rd.sig_responses.lle(2) == 0
        text(-1,y_max+y_max*.25,'X','FontSize',6,'Color',[0 1 1])
    else
        plot([glr_rd.sig_responses.lle(1,10) glr_rd.sig_responses.lle(1,11)], [y_max+y_max*.25 y_max+y_max*.25],'c','LineWidth',2)
    end
    set(gca,'Box','off','xlim',[0 size(glr_rd.gauss_traces.one_ms,1)/10],'ylim',[0 max_gauss_rd + .5*max_gauss_rd])
elseif gauss_width == 5
    subplot(sp1,sp2,sp3)
    max_gauss_rd = max(glr_rd.gauss_traces.five_ms(:,1));
    plot(glr_rd.gauss_traces.five_ms(:,2),glr_rd.gauss_traces.five_ms(:,1),'.','MarkerSize',1); hold on
    plot([glr_rd.periods.reference(1) glr_rd.periods.reference(2)],[glr_rd.std.five_ms.plus glr_rd.std.five_ms.plus],'r:')
    plot([glr_rd.periods.response(1) glr_rd.periods.response(2)],[glr_rd.std.five_ms.plus glr_rd.std.five_ms.plus],'r')
    plot([glr_rd.periods.reference(1) glr_rd.periods.reference(2)],[glr_rd.std.five_ms.minus glr_rd.std.five_ms.minus],'r:')
    plot([glr_rd.periods.response(1) glr_rd.periods.response(2)],[glr_rd.std.five_ms.minus glr_rd.std.five_ms.minus],'r')


    if max_gauss_rd+.05*max_gauss_rd > glr_rd.std.five_ms.plus
        y_max = max_gauss_rd+.05*max_gauss_rd;
    elseif glr_rd.std.five_ms.plus > max_gauss_rd+.05*max_gauss_rd 
        y_max = glr_rd.std.five_ms.plus+.1*glr_rd.std.five_ms.plus;
    end

    if glr_rd.sig_responses.sle(2) == 0
        text(-1,y_max+y_max*.05,'X','FontSize',6,'Color',[1 0 1])
    else
        plot([glr_rd.sig_responses.sle(1,10) glr_rd.sig_responses.sle(1,11)], [y_max+y_max*.05 y_max+y_max*.05],'m','LineWidth',2)
    end

    if glr_rd.sig_responses.inh(2) == 0
        text(-1,y_max+y_max*.15,'X','FontSize',6,'Color',[0 1 0])
    else
        plot([glr_rd.sig_responses.inh(1,10) glr_rd.sig_responses.inh(1,11)], [y_max+y_max*.15 y_max+y_max*.15],'g','LineWidth',2)
    end

    if glr_rd.sig_responses.lle(2) == 0
        text(-1,y_max+y_max*.25,'X','FontSize',6,'Color',[0 1 1])
    else
        plot([glr_rd.sig_responses.lle(1,10) glr_rd.sig_responses.lle(1,11)], [y_max+y_max*.25 y_max+y_max*.25],'c','LineWidth',2)
    end
    set(gca,'Box','off','xlim',[0 size(glr_rd.gauss_traces.five_ms,1)/10],'ylim',[0 max_gauss_rd + .5*max_gauss_rd])
elseif gauss_width == 20
    subplot(sp1,sp2,sp3)
    max_gauss_rd = max(glr_rd.gauss_traces.twenty_ms(:,1));
    plot(glr_rd.gauss_traces.twenty_ms(:,2),glr_rd.gauss_traces.twenty_ms(:,1),'.','MarkerSize',1); hold on
    plot([glr_rd.periods.reference(1) glr_rd.periods.reference(2)],[glr_rd.std.twenty_ms.plus glr_rd.std.twenty_ms.plus],'r:')
    plot([glr_rd.periods.response(1) glr_rd.periods.response(2)],[glr_rd.std.twenty_ms.plus glr_rd.std.twenty_ms.plus],'r')
    plot([glr_rd.periods.reference(1) glr_rd.periods.reference(2)],[glr_rd.std.twenty_ms.minus glr_rd.std.twenty_ms.minus],'r:')
    plot([glr_rd.periods.response(1) glr_rd.periods.response(2)],[glr_rd.std.twenty_ms.minus glr_rd.std.twenty_ms.minus],'r')


    if max_gauss_rd+.05*max_gauss_rd > glr_rd.std.twenty_ms.plus
        y_max = max_gauss_rd+.05*max_gauss_rd;
    elseif glr_rd.std.twenty_ms.plus > max_gauss_rd+.05*max_gauss_rd 
        y_max = glr_rd.std.twenty_ms.plus+.1*glr_rd.std.twenty_ms.plus;
    end

    if glr_rd.sig_responses.sle(2) == 0
        text(-1,y_max+y_max*.05,'X','FontSize',6,'Color',[1 0 1])
    else
        plot([glr_rd.sig_responses.sle(1,10) glr_rd.sig_responses.sle(1,11)], [y_max+y_max*.05 y_max+y_max*.05],'m','LineWidth',2)
    end

    if glr_rd.sig_responses.inh(2) == 0
        text(-1,y_max+y_max*.15,'X','FontSize',6,'Color',[0 1 0])
    else
        plot([glr_rd.sig_responses.inh(1,10) glr_rd.sig_responses.inh(1,11)], [y_max+y_max*.15 y_max+y_max*.15],'g','LineWidth',2)
    end

    if glr_rd.sig_responses.lle(2) == 0
        text(-1,y_max+y_max*.25,'X','FontSize',6,'Color',[0 1 1])
    else
        plot([glr_rd.sig_responses.lle(1,10) glr_rd.sig_responses.lle(1,11)], [y_max+y_max*.25 y_max+y_max*.25],'c','LineWidth',2)
    end
    set(gca,'Box','off','xlim',[0 size(glr_rd.gauss_traces.twenty_ms,1)/10],'ylim',[0 max_gauss_rd + .5*max_gauss_rd])
end