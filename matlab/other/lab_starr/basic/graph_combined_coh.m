%% Combining data and comparing across disease states
%REST
% figure;
% subplot(2,1,1)
%     plot(adam_trans_coh.freq, adam_trans_coh.rest(:,5), 'DisplayName', 'Adam','color','b') 
%     hold on;
% %     plot(andr_trans_coh.freq,andr_trans_coh.rest_M1,'DisplayName', 'Andr','color', 'b')
% %     hold on;
% %     plot(hase_trans_coh.freq,hase_trans_coh.rest_M1,'DisplayName', 'Hase','color', 'b')
% %     hold on;
% %     plot(lave_trans_coh.freq,lave_trans_coh.rest_M1,'DisplayName', 'Lave','color', 'b')
% %     hold on;
%     plot(espi_trans_coh.freq, espi_trans_coh.rest(:,5),'DisplayName', 'Espi','color','b')
%     hold on;
% %     plot(beat_trans_coh.freq,beat_trans_coh.rest_M1,'DisplayName', 'Beat','color', 'r')
% %     hold on;
%     plot(paol_trans_coh.freq,paol_trans_coh.rest_M1,'DisplayName', 'Paol','color', 'r')
%     hold on;
%     plot(mill_trans_coh.freq,mill_trans_coh.rest(:,5), 'DisplayName', 'Mill','color', 'r')
%     axis([0 150 0 1])
%     set (gca,'xtick',[0:20:150]);
%     set (gca,'ytick',[0:0.2:1]);
%     
% subplot(2,1,2)
%     plot(adam_trans_coh.freq, adam_trans_coh.rest(:,5), 'DisplayName', 'Adam','color','b') 
%     hold on;
% %     plot(andr_trans_coh.freq,andr_trans_coh.rest_M1,'DisplayName', 'Andr','color', 'b')
% %     hold on;
% %     plot(hase_trans_coh.freq,hase_trans_coh.rest_M1,'DisplayName', 'Hase','color', 'b')
% %     hold on;
% %     plot(lave_trans_coh.freq,lave_trans_coh.rest_M1,'DisplayName', 'Lave','color', 'b')
% %     hold on;
%     plot(espi_trans_coh.freq, espi_trans_coh.rest(:,5),'DisplayName', 'Espi','color','b')
%     hold on;
% %     plot(beat_trans_coh.freq,beat_trans_coh.rest_M1,'DisplayName', 'Beat','color', 'r')
% %     hold on;
%     plot(paol_trans_coh.freq,paol_trans_coh.rest_M1,'DisplayName', 'Paol','color', 'r')
%     hold on;
%     plot(mill_trans_coh.freq,mill_trans_coh.rest(:,5), 'DisplayName', 'Mill','color', 'r')
%     axis([0 50 0 1])
%     set (gca,'xtick',[0:10:50]);
%     set (gca,'ytick',[0:0.2:1]);

% ACTIVE
figure;
subplot(2,1,1)
    plot(adam_trans_coh.freq, adam_trans_coh.active(:,5), 'DisplayName', 'Adam','color','b') 
    hold on;
%     plot(andr_trans_coh.freq,andr_trans_coh.rest_M1,'DisplayName', 'Andr','color', 'b')
%     hold on;
%     plot(hase_trans_coh.freq,hase_trans_coh.rest_M1,'DisplayName', 'Hase','color', 'b')
%     hold on;
%     plot(lave_trans_coh.freq,lave_trans_coh.rest_M1,'DisplayName', 'Lave','color', 'b')
%     hold on;
    plot(espi_trans_coh.freq, espi_trans_coh.active(:,5),'DisplayName', 'Espi','color','b')
    hold on;
%     plot(beat_trans_coh.freq,beat_trans_coh.rest_M1,'DisplayName', 'Beat','color', 'r')
%     hold on;
    plot(paol_trans_coh.freq,paol_trans_coh.active_M1,'DisplayName', 'Paol','color', 'r')
    hold on;
    plot(mill_trans_coh.freq,mill_trans_coh.rest(:,5), 'DisplayName', 'Mill','color', 'r')
    axis([0 150 0 1])
    set (gca,'xtick',[0:20:150]);
    set (gca,'ytick',[0:0.2:1]);
    
subplot(2,1,2)
    plot(adam_trans_coh.freq, adam_trans_coh.active(:,5), 'DisplayName', 'Adam','color','b') 
    hold on;
%     plot(andr_trans_coh.freq,andr_trans_coh.rest_M1,'DisplayName', 'Andr','color', 'b')
%     hold on;
%     plot(hase_trans_coh.freq,hase_trans_coh.rest_M1,'DisplayName', 'Hase','color', 'b')
%     hold on;
%     plot(lave_trans_coh.freq,lave_trans_coh.rest_M1,'DisplayName', 'Lave','color', 'b')
%     hold on;
    plot(espi_trans_coh.freq, espi_trans_coh.active(:,5),'DisplayName', 'Espi','color','b')
    hold on;
%     plot(beat_trans_coh.freq,beat_trans_coh.rest_M1,'DisplayName', 'Beat','color', 'r')
%     hold on;
    plot(paol_trans_coh.freq,paol_trans_coh.active_M1,'DisplayName', 'Paol','color', 'r')
    hold on;
    plot(mill_trans_coh.freq,mill_trans_coh.active(:,5), 'DisplayName', 'Mill','color', 'r')
    axis([0 50 0 1])
    set (gca,'xtick',[0:10:50]);
    set (gca,'ytick',[0:0.2:1]);
    