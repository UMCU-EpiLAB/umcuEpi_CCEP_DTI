function H = plot_ccep(dataBase,subj,run,tt,chan,stimp)
% figure with left the epoch, and right zoomed in
H=figure(1);
H.Units = 'normalized';
H.Position = [0.01 0.35 0.77 0.5];

usedSignal = squeeze(dataBase(subj).metadata(run).cc_epoch_sorted_reref(chan,:,stimp,:));

low_ci = quantile(usedSignal,.16,1); % find "mean - SD"
high_ci = quantile(usedSignal,.84,1); % find "mean + SD"

tt0 = find(tt>= -0.1,1,'first');
tt1 = find(tt>= 0.5,1,'first');

% /// SUBPLOT 1
subplot(1,2,1)
fill([tt(tt0:tt1) tt(tt1:-1:tt0)],[low_ci(tt0:tt1) high_ci(tt1:-1:tt0)],'r','LineStyle','none','FaceAlpha',0.3);
hold on
plot(tt,usedSignal,':'); % individual responses
plot(tt,low_ci,'-','Color',[0.5 0.5 0.5],'LineWidth',1) % plot confidence interval
plot(tt,high_ci,'-','Color',[0.5 0.5 0.5],'LineWidth',1)

% plot the part with stimulus artefact
fill([tt(tt>=0 & tt<0.01) flip(tt(tt>=0 & tt<0.01))],...
    [-2000*ones(1,size(tt(tt>=0 & tt<0.01),2)) 2000*ones(1,size(tt(tt>=0 & tt<0.01),2))],'k','LineStyle','none','FaceAlpha',0.3)

% plot the window in which the n1-CCEP should be scored
plot(tt(find(tt<0.1,1,'last'))*ones(2,1),[-2000 2000],'k-');

% plot average signal and re-ref signal
plot(tt,squeeze(dataBase(subj).metadata(run).cc_epoch_sorted_avg(chan,stimp,:)),'k-.','linewidth',1);
plot(tt,squeeze(dataBase(subj).metadata(run).cc_epoch_sorted_reref_avg(chan,stimp,:)),'k','linewidth',2);

set(gca,'Fontsize',14)

hold off

xlim([-2 2])
ylim([-2000 2000])
xlabel('Time (s)')
ylabel('Amplitude (uV)')
title(sprintf('Electrode %s, stimulating %s-%s',...
    dataBase(subj).metadata(run).ch{chan},...
    dataBase(subj).metadata(run).cc_stimchans{stimp,1},...
    dataBase(subj).metadata(run).cc_stimchans{stimp,2}))

% /// SUBPLOT 2
subplot(1,2,2)
h1 = fill([tt(tt0:tt1) tt(tt1:-1:tt0)],[low_ci(tt0:tt1) high_ci(tt1:-1:tt0)],'r','LineStyle','none','FaceAlpha',0.2);
hold on
h2 = plot(tt,usedSignal,':','LineWidth',1.5); % individual responses/

plot(tt,low_ci,'-','Color',[0.5 0.5 0.5],'LineWidth',1) % plot confidence interval
plot(tt,high_ci,'-','Color',[0.5 0.5 0.5],'LineWidth',1)

% plot the part with stimulus artefact
fill([tt(tt>=0 & tt<0.01) flip(tt(tt>=0 & tt<0.01))],...
    [-750*ones(1,size(tt(tt>=0 & tt<0.01),2)) 750*ones(1,size(tt(tt>=0 & tt<0.01),2))],'k','LineStyle','none','FaceAlpha',0.3)

% plot the window in which the n1-CCEP should be scored
h6 = plot(tt(find(tt<0.1,1,'last'))*ones(2,1),[-750 750],'k-');

% plot average signal and re-ref signal
h3 = plot(tt,squeeze(dataBase(subj).metadata(run).cc_epoch_sorted_avg(chan,stimp,:)),'k-.','linewidth',1);
h4 = plot(tt,squeeze(dataBase(subj).metadata(run).cc_epoch_sorted_reref_avg(chan,stimp,:)),'k','linewidth',2);

% plot detected N1 (if detected)
h5 = plot(tt(dataBase(subj).metadata(run).ccep.n1_peak_sample(chan,stimp)),dataBase(subj).metadata(run).ccep.n1_peak_amplitude(chan,stimp),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4);

hold off

% show legend
legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'CI','indiv responses','average','average reref','detected ER','time window ERP')
set(gca,'Fontsize',14)
xlim([-0.2 0.4])
ylim([-750 750])
xlabel('Time (s)')
ylabel('Amplitude (uV)')
title('Zoomed average signal')
end