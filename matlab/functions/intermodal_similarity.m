function [J,R] = intermodal_similarity(JI, JI_expected, ratio, modality)
% INPUT:
% cfg.sub_label = {'RESPXXXX', 'RESPXXXX'}
% vector_v1 = topography vector modality 1
% vector_v2 = topography vector modality 2
% v1_soz = topography vector soz channels
% v1_nsoz = topography vector outside soz channels
% i,j = place in the figure/subject

% OUTPUT:
% V: figure with 1 subplots per patient
% black dots: the value of the topography measure outside the soz
% purple dots: the value of the topography measure inside the soz
% orange line: the linear fit trough the data points
%
% Visualize the inter-modal similarity as determined by the Jaccard index and the set of difference connections from structural and effective networks

% Jaccard
J = figure();
pt = 1:size(JI,1);
s = scatter(pt,JI_expected,400,'LineWidth', 1.8); % plot expected jaccard
color = [0.99 0.38 0.22];
s.MarkerEdgeColor = color;
s.Marker = '_';

hold on
s2 = scatter(pt,JI,400,'LineWidth', 1.8); % plot observed jaccard
color = [0.99 0.38 0.22]; %orange
s2.MarkerEdgeColor = color;
s2.Marker = '.';

xlim([0,size(JI,1)])
ylim([0,1])

set(gca, 'XTick',1:13)
set(gca, 'YTick',[0,0.5,1])

set(gca, 'FontSize',16)

ylabel('Jaccard Index')
xlabel('Patient')
legend(' Expected Jaccard Index', ' Observed Jaccard Index','FontSize',16)

% Ratio
R = figure('Renderer', 'painters', 'Position', [600 300 600 400]);
ecog = contains(modality,'ecog');
w = rand(size(ratio(ecog),1),1)+2; % make a XJitter such that the individual points do not overlap
scatter(w,ratio(ecog),40,'o','k') % scatter the ecog patients
hold on
seeg = contains(modality,'seeg');
v = rand(size(ratio(seeg),1),1); % make a XJitter such that the individual points do not overlap
scatter(v,ratio(seeg),400,'.','k','LineWidth', 1.8,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1) % scatter the seeg patients
hold on
yline(0)
ylim([-1.1 1.1])
xlim([-1 4])

end



