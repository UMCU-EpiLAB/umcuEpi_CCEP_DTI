function V = visual_topography(vector_v1, vector_v2, v1_soz, v1_nsoz, v2_soz, v2_nsoz,subj)
% INPUT:
% cfg.sub_label = {'RESPXXXX', 'RESPXXXX'}
% vector_v1 = topography vector modality 1
% vector_v2 = topography vector modality 2
% v1_soz = topography vector soz channels
% v1_nsoz = topography vector outside soz channels

% OUTPUT:
% V: figure with 1 subplots per patient
% black dots: the value of the topography measure outside the soz
% purple dots: the value of the topography measure inside the soz
% orange line: the linear fit trough the data points
% 
% Visualize the correlation between a topography measure from structural and effective networks

markers = ['o';'+';'*';'.';'x';'s';'d';'p';'^';'v';'>';'<';'_'];
i=0;j=0; 
set(gcf,'renderer','Painters')
V = figure(2);
V.WindowState = 'maximized';
i=i+1;j=j+1;

subplot(4,2,i)

ns = scatter(v1_nsoz,v2_nsoz,60,'LineWidth', 1.8); 
color = 'k'; %[0.7,0.7,0.7];
ns.MarkerEdgeColor = color;
ns.Marker = markers(subj);
hold on
s = scatter(v1_soz,v2_soz,60,'LineWidth', 1.8); 
color = [0.4940 0.1840 0.5560];
s.MarkerEdgeColor = color;
s.Marker = markers(subj);

yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
xlim([0,round(max(vector_v1)/5)*5])
set(gca, 'XTick',[min(xt),round(max(vector_v1/2)/5)*5,round(max(vector_v1)/5)*5])%numelxt
set(gca, 'XTickLabel',[min(xt),round(max(vector_v1/2)/5)*5,round(max(vector_v1)/5)*5])%numelxt

set(gca, 'YTick',[min(yt),round(max(vector_v2/2)/5)*5,round(max(vector_v2)/5)*5])
set(gca, 'YTickLabel',[min(yt),round(max(vector_v2/2)/5)*5,round(max(vector_v2)/5)*5])

[P,S] = polyfit(vector_v1,vector_v2,1);
[y_fit, delta] = polyval(P,vector_v1,S);

[RHO,PVAL] = corr(vector_v2,vector_v1,'Type','Spearman'); % compute the correlation between the degrees
% make the p-values presentable in a figure
if PVAL < 0.0001
    PVAL = '< 0.0001';
elseif PVAL < 0.001
    PVAL = '< 0.001';
elseif PVAL < 0.01
    PVAL = '< 0.01';
elseif PVAL < 0.05
    PVAL = '< 0.05';
else 
    PVAL = sprintf('= %g',round(PVAL,2));
end
            
hold on
color = [0.99 0.38 0.22];

% plot polyfit throught data points
plot(vector_v1,y_fit,'Color',color,'LineWidth',2)
set(gca, 'FontSize',12)

if ismember(subj,[1,6,8,9,10])
    title(sprintf('%s Degree patient %g, p %s , r_{s}= %g','ECoG',j,PVAL,round(RHO,2)))
else
title(sprintf('%s Degree patient %g, p %s , r_{s}= %g','sEEG',j,PVAL,round(RHO,2)))
end
