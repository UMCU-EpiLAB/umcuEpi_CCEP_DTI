function T = visual_topology(vector_v1, vector_v2, elec_info, elec_include, i)
% INPUT:
% vector_v1 = topology vector modality 1
% vector_v2 = topology vector modality 2
% v1_soz = topology vector soz channels
% v1_nsoz = topology vector outside soz channels
% i = place in the figure/subject

% OUTPUT:
% T: figure with 1 subplots per patient
% black dots: the value of the topology measure outside the soz
% red dots: the value of the topology measure inside the soz
% orange line: the linear fit trough the data points
% 
% Visualize the correlation between a topology measure from structural and effective networks
soz_all = strcmpi(elec_info.soz,'yes'); % electrodes in soz
soz = soz_all(elec_include); 

marker = '.'; 

set(gcf,'renderer','Painters')
T = figure(1);
T.WindowState = 'maximized';
 
subplot(2,3,i)
% plot the topology values outside the soz
ns = scatter(vector_v1(~soz),vector_v2(~soz),400,'LineWidth', 1.8); 
color = 'k'; 
ns.MarkerEdgeColor = color;
ns.Marker = marker; 
hold on

% plot the topology values inside the soz
s = scatter(vector_v1(soz),vector_v2(soz),400,'LineWidth', 1.8); 
color = [204/250 37/250 41/250]; %red
s.MarkerEdgeColor = color;
s.Marker = marker; 

% set limits and ticks
xlim([0,round(max([vector_v1;vector_v2])/5)*5])
ylim([0,round(max([vector_v1;vector_v2])/5)*5])

set(gca, 'XTick',[0,round(max([vector_v1;vector_v2]/2)/5)*5,round(max([vector_v1;vector_v2])/5)*5])%numelxt
set(gca, 'XTickLabel',[0,round(max([vector_v1;vector_v2]/2)/5)*5,round(max([vector_v1;vector_v2])/5)*5])%numelxt

set(gca, 'YTick',[0,round(max([vector_v1;vector_v2]/2)/5)*5,round(max([vector_v1;vector_v2])/5)*5])
set(gca, 'YTickLabel',[0,round(max([vector_v1;vector_v2]/2)/5)*5,round(max([vector_v1;vector_v2])/5)*5])

axis square

% do a linear fit
[P,S] = polyfit(vector_v1,vector_v2,1);
[y_fit,~] = polyval(P,vector_v1,S);

%  compute the correlation between the topology measure
[RHO,PVAL] = corr(vector_v2,vector_v1,'Type','Spearman'); 

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
plot(vector_v1,y_fit,'Color',color,'LineWidth',2.8)
set(gca, 'FontSize',16) %12

title(sprintf('p %s , r_{s}= %g', PVAL,round(RHO,2)))

end
