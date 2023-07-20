% STReEF06_figures
% make figures to present the inter-modal similarity and network topography comparison between structural (derived from DWI) and effective (derived from SPES) networks

% author: Susanne Jelsma & Dorien van Blooijs
% date: June 2022

% visualize the inter-modal similarity with connectivity matrices and visualize the network topography per patient

%% visualize the inter-modal similarity with connectivity matrices

for subj = 1:size(dataBase,2)
V = visual_networks(dataBase(subj).dwi.SC_matrix,dataBase(subj).ccep.EC_matrix,dataBase(subj).sub_label) % function to plot the connectivity matrices
V.WindowState = 'maximized';
print('-vector','-depsc',V,sprintf('visual_networks_symetric%s',dataBase(subj).sub_label))% save the figure for further processing with Adobe Illustrator
end
%%  visualize the network topography per patient (must become a function)
% example of a plot with degree
markers = ['o';'+';'*';'.';'x';'s';'d';'p';'^';'v';'>';'<';'_'];
d_SC = [];d_EC = [];bc_SC = [];bc_EC = [];cc_SC = [];cc_EC = [];
i=0;j=8; % j=8
set(gcf,'renderer','Painters')
f4 = figure(5);
f4.WindowState = 'maximized';

for subj = 1:size(dataBase,2)
    i=i+1;j=j+1;
sub_label = erase(dataBase(subj).sub_label,'sub-');
SC = dataBase(subj).dwi.SC_matrix;
EC =  dataBase(subj).ccep.EC_matrix;  

ch = dataBase(subj).ccep.ch;
elec_include = dataBase(subj).dwi.elec_include;
ch_include = ch(elec_include);
soz_ind = dataBase(subj).soz_ind;
no_soz_ind = dataBase(subj).no_soz_ind;

% compute the DEGREE
G_SC = graph(SC,ch_include); % graph structure
G_EC = graph(EC,ch_include);

degree_SC = degree(G_SC);
degree_EC = degree(G_EC);
degree_SC_soz = degree_SC(soz_ind);
degree_EC_soz = degree_EC(soz_ind);
degree_SC_nsoz = degree_SC(no_soz_ind);
degree_EC_nsoz = degree_EC(no_soz_ind);

 
nc_EC = degree_EC_soz;nc_SC = degree_SC_soz;
nnc_EC = degree_EC_nsoz;nnc_SC = degree_SC_nsoz;

subplot(4,2,i)

ns = scatter(nnc_EC,nnc_SC,60,'LineWidth', 1.8); %200 80
color = 'k'; %[0.7,0.7,0.7];
ns.MarkerEdgeColor = color;
ns.Marker = markers(subj);
hold on
s = scatter(nc_EC,nc_SC,60,'LineWidth', 1.8); 
color = [0.4940 0.1840 0.5560];
s.MarkerEdgeColor = color;
s.Marker = markers(subj);

yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
xlim([0,round(max(degree_EC)/5)*5])
set(gca, 'XTick',[min(xt),round(max(degree_EC/2)/5)*5,round(max(degree_EC)/5)*5])%numelxt
set(gca, 'XTickLabel',[min(xt),round(max(degree_EC/2)/5)*5,round(max(degree_EC)/5)*5])%numelxt

set(gca, 'YTick',[min(yt),round(max(degree_SC/2)/5)*5,round(max(degree_SC)/5)*5])
set(gca, 'YTickLabel',[min(yt),round(max(degree_SC/2)/5)*5,round(max(degree_SC)/5)*5])

[P,S] = polyfit(degree_EC,degree_SC,1);
[y_fit, delta] = polyval(P,degree_EC,S);

[RHO,PVAL] = corr(degree_SC,degree_EC,'Type','Spearman'); % compute the correlation between the degrees
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
plot(degree_EC,y_fit,'Color',color,'LineWidth',2)
set(gca, 'FontSize',12)

if ismember(subj,[1,6,8,9,10])
    title(sprintf('%s Degree patient %g, p %s , r_{s}= %g','ECoG',j,PVAL,round(RHO,2)))
else
title(sprintf('%s Degree patient %g, p %s , r_{s}= %g','sEEG',j,PVAL,round(RHO,2)))
end
end

saveas(f4,'correlation degree all','epsc') % save the figure for further processing with Adobe Illustrator