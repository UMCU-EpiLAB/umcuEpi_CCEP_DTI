function V = visual_VEA(vector_VEA, soz, count)
% INPUT:
% vector_VEA = vector with volume of electrode contact areas (VEA)
% soz =  containing soz information
% count= place in the figure/subject

% OUTPUT:
% T: figure with 1 subplots per patient
% black bars: the number of electrode contact areas with a certain volume outside the soz
% red bars: the number of electrode contact areas with a certain volume inside the soz
% 
% Make a histogram of the volume of electrode contact areas

set(gcf,'renderer','Painters')
V = figure(1);
V.WindowState = 'maximized';
 
subplot(2,3,count)
histogram(vector_VEA(~soz),'BinWidth',16,'BinLimits',[0,64],'FaceColor','k','FaceAlpha',1,"LineStyle","-") % plot volume electrode contact areas
axis square
yt = get(gca, 'YTick');
hold on
if ~isempty(vector_VEA(soz))
histogram(vector_VEA(soz),'BinWidth',16,'BinLimits',[0,64],'FaceColor',[204/250 37/250 41/250],'FaceAlpha',1,'EdgeColor','none') % plot volume electrode contact areas in the soz
end
histogram(vector_VEA(~soz),'BinWidth',16,'BinLimits',[0,64],'EdgeColor','k','DisplayStyle','stairs') % plot line histogram of non soz for when soz counts overlap the non soz counts

xlim([0,64])

set(gca, 'XTick',[0,16,32,48,64])
set(gca, 'XTickLabel',[0,16,32,48,64])
set(gca, 'YTick',[0,max(yt)/2,max(yt)])     
set(gca, 'YTickLabel',[0,max(yt)/2,max(yt)])
set(gca, 'FontSize',16) 
end

