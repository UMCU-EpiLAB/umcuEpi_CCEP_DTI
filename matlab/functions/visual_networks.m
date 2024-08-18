function V = visual_networks(matrix_m1, matrix_m2, subj)
% INPUT:
% sub_label = {'STREEFXXXX', 'STREEFXXXX'}
% matrix_m1 = connectivity matrix modality 1
% matrix_m2 = connectivity matrix modality 2

% OUTPUT:
% V: figure with 3 subplots
% 1: the visual network of modality one (blue)
% 2: the visual network of modality two (green)
% 3: the union (orange) and intersection (blue & green) of the two visual
% 
% Visualize the inter-modal similarity between structural and effective networks
V = figure ();
V.WindowState = 'maximized';
f =17;
ax1 = subplot(1,3,1);
imagesc(matrix_m1) 
set(gca, 'YDir','reverse')

map_umc = [1 1 1;0.07 0.57 0.98];
colormap(ax1,map_umc)

daspect([1 1 1]) 
for x=1:size(matrix_m1,1)
    xline(x+0.5)
    yline(x+0.5)
end
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel(sprintf('Electrode area 1-%.2d',size(matrix_m2,1)),'Fontsize',f)
ylabel(sprintf('Electrode area 1-%.2d',size(matrix_m2,1)),'Fontsize',f)


ax2 = subplot(1,3,2);
imagesc(matrix_m2)
set(gca, 'YDir','reverse')
map_umc2 = [1 1 1;0 0.78 0.47];
colormap(ax2,map_umc2)
daspect([1 1 1])

siz = size(matrix_m2,1)/100;
t = text(ax2,25*siz,-15*siz,sprintf('Patient %s',subj));
t.FontSize = 16;
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel(sprintf('Electrode area 1-%.2d',size(matrix_m2,1)),'Fontsize',f)
ylabel(sprintf('Electrode area 1-%.2d',size(matrix_m2,1)),'Fontsize',f)

for x=1:size(matrix_m1,1)
    xline(x+0.5)
    yline(x+0.5)
end

ax3 = subplot(1,3,3);
EC_deel = matrix_m2;
EC_deel(matrix_m2>0) = 0.5;
overlap = EC_deel + matrix_m1; % calculate the intersection and union (intersection 1.5, matrix_m2 0.5, matrix_m1 1)
imagesc(overlap)
set(gca, 'YDir','reverse')
map2 = [1 1 1;0 0.78 0.47; 0.07 0.57 0.98 ;0.99 0.38 0.22];
colormap(ax3,map2)
daspect([1 1 1])

for x=1:size(matrix_m1,1)
    xline(x+0.5)
    yline(x+0.5)
end

set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel(sprintf('Electrode area 1-%.2d',size(matrix_m2,1)),'Fontsize',f)
ylabel(sprintf('Electrode area 1-%.2d',size(matrix_m2,1)),'Fontsize',f)
title(ax2,'Effective','Fontsize',f);
title(ax1,'Structural','Fontsize',f);
title(ax3,'Union and Intersection','Fontsize',f) 
end


