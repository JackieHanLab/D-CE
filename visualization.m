%visualization of spatial reconstruction of geo-seq dataset
coords = readtable("3Dcoordinate.txt");
coords = table2array(coords);
d=size(coords);
[x y z]=sphere(7);
load('info.mat');
rad=0.05
%visualization the reconstructed structure colored by position
for i = 1:1:d(1)
    surf(rad*x+coords(i,1), rad*y+coords(i,3), rad*z+coords(i,2),'EdgeColor',color(i,:),'FaceColor',color(i,:),'linestyle','none','EdgeColor','none')
    hold on
end
lighting phong;
% change the orintation of the reconstructed structure to get better visualization
h=allchild(gca)
direction = [0 0 1];
rotate(h,direction,75)
direction = [0 1 0];
rotate(h,direction,15)
direction = [1 0 0];
rotate(h,direction,-15)
direction = [0 0 1];
rotate(h,direction,-21)
set(gca, 'XDir','reverse');
set(gca, 'ZDir','reverse');
light('Position',[0.58674 -0.05336 -0.80801]);
savefig('E7.5_all_germ_layer.fig')
close all
%visualization of each germ layer
lab = unique(germlayerinfo)
for j = 1:1:3
%j=1,2,3 for ectoderm, endoderm and mesoderm 
layerinfo2=layerinfo(germlayerinfo==lab(j),1);
coordstolink = coords(germlayerinfo==lab(j),:);
d=size(coordstolink)
if(mod(d(1),2)==1)
    layerinfo2=layerinfo2(2:d(1))
    coordstolink=coordstolink(2:d(1),:)
end
for i = min(layerinfo2):max(layerinfo2)
        try
        APnode = coordstolink(layerinfo2==i,:);
        APnode = APnode(1:2,:);
        colorlink = [1-0.09*i 1-0.09*i 1-0.09*i]
        colorlink(j)=1
        surf(rad*x+APnode(1,1), rad*y+APnode(1,3), rad*z+APnode(1,2),'EdgeColor',[0.9 0.9 0.9],'FaceColor',[0.9 0.9 0.9],'linestyle','none','EdgeColor','none')
        hold on
        surf(rad*x+APnode(2,1), rad*y+APnode(2,3), rad*z+APnode(2,2),'EdgeColor',[0.9 0.9 0.9],'FaceColor',[0.9 0.9 0.9],'linestyle','none','EdgeColor','none')
        hold on
        plot3([APnode(1,1) APnode(2,1)],[APnode(1,3) APnode(2,3)], [APnode(1,2) APnode(2,2)], ...
            'Color', colorlink, 'LineWidth', 5)
        hold on
        end
end
%link L/R in ectoderm
if(j==1)
    layerinfo2=layerinfo(germlayerinfo==lab(j),1);
    coordstolink = coords(germlayerinfo==lab(j),:);
    i=2
    APnode = coordstolink(layerinfo2==i,:);
    colorlink = [1-0.09*i 1-0.09*i 1-0.09*i]
    colorlink(1)=1
    surf(rad*x+APnode(3,1), rad*y+APnode(3,3), rad*z+APnode(3,2),'EdgeColor',[0.9 0.9 0.9],'FaceColor',[0.9 0.9 0.9],'linestyle','none','EdgeColor','none')
        hold on
    surf(rad*x+APnode(4,1), rad*y+APnode(4,3), rad*z+APnode(4,2),'EdgeColor',[0.9 0.9 0.9],'FaceColor',[0.9 0.9 0.9],'linestyle','none','EdgeColor','none')
        hold on
    plot3([APnode(3,1) APnode(4,1)],[APnode(3,3) APnode(4,3)], [APnode(3,2) APnode(4,2)], ...
       'Color', colorlink, 'LineWidth', 5)
    hold on
    plot3([APnode(3,1) APnode(1,1)],[APnode(3,3) APnode(1,3)], [APnode(3,2) APnode(1,2)], ...
       'Color', colorlink, 'LineWidth', 5)
    hold on
    plot3([APnode(2,1) APnode(4,1)],[APnode(2,3) APnode(4,3)], [APnode(2,2) APnode(4,2)], ...
       'Color', colorlink, 'LineWidth', 5)
    hold on
    i=9
    APnode = coordstolink(layerinfo2==i,:);
    colorlink = [1-0.09*i 1-0.09*i 1-0.09*i]
    colorlink(1)=1
    surf(rad*x+APnode(3,1), rad*y+APnode(3,3), rad*z+APnode(3,2),'EdgeColor',[0.9 0.9 0.9],'FaceColor',[0.9 0.9 0.9],'linestyle','none','EdgeColor','none')
        hold on
    surf(rad*x+APnode(4,1), rad*y+APnode(4,3), rad*z+APnode(4,2),'EdgeColor',[0.9 0.9 0.9],'FaceColor',[0.9 0.9 0.9],'linestyle','none','EdgeColor','none')
        hold on
    plot3([APnode(3,1) APnode(4,1)],[APnode(3,3) APnode(4,3)], [APnode(3,2) APnode(4,2)], ...
       'Color', colorlink, 'LineWidth', 5)
    hold on
    plot3([APnode(3,1) APnode(2,1)],[APnode(3,3) APnode(2,3)], [APnode(3,2) APnode(2,2)], ...
       'Color', colorlink, 'LineWidth', 5)
    hold on
    plot3([APnode(1,1) APnode(4,1)],[APnode(1,3) APnode(4,3)], [APnode(1,2) APnode(4,2)], ...
       'Color', colorlink, 'LineWidth', 5)
    hold on
    for i = 3:max(layerinfo2)
        try
        APnode = coordstolink(layerinfo2==i,:);
        distAPLR = pdist2(APnode,APnode,'euclidean')
        colorlink = [1-0.09*i 1-0.09*i 1-0.09*i]
        colorlink(1)=1
        surf(rad*x+APnode(3,1), rad*y+APnode(3,3), rad*z+APnode(3,2),'EdgeColor',[0.9 0.9 0.9],'FaceColor',[0.9 0.9 0.9],'linestyle','none','EdgeColor','none')
        hold on
        surf(rad*x+APnode(4,1), rad*y+APnode(4,3), rad*z+APnode(4,2),'EdgeColor',[0.9 0.9 0.9],'FaceColor',[0.9 0.9 0.9],'linestyle','none','EdgeColor','none')
        hold on
        surf(rad*x+APnode(5,1), rad*y+APnode(5,3), rad*z+APnode(5,2),'EdgeColor',[0.9 0.9 0.9],'FaceColor',[0.9 0.9 0.9],'linestyle','none','EdgeColor','none')
        hold on
        surf(rad*x+APnode(6,1), rad*y+APnode(6,3), rad*z+APnode(6,2),'EdgeColor',[0.9 0.9 0.9],'FaceColor',[0.9 0.9 0.9],'linestyle','none','EdgeColor','none')
        hold on
        plot3([APnode(3,1) APnode(5,1)],[APnode(3,3) APnode(5,3)], [APnode(3,2) APnode(5,2)], ...
        'Color', colorlink, 'LineWidth', 5)
        hold on
        plot3([APnode(4,1) APnode(6,1)],[APnode(4,3) APnode(6,3)], [APnode(4,2) APnode(6,2)], ...
        'Color', colorlink, 'LineWidth', 5)
        hold on
        N=max(min(distAPLR([3 5],1)))    %或者N=max(A(:))
        [r,c]=find(N==distAPLR)
        plot3([APnode(r(1),1) APnode(c(1),1)],[APnode(r(1),3) APnode(c(1),3)], [APnode(r(1),2) APnode(c(1),2)], ...
        'Color', colorlink, 'LineWidth', 5)
        hold on
        N=max(min(distAPLR([4 6],2)))    %或者N=max(A(:))
        [r,c]=find(N==distAPLR)
        plot3([APnode(r(1),1) APnode(c(1),1)],[APnode(r(1),3) APnode(c(1),3)], [APnode(r(1),2) APnode(c(1),2)], ...
        'Color', colorlink, 'LineWidth', 5)
        hold on
        N=min(min(distAPLR([3 5],[4 6])))    %或者N=max(A(:))
        [r,c]=find(N==distAPLR)
        plot3([APnode(r(1),1) APnode(c(1),1)],[APnode(r(1),3) APnode(c(1),3)], [APnode(r(1),2) APnode(c(1),2)], ...
       'Color', colorlink, 'LineWidth', 5)
        hold on
    end
end

end
h=allchild(gca)
direction = [0 0 1];
rotate(h,direction,60)
direction = [0 1 0];
rotate(h,direction,15)
set(gca, 'XDir','reverse');
set(gca, 'ZDir','reverse');
figurename = join(['E7.5_',lab(j),'_only.fig'],'')
savefig(char(figurename))
close all;
end
