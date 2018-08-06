function [  ] = plotImage( i, t, frame, pixel_size, a, box_x, box_y, ms, axis_params, time_vec )
%PLOTIMAGE Summary of this function goes here
%   Detailed explanation goes here

if box_x>0 || box_y>0 % if periodic box, then scale markersize with radius
%     hscat=scatter(mod(pos_COMs(t,:,2),box_y)/a,mod(pos_COMs(t,:,1),box_x)/a,ms,pos_COMs(t,:,3)/a,'.');
    dims = size(frame);
    imshow(frame==1, 'XData', [0 dims(2)]*pixel_size/a, 'YData', [0 dims(1)]*pixel_size/a, 'Colormap', gray)
    set(gca,'Ydir','Normal')
    axis on
    axis equal
    axis(axis_params)
    curunits = get(gca, 'Units');
    set(gca, 'Units', 'Points');
    cursize = get(gca, 'Position');
    set(gca, 'Units', curunits);             
    Lybox = cursize(3); 
    LyboxInRadii = box_y/a;   
    ms = 10;
    ms = pi*(Lybox/ LyboxInRadii)^2*ms; %markersize
%     hscat=scatter(mod(pos_COMs(t,:,2),box_y)/a,mod(pos_COMs(t,:,1),box_x)/a,ms,pos_COMs(t,:,3)/a,'.');
else
    ms = 10;
%     hscat=scatter(mod(pos_COMs(t,:,2),box_y)/a,mod(pos_COMs(t,:,1),box_x)/a,ms,pos_COMs(t,:,3)/a,'.');
end
colormap(parula)
caxis([0 0])


end

