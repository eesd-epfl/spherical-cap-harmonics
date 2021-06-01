function fig = paraview_patch(v, f, maps)
% Plot STL meshes with Paraview maps.
% Author: Mahmoud Shaqfa (EPFL)

fig = figure; 
% set(fig,'renderer','painters'); % To save vectorized SVG and PDF formats

if nargin < 3
    patch('Faces',f,'Vertices',v,'FaceColor',[42,63,188]./255,'LineWidth',0.5);
else
    patch('Faces',f,'Vertices',v,'FaceColor','interp','FaceVertexCData',maps,...
        'EdgeColor','k', 'LineWidth', 0.1, 'LineStyle', '-');
%     patch('Faces',f,'Vertices',v,'FaceColor','interp','FaceVertexCData',maps,...
%         'EdgeColor','k', 'LineStyle', 'none');
     
    % Define a color map similar to Paraview's shades
    red_color = zeros(1, 255); green_color = zeros(1, 255); blue_color = zeros(1, 255);
    red_color(1:127) = linspace(42,220,127); red_color(128:255) = linspace(220,174,127+1);
    green_color(1:127) = linspace(63,220,127); green_color(128:255) = linspace(220,0,127+1);
    blue_color(1:127) = linspace(181,220,127); blue_color(128:255) = linspace(220,22,127+1);
    
    ParaviewMap = [red_color', green_color', blue_color']./255;
    colormap(ParaviewMap);
    set(gcf,'color','w');
    cb = colorbar;
    set(cb,'position',[.15 .1 .05 .2])
end
axis equal tight off
view(45, 45);
gca.Clipping = 'off';