function plotSol(u,mesh)

% Construct a Delaunay triangulation for visualisation
tri = delaunay(mesh.X(:,1),mesh.X(:,2)); 
figure; 
trisurf(tri,mesh.X(:,1),mesh.X(:,2),u,'EdgeColor','none','FaceCOlor','interp')
camlight
colorbar
set(gca,'FontSize',16)