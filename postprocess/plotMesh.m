function plotMesh(mesh,refElem,approxCG)

deg = refElem.degree; 
aux = 1:deg-1; 
if refElem.type == 1
    faceNodesOrder = [1, aux+4, 2, aux+4+deg-1, 3,aux+4+2*(deg-1), 4, aux+4+3*(deg-1), 1]; 
elseif refElem.type == 2
    faceNodesOrder = [1, aux+3, 2, aux+3+deg-1,3, aux+3+2*(deg-1),1];
end

figure
hold on
plot(mesh.X(:,1), mesh.X(:,2),'bo');
for ielem = 1:mesh.nOfElem
    Te = mesh.T(ielem,:); 
    Xe = mesh.X(Te,:); 
    plot(Xe(faceNodesOrder,1), Xe(faceNodesOrder,2),'k');
end
axis equal; axis tight

if approxCG.isSC == 1
    figure
    hold on
    % TO BE COMPLETED: Plot the active nodes of the statically condensed 
    % problem
    %
    % plot(,'rd');
    for ielem = 1:mesh.nOfElem
        Te = mesh.T(ielem,:); 
        Xe = mesh.X(Te,:); 
        plot(Xe(faceNodesOrder,1), Xe(faceNodesOrder,2),'k');
    end
    axis equal; axis tight
end