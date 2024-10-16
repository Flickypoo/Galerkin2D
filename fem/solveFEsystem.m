function u = solveFEsystem(A_red,b_red,elemMat,mesh,approxCG,refElem,solData)

% Solution of the reduced system
sol = A_red \ b_red;

% Solution accounting for Dirichlet BC
u = zeros(mesh.nOfNodes,1); 

if approxCG.isSC == 1
    u_tmp = zeros(mesh.nOfActiveNodes,1); 
    u_tmp(mesh.dirNodes) = solData.uD; 
    u_tmp(mesh.solNodes) = sol; 

    for iElem = 1:mesh.nOfElem
        Te = mesh.T(iElem,:);
        TeActive = mesh.Tact(iElem,:);

        % TO BE COMPLETED: determine the solution in all the nodes, both 
        % active and condensed
        %
        %u(Te(1:refElem.nOfActiveNodes)) = ;
        %u(Te(refElem.nOfActiveNodes+1:refElem.nOfElemNodes)) = ;
    end
else
    u(mesh.dirNodes) = solData.uD; 
    u(mesh.solNodes) = sol; 
end

