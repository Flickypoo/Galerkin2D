function [A,b,elemMat] = createFEsystem(mesh,approxCG,refElem,problemParams,solData)

% Indices for assembly
n = mesh.nOfElem*refElem.nOfActiveNodes; 
ind_i  = zeros(1,n);
ind_j  = zeros(1,n);
coef_A = zeros(1,n);
b = zeros(mesh.nOfActiveNodes,1);

% Structure to store information of static condensation
elemMat(mesh.nOfElem).Ae = [];
elemMat(mesh.nOfElem).be = [];
elemMat(mesh.nOfElem).matPost = [];
elemMat(mesh.nOfElem).vecPost = [];


% Compute local contributions
ind = 1;
for iElem = 1:mesh.nOfElem
    Te = mesh.T(iElem,:);
    TeActive = mesh.Tact(iElem,:);
    Xe = mesh.X(Te,:); 
    elemMat(iElem) = computeElemMatrices(Xe,approxCG,refElem,problemParams,solData); 
    
    for irow = 1:refElem.nOfActiveNodes
        for icol = 1:refElem.nOfActiveNodes
            ind_i(ind)  = TeActive(irow);
            ind_j(ind)  = TeActive(icol);
            coef_A(ind) = elemMat(iElem).Ae(irow,icol);
            ind = ind+1;
        end
    end
    % Assembly of the right-hand side of the global system
    b(TeActive) = b(TeActive) + elemMat(iElem).be; 
end

% Assembly of the matrix of the global system
ind_i  = ind_i(1:ind-1);
ind_j  = ind_j(1:ind-1);
coef_A = coef_A(1:ind-1);
A = sparse(ind_i,ind_j,coef_A);