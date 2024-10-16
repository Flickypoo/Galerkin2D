function elemMat = computeElemMatrices(Xe,approxCG,refElem,problemParams,solData)

% Allocate local matrices and vector
Ke = zeros(refElem.nOfElemNodes);
Me = zeros(refElem.nOfElemNodes);
be_tmp = zeros(refElem.nOfElemNodes,1);

% Compute elemental contributions
for ig = 1:refElem.nIP
    N_ig = refElem.N(ig,:); 
    Nxi_ig = refElem.Nxi(ig,:); 
    Neta_ig = refElem.Neta(ig,:); 
    
	J = [
        Nxi_ig*Xe(:,1)    Nxi_ig*Xe(:,2)
        Neta_ig*Xe(:,1)   Neta_ig*Xe(:,2)
        ];
	dvolu = refElem.wIP(ig)*det(J); 
    grad = J\[Nxi_ig; Neta_ig]; 
    Nx_ig = grad(1,:); 
    Ny_ig = grad(2,:); 
    
    % Stiffness and mass matrices
    Ke = Ke + (Nx_ig'*Nx_ig + Ny_ig'*Ny_ig)*dvolu; 
    Me = Me + N_ig'*N_ig*dvolu; 
    x_ig = N_ig*Xe; 
    be_tmp = be_tmp + computeSourceTerm(x_ig,problemParams,solData)*N_ig'*dvolu; 
end

% Problem matrix
Ae_tmp = problemParams.nu*Ke + problemParams.sigma*Me; 

% Reduce the matrix if static condensation is active
if approxCG.isSC == 1
    % TO BE COMPLETED: compute the matrices required for static
    % condensation
    % a: active; c: condensed
    %
    %A_aa = ;
    %A_cc = ;
    %A_ac = ;
    %A_ca = ;
    %f_a  = ;
    %f_c  = ;
    %
    %elemMat.matPost = ;
    %elemMat.vecPost = ;
    %elemMat.Ae = ;
    %elemMat.be = ;
else
    elemMat.Ae = Ae_tmp;
    elemMat.be = be_tmp;
    elemMat.matPost = [];
    elemMat.vecPost = [];
end