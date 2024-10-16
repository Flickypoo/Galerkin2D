function err = computeErrors(u,mesh,refElem,solData)

% Initialise error
uErrL2 = 0; 
uErrH1semi = 0; 

% Initialise norm of the analytical solution
uL2norm = 0;
uH1semiNorm = 0;

% Element-by-element computation
for ielem = 1:mesh.nOfElem
    Te = mesh.T(ielem,:); 
    Xe = mesh.X(Te,:); 
    ue = u(Te); 
    
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

        x_ig = N_ig*Xe; 
        [u_ex,ux_ex,uy_ex] = computeExactSol(x_ig,solData); 
        
        u_num = N_ig*ue; 
        ux_num = Nx_ig*ue;
        uy_num = Ny_ig*ue;
        
        uErrL2 = uErrL2 + (u_num - u_ex)^2*dvolu;
        uErrH1semi = uErrH1semi + ((ux_num - ux_ex)^2 + (uy_num - uy_ex)^2)*dvolu;

        uL2norm = uL2norm + u_ex^2*dvolu;
        uH1semiNorm = uH1semiNorm + (ux_ex^2 + uy_ex^2)*dvolu;
    end
end

% H1 norm computation
uErrH1 = uErrL2 + uErrH1semi;
uH1norm = uL2norm + uH1semiNorm;

% Relative error
if uL2norm > 1e-9
    uErrL2 = sqrt(uErrL2)/sqrt(uL2norm);
else
    disp('Warning: the error in L2 norm is absolute.')
end

if uH1norm > 1e-9
    uErrH1 = sqrt(uErrH1)/sqrt(uH1norm);
else
    disp('Warning: the error in H1 norm is absolute.')
end

% Create structure to store the error information
err.L2 = uErrL2;
err.H1 = uErrH1;