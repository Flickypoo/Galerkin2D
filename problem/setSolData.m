function solData = setSolData(mesh)

%% Analytical solution 
% Integer coefficients in the definition of solution (and source term)
% u = sin(alpha*pi*x).*sin(beta*pi*y); 
solData.alpha = 1;
solData.beta = 1;

%% Value of the Dirichlet boundary datum
solData.uD = zeros(size(mesh.dirNodes)); 