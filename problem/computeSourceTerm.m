function s = computeSourceTerm(pts,problemParams,solData)

% Data
nu = problemParams.nu;
sigma = problemParams.sigma;
alpha = solData.alpha;
beta = solData.beta;


% Coordinates
x = pts(:,1);
y = pts(:,2); 

s = sin(alpha*pi*x)*sin(beta*pi*y)*(nu*pi^2*(alpha^2 + beta^2) + sigma); 