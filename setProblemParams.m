function problemParams = setProblemParams

%% Domain
% Rectangular domain: [x1,x2,y1,y2]
problemParams.dom = [0,1,0,1]; 

%% Mesh properties
% Number of mesh subdivision in x and y direction
problemParams.nx = 5;
problemParams.ny = problemParams.nx; 

%% Problem setting
% Diffusion and reaction coefficients
problemParams.nu = 1;
problemParams.sigma = 1;