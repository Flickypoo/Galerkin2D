function approxCG = setMethodOptions

%% Nodal distribution in the reference element
approxCG.nodalDist = 0; % 0: Fekete, 1: Uniform
                           
%% Spatial and functional approximation
approxCG.elemType = 1; % 1: QUA, 2: TRI 
approxCG.degree = 4; % Polynomial degree of approximation

%% Static condensation
% Boolean variable for static condensation (Default: FALSE)
approxCG.isSC = 0;