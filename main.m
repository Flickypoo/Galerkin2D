% This codes uses continuous Galerkin FEM to solve the 2D
% diffusion-reaction equation with Dirichlet BC on the entire boundary
%
% Options:
%  - Triangular and quadrilateral mesh elements 
%  - Polynomial degree from 1 to 5
%  - Allows static condensation 

clearvars;
close all;
clc

setPath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
% Problem domain and parameters
problemParams = setProblemParams;

% Method options
approxCG = setMethodOptions;

% Reference element
refElem = setReferenceElement(approxCG);

% Computational mesh
mesh = createMeshRect(problemParams,refElem,approxCG); 
plotMesh(mesh,refElem,approxCG);

% Data for analytical solution and BC
solData = setSolData(mesh);

%% Construct the FE system
[A,b,elemMat] = createFEsystem(mesh,approxCG,refElem,problemParams,solData); 

%% Impose Dirichlet boundary conditions (System reduction)  
[A_red,b_red] = imposeDirBC(A,b,mesh,solData);

%% Solve the system
u = solveFEsystem(A_red,b_red,elemMat,mesh,approxCG,refElem,solData);

%% Compute errors
err = computeErrors(u,mesh,refElem,solData); 

%% Postprocess
plotSol(u,mesh); % not accurate, just for visualisation purposes

%% Print log
printLog(mesh,approxCG,err);