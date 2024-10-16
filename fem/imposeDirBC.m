function [A_red,b_red] = imposeDirBC(A,b,mesh,solData)

b_tmp = b - A(:,mesh.dirNodes)*solData.uD; 
A_red = A(mesh.solNodes,mesh.solNodes); 
b_red = b_tmp(mesh.solNodes);