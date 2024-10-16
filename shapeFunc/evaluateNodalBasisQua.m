function [N,dNdxi,dNdeta]=evaluateNodalBasisQua(XI,XInodes,degree)
% [N,dNdxi,dNdeta]evaluateNodalBasisQua(XIs,XInodes,degree)
% Evaluates at XI the nodal basis of polynomials for the given degree
% and the nodes in XInodes.



p = degree; nen = (p+1)^2;
points = XI;
npt = size(points,1);

coord = XInodes([1, 5:p+3, 2],1);
[nx,dnx] = evaluateNodalBasis1D(points(:,1),coord,degree); 
[ny,dny] = evaluateNodalBasis1D(points(:,2),coord,degree); 

nx = reshape(nx,npt*(p+1),1); nx = repmat(nx,1,p+1);
ny = repmat(ny,p+1,1);
dnx = reshape(dnx,npt*(p+1),1); dnx = repmat(dnx,1,p+1);
dny = repmat(dny,p+1,1);

nodes = 1:(p+1)^2;
nodes = reshape(nodes,p+1,p+1)';
nodes(1,:)=[]; nodes(end,:)=[]; nodes(:,1)=[]; nodes(:,end)=[];
ordre_int = reshape(nodes',1,(p-1)^2);
ordre = [1,p+1,(p+1)^2,p*(p+1)+1,...
    2:p, 2*(p+1):p+1:p*(p+1), ...
    (p+1)^2-1:-1:p*(p+1)+2, ...
    (p-1)*(p+1)+1:-(p+1):p+2, ordre_int];
N  = nx.*ny;  N  = reshape(N,npt,nen);  N  = N(:, ordre);
Nx = dnx.*ny; Nx = reshape(Nx,npt,nen); Nx = Nx(:, ordre);
Ny = nx.*dny; Ny = reshape(Ny,npt,nen); Ny = Ny(:, ordre);

dNdxi = Nx; dNdeta = Ny;



