function mesh = createMeshRect(problemParams,refElem,approxCG)

% Discretisation details
nx = problemParams.nx;
ny = problemParams.ny;
elemType = approxCG.elemType;
degree = approxCG.degree;

% Domain
x1 = problemParams.dom(1); 
x2 = problemParams.dom(2);
y1 = problemParams.dom(3); 
y2 = problemParams.dom(4);
domSize = abs(x2-x1)*abs(y2-y1);

% Nodal coordinates
npx = degree*nx + 1;
npy = degree*ny + 1;
npt = npx*npy;
x = linspace(x1,x2,npx);
y = linspace(y1,y2,npy);
[x,y] = meshgrid(x,y);
X = [reshape(x',npt,1), reshape(y',npt,1)];

% Mesh connectivity
if elemType == 1
    if degree == 1
        T = zeros(nx*ny,4);
        for i=1:ny
            for j=1:nx
                iElem = (i-1)*nx+j;
                inode = (i-1)*(npx)+j;
                T(iElem,:) = [inode   inode+1   inode+npx+1   inode+npx];
            end
        end
    elseif degree == 2
        T = zeros(nx*ny,9);
        for i=1:ny
            for j=1:nx
                iElem = (i-1)*nx + j;
                inode = (i-1)*2*npx + 2*(j-1) + 1;
                nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
                T(iElem,:) = nodes_aux([1  3  9  7  2  6  8  4  5]);
            end
        end
    elseif degree == 3
        T = zeros(nx*ny,16);
        for i=1:ny
            for j=1:nx
                iElem = (i-1)*nx + j;
                inode = (i-1)*3*npx + 3*(j-1) + 1;
                nodes_aux = [inode+(0:3)  inode+npx+(0:3)  ...
                    inode+2*npx+(0:3)  inode+3*npx+(0:3)];
                T(iElem,:) = nodes_aux([1  4  16  13  2  3  8  12  ...
                    15  14  9  5  6  7  10  11]);
            end
        end
    elseif degree == 4
        T = zeros(nx*ny,25);
        for i=1:ny
            for j=1:nx
                iElem = (i-1)*nx + j;
                inode = (i-1)*4*npx + 4*(j-1) + 1;
                nodes_aux = [inode+(0:4)  inode+npx+(0:4)  ...
                    inode+2*npx+(0:4)  inode+3*npx+(0:4)  inode+4*npx+(0:4)];
                T(iElem,:) = nodes_aux([1  5  25  21  2  3  4  10  15  20 ...
                    24  23  22  16  11  6  7  8  9  12  13  14  17  18  19]);
            end
        end
    elseif degree == 5
        T = zeros(nx*ny,36);
        for i=1:ny
            for j=1:nx
                iElem = (i-1)*nx + j;
                inode = (i-1)*5*npx + 5*(j-1) + 1;
                nodes_aux = [inode+(0:5)  inode+npx+(0:5)  inode+2*npx+(0:5) ...
                    inode+3*npx+(0:5)  inode+4*npx+(0:5)  inode+5*npx+(0:5)];
                T(iElem,:) = nodes_aux([1  6  36  31  ...
                    2  3  4  5  12  18  24  30  35  34  33  32  25  19  13  7 ...
                    8  9  10  11  14  15  16  17  20  21  22  23  26  27  28  29]);
            end
        end
    elseif degree == 6
        T = zeros(nx*ny,49);
        for i=1:ny
            for j=1:nx
                iElem = (i-1)*nx + j;
                inode = (i-1)*6*npx + 6*(j-1) + 1;
                nodes_aux = [inode+(0:6)  inode+npx+(0:6)  inode+2*npx+(0:6) ...
                    inode+3*npx+(0:6)  inode+4*npx+(0:6)  inode+5*npx+(0:6) ...
                    inode+6*npx+(0:6)];
                T(iElem,:) = nodes_aux([1  7  49  43  ...
                    2  3  4  5  6  14  21  28  35  42  48  47  46  45  44  36  29  22 15  8 ...
                    9:13   16:20   23:27   30:34   37:41]);
            end
        end
    elseif degree == 7
        T = zeros(nx*ny,64);
        for i=1:ny
            for j=1:nx
                iElem = (i-1)*nx + j;
                inode = (i-1)*7*npx + 7*(j-1) + 1;
                nodes_aux = [inode+(0:7)  inode+npx+(0:7)  inode+2*npx+(0:7) ...
                    inode+3*npx+(0:7)  inode+4*npx+(0:7)  inode+5*npx+(0:7) ...
                    inode+6*npx+(0:7)  inode+7*npx+(0:7)];
                T(iElem,:) = nodes_aux([1  8  64  57  ...
                    2:7   16:8:56   63:-1:58   49:-8:9 ...
                    10:15   18:23   26:31   34:39   42:47   50:55]);
            end
        end
    else
        error('Degree not implemented for quadrilaterals in createMeshRect')
    end
elseif elemType == 2
    if degree == 1
        nx_2 = round(nx/2); ny_2 = round(ny/2);
        T = zeros(2*nx*ny,3);
        for i=1:ny
            for j=1:nx
                iElem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                nodes = [inode   inode+1   inode+npx+1    inode+npx];
                if (i<=ny_2 && j<=nx_2) || (i>ny_2 && j>nx_2)
                    T(iElem,:) = nodes([1,2,3]);
                    T(iElem+1,:) = nodes([1,3,4]);
                else
                    T(iElem,:) = nodes([1,2,4]);
                    T(iElem+1,:) = nodes([2,3,4]);
                end
            end
        end
    elseif degree == 2
        nx_2 = round(nx/2); ny_2 = round(ny/2);
        T = zeros(2*nx*ny,6);
        for i=1:ny
            for j=1:nx
                iElem=2*((i-1)*nx+j)-1;
                inode=(i-1)*2*(npx)+2*(j-1)+1;
                nodes = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
                if (i<=ny_2 && j<=nx_2) || (i>ny_2 && j>nx_2)
                    T(iElem,:)   = nodes([1  9  7  5  8  4]);
                    T(iElem+1,:) = nodes([1  3  9  2  6  5]);
                else
                    T(iElem,:)   = nodes([1  3  7  2  5  4]);
                    T(iElem+1,:) = nodes([3  9  7  6  8  5]);
                end
            end
        end
    elseif degree == 3
        nx_2 = round(nx/2); ny_2 = round(ny/2);
        T = zeros(2*nx*ny,10);
        for i=1:ny
            for j=1:nx
                iElem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*3*npx + 3*(j-1) + 1;
                nodes = [inode+(0:3)   inode+npx+(0:3)   inode+2*npx+(0:3)   inode+3*npx+(0:3)];
                if (i<=ny_2 && j<=nx_2) || (i>ny_2 && j>nx_2)
                    T(iElem,:)=nodes([1  16  13  6  11  15  14  9  5  10]);
                    T(iElem+1,:)=nodes([1  4  16  2  3  8  12  11  6  7]);
                else
                    T(iElem,:)=nodes([1  4  13  2  3  7  10  9  5  6]);
                    T(iElem+1,:)=nodes([4  16  13  8  12  15  14  10  7  11]);
                end
            end
        end
    elseif degree == 4
        nx_2 = round(nx/2); ny_2 = round(ny/2);
        T = zeros(2*nx*ny,15);
        for i=1:ny
            for j=1:nx
                iElem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*4*npx + 4*(j-1) + 1;
                nodes = [inode+(0:4)  inode+npx+(0:4)  ...
                    inode+2*npx+(0:4)  inode+3*npx+(0:4)  inode+4*npx+(0:4)];
                if (i<=ny_2 && j<=nx_2) || (i>ny_2 && j>nx_2)
                    T(iElem,:)   = nodes([1  25  21  7  13  19  24  23  22  16  11  6  12  18  17]);
                    T(iElem+1,:) = nodes([1  5  25  2  3  4  10  15  20  19  13  7  8  9  14]);
                else
                    T(iElem,:)   = nodes([1  5  21  2  3  4  9  13 ...
                        17  16  11  6  7  8  12]);
                    T(iElem+1,:) = nodes([5  25  21  10  15  20  24  23 ...
                        22  17  13  9  14  19  18]);
                end
            end
        end
    elseif degree == 5
        nx_2 = round(nx/2); ny_2 = round(ny/2);
        T = zeros(2*nx*ny,21);
        for i=1:ny
            for j=1:nx
                iElem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*5*npx + 5*(j-1) + 1;
                nodes = [inode+(0:5)  inode+npx+(0:5)  inode+2*npx+(0:5) ...
                    inode+3*npx+(0:5)  inode+4*npx+(0:5)  inode+5*npx+(0:5)];
                if (i<=ny_2 && j<=nx_2) || (i>ny_2 && j>nx_2)
                    T(iElem,:)   = nodes([1  36  31  8  15  22  29  35  34  33  32  25  19  13  7  14  21  28  20  27  26]);
                    T(iElem+1,:) = nodes([1  6  36  2  3  4  5  12  18  24  30  29  22  15  8  9  10  11  16  17  23]);
                else
                    T(iElem,:) = nodes([1  6  31  ...
                        2  3  4  5  11  16  21  26  25  19  13  7  ...
                        8  9  10  14  15  20]);
                    T(iElem+1,:) = nodes([6  36  31 ...
                        12  18  24  30  35  34  33  32  26  21  16  11 ...
                        17  23  29  22  28  27]);
                end
            end
        end
    else
        error('Degree not implemented for triangles in createMeshRect')
    end
else
    error('Element type not implemented in createMeshRect')
end

% If Fekete points are employed, modify coordinates
if approxCG.nodalDist == 0
    Xe_ref = refElem.NodesCoord;
    if elemType == 1 %QUA
        nOfElemVertices = 4;
        N1 = evaluateNodalBasisQua(Xe_ref,Xe_ref(1:nOfElemVertices,:),1);
    elseif elemType == 2 %TRI
        nOfElemVertices = 3;
        N1 = evaluateNodalBasisTri(Xe_ref,Xe_ref(1:nOfElemVertices,:),1);
    end
    for iElem = 1:size(T,1)
        Te = T(iElem,:);
        Xe = X(Te,:);
        Xe_mod = N1*Xe(1:nOfElemVertices,:);
        X(Te,:) = Xe_mod;
    end
end


% Store coordinates and connectivity
mesh.X = X;
mesh.T = T;

% Additional mesh information
mesh.nOfElem = size(T,1);
mesh.nOfNodes = size(mesh.X,1);
mesh.h = sqrt(domSize/mesh.nOfElem); % Estimate of the  mesh size for a uniform mesh

% Identify Dirichlet nodes and solution nodes
if approxCG.isSC == 1
    idxActiveNodes = unique(reshape(mesh.T(:,1:refElem.nOfActiveNodes),mesh.nOfElem*refElem.nOfActiveNodes,1));
    idxCancelledNodes = unique(reshape(mesh.T(:,refElem.nOfActiveNodes+1:refElem.nOfElemNodes),mesh.nOfElem*(refElem.nOfElemNodes-refElem.nOfActiveNodes),1));
    % Coordinates of the active nodes
    mesh.nOfActiveNodes = length(idxActiveNodes);
    mesh.Xact = mesh.X(idxActiveNodes,:);
    mesh.dirNodes = find(abs(mesh.Xact(:,1)-x1)<1e-6 | abs(mesh.Xact(:,1)-x2)<1e-6 | abs(mesh.Xact(:,2)-y1)<1e-6 | abs(mesh.Xact(:,2)-y2)<1e-6); 
    % Connectivity matrix of the active nodes
    Tact = zeros(mesh.nOfElem,refElem.nOfActiveNodes);
    for iElem = 1: mesh.nOfElem
        Tact(iElem,:) = mesh.T(iElem,1:refElem.nOfActiveNodes) - sum(mesh.T(iElem,1:refElem.nOfActiveNodes) > idxCancelledNodes);
    end
    mesh.Tact = Tact;
else
    mesh.nOfActiveNodes = mesh.nOfNodes;
    mesh.dirNodes = find(abs(mesh.X(:,1)-x1)<1e-6 | abs(mesh.X(:,1)-x2)<1e-6 | abs(mesh.X(:,2)-y1)<1e-6 | abs(mesh.X(:,2)-y2)<1e-6); 
    mesh.Tact = mesh.T;
end
%
mesh.solNodes = setdiff([1:mesh.nOfActiveNodes]',mesh.dirNodes);
