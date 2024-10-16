function refElem = setReferenceElement(approxCG)

elemType = approxCG.elemType;
degree = approxCG.degree;
nodalDist = approxCG.nodalDist;

if degree > 5
     error('Error: polynomial degree not available for simulation')
end


if elemType == 1
    if nodalDist == 1
        coord1d = linspace(-1,1,degree+1)';
    else
        coord1d = feketeNodes1D(degree,1:degree+1);
    end
	[x,y] = meshgrid(coord1d);
    coord2d = [reshape(x',(degree+1)^2,1)  reshape(y',(degree+1)^2,1)];
    nodes = 1:(degree+1)^2;
    nodes = reshape(nodes,degree+1,degree+1)';
    nodes(1,:)=[]; nodes(end,:)=[]; nodes(:,1)=[]; nodes(:,end)=[];
    ordre_int = reshape(nodes',1,(degree-1)^2);
    ordre = [1,degree+1,(degree+1)^2,degree*(degree+1)+1,...
        2:degree, 2*(degree+1):degree+1:degree*(degree+1), ...
        (degree+1)^2-1:-1:degree*(degree+1)+2, ...
        (degree-1)*(degree+1)+1:-(degree+1):degree+2, ordre_int];
    coord2d = coord2d(ordre,:);
    
    

    nOfGaussPoints1d = 2*degree+1; %number of integration points 1D
    [gp1d,gw1d] = gaussLegendre(nOfGaussPoints1d);
    nOfGaussPoints = nOfGaussPoints1d^2;
    gp = zeros(nOfGaussPoints,2);
    gw = zeros(1,nOfGaussPoints);
    %position and weights of gauss points
    iGauss = 1;
    for i = 1:nOfGaussPoints1d
        for j = 1:nOfGaussPoints1d
            gp(iGauss,:) = [gp1d(i),gp1d(j)];
            gw(iGauss) = gw1d(i)*gw1d(j);
            iGauss = iGauss + 1;
        end
    end    
    
    [N,Nxi,Neta]=evaluateNodalBasisQua(gp,coord2d,degree);

elseif elemType == 2
    % Nodal coordinates in the reference element
    if degree == 1
        coord2d = [-1 -1; 1 -1; -1 1];
    elseif degree == 2
        coord2d = [-1 -1; 1 -1; -1 1; 0 -1; 0 0; -1 0];
    else
        if nodalDist == 1
            uniformFile = load('positionUniformNodesTri2d.mat'); 
            coord2d = uniformFile.uniformNodesPosition.(['P' num2str(degree)]);
        else
            feketeFile = load('positionFeketeNodesTri2D_EZ4U.mat'); % EZ4U reference element
            coord2d = feketeFile.feketeNodesPosition.(['P' num2str(degree)]);
        end
    end
    
    % Quadrature (integration points and weigths)
    switch degree
        case 1
            OrderCubature = 5;
        case 2
            OrderCubature = 10;
        case 3
            OrderCubature = 10;
        case 4
            OrderCubature = 15;
        case 5
            OrderCubature = 15;
        otherwise
            OrderCubature = 25; 
    end
    [z,w] = GaussLegendreCubature2D(OrderCubature);
    gw = 2*w'; gp = 2*z -1; % mapping onto the normal reference triangle
    
    % Shape functions evaluated at the integration points
    [N,Nxi,Neta]=evaluateNodalBasisTri(gp,coord2d,degree);
end

nIP = length(gw); 
nOfElemNodes = size(coord2d,1);

%Creating reference element structure
refElem = struct('type',elemType,'degree',degree,...
    'NodalDistrib', nodalDist, ...
    'nOfElemNodes',nOfElemNodes,'NodesCoord',coord2d, ...
    'nIP', nIP, 'xIP', gp, 'wIP', gw,...
    'N',N, 'Nxi', Nxi, 'Neta',Neta);

%
if approxCG.isSC == 1
    % TO BE COMPLETED: specify the number of nodes interior to the element
    % according to its type (quadrilateral/triangle) and its degree of 
    % approximation. The vector nOfIntNodes has dimension 2x5, specifically:
    %  - Row 1 is associated with quadrilateral elements.
    %  - Row 2 is associated with triangular elements.
    %  - Each column (from 1 to 5) is associated with a degree of approximation.
    %
    %nOfIntNodes = ;

    refElem.nOfCondensedNodes = nOfIntNodes(elemType,degree);
else
    refElem.nOfCondensedNodes = 0;
end
%
refElem.nOfActiveNodes = nOfElemNodes - refElem.nOfCondensedNodes;