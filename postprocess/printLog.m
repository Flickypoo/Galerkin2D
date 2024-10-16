function printLog(mesh,approxCG,err)

disp('============= Log =============');
%
fprintf(' >Number of mesh elements: %d\n', mesh.nOfElem);
%
switch approxCG.elemType
    case 1
        disp(' >Element type: Quadrilaterals')
    case 2
        disp(' >Element type: Triangles')
end
%
switch approxCG.nodalDist
    case 0
        disp(' >Nodal distribution: Fekete points')
    case 1
        disp(' >Nodal distribution: Uniform points')
end

%% Print mesh information
fprintf('\n >Characteristic mesh size: %5.2e\n', mesh.h);
fprintf(' >Degree of approximation: %d\n\n', approxCG.degree);

%% Static condensation
switch approxCG.isSC
    case 0
        disp(' >Static condensation: OFF')
    case 1
        disp(' >Static condensation: ON')
end

%% Print error information
fprintf('\n >Relative L2 error = %5.2e\n', err.L2);
fprintf(' >Relative H1 error = %5.2e\n', err.H1);
