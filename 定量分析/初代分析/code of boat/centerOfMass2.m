function COM = centerOfMass2(masses,mesh)
    % centerOfMass2: computes center of mass in 2D
    % masses: matrix of masses
    % mesh: structure containing ygrid and zgrid
    % returns: Vector [ycom,zcom]
    M = matrixSum(masses);
    ycom = matrixSum(masses .* mesh.ygrid) / M;
    zcom = matrixSum(masses .* mesh.zgrid) / M;
    COM = [ycom,zcom];
end

