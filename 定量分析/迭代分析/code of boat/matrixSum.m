function M = matrixSum(masses)
    % matrixSum: returns total of all elements in the matrix
    
    % normally sum(m) computes the sums of the columns
    % selecting m(:) flattens the matrix and computes the sum of all elements
    % see https://stackoverflow.com/questions/1721987/what-are-the-ways-to-sum-matrix-elements-in-matlab
    M = sum(masses(:));
end

