function [ matrix ] = applyMatrixBoundaryConditions( matrix, bc )
%APPLYMATRIXBOUNDARYCONDITIONS Removes bc lines and columns from a matrix
%   Detailed explanation goes here
    matrix(bc, :) = [];
    matrix(:, bc) = [];

end

