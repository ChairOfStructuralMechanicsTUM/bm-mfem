function [ len ] = computeLength( coords1, coords2 )
%COMPUTELENGTH Computes the length between two 2-d or 3-d points
%   Detailed explanation goes here
if (length(coords1) == 2)
    len = sqrt((coords2(1)-coords1(1))^2 ...
        + (coords2(2)-coords1(2))^2);
else
    len = sqrt((coords2(1)-coords1(1))^2 ...
        + (coords2(2)-coords1(2))^2 ...
        + (coords2(3)-coords1(3))^2);
end

end

