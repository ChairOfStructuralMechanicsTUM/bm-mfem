function [ check ] = isOnLineBetweenTwoPoints( c1, c2, c3 )
%ISONLINEBETWEENTWOPOINTS Checks, if C3 is on a line between C1 and C2
%   C = isOnLineBetweenTwoPoints(C1,C2,C3)
check = false;

if all(abs(cross(c1-c2,c1-c3)) < 1e-8)
    d = c2 - c1;
    [~, ind] = max(abs(d));
    if d(ind) > 0
        check = c1(ind) <= c3(ind) && c3(ind) <= c2(ind);
    else
        check = c2(ind) <= c3(ind) && c3(ind) <= c1(ind);
    end
end


end

