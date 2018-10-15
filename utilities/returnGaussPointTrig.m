function [ w,g ] = returnGaussPointTrig(number, i)
%returnGaussPoint Returns evaluation point g and weighting w depending on
%number of GaussPoints
%   Detailed explanation goes here
if number == 1
    g=[1/3,1/3,1/3];
    w=1;
elseif number == -3 % midpoint rule
    g=[1/2,1/2,1/2];
    g(i)=0;
    w=1/3;
elseif number == 3
    g=[1/6,1/6,1/6];
    g(i)=2/3;
    w=1/3;
elseif number == 7
    if i == 1
        g=[1/3,1/3,1/3];
        w=9/40;
    elseif i <= 4
        g=[(6-sqrt(15))/21,(6-sqrt(15))/21,(6-sqrt(15))/21];
        g(i-1)=(9+2*sqrt(15))/21;
        w=(155-sqrt(15))/1200;
    else
        g=[(6+sqrt(15))/21,(6+sqrt(15))/21,(6+sqrt(15))/21];
        g(i-4)=(9-2*sqrt(15))/21;
        w=(31/120)-(155-sqrt(15))/1200;
    end
else
    error('Number of Gauss Points not possible for Triangle Element or not yet implemented.')
end

end