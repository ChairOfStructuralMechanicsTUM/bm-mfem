function [y1,y2] = computeEllipse (x,a,b)
y1 = sqrt(b^2*(1-x.^2/a^2));
y2 = -y1;
end
