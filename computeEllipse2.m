function [y1,y2] = computeEllipse2(x)

a=20/(5*pi())/100; % horizontal radius
b=5/100; % vertical radius

x = x-0.1;
y1 = sqrt(b^2*(1-x^2/a^2));
y2 = -y1;
y1 = y1+ 0.075;
y2 = y2+ 0.075;

end