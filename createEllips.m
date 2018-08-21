function [posPoints,negPoints] = createEllips ()

a=20/(5*pi())/100; % horizontal radius
b=5/100; % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
cE = @computeEllipse;
x = linspace(-a,a,9);
[y1,y2] = cE(x,a,b);

x = x + 0.1;
y1 = y1 + 0.075;
y2 = y2 + 0.075;

posPoints = [x',y1'];
negPoints = [x',y2'];



end