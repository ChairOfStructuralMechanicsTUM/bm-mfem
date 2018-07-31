function [left, right]=Matrixsplit(Matrix)
%input: nodematrix
%output: vektor: erster und letzter eintrag sind die coners alle anderen
%einträge boundry remainders
%matrix wird in 2 gleich große matrizen geteilt, Teilung findet über median
%statt
[n,m]=size(Matrix);
%n:Zeilenanzahl  m:Spaltenanzahl
left=Matrix(1:n,1:ceil(m/2));
right=Matrix(1:n, ceil(m/2):m);
end