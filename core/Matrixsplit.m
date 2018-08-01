function [left, right,bc,br]=Matrixsplit(Matrix)
%input: nodematrix
%output: vektor: erster und letzter eintrag sind die coners alle anderen
%eintr‰ge boundry remainders
%matrix wird in 2 gleich groﬂe matrizen geteilt, Teilung findetbei der
%h‰lfte statt
[n,m]=size(Matrix);
%n:Zeilenanzahl  m:Spaltenanzahl
left=Matrix(1:n,1:ceil(m/2));
right=Matrix(1:n, ceil(m/2):m);
bc=[Matrix(1,ceil(m/2));Matrix(n,ceil(m/2))];
br=Matrix(2:n-1,ceil(m/2));
end