function [left, right,bc,br]=Matrixsplit(Matrix)
%input: nodematrix
%output: vektor: erster und letzter eintrag sind die coners alle anderen
%eintr‰ge boundry remainders
%matrix wird in 2 gleich groﬂe matrizen geteilt, Teilung findetbei der
%h‰lfte statt
[n,m]=size(Matrix);
%n:Zeilenanzahl  m:Spaltenanzahl
left=Matrix(1:n,1:floor(m/2)+1);
right=Matrix(1:n, floor(m/2)+1:m);
bc=[Matrix(1,floor(m/2)+1);Matrix(n,floor(m/2)+1)];
br=Matrix(2:n-1,floor(m/2)+1);
end