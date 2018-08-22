function [Kmatrix,Ksort]=splitKmatrix (n)

%Kmatrix= globale Systemsteifigkeitsmatrix, lokale Elemente bereits
%assembliert
%System: Fachwerk mit Grundelement 1; 4 Elemente 4 Knoten 8 globalen Dofs
%enspricht linker oberer Ecke einer minimalen Substruktur
%% setup 
%Kmatrix
Kmatrix=string(ones(n));
for i=1:n
    for j=1:n
        id=[i,j];
        Kmatrix(i,j)=['K' num2str(id(1)) num2str(id(2))];
    end
end
%bc: Freiheitsgrade der Eckknoten
bc=[1;2];
%br
br=[3;4];
%i
i=[5;6;7;8];
%Vektor der Verschiebungsfreiheitsgrade
u=[i;br;bc];
Ksort=string(zeros(n,n));
k=1;
for i=1:n
    for j=1:n
        Ksort(i,j)=Kmatrix(u(i),u(j));
    end
end
end