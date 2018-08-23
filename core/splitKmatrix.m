function [Ksort,Krr,Kcc,Krc,Kcr]=splitKmatrix (Kmatrix,n)

%Kmatrix= globale Systemsteifigkeitsmatrix, lokale Elemente bereits
%assembliert
%System: Fachwerk mit Grundelement 1; 4 Elemente 4 Knoten 8 globalen Dofs
%enspricht linker oberer Ecke einer minimalen Substruktur
%% setup 
%Kmatrix
% Kmatrix=string(ones(n));
% for i=1:n
%     for j=1:n
%         id=[i,j];
%         Kmatrix(i,j)=['K' num2str(id(1)) num2str(id(2))];
%     end
% end
%bc: Freiheitsgrade der Eckknoten
bc=[1;2];
%br
br=[3;4];
%i
i=[5;6;7;8];
%r
r=[i;br];
%Vektor der Verschiebungsfreiheitsgrade
u=[i;br;bc];
%Umsortierte Stiefigkeitsmatrix
Ksort=string(zeros(n,n));
for i=1:n
    for j=1:n
        Ksort(i,j)=Kmatrix(u(i),u(j));
    end
end
%Steifigkeitsmatrix der remainder (br und i): Krr
%Krr=string(zeros(size(r,1)));
Krr=Ksort(1:size(r,1),1:size(r,1));
%Steifigkeitsmatrix der corner Freiheitsgrade (bc):Kcc
%Kcc=string(zeros(size(bc,1)));
Kcc=Ksort(size(r,1)+1:size(r,1)+size(bc,1),size(r,1)+1:size(r,1)+size(bc,1));
%Steifigkeitsmatrizen der Kombinierten Freiheitsgrade rbc, bcr: Krc, Kcr
Krc=Ksort(1:size(r,1),size(r,1)+1:size(Ksort,1));
Kcr=Ksort(size(r,1)+1:size(Ksort,1),1:size(r,1));
end