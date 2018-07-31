function [Kmatrix]=substructuretest1 (N)
%Test-steifigkeitsmatrix wird eingeführt um substrrukturierung testen zu
%können
%dof: [D1 D2 D3]^T werden an K ranmultipliziert: erste Spalte von K bezieht
%sich auf D1... 
%Benennung von K: Kik: Festhaltekraft am Freiheitsgrad(dof)
%i infolge Einheitsverschiebung Dk=1


%Kmatrix= globale Systemsteifigkeitsmatrix, lokale Elemente bereits
%assembliert
%System: Fachwerk mit Grundelement 1; 7 Elemente 6 Knoten 12 globalen Dofs
Kmatrix=string(ones(12));
%id=string();
for i=1:12
    for j=1:12
        id=[i,j]
        Kmatrix(i,j)=['K' num2str(id(1)) num2str(id(2))]
    end
end


end