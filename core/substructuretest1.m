function [Kmatrix]=substructuretest1 (N)
%Test-steifigkeitsmatrix wird eingeführt um substrrukturierung testen zu
%können
%dof: [D1 D2 D3]^T werden an K ranmultipliziert: erste Spalte von K bezieht
%sich auf D1... 
%Benennung von K: Kik: Festhaltekraft am Freiheitsgrad(dof)
%i infolge Einheitsverschiebung Dk=1

Kmatrix=[K11 K12  K13;
         K21 K22  K23
         K31 K32  K33];

end