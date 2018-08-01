function [gbc,gbr]=substructure1(N,dim,nodearray)
%N: Anzahl der gew�nschten substructures
%dim= dimension des Fachwerks z.B. 3x5, muss gleiche Anzahl Knoten wie
%nodeaaray habe!
%Test-steifigkeitsmatrix wird eingef�hrt um substrrukturierung testen zusubstructure1 (N, nodearray, dofarray, Kmatrix, dim)
%k�nnen
%dof: [D1 D2 D3]^T werden an K ranmultipliziert: erste Spalte von K bezieht
%sich auf D1... 
%Benennung von K: Kik: Festhaltekraft am Freiheitsgrad(dof)
%i infolge Einheitsverschiebung Dk=1


%Kmatrix= globale Systemsteifigkeitsmatrix, lokale Elemente bereits
%assembliert
%System: Fachwerk mit Grundelement 1; 12 Elemente 9 Knoten 18 globalen Dofs
%%
%Kmatrix=string(ones(18));
% for i=1:18
%     for j=1:18
%         id=[i,j]
%         Kmatrix(i,j)=['K' num2str(id(1)) num2str(id(2))]
%     end
% end

%% nodearray-->nodematrix
nodematrix=zeros(dim);
k=1;
for i=1:dim(2)
    for j=1:dim(1)
        nodematrix(j,i)=nodearray(k);
        k=k+1;
    end
end
%% nodematrix split in N substructures
%gbc: globaler Vektor der corner nodes  bc:lokaler Vektor der corner nodes
%einer subdomain in einer Iteration
%gbr: globaler Vektor der corner remainders  br: lokaler Vektor der boundry
%remainders in einer iteration

[left,right,gbc,gbr]=Matrixsplit(nodematrix);

l=1; %Schleifenz�hler
b=3;%Z�hlervariable f�r gbc
r=dim(1)-1; %Z�hlervariable f�r gbr
while l<N-1
    if size(left,2)<size(right,2)
        [left,right,bc,br]=Matrixsplit(right)
       
    else
        [left,right,bc,br]=Matrixsplit(left)
    end
    gbc(b:b+1,1)=bc
    gbr(r:r+dim(1)-3,1)=br
    l=l+1
    b=b+2
    r=r+dim(1)-2
end
end


