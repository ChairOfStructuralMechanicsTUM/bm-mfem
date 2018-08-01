function [gbc,gbr,ufinal]=substructure1(N,dim,nodearray)
%N: Anzahl der gewünschten substructures
%dim= dimension des Fachwerks z.B. 3x5, muss gleiche Anzahl Knoten wie
%nodeaaray habe!
%Test-steifigkeitsmatrix wird eingeführt um substrrukturierung testen zusubstructure1 (N, nodearray, dofarray, Kmatrix, dim)
%können
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
%% nodematrix split in N substructures, gbc and gbr filled
%gbc: globaler Vektor der corner nodes  bc:lokaler Vektor der corner nodes
%einer subdomain in einer Iteration
%gbr: globaler Vektor der corner remainders  br: lokaler Vektor der boundry
%remainders in einer iteration

[left,right,gbc,gbr]=Matrixsplit(nodematrix);

l=1; %Schleifenzähler
b=3;%Zählervariable für gbc
r=dim(1)-1; %Zählervariable für gbr
while l<N-1
    if size(left,2)<size(right,2)
        [left,right,bc,br]=Matrixsplit(right);
       
    else
        [left,right,bc,br]=Matrixsplit(left);
    end
    gbc(b:b+1,1)=bc;
    gbr(r:r+dim(1)-3,1)=br;
    l=l+1;
    b=b+2;
    r=r+dim(1)-2;
end
%% merge gbc and gbr to global u vector
%gi: vector of internal nodes
u1=[gbr;gbc]; %combined array, helps to find internal nodes
ufinal=zeros(size(nodearray,1),1);
%ufinal=sorted node vector: i;br;bc
for i=1:size(nodearray,1)
        if any(u1==nodearray(i))==false
            ufinal(i,1)=nodearray(i,1);
        else
          %%fullfill: i zält weiter, aber erst eintrag weiter hinten kommt
          %%an die stelle, id marker mitlaufen lassen
        end
end

ufinal=[ufinal;gbr;gbc];            



end


