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
if N>size(nodematrix,2)-1
    fprintf ('warning: too many subdomains!')
    return
    
    elseif N==1
        fprintf ('Trivialfall')
        return

    elseif N==2  
        [~,~,gbc,gbr]=Matrixsplit(nodematrix);

    elseif N==3
        [left,~,bc,br]=Matrixsplit(nodematrix);
        [~,~,bc2,br2]=Matrixsplit(left);
        gbc=[bc2;bc];
        gbr=[br2;br];
else 
[left,right,gbc,gbr]=Matrixsplit(nodematrix);%Initialisierung, left und right bleiben als Ausgangsbasis erhalten
[left2,right2,gbc,gbr]=Matrixsplit(left);
%left2 und right 2 sind die zweite ebene der Ausgangsbasis die mit left und right verglichen werden...
%als erster Schritt wird left nochmals geteilt da aufgrund der
%Implementierung von Matrixsplit leftimmer >=right ist
l=1; %Schleifenzähler 
b=3;%Zählervariable für gbc
r=dim(1)-1; %Zählervariable für gbr
while l<N-2 %2 Schritte bei Initialisierung
    if size(left2,2)<= size(right,2)
        %jump to other branch
        [left2,right2,bc,br]=Matrixsplit(right);
        left=left2;
    else 
        [left2,right2,bc,br]=Matrixsplit(left2);
    end
    gbc(b:b+1,1)=bc;
    gbr(r:r+dim(1)-3,1)=br;
    l=l+1;
    b=b+2;
    r=r+dim(1)-2;
end
end
    
%% merge gbc and gbr to global ufinal vector

%ufinal=sorted node vector: [i;br;bc]
%gi: vector of internal nodes
u1=[gbr;gbc]; %combined array, helps to find internal nodes
ufinal=zeros(size(nodearray,1)-size(u1,1),1);

v=1;
for i=1:size(nodearray,1)
        if any(u1==nodearray(i))==false
            ufinal(v,1)=nodearray(i,1);
        else
            v=v-1;
          %%fullfill: i zält weiter, aber erst eintrag weiter hinten kommt
          %%an die stelle, id marker mitlaufen lassen
        end
        v=v+1;
end

ufinal=[ufinal;gbr;gbc];            

%% multiply Kmatrix and ufinal(converted to doffinal) to sort Kmatrix



end


