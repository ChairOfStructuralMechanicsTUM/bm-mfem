function [nodematrix,K]=substructure1(Ns,hz,v,nodearray,dim)
clc
%Ns: Anzahl der gew�nschten substructures
%hz: Anzahl der Substrukturen in horizontale richtung
%v: Anzahl der Substrukturen in vertikale Richtung
%hz und v sollten gleich sein oder zumindest so nah wie m�glich beieinander
%liegen um eine sinnvolle Substrukturierung zu erhalten also z.B. 10x10 und
%nicht 50x2
%dim(~,~)= dimension der Knotenmatrix= Zeilen und Spaltenanzahl
if hz*v~=Ns
    fprintf('unzul�ssige Substrukturierung, bitte Anzahl und Unterteilung der Substrukturen �berpr�fen')
   return
else
%% nodearray-->nodematrix

nodematrix=zeros(dim);
k=1;
for i=1:dim(2)
    for j=1:dim(1)
        nodematrix(j,i)=nodearray(k);
        k=k+1;
    end
end

%% substructuring of nodematrix

K=cell(v,hz); %cell um die verschiedenen Substructures als arrays darin zu speichern
a=size(nodematrix,1)/v; %Anzahl Knoten einer Spalte einer Substtruktur, Zeilenanzahl
b=size(nodematrix,2)/hz; %%Anzahl Knoten einer Zeile einer Substtruktur, Spaltenanzahl
if round(a)==a
    a=a+1
else
    a=round(a)
end
if round(b)==b
    b=b+1
else
    b=round(b)
end    
%testen auf minimalsubstruktur: 3x3 Fachwerk
if or(a<3,b<3)
    fprintf('Substrukturierung erzeugt zu kleine Substrukturen, bitte kleinere Anzahl an Substrukturen w�hlen!');
    return
else
    
for i=1:hz
    for j=1:v
        if i==hz && j~=v
            K{j,i}={nodematrix((j-1)*a+1:j*(a-1),(i-1)*b:i*(b-1))};
        elseif j==v && i~=hz
            K{j,i}={nodematrix((j-1)*a:j*(a-1),(i-1)*b+1:i*(b-1))};
        elseif j==v && i==hz
            K{j,i}={nodematrix((j-1)*a:j*(a-1),(i-1)*b:i*(b-1))};
        else
            K{j,i}={nodematrix((j-1)*a+1:j*a,(i-1)*b+1:i*b)};
        end
    end
end
% for i=(hz-1):hz
%     for j=(v-1):v
%     K{j,i}={nodematrix(j*a:dim(1),i*b:dim(2))};
%     end
% end


%% substructure first aproach
%% 
%dim= dimension des Fachwerks, bezieht sich auf Knoten!!! z.B. 3x5, muss gleiche Anzahl Knoten
%wie nodeaaray haben

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

%% nodematrix split in N substructures, gbc and gbr filled

%gbc: globaler Vektor der corner nodes  bc:lokaler Vektor der corner nodes
%einer subdomain in einer Iteration
%gbr: globaler Vektor der corner remainders  br: lokaler Vektor der boundry
%remainders in einer iteration
% if N>size(nodematrix,2)-1
%     fprintf ('warning: too many subdomains!')
%     return
%     
%     elseif N==1
%         fprintf ('Trivialfall')
%         return
% 
%     elseif N==2  
%         [~,~,gbc,gbr]=Matrixsplit(nodematrix);
% 
%     elseif N==3
%         [left,~,bc,br]=Matrixsplit(nodematrix);
%         [~,~,bc2,br2]=Matrixsplit(left);
%         gbc=[bc2;bc];
%         gbr=[br2;br];
% else 
% [left,right,bc,br]=Matrixsplit(nodematrix);%Initialisierung, left und right bleiben als Ausgangsbasis erhalten
% [left2,right2,bc2,br2]=Matrixsplit(left);
% gbc=[bc2;bc];
% gbr=[br2;br];
%left2 und right 2 sind die zweite ebene der Ausgangsbasis die mit left und right verglichen werden...
%als erster Schritt wird left nochmals geteilt da aufgrund der
%Implementierung von Matrixsplit leftimmer >=right ist
% l=1; %Schleifenz�hler 
% b=5;%Z�hlervariable f�r gbc
% r=size(gbr,1); %Z�hlervariable f�r gbr
% while l<N-2 %2 Schritte bei Initialisierung
%     if size(left2,2)<= size(right,2)
%         %jump to other branch
%         %left=left2;
%         [right,right2,bc,br]=Matrixsplit(right);
%         gbc(b:b+1,1)=bc;
%         gbr(r+1:r+dim(1)-2,1)=br;
%         b=b+2;
%         r=r+dim(1)-2;
%     else 
%         [left2,right2,bc,br]=Matrixsplit(left2);
%         gbc(b:b+1,1)=bc;
%         gbr(r+1:r+dim(1)-2,1)=br;
%         b=b+2;
%         r=r+dim(1)-2;
%     end
%     %gbc(b:b+1,1)=bc;
%     %gbr(r:r+dim(1)-3,1)=br;
%     l=l+1;
%     %b=b+2;
%     %r=r+dim(1)-2;
%     gbc=sort(gbc);
%     gbr=sort(gbr);
% end
% end
%     
% %% merge gbc and gbr to global ufinal vector
% 
% %ufinal=sorted node vector: [i;br;bc]
% %gi: vector of internal nodes
% u1=[gbr;gbc]; %combined array, helps to find internal nodes
% ufinal=zeros(size(nodearray,1)-size(u1,1),1);
% 
% v=1;
% for i=1:size(nodearray,1)
%         if any(u1==nodearray(i))==false
%             ufinal(v,1)=nodearray(i,1);
%         else
%             v=v-1;
%           %%fullfill: i z�lt weiter, aber erst eintrag weiter hinten kommt
%           %%an die stelle, id marker mitlaufen lassen
%         end
%         v=v+1;
% end
% 
% ufinal=[ufinal;gbr;gbc];            
% 
% %% multiply Kmatrix and ufinal(converted to doffinal) to sort Kmatrix
% %boolean Matrizen B aufstellen um K zu sortieren
% %B1: bezieht sich auf gbc
% Kfinal=Kmatrix(:,[ufinal.']);
% 
end
end
end

