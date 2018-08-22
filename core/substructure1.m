function [nodematrix,K,bc,br,i]=substructure1(Ns,v,hz,nodearray,dim)
clc
%Ns: Anzahl der gewünschten substructures
%hz: Anzahl der Substrukturen in horizontale richtung
%v: Anzahl der Substrukturen in vertikale Richtung
%hz und v sollten gleich sein oder zumindest so nah wie möglich beieinander
%liegen um eine sinnvolle Substrukturierung zu erhalten also z.B. 10x10 und
%nicht 50x2
%dim(~,~)= dimension der Knotenmatrix= Zeilen und Spaltenanzahl

%% Abbruchkriterien

if hz*v~=Ns
    fprintf('unzulässige Substrukturierung, bitte Anzahl und Unterteilung der Substrukturen überprüfen')
   return
elseif Ns==1
    fprintf('keine Substrukturierung ausgewählt bitte Ns ungleich 1 wählen')
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
a=floor(size(nodematrix,1)/v); %Anzahl Knoten einer Spalte einer Substtruktur, Zeilenanzahl
b=floor(size(nodematrix,2)/hz); %%Anzahl Knoten einer Zeile einer Substtruktur, Spaltenanzahl
bc=cell(v,hz);
br=cell(v,hz);
in=cell(v,hz);

%% Testen auf minimalsubstruktur: 3x3 Fachwerk mit 9 Knoten

if or(a<2,b<2)
    fprintf('Substrukturierung erzeugt zu kleine Substrukturen, bitte kleinere Anzahl an Substrukturen wählen!');
    return
else
%% Unterteilung der nodematrix in Substrukturen Kij 
% Sortierung der Vektoren bc, br, in

for i=1:hz
    for j=1:v
        if i==1 && j==1 %Fall 1: linke obere Ecke
            if v==1 %falls nur eine Substruktur in vertikale Richtung vorhanden ist
                bc(j,i)={[nodematrix(1,b+1);nodematrix(a,b+1)]};
                br(j,i)={nodematrix(2:a-1,b+1)};
                K(j,i)={nodematrix((j-1)*a+1:j*a,(i-1)*b+1:i*b+1)};
            else
                bc(j,i)={[nodematrix(a+1,1);nodematrix(1,b+1);nodematrix(a+1,b+1)]};
                br(j,i)={[nodematrix(a+1,2:b);nodematrix(2:a,b+1)]};
                K(j,i)={nodematrix((j-1)*a+1:j*a+1,(i-1)*b+1:i*b+1)};
            end
        elseif i==1 && j~=1 && j~=v %Fall 2: erste Spalte Mitte
            bc(j,i)={[nodematrix((j-1)*a+1,1);nodematrix((j-1)*a+1,b+1);nodematrix(j*a+1,1);nodematrix(j*a+1,b+1)]};
            br(j,i)={[nodematrix((j-1)*a+1,2:b).';nodematrix(j*a+1,2:b).';nodematrix((j-1)*a+2:j*a,b+1)]};
            K(j,i)={nodematrix((j-1)*a+1:j*a+1,(i-1)*b+1:i*b+1)};
        elseif i==1 && j==v %Fall 3 linke untere Ecke
            bc(j,i)={[nodematrix(dim(1)-(a),1);nodematrix(dim(1)-a,b+1);nodematrix(dim(1),b+1)]};
            br(j,i)={[nodematrix(dim(1)-a,2:b);nodematrix(dim(1)-a+1:dim(1)-1,b+1)]};
            K(j,i)={nodematrix((j-1)*a+1:dim(1),(i-1)*b+1:i*b+1)};  
        elseif j==1 && i~=1 && i~= hz    %Fall 4 erste Zeile Mitte
            if v==1 %falls nur eine Substruktur in vertikale Richtung vorhanden ist
                bc(j,i)={[nodematrix(1,(i-1)*b+1);nodematrix(a,(i-1)*b+1);nodematrix(1,i*b+1);nodematrix(a,i*b+1)]};
                br(j,i)={[nodematrix(2:j*a-1,(i-1)*b+1);nodematrix(2:j*a-1,i*b+1)]};
                K(j,i)={nodematrix((j-1)*a+1:j*a,(i-1)*b+1:i*b+1)};
            else
                bc(j,i)={[nodematrix((j-1)*a+1,(i-1)*b+1);nodematrix(j*a+1,(i-1)*b+1);nodematrix((j-1)*a+1,i*b+1);nodematrix(j*a+1,i*b+1)]};
                br(j,i)={[nodematrix(2:j*a,(i-1)*b+1);nodematrix(a+1,(i-1)*b+2:i*b).';nodematrix(2:j*a,i*b+1)]};
                K(j,i)={nodematrix((j-1)*a+1:j*a+1,(i-1)*b+1:i*b+1)};
            end
        elseif i>1 && i<hz && j>1 && j<v %Fall 5 Mitte Mitte
            bc(j,i)={[nodematrix((j-1)*a+1,(i-1)*b+1);nodematrix(j*a+1,(i-1)*b+1);nodematrix((j-1)*a+1,i*b+1);nodematrix(j*a+1,i*b+1)]};
            br(j,i)={[nodematrix((j-1)*a+2:j*a,(i-1)*b+1);nodematrix((j-1)*a+1,(i-1)*b+2:i*b).';nodematrix(j*a+1,(i-1)*b+2:i*b).';nodematrix((j-1)*a+2:j*a,i*b+1)]};
            K(j,i)={nodematrix((j-1)*a+1:j*a+1,(i-1)*b+1:i*b+1)};
        elseif j==v && i~=1 && i~= hz %Fall 6 unterste Zeile Mitte
            bc(j,i)={[nodematrix((j-1)*a+1,(i-1)*b+1);nodematrix(dim(1),(i-1)*b+1);nodematrix((j-1)*a+1,i*b+1);nodematrix(dim(1),i*b+1)]};
            br(j,i)={[nodematrix((j-1)*a+2:dim(1)-1,(i-1)*b+1);nodematrix((j-1)*a+1,(i-1)*b+2:i*b).';nodematrix((j-1)*a+2:dim(1)-1,i*b+1)]};
            K(j,i)={nodematrix((j-1)*a+1:dim(1),(i-1)*b+1:i*b+1)};
        elseif j==1 && i==hz %Fall 7 rechte obere Ecke
            if v==1 %falls nur eine Substruktur in vertikale Richtung vorhanden ist
                bc(j,i)={[nodematrix(1,dim(2)-b);nodematrix(a,dim(2)-b)]};
                br(j,i)={[nodematrix(2:a-1,dim(2)-b)]};
                K(j,i)={nodematrix((j-1)*a+1:j*a,(i-1)*b+1:dim(2))};
            else
                bc(j,i)={[nodematrix(1,dim(2)-b);nodematrix(a,dim(2)-b);nodematrix(a,dim(2))]};
                br(j,i)={[nodematrix(2:a,dim(2)-b);nodematrix(a+1,dim(2)-b+1:dim(2)-1).']};
                K(j,i)={nodematrix((j-1)*a+1:j*a+1,(i-1)*b+1:dim(2))};
            end
        elseif i==hz && j~=1 && j~=v %Fall 8 letzte Spalte Mitte
            bc(j,i)={[nodematrix((j-1)*a+1,(i-1)*b+1);nodematrix(j*a+1,(i-1)*b+1);nodematrix((j-1)*a+1,dim(2));nodematrix(j*a+1,dim(2))]};
            
            K(j,i)={nodematrix((j-1)*a+1:j*a+1,(i-1)*b+1:dim(2))};
        else %Fall 9 rechte untere Ecke
            bc(j,i)={[nodematrix((j-1)*a+1,(i-1)*b+1);nodematrix(dim(1),(i-1)*b+1);nodematrix((j-1)*a+1,dim(2))]};
            K(j,i)={nodematrix((j-1)*a+1:dim(1),(i-1)*b+1:dim(2))};
        end
    end
end

%% Identifizierung der Eckknoten bc, Interface Knoten br und internen Knoten i

%K1=cell2mat(K{1,1});



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
% l=1; %Schleifenzähler 
% b=5;%Zählervariable für gbc
% r=size(gbr,1); %Zählervariable für gbr
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
%           %%fullfill: i zält weiter, aber erst eintrag weiter hinten kommt
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


