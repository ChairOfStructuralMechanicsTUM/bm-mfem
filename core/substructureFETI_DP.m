classdef substructureFETI_DP < handle
    
   %% properties
   properties (Access=private)
       nodematrix=Node.empty  %Matrix der Knoten ids
       nodearray=Node.empty
       K        % cell array, speichert die Knoten ids jeder Substruktur
       bc
       br
       in
       gbc
       gbr
       gin
   end
   
   
   
   %% test and constructor
   methods
       
       function substructuring= substructureFETI_DP(femModel)
             if (nargin > 0)
                
            else
                error('input model is missing');
            end
       end
       
   end
   
   
  %% methods
   methods (Static)
       
       %% nodematrix
       function [nodematrix]=setupNodeMatrix(femModel,dim)
       nodearray=femModel.getAllNodes;
       nodeIdArray=zeros(length(nodearray));
       for itnod = 1:length(nodearray)
           nodeIdArray(itnod)=nodearray(itnod).getId;
       end
       
       nodematrix=zeros(dim);
       k=1;
       for i=1:dim(2)
           for j=1:dim(1)
           nodematrix(j,i)=nodeIdArray(k);
           k=k+1;
           end
       end
       end
       
       
       %% Unterteilung der nodematrix in Substrukturen Kij
       % Sortierung der Vektoren bc, br, in
       
       function [K,bc,br,in,gbc,gbr,gin]= substructureNodeMatrix(femModel,nodematrix,Ns,v,hz,dim)
           
           if hz*v~=Ns
                fprintf('unzul�ssige Substrukturierung, bitte Anzahl und Unterteilung der Substrukturen �berpr�fen')
                return
           elseif Ns==1
                fprintf('keine Substrukturierung ausgew�hlt bitte Ns ungleich 1 w�hlen')
                return
           else
               
            K=cell(v,hz); %cell um die verschiedenen Substructures als arrays darin zu speichern
            a=floor(size(nodematrix,1)/v); %Anzahl Knoten einer Spalte einer Substtruktur, Zeilenanzahl
            b=floor(size(nodematrix,2)/hz); %%Anzahl Knoten einer Zeile einer Substtruktur, Spaltenanzahl
            bc=cell(v,hz);
            br=cell(v,hz);
            in=cell(v,hz);
            gbc=[];
            gbr=[];
            gin=[];
            
            if or(a<2,b<2)
                fprintf('Substrukturierung erzeugt zu kleine Substrukturen, bitte kleinere Anzahl an Substrukturen w�hlen!');
            return
            else
            n=1;
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
                                br(j,i)={nodematrix(2:a-1,dim(2)-b)};
                                K(j,i)={nodematrix((j-1)*a+1:j*a,(i-1)*b+1:dim(2))};
                            else
                                bc(j,i)={[nodematrix(1,dim(2)-b);nodematrix(a+1,dim(2)-b);nodematrix(a+1,dim(2))]};
                                br(j,i)={[nodematrix(2:a,dim(2)-b);nodematrix(a+1,dim(2)-b+1:dim(2)-1).']};
                                K(j,i)={nodematrix((j-1)*a+1:j*a+1,(i-1)*b+1:dim(2))};
                            end
                        elseif i==hz && j~=1 && j~=v %Fall 8 letzte Spalte Mitte
                            bc(j,i)={[nodematrix((j-1)*a+1,(i-1)*b+1);nodematrix(j*a+1,(i-1)*b+1);nodematrix((j-1)*a+1,dim(2));nodematrix(j*a+1,dim(2))]};
                            br(j,i)={[nodematrix((j-1)*a+2:j*a,(i-1)*b+1);nodematrix((j-1)*a+1,dim(2)-b+1:dim(2)-1).';nodematrix(j*a+1,dim(2)-b+1:dim(2)-1).']};
                            K(j,i)={nodematrix((j-1)*a+1:j*a+1,(i-1)*b+1:dim(2))};
                        else %Fall 9 rechte untere Ecke
                            bc(j,i)={[nodematrix((j-1)*a+1,(i-1)*b+1);nodematrix(dim(1),(i-1)*b+1);nodematrix((j-1)*a+1,dim(2))]};
                            br(j,i)={[nodematrix((j-1)*a+2:dim(1)-1,dim(2)-b);nodematrix(dim(1)-a,(dim(2)-b+1:dim(2)-1)).']};
                            K(j,i)={nodematrix((j-1)*a+1:dim(1),(i-1)*b+1:dim(2))};
                        end
                        in(j,i)={setdiff(cell2mat(K(j,i)),union(cell2mat(bc(j,i)),cell2mat(br(j,i))))};
                        
                        gbc1=cell2mat(bc(j,i));
                        gbr1=cell2mat(br(j,i));
                        gin1=cell2mat(in(j,i));
                        
                        %size (gbc1...) l�uft als Variable mit
                        
                        gbc(size(gbc,1)+1:size(gbc,1)+size(gbc1,1),1)=gbc1;  %globaler Vektor der Eckknoten
                        gbr(size(gbr,1)+1:size(gbr,1)+size(gbr1,1),1)=gbr1;  %globaler Vektor der interface Knoten
                        gin(size(gin,1)+1:size(gin,1)+size(gin1,1),1)=gin1;  % globaler Vektor der internen Knoten
                        n=n+1;
                 end
            end
            end  
            end
            end      
       
       
       %% Verdopplung und Neubenneung der br Knoten und interface Elemente, speichern der neuen infos im femModel
       
%        function [doubleNodes]= getDoubleNodes(femModel,gbr)
%            k=1;
%            for i=1:length(gbr)
%            for j=i:length(gbr)
%                if gbr(i)==gbr(j)
%                    doubleNodes(k)=gbr(i);
%                    %idVector=[doubleNodes.getId];
%                    k=k+1;
%                end
%            end
%            end
%        end
%        
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
   end
end