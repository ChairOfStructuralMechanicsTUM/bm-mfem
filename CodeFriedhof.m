
      %             dofIds = Dof.getId;      %%%oder femModel.getId? erst mit femModel "objekt erstellen" zb femmodel.getarray
%             %nodeIds = Node.getId;    nur von 1 Knoten
%             %nodeCoords = node.getCoords   nur von 1 Knoten...
%             x = Node.getX;
%             y = Node.getY;
%             ysorted=sort(y);
%             xsorted=sort(x);
%             b=1;
%             while xsorted(1)==xsorted(b)
%                 b=b+1; %count the amount of same x-Coordinates/number of rows
%             end
%             a=zeros(1,b);
%             for i=1:b
%                 a(i)=1;
%                 x=1;
%                 while ysorted((i-1)*x+1) == ysorted(a(i))
%                     a(i)=a(i)+1;   
%                 end
%                 %a(i)= Anzahl gleicher y-Werte mit Koordinate x(i)
%                 x = a(1); 
%                 if a(1)~=a(i) %for every x=f, there is the same amount of nodes 
%                                 %with the same y-Coordinates
%                     disp('y-Coordinates in every row have to be equal')
%                 end
%                     
%             end
%             
%             
            
            %             nodes = femModel.getAllNodes;
%             a=2;
%             while nodes(3,1) == nodes(3,a) %Compare y-Coordinates of first row
%                 a=a+1;
%             end
%             n=a-1; %length of the beam
% 
%             a=1;
%             while nodes(2,1) == nodes(2,1+a*n) 
%                 %Compare x-Coordinates of first column
% 
%                 a=a+1;
%                 if size(nodes,2) < 1+a*n
%                     break
%                 end
%                 
%             end
%             m=a; %length (in y-direction) of the beam
% 
%             
%        
%             for i=1:m            
%                 x=1+n*(i-1);    %left node ID
%                 for j=1:n
% 
%                     if nodes(3,x)==nodes(3,x+j-1)
%                     else
%                         disp('y-Coordinates in the same row have to be equal')
%                     end 
%                 end
%           end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Aus function [...] = eliminateFixedDofs(obj,leftDofs,innerDofs,rightDofs)

% % %            %innerDofs(fixedDofIds==innerDofs)=[];  logical Indexing
% % %            a=0;
% % %            for i = 1:length(innerDofs)           
% % %                for j = 1:length(fixedDofIds)                                
% % %                    if innerDofs(i)+ a == fixedDofIds(j)
% % %                         for k = length(innerDofs):-1:i
% % %                             innerDofs(k) = innerDofs(k)-1;
% % %                         end
% % %                         a=a+1;
% % %                     end
% % %                end
% % %            end
% % %            a=0;
% % %            for i = 1:length(leftDofs)
% % %                for j = 1:length(fixedDofIds)              
% % %                     if leftDofs(i)+a == fixedDofIds(j)
% % %                         for k = length(leftDofs):-1:i
% % %                             leftDofs(k) = leftDofs(k)-1;
% % %                         end
% % %                         a = a+1;
% % %                     end
% % %                end
% % %            end
% % %            a=0;
% % %            for i = 1:length(rightDofs)
% % %                for j = 1:length(fixedDofIds)              
% % %                     if rightDofs(i)+a == fixedDofIds(j)
% % %                         for k = length(rightDofs):-1:i
% % %                             rightDofs(k) = rightDofs(k)-1;
% % %                         end
% % %                         a=a+1;
% % %                     end
% % %                end
% % %            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%        function [kx,miu] = propConst(obj,numberOfWaveNumbers) %obj kann auch entfernt werden -> dann muss die Funktion aber statisch sein!
%            kx = linspace(0,pi,numberOfWaveNumbers);    %15 ist viel zu wenig
%            miu = exp(i*kx); %kx:Phase
%        end
%           
%        function R = transformationMatrix(obj,miu,index)
%            nodeIdsRight = obj.rightNodes;
%            nodeIdsLeft = obj.leftNodes; 
%            leftDofs = obj.leftDofs;
%            rightDofs = obj.rightDofs;
%            M = obj.massMatrix;
%            R = [eye(length(leftDofs)), zeros(length(leftDofs),length(M)-2*length(leftDofs)); ...
%                zeros(length(M)-2*length(leftDofs),length(leftDofs)), eye(length(M)-2*length(leftDofs));...
%                miu(index)*eye(length(leftDofs)), zeros(length(leftDofs),length(M)-2*length(leftDofs))];
% 
%            if miu == exp(i*pi/2)
%                disp('R_pi/2 = ')
%                disp(R)
%            end
%        end          
%         
%        function [Ksorted,Msorted] = sortKandM(obj,K,M) 
%            leftDofs = obj.leftDofs;
%            rightDofs = obj.rightDofs;
%            innerDofs = getInnerDofIds(obj);         
%            vecdofs = [leftDofs,innerDofs,rightDofs];
%            Ksorted = K(vecdofs,vecdofs);
%            Msorted = M(vecdofs,vecdofs);
%        end
%        
%        
%        
%        function [Kred,Mred] = reducedStiffnesAndMass (obj,K,M)
%            numberOfWaveNumbers = 10;
%             
%            [kx,miu] = propConst(obj,numberOfWaveNumbers);
%             Kred = cell(numberOfWaveNumbers,1);
%             Mred = cell(numberOfWaveNumbers,1);
%             
%             
%             for i=1:numberOfWaveNumbers
%             
%                 
%             R = transformationMatrix(obj,miu,i);    
% %             Kred{i,1} = conj(R)'*K*R;
% %             Mred{i,1} = conj(R)'*M*R;
%               Mred{i,1} = R'*M*R;
%               Kred{i,1} = R'*K*R;
%             end
%        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%