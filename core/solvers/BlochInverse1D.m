classdef BlochInverse1D < Solver
    
    properties (Access = private)
        femModel
        assembler
        isInitialized
        
        stiffnessMatrix
        massMatrix    
        
        leftDofs
        rightDofs
        leftNodes
        rightNodes     
    end
    
    methods
        function obj = BlochInverse1D(femModel,leftNodes,rightNodes)
            if nargin > 0 
                obj.femModel = femModel;
                obj.assembler = SimpleAssembler(femModel);
                obj.isInitialized = false;
            else
                error('Error (BlochInverse1D): no fem model defined!')
            end        
            obj.leftNodes = leftNodes;
            obj.rightNodes = rightNodes;
        end %blochInverse1D
        
        function solve(obj, ~)
            if ~ obj.isInitialized
                obj.initialize();
            end
            
            %rest von solve fehlt (mit Blochtheorem)
        end
       
        
        function initialize(obj)
            if ~ obj.femModel.isInitialized()
                obj.femModel.initialize;
            end
            
            % assemble and reduce matrices
            obj.massMatrix = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            obj.stiffnessMatrix = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            
            %%% suche nach weiteren Fehlern einfügen
            
            
            
            
%             dofIds = Dof.getId;      %%%oder femModel.getId? ersrt mit femModel "objekt erstellen" zb femmodel.getarray
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

        end %end initialize
        
        function [leftNodes,rightNodes]=getLeftRightNodes(NodeArray)
            
        
        end
        
        function [leftDofs,rightDofs]=getLeftRightDofs(NodeArray)
            
        
        end
        
      
        
      
        
    end %end methods
    
end