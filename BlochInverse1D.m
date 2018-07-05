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
        function obj = BlochInverse1D(femModel)
            if nargin > 0 
                obj.femModel = femModel;
                obj.assembler = SimpleAssembler(femModel);
                obj.isInitialized = false;
            else
                error('Error (BlochInverse1D): no fem model defined!')
            end        
            obj.leftNodes = obj.findLeftNodes();
            
            obj.rightNodes = obj.findRightNodes();
        end %blochInverse1D
        
        function solve(obj, ~)
            if ~ obj.isInitialized
                obj.initialize();
            end
            
            %rest von solve fehlt (mit Blochtheorem)
        end
       
        
        function [nodeIdsLeft] = findLeftNodes(obj)
            %function [leftNodes] = findLeftNodes(obj)
            nodeArray = obj.femModel.getAllNodes;
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            nodeXcoords = arrayfun(@(node) node.getX, nodeArray);
            %nodeYcoords = arrayfun(@(node) node.gety, nodeArray);
            sortedX = sort(nodeXcoords);
            minX = sortedX(1);
            
            
            n=0;         
            for i=1:length(nodeXcoords)
                if nodeXcoords(i) == minX
                    n = n+1;
                    nodeIdsLeft(n) = nodeIds(i);        %Iterative Vergrößerung schlimm?
                    %%nodeXcoordsLeft(n) = nodeXcoords(i);
                    %[leftNodes] = NodeArray(nodeIds(i));  kann ich davon
                    %ausgehen, dass alle Knoten in nodeArray
                    %aufeinanderfolgende Ids haben? ->in NodeArray hat
                    %NodeArray(1) die ID 1
                end
            end
            fprintf('Number of left boundary nodes is %s. \n', num2str(n))
            
        end
            
       function [nodeIdsRight] = findRightNodes(obj) 
           %function [rightNodes] = findRightNodes(obj)    
            nodeArray = obj.femModel.getAllNodes;
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            nodeXcoords = arrayfun(@(node) node.getX, nodeArray);
            sortedX = sort(nodeXcoords);
 
            maxX = sortedX(length(sortedX));
            n=0;         
            for i=1:length(nodeXcoords)
                if nodeXcoords(i) == maxX
                    n = n+1;
                    nodeIdsRight(n) = nodeIds(i);
                    %nodeXcoordsRight(n) = nodeXcoords(i); für function
                    %[nodeXcoordsRight,nodeIdsRight]=...                                     
                end
            end
%             rightNodes = [nodeIds;nodeXcoordsRight];
            fprintf('Number of right boundary nodes is %s. \n', num2str(n))
       end
             
            
            
        
        function initialize(obj)
            if ~ obj.femModel.isInitialized()
                obj.femModel.initialize;
            end
            
            % assemble and reduce matrices
            obj.massMatrix = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            obj.stiffnessMatrix = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            
            %%% suche nach weiteren Fehlern einfügen
            nodeArray = obj.femModel.getAllNodes;
            nodeIdsLeft = obj.findLeftNodes;
            
            %nodeIdsRight = findRightNodes(obj)
            nodeIdsRight = obj.findRightNodes;
            
            
            leftNodes = getNodes(obj.femModel, nodeIdsLeft); 
            rightNodes = getNodes(obj.femModel, nodeIdsRight);
            
           
            
            
            leftNodeX = getX(leftNodes);
            leftNodeY = getY(leftNodes);
            rightNodeX = getX(rightNodes);
            rightNodeY = getY(rightNodes);
            
            disp('left Nodes: [id,x,y]')
            X=[nodeIdsLeft.' leftNodeX.' leftNodeY.'];
            disp(X)
            
            disp('right Nodes: [id,x,y]')
            Y=[nodeIdsRight.' rightNodeX.' rightNodeY.'];
            disp(Y)
            
            if length(nodeIdsLeft) ~= length(nodeIdsRight)
                disp('Same amount of left and right boundary nodes are requiered')
            end
            
            for i = 1:length(leftNodeY)
                if leftNodeY(i) ~= rightNodeY(i)
                    disp('corresponding boundary nodes must have the same y-coordinates')
                end
            end
            
            for i = 1:length(leftNodeX)
                if leftNodeX(1) ~= leftNodeX(i)
                    disp('All left boundary nodes must have the same x-coordinates')
                end
                if rightNodeX(1) ~= rightNodeX(i)
                    disp('All right boundary nodes must have the same x-coordinates')
                end
            end
                
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

        end %end initialize
        
        
        function [leftDofs,rightDofs] = getLeftRightDofs(NodeArray)
            
        
        end
        
      
        
      
        
    end %end methods
    
end