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
%             obj.leftNodes = obj.findLeftNodes();
%             
%             obj.rightNodes = obj.findRightNodes();
        end 
        
        function solve(obj, ~)
            if ~ obj.isInitialized
                obj.initialize();
            end            
            %rest von solve fehlt (mit Blochtheorem)
            
            
        end
       
        
        function [nodeIdsLeft] = findLeftNodes(obj)
           
            nodeArray = obj.femModel.getAllNodes;
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            nodeXcoords = arrayfun(@(node) node.getX, nodeArray);
           
            sortedX = sort(nodeXcoords);
            minX = sortedX(1);       
            n=0;         
            for i=1:length(nodeXcoords)
                if nodeXcoords(i) == minX
                    n = n+1;
                    nodeIdsLeft(n) = nodeIds(i);                     
                end
            end
            fprintf('Number of left boundary nodes is %s. \n', num2str(n))
            %obj.leftNodes = nodeIdsLeft;    %anstatt  obj.leftNodes(n) = nodeIds(i);
        end
            
       function [nodeIdsRight] = findRightNodes(obj) 
             
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
                end
            end
%          
            fprintf('Number of right boundary nodes is %s. \n', num2str(n))
       end
             
            
       function [leftDofs,rightDofs] = getLeftRightDofIds(obj)  %right and left nodes are requiered
            nodeArray = obj.femModel.getAllNodes;
            nodeIdsRight = obj.rightNodes;
            nodeIdsLeft = obj.leftNodes;      

            dofArray1 = arrayfun(@(node) node.getDofArray, nodeArray, 'UniformOutput', false)';  %dof Array1=[(2Spaltex1Zeile);(2x1);(2x1)...],  2=x und y FG
            %a=dofArray1{2,1}(1,2) works here        
             dofArray = [dofArray1{:}];   %dofArray=[1x1 1x1 1x1] (mit doppelter Länge, da x und y jew eigene Einträge haben)
%             femModel.dofArray = reshape(femModel.dofArray,1,size(femModel.dofArray,1)*size(femModel.dofArray,2));
            for ii = 1:length(dofArray)
                dofArray(ii).setId(ii);
            end           

            n = length(dofArray1{1,1});
            fprintf('%s degrees of freedom',num2str(n))
            leftDofs = zeros(1,length(nodeIdsLeft)*n);
            rightDofs = zeros(1,length(nodeIdsRight)*n);
            x=0;
            for i=1:(length(nodeIdsLeft))
                for j=1:n
                    x=x+1;
                    leftDofs(x) = getId(dofArray(nodeIdsLeft(i)*n-n+j));
                    rightDofs(x) = getId(dofArray(nodeIdsRight(i)*n-n+j));
                end
            end
                                        
       end
       
       function innerDofs = getInnerDofIds(obj)
            nodeArray = obj.femModel.getAllNodes;   
            leftDofs = obj.leftDofs;
            rightDofs = obj.rightDofs;
            dofArray1 = arrayfun(@(node) node.getDofArray, nodeArray, 'UniformOutput', false)';  
            dofArray = [dofArray1{:}];
            for ii = 1:length(dofArray)
                dofArray(ii).setId(ii);
            end      
            
            innerDofs = zeros(1,length(dofArray)-length(leftDofs)*2);
            x=0;
            for i = 1:length(dofArray)
                a=0;
                for j = 1:length(leftDofs)
                    
                    if getId(dofArray(i)) == leftDofs(j) || getId(dofArray(i)) == rightDofs(j)
                        break
                    else 
                        a=a+1;
                    end
                    
                end
                if a==length(leftDofs)
                    x = x+1;
                    innerDofs(x) = getId(dofArray(i));
                   
                end
            end
%             disp(innerDofs)%not necessary
       end
       
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

       function [kx,miu] = propConst(obj,numberOfPhases) %obj kann auch entfernt werden -> dann muss die Funktion aber statisch sein!
           kx = linspace(0,pi,numberOfPhases);    
           miu = exp(i*kx); %kx:Phase
       end
           
        
       function R = transformationMatrix(obj,miu,index)
           nodeIdsRight = obj.rightNodes;
           nodeIdsLeft = obj.leftNodes; 
           leftDofs = obj.leftDofs;
           rightDofs = obj.rightDofs;
           M = obj.massMatrix;
           R = [eye(length(leftDofs)), zeros(length(leftDofs),length(M)-2*length(leftDofs)); ...
               zeros(length(M)-2*length(leftDofs),length(leftDofs)), eye(length(M)-2*length(leftDofs));...
               miu(index)*eye(length(leftDofs)), zeros(length(leftDofs),length(M)-2*length(leftDofs))];


       end          
        
       function [Ksorted,Msorted] = sortKandM(obj,K,M) 
           leftDofs = obj.leftDofs;
           rightDofs = obj.rightDofs;
           innerDofs = getInnerDofIds(obj);         
           vecdofs = [leftDofs,innerDofs,rightDofs];
           Ksorted = K(vecdofs,vecdofs);
           Msorted = M(vecdofs,vecdofs);
       end
       
       
       
       function [Kred,Mred] = reducedStiffnesAndMass (obj,K,M,numberOfPhases)
          
            
           [kx,miu] = propConst(obj,numberOfPhases);
            Kred = cell(numberOfPhases,1);
            Mred = cell(numberOfPhases,1);
            
            
            for i=1:numberOfPhases
            
                
            R = transformationMatrix(obj,miu,i);    
%             Kred{i,1} = conj(R)'*K*R;
%             Mred{i,1} = conj(R)'*M*R;
              Mred{i,1} = R'*M*R;
              Kred{i,1} = R'*K*R;
            end
       end

       
        function initialize(obj)
            if ~ obj.femModel.isInitialized()
                obj.femModel.initialize;
            end
                       
            % assemble and reduce matrices
            obj.massMatrix = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            obj.stiffnessMatrix = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            
            
            %nodeArray = obj.femModel.getAllNodes; 
            %momentan auch in"findright/left Nodes/Dofs"
            
            nodeIdsLeft = obj.findLeftNodes;            
            %nodeIdsRight = findRightNodes(obj);
            nodeIdsRight = obj.findRightNodes;
            
            obj.leftNodes = nodeIdsLeft; %auch in Funktionen möglich
            obj.rightNodes = nodeIdsRight;
            
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
            
            
            [leftDofs,rightDofs] = getLeftRightDofIds(obj);
%             [leftDofs2,rightDofs2] = obj.getLeftRightDofIds
                        
            obj.leftDofs = leftDofs;
            obj.rightDofs = rightDofs;
            
            
            
            if length(nodeIdsLeft) ~= length(nodeIdsRight)
                error('Same amount of left and right boundary nodes are requiered')
            end
            
            for i = 1:length(leftNodeY)
                if leftNodeY(i) ~= rightNodeY(i)
                    error('corresponding boundary nodes must have the same y-coordinates')
                end
            end
            
            for i = 1:length(leftNodeX)
                if leftNodeX(1) ~= leftNodeX(i)
                    error('All left boundary nodes must have the same x-coordinates')
                end
                if rightNodeX(1) ~= rightNodeX(i)
                    error('All right boundary nodes must have the same x-coordinates')
                end
            end
                


        end %end initialize
        

        
    end %end methods
    
    methods (Static)
        
%        function omega = calcOmega(Kred,Mred)
%            omega2 = eigs(Kred,Mred,4,'sm');
%            omega = sqrt(abs(omega2));          
%        end
       function omega = calcOmega(Kred,Mred,numberOfBands)
           omega2 = eigs(Kred,Mred,numberOfBands,'sm');
           omega = sqrt(abs(omega2));          
       end
    end
    
end