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
        end
        
        function solve(obj, ~)
            if ~ obj.isInitialized
                obj.initialize();
            end
            
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
            
            fprintf('Number of right boundary nodes is %s. \n', num2str(n))
        end
        
        
        function [leftDofs,rightDofs] = getLeftRightDofIds(obj)
            nodeArray = obj.femModel.getAllNodes;
            nodeIdsRight = obj.rightNodes;
            nodeIdsLeft = obj.leftNodes;
            
            dofArray1 = arrayfun(@(node) node.getDofArray, nodeArray, 'UniformOutput', false)';
            
            dofArray = [dofArray1{:}];
            
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
            
        end
        
        function [kx,miu] = propConst(obj,numberOfPhases)
            kx = linspace(0,pi,numberOfPhases);
            miu = exp(1i*kx);
        end
        
        function R = transformationMatrix(obj,miu,index)
            nodeIdsRight = obj.rightNodes;
            nodeIdsLeft = obj.leftNodes;
            leftDofs = obj.leftDofs;
            rightDofs = obj.rightDofs;
            [freeDofs, fixedDofs] = getDofConstraints(obj.femModel);
            fixedDofIds = getId(fixedDofs);
            
            indices1=find(ismember(leftDofs,fixedDofIds));
            leftDofs(indices1)=[];
            M = obj.massMatrix;
            X = length(M)-length(fixedDofs);
            R = [eye(length(leftDofs)), zeros(length(leftDofs),X-2*length(leftDofs)); ...
                zeros(X-2*length(leftDofs),length(leftDofs)), eye(X-2*length(leftDofs));...
                miu(index)*eye(length(leftDofs)), zeros(length(leftDofs),X-2*length(leftDofs))];
        end
        
        function [Ksorted,Msorted] = sortKandM(obj,K,M)
            leftDofs = obj.leftDofs;
            rightDofs = obj.rightDofs;
            innerDofs = getInnerDofIds(obj);
            
            [reducedLeftDofs,reducedInnerDofs, reducedRightDofs,vecdofsAll] = eliminateFixedDofs(obj,leftDofs,innerDofs,rightDofs);
            
            Ksorted = K(vecdofsAll,vecdofsAll);
            Msorted = M(vecdofsAll,vecdofsAll);
        end
        
        function [reducedLeftDofs,reducedInnerDofs, reducedRightDofs,allDofs] = eliminateFixedDofs(obj,leftDofs,innerDofs,rightDofs)
            
            
            [freeDofs, fixedDofs] = getDofConstraints(obj.femModel);
            fixedDofIds = getId(fixedDofs);
            
            allDofs=[leftDofs,innerDofs,rightDofs];
            indicesAll =find(ismember(allDofs,fixedDofIds));
            
            a=0;
            for j = 1:length(fixedDofIds)
                for i=1:length(allDofs)
                    if allDofs(i) + a == fixedDofIds(j)
                        for k=1:length(allDofs)
                            if allDofs(k)>=fixedDofIds(j)
                                allDofs(k)=allDofs(k)-1;
                            end
                        end
                        a=a+1;
                    end
                end
            end
            
            allDofs(indicesAll)=[];
            
            indices1=find(ismember(leftDofs,fixedDofIds));
            indices2=find(ismember(rightDofs,fixedDofIds));
            indices3=find(ismember(innerDofs,fixedDofIds));
            
            leftDofs(indices1)=[];
            rightDofs(indices2)=[];
            innerDofs(indices3)=[];
            reducedLeftDofs = leftDofs;
            reducedRightDofs = rightDofs;
            reducedInnerDofs = innerDofs;
        end
        
        function [Kred,Mred] = reducedStiffnesAndMass (obj,K,M,numberOfPhases)
            
            [kx,miu] = propConst(obj,numberOfPhases);
            Kred = cell(numberOfPhases,1);
            Mred = cell(numberOfPhases,1);
            
            
            for i=1:numberOfPhases
                
                
                R = transformationMatrix(obj,miu,i);
                
                Mred{i,1} = R'*M*R;
                Kred{i,1} = R'*K*R;
            end
        end
        
        
        function initialize(obj)
            if ~ obj.femModel.isInitialized()
                obj.femModel.initialize;
            end
            
            obj.massMatrix = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            obj.stiffnessMatrix = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            
            nodeIdsLeft = obj.findLeftNodes;
            nodeIdsRight = obj.findRightNodes;
            
            obj.leftNodes = nodeIdsLeft;
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
            
        end
        
    end
    
    methods (Static)
        
        function omega = calcOmega(Kred,Mred,numberOfBands)
            omega2 = eigs(Kred,Mred,numberOfBands,'sm');
            omega = sqrt(abs(omega2));
        end
    end
    
end