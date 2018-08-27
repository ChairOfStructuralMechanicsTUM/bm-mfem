classdef BlochInverse1D_mm < Solver
    
    properties (Access = private)
        femModel
        assembler
        isInitialized
        
        stiffnessMatrix
        massMatrix
        
        Ksorted
        Msorted
        
        leftDofs
        rightDofs
        interiorDofs
        leftNodes
        rightNodes
    end
    
    methods
        function obj = BlochInverse1D_mm(femModel)
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
        
        function [phases,frequencies] = solve(obj,numberOfPhases,numberOfBands)
            if ~ obj.isInitialized
                obj.initialize();
            end
            
            % calculate the transformed matrices
            [phases,mius] = obj.propConst(numberOfPhases);
            
            % initialize empty solution vector for frequencies
            % each row of represents one band
            frequencies=zeros(numberOfBands,length(phases));
            
            % iterate through all mius/phases
            for i=1:length(mius)
                miu=mius(i);
                
                % compute reduced matrices
                [Kred,Mred] = reducedStiffnesAndMass(obj,miu);
                
                % solve eigenvalue problem
                omega2 = eigs(Kred,Mred,numberOfBands,'sm');
                freq = sqrt(abs(omega2))/(2*pi);
                
                % store the solution
                frequencies(:,i)=freq;
            end
        end
        
        
        function initialize(obj)
            if ~ obj.femModel.isInitialized()
                obj.femModel.initialize;
            end
            
            % find left and right Nodes
            nodeIdsLeft = obj.findLeftNodes();
            nodeIdsRight = obj.findRightNodes();
            
            % save it as property of the object
            obj.leftNodes = nodeIdsLeft;
            obj.rightNodes = nodeIdsRight;
            
            % find dof Ids of repective nodes
            [leftDofIds,rightDofIds,interiorDofIds] = getLeftRightInteriorDofIds(obj);
            
            % save it as property of the object
            obj.leftDofs = leftDofIds;
            obj.rightDofs = rightDofIds;
            obj.interiorDofs = interiorDofIds;
            
            % test the assignement of the nodes
            obj.testAssignmentOfNodes()
            
            % get FE-Matrices from femModel
            obj.massMatrix = obj.assembler.assembleGlobalMassMatrix(obj.femModel);
            obj.stiffnessMatrix = obj.assembler.assembleGlobalStiffnessMatrix(obj.femModel);
            
            % sort matrices with respect to left interior and right dof
            [obj.Ksorted,obj.Msorted] = obj.sortKandM();
             
        end %end initialize
        
        function testAssignmentOfNodes(obj)
            leftNodes = getNodes(obj.femModel, obj.leftNodes);
            rightNodes = getNodes(obj.femModel, obj.rightNodes);
            
            leftNodeX = getX(leftNodes);
            leftNodeY = getY(leftNodes);
            rightNodeX = getX(rightNodes);
            rightNodeY = getY(rightNodes);
                        
            disp('left Nodes: [id,x,y]')
            X=[obj.leftNodes.' leftNodeX.' leftNodeY.'];
            disp(X)
            
            disp('right Nodes: [id,x,y]')
            Y=[obj.rightNodes.' rightNodeX.' rightNodeY.'];
            disp(Y) 
            
            if length(obj.leftNodes) ~= length(obj.rightNodes)
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
            obj.leftNodes = nodeIdsLeft;    %anstatt  obj.leftNodes(n) = nodeIds(i);
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
            obj.rightNodes = nodeIdsRight;
            fprintf('Number of right boundary nodes is %s. \n', num2str(n))
        end
        
        
        function [leftDofIds,rightDofIds,interiorDofIds] = getLeftRightInteriorDofIds(obj)  %right and left nodes are requiered
            nodeArray = obj.femModel.getAllNodes;
            nodeIdsRight = obj.rightNodes;
            nodeIdsLeft = obj.leftNodes;
            
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            dofArray1 = arrayfun(@(node) node.getDofArray, nodeArray, 'UniformOutput', false)';  %dof Array1=[(2Spaltex1Zeile);(2x1);(2x1)...],  2=x und y FG
            %a=dofArray1{2,1}(1,2) works here
            % dofArray = [dofArray1{:}];   %dofArray=[1x1 1x1 1x1] (mit doppelter Länge, da x und y jew eigene Einträge haben)
            %             femModel.dofArray = reshape(femModel.dofArray,1,size(femModel.dofArray,1)*size(femModel.dofArray,2));
            %             for ii = 1:length(dofArray)
            %                 dofArray(ii).setId(ii);
            %             end
            
            indicesRightNodes=find(ismember(nodeIds,nodeIdsRight));
            indicesLeftNodes=find(ismember(nodeIds,nodeIdsLeft));
            
            dofArrayRight=dofArray1(indicesRightNodes);  %%dofArray
            dofArrayRight = [dofArrayRight{:}];
            
            dofArrayLeft=dofArray1(indicesLeftNodes);
            dofArrayLeft = [dofArrayLeft{:}];
            
            dofArrayInterior=dofArray1;
            dofArrayInterior([indicesRightNodes,indicesLeftNodes])=[];
            dofArrayInterior = [dofArrayInterior{:}];
            
            rightDofIds = arrayfun(@(dof) dof.getId, dofArrayRight);
            leftDofIds = arrayfun(@(dof) dof.getId, dofArrayLeft);
            interiorDofIds = arrayfun(@(dof) dof.getId, dofArrayInterior);
            
            %%Eliminate FixedDofs
            [freeDofs, fixedDofs] = getDofConstraints(obj.femModel);
            fixedDofIds = getId(fixedDofs);
            %             fixedDofIds = arrayfun(@(dof) dof.getId, fixedDofs);
            
            indicesFixedRightIds = find(ismember(rightDofIds,fixedDofIds));
            rightDofIds(indicesFixedRightIds)=[];
            
            indicesFixedLeftIds = find(ismember(leftDofIds,fixedDofIds));
            leftDofIds(indicesFixedLeftIds)=[];
            
            indicesFixedInteriorIds = find(ismember(interiorDofIds,fixedDofIds));
            interiorDofIds(indicesFixedInteriorIds)=[];
            
            % save dofIds in object properties
            obj.leftDofs=leftDofIds;
            obj.rightDofs=rightDofIds;
            obj.interiorDofs=interiorDofIds;
            
        end
        
        
        %        function innerDofs = getInnerDofIds(obj)
        %             nodeArray = obj.femModel.getAllNodes;
        %             leftDofs = obj.leftDofs;
        %             rightDofs = obj.rightDofs;
        %             dofArray1 = arrayfun(@(node) node.getDofArray, nodeArray, 'UniformOutput', false)';
        %             dofArray = [dofArray1{:}];
        % %             for ii = 1:length(dofArray)
        % %                 dofArray(ii).setId(ii);
        % %             end
        %
        %             dofIDs = arrayfun(@(Dof) Dof.getId , dofArray);
        %
        %             innerDofs = zeros(1,length(dofArray)-length(leftDofs)*2);
        %             x=0;
        %             for i = 1:length(dofArray)
        %                 a=0;
        %                 for j = 1:length(leftDofs)
        %
        %                     if getId(dofArray(i)) == leftDofs(j) || getId(dofArray(i)) == rightDofs(j)
        %                         break
        %                     else
        %                         a=a+1;
        %                     end
        %
        %                 end
        %                 if a==length(leftDofs)
        %                     x = x+1;
        %                     innerDofs(x) = getId(dofArray(i));
        %
        %                 end
        %             end
        % %             disp(innerDofs)%not necessary
        %        end
        
        
        
        
        
        
        function R = transformationMatrix(obj,miu)   %%funktioniert, da nur Dimensionen geändert werden mussten
            %            nodeIdsRight = obj.rightNodes;
            %            nodeIdsLeft = obj.leftNodes;
            %            leftDofs = obj.leftDofs;
            %            rightDofs = obj.rightDofs;
            %            [freeDofs, fixedDofs] = getDofConstraints(obj.femModel);
            %            fixedDofIds = getId(fixedDofs);
            %
            %            indices1=find(ismember(leftDofs,fixedDofIds));
            %            leftDofs(indices1)=[];
            %            M = obj.massMatrix;
            %            X = length(M)-length(fixedDofs);
            %            R = [eye(length(leftDofs)), zeros(length(leftDofs),X-2*length(leftDofs)); ...
            %                zeros(X-2*length(leftDofs),length(leftDofs)), eye(X-2*length(leftDofs));...
            %                miu(index)*eye(length(leftDofs)), zeros(length(leftDofs),X-2*length(leftDofs))];
            
            
            numL = length(obj.leftDofs);
            numR = length(obj.rightDofs);
            numI = length(obj.interiorDofs);
            
            
            R = [eye(numL)       , zeros(numL,numI) ; ...
                zeros(numI,numL), eye(numI)        ; ...
                miu*eye(numR)   , zeros(numR,numI)];
        end
        
        function [Ksorted,Msorted] = sortKandM(obj)
            vecdofsAll = [obj.leftDofs,obj.interiorDofs,obj.rightDofs];
            Ksorted = obj.stiffnessMatrix(vecdofsAll,vecdofsAll);
            Msorted = obj.massMatrix(vecdofsAll,vecdofsAll);
        end
        
        
        function [Kred,Mred] = reducedStiffnesAndMass(obj,miu)
            
            R = transformationMatrix(obj,miu);
            
            Mred = R'*obj.Msorted*R;
            Kred = R'*obj.Ksorted*R;
        end
        
    end %end methods
    
    methods(Static)
        function [phase,miu] = propConst(numberOfPhases)
            phase = linspace(10e-6,pi,numberOfPhases);
            miu = exp(1i*phase); %kx:Phase
        end
        
    end %end static methods
    
    
    
    
end %end classdef


