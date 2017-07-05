classdef SimpleAssembler < Assembler
    %SIMPLEASSEMBLER Simple elimination-based assembler
    %   Detailed explanation goes here
    
    properties
        stiffnessMatrix
        reducedStiffnessMatrix
        forceVector
        reducedForceVector
    end
    
    methods %test
        
        function assembling = SimpleAssembler(femModel)
            if (nargin > 0)
                
                [assembling.stiffnessMatrix, assembling.reducedStiffnessMatrix] = SimpleAssembler.assembleGlobalStiffnessMatrix(femModel);
                [assembling.forceVector, assembling.reducedForceVector] = SimpleAssembler.applyExternalForces(femModel);
                
            else
                error('input model is missing');
            end
        end
    
      
    end
    
    methods (Static)
        
        function [stiffnessMatrix, reducedStiffnessMatrix] = assembleGlobalStiffnessMatrix(femModel)
            nDofs = length(femModel.getDofArray);
            nNodalDofs = nDofs / length(femModel.getAllNodes);
            stiffnessMatrix = zeros(nDofs);
            
            for itEle = 1:length(femModel.getAllElements)
                elements = femModel.getAllElements;
                currentElement = elements(itEle);
                elementFreedomTable = {};
                
                for itNode = 1:length(currentElement.getNodes)
                    nodes = currentElement.getNodes;
                    currentNode = nodes(itNode);
                    globalDofArray = zeros(1,nNodalDofs);   %ids of the global dofs
                    globalDofArray(nNodalDofs) = nNodalDofs * currentNode.getId;
                   
                    for i = (nNodalDofs - 1) : -1 : 1
                        globalDofArray(i) = globalDofArray(i+1) - 1;
                    end
                    
                    elementFreedomTable = [elementFreedomTable globalDofArray];
                    
                end
                
                elementFreedomTable = [elementFreedomTable{:}];
                elementStiffnessMatrix = currentElement.computeLocalStiffnessMatrix;
                
                stiffnessMatrix(elementFreedomTable, elementFreedomTable) ...
                    = stiffnessMatrix(elementFreedomTable, elementFreedomTable) ...
                    + elementStiffnessMatrix;
                
            end
            
            reducedStiffnessMatrix = SimpleAssembler.applyBoundaryConditions(femModel, stiffnessMatrix);
            
        end
        
       
        
        function [forceVector, reducedForceVector] = applyExternalForces(femModel)
            dofs = femModel.getDofArray;
            nDofs = length(dofs);
            forceVector = zeros(1,nDofs);
            fixedDofs = [];
            
            % get the point load on the dofs
            for itDof = 1:nDofs
                if (dofs(itDof).isFixed)
                    fixedDofs = [fixedDofs itDof];   % array of fixed dofs and their location
                else
                    forceVector(itDof) = dofs(itDof).getDofLoad;
                end
            end
            
            % get the external forces from every element
            elements = femModel.getAllElements;
            for itEle = 1:length(elements)
                elementalForceVector = elements(itEle).computeLocalForceVector;
                elementalDofs = elements(itEle).getDofs;
                forceVector(elementalDofs.getId) = forceVector(elementalDofs.getId) + elementalForceVector;
            end
            
            reducedForceVector = forceVector;
            reducedForceVector(fixedDofs) = [];
            
            
        end
        
        
        
        
        
        function assignResultsToDofs(femModel, resultVector)
            dofs = femModel.getDofArray;
            nDofs = length(dofs);
            itResult = 1;
            nResultSteps = size(resultVector,2);
            
            for itDof = 1:nDofs
                if dofs(itDof).isFixed
                    values(1:nResultSteps) = dofs(itDof).getValue(1);
                    dofs(itDof).setValue(values);
                else
                    dofs(itDof).setValue(resultVector(itResult,:));
                    itResult = itResult + 1;
                end
                
            end
            
        end
        
        function appendValuesToDofs(femModel, values)
            [freeDofs, fixedDofs] = femModel.getDofConstraints();
            if length(freeDofs) ~= length(values)
                error('the arrays of dofs and values are not of the same size')
            end
            
            freeDofs.appendValue(values);
            fixedDofs.appendValue(zeros(1,length(fixedDofs)));
        end
        
        function appendValuesToNodes(femModel, valueName, values)
            [freeDofs, fixedDofs] = femModel.getDofConstraints();
            if length(freeDofs) ~= length(values)
                error('the arrays of dofs and values are not of the same size')
            end
            
            %append the values
            for ii = 1:length(freeDofs)
                dof = freeDofs(ii);
                dofName = dof.getValueType;
                dofDirection = dofName(end-1:end);
                node = dof.getNode;
                node.appendStepValue(strcat(valueName, dofDirection), values(ii));
            end
            
            %append 0 to the fixed dofs
            for ii = 1:length(fixedDofs)
                dof = fixedDofs(ii);
                dofName = dof.getValueType;
                dofDirection = dofName(end-1:end);
                node = dof.getNode;
                node.appendStepValue(strcat(valueName, dofDirection), 0);
            end
        end
        
        function assignValuesToNodes(femModel, valueName, values, step)
            [freeDofs, fixedDofs] = femModel.getDofConstraints();
            if length(freeDofs) ~= length(values)
                error('the arrays of dofs and values are not of the same size')
            end
            
            %append the values
            for ii = 1:length(freeDofs)
                dof = freeDofs(ii);
                dofName = dof.getValueType;
                dofDirection = dofName(end-1:end);
                node = dof.getNode;
                node.setStepValue(strcat(valueName, dofDirection), values(ii), step);
            end
            
            %append 0 to the fixed dofs
            for ii = 1:length(fixedDofs)
                dof = fixedDofs(ii);
                dofName = dof.getValueType;
                dofDirection = dofName(end-1:end);
                node = dof.getNode;
                node.setStepValue(strcat(valueName, dofDirection), 0, step);
            end
        end
        
        function massMatrix = assembleGlobalMassMatrix(femModel)
            elements = femModel.getAllElements;
            ndofs = length(femModel.getDofArray);
            massMatrix = zeros(ndofs);
            
            for itEle = 1:length(elements)
               elementalMassMatrix = elements(itEle).computeLocalMassMatrix;
               elementalDofs = elements(itEle).getDofs;
               massMatrix(elementalDofs.getId, elementalDofs.getId) = ...
                   massMatrix(elementalDofs.getId, elementalDofs.getId) + elementalMassMatrix;
            end
        end
        
        function dampingMatrix = assembleGlobalDampingMatrix(femModel)
            elements = femModel.getAllElements;
            ndofs = length(femModel.getDofArray);
            dampingMatrix = zeros(ndofs);
            
            for itEle = 1:length(elements)
               elementalDampingMatrix = elements(itEle).computeLocalDampingMatrix;
               elementalDofs = elements(itEle).getDofs;
               dampingMatrix(elementalDofs.getId, elementalDofs.getId) = ...
                   dampingMatrix(elementalDofs.getId, elementalDofs.getId) + elementalDampingMatrix;
            end
        end
        
        function vals = assemble3dDofVector(femModel, dofName)
            dofs = femModel.getDofArray;
            vals = zeros(1,length(dofs));
            for itDof = 1:length(dofs)
                dof = dofs(itDof);
                itDofName = dof.getValueType();
                itDofName = itDofName(1:end-2);
                if strcmp(dofName, itDofName)
                    val = dof.getValue();
                    vals(dof.getId) = val(end);
                end
            end
        end
        
        function vals = assemble3dValueVector(femModel, valueName)
            %ASSEMBLE3DVALUEVECTOR returns a vector with all values from
            %VALUENAME corresponding to the global dofs
            dofs = femModel.getDofArray;
            vals = zeros(1,length(dofs));
            for itDof = 1:length(dofs)
                dof = dofs(itDof);
                dofName = dof.getValueType;
                dofDirection = dofName(end-1:end);
                node = dof.getNode;
                val = node.getValue(strcat(valueName, dofDirection));
                val = val(end);
                vals(dof.getId) = val;
            end
        end
        
        function stiffnessMatrix = getStiffnessMatrix(assembler)
            stiffnessMatrix = assembler.stiffnessMatrix;
        end
        
        function reducedStiffnessMatrix = getReducedStiffnessMatrix(assembler)
            reducedStiffnessMatrix = assembler.reducedStiffnessMatrix;
        end
        
        function forceVector = getForceVector(assembler)
            forceVector = assembler.forceVector;
        end
        
        function reducedForceVector = getReducedForceVector(assembler)
            reducedForceVector = assembler.reducedForceVector;
        end
        
    end
    
    methods (Static, Access = private)
        
        function stiffnessMatrix = applyBoundaryConditions(femModel, stiffnessMatrix)
            dofs = femModel.getDofArray;
            fixedDofs = [];
            
            for itDof = 1:length(dofs)
               if dofs(itDof).isFixed
                   fixedDofs = [fixedDofs itDof];
               end
            end
            
            stiffnessMatrix(fixedDofs,:) = [];
            stiffnessMatrix(:,fixedDofs) = [];
            
        end
        
    end
    
end

