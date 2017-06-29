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
               dampingMatrix(elementalDofs.getId) = dampingMatrix(elementalDofs.getId) + elementalDampingMatrix;
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

