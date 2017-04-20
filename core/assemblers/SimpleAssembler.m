classdef SimpleAssembler < Assembler
    %SIMPLEASSEMBLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function stiffnessMatrix = assembleGlobalStiffnessMatrix(femModel)
            nDofs = length(femModel.getDofArray);
            nNodalDofs = nDofs / length(femModel.getNodeArray);
            stiffnessMatrix = zeros(nDofs);
            
            for itEle = 1:length(femModel.getElementArray)
                elements = femModel.getElementArray;
                currentElement = elements(itEle);
                elementFreedomTable = {};
                
                for itNode = 1:length(currentElement.getNodeArray)
                    nodes = currentElement.getNodeArray;
                    currentNode = nodes(itNode);
                    globalDofArray = zeros(1,nNodalDofs);
                    globalDofArray(nNodalDofs) = nNodalDofs * currentNode.getId;
                   
                    for i = (nNodalDofs - 1) : -1 : 1
                        globalDofArray(i) = globalDofArray(i+1) - 1;
                    end
                    % das ist super hässlich, geht das auch anders?
                    elementFreedomTable =[elementFreedomTable globalDofArray];
                    
                end
                
                elementFreedomTable = [elementFreedomTable{:}];
                elementStiffnessMatrix = currentElement.computeLocalStiffnessMatrix;
                
                stiffnessMatrix(elementFreedomTable, elementFreedomTable) ...
                    = stiffnessMatrix(elementFreedomTable, elementFreedomTable) ...
                    + elementStiffnessMatrix;
                
            end
        end
        
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
        
        function forceVector = applyExternalForces(femModel)
            dofs = femModel.getDofArray;
            nDofs = length(dofs);
            forceVector = zeros(1,nDofs);
            fixedDofs = [];
            
            for itDof = 1:nDofs
                if (dofs(itDof).isFixed)
                    fixedDofs = [fixedDofs itDof];
                else
                    forceVector(itDof) = dofs(itDof).getValue;
                end
            end
            
            forceVector(fixedDofs) = [];
            
        end
        
        function assignResultsToDofs(femModel, resultVector)
            dofs = femModel.getDofArray;
            nDofs = length(dofs);
            itResult = 1;
            
            for itDof = 1:nDofs
               if (~dofs(itDof).isFixed)
                  dofs(itDof).setValue(resultVector(itResult));
                  itResult = itResult + 1;
               end
            end
            
        end
        
    end
    
end

