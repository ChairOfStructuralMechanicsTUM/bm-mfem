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
            
            %PROBLEM: If the node Ids in a substructure are not in
            %increasing order one can't just simply subtract the smallest
            %Id from all Ids, there will be an inconsistency. Maybe add
            %local Ids in a substructure?
            
            nDofs = length(femModel.getDofArray);
            nNodalDofs = nDofs / length(femModel.getAllNodes);
            stiffnessMatrix = zeros(nDofs);
            
            %find smallest nodeId
            nodes = femModel.getAllNodes;
            mini = min(nodes.getId)-1;
            
            for itEle = 1:length(femModel.getAllElements)
                elements = femModel.getAllElements;
                currentElement = elements(itEle);
                elementFreedomTable = {};
                
                
                for itNode = 1:length(currentElement.getNodes)
                    nodes = currentElement.getNodes;
                    currentNode = nodes(itNode);
                    globalDofArray = zeros(1,nNodalDofs);
                    
                    %%%CHANGED: currentNode.getId is changed. The minimum
                    %nodeId is always substracted. This makes the nodeId
                    %always start from 1.
                    %globalDofArray(nNodalDofs) = nNodalDofs * currentNode.getId;
                    globalDofArray(nNodalDofs) = nNodalDofs * (currentNode.getId-mini);
                   
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
        

%         BEFORE MODIFICATION --> funktioniert
%         function reducedForceVector = applyExternalForces(femModel)
%         
%             dofs = femModel.getDofArray;
%             nDofs = length(dofs);
%             reducedForceVector = zeros(1,nDofs);
%             fixedDofs = [];
%             
%             for itDof = 1:nDofs
%                 if (dofs(itDof).isFixed)
%                     fixedDofs = [fixedDofs itDof];                          % array of fixed dofs and their location
%                 else
%                     reducedForceVector(itDof) = dofs(itDof).getValue;
%                 end
%             end
%             
%             reducedForceVector(fixedDofs) = [];
%             
%             
%         end
        
        
    function [forceVector, reducedForceVector] = applyExternalForces(femModel)
        dofs = femModel.getDofArray;
        nDofs = length(dofs);
        forceVector = zeros(1,nDofs);
        fixedDofs = [];
    
        for itDof = 1:nDofs
            if (dofs(itDof).isFixed)
                fixedDofs = [fixedDofs itDof];                          % array of fixed dofs and their location
            else
                %forceVector(itDof) = dofs(itDof).getValue;
                %forceVector(itDof)
                %x = dofs(itDof).getDofLoad
                
                forceVector(itDof) = dofs(itDof).getDofLoad;
                
            end
        end
        
        reducedForceVector = forceVector;
        reducedForceVector(fixedDofs) = [];
            
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

