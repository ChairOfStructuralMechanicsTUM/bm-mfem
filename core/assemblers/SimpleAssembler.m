classdef SimpleAssembler < Assembler
    %SIMPLEASSEMBLER Simple elimination-based assembler
    %   Detailed explanation goes here
    
    properties
%         stiffnessMatrix
%         reducedStiffnessMatrix
%         forceVector
%         reducedForceVector
    end
    
    methods %test
        
        function assembling = SimpleAssembler(femModel)
            if (nargin > 0)
                
%                 [assembling.stiffnessMatrix, assembling.reducedStiffnessMatrix] = SimpleAssembler.assembleGlobalStiffnessMatrix(femModel);
%                 [assembling.forceVector, assembling.reducedForceVector] = SimpleAssembler.applyExternalForces(femModel);
                
            else
                error('input model is missing');
            end
        end
    end
    
    methods (Static)
        
        function [stiffnessMatrix, reducedStiffnessMatrix] = assembleGlobalStiffnessMatrix(femModel)
            elements = femModel.getAllElements;
            ndofs = length(femModel.getDofArray);
            stiffnessMatrix = zeros(ndofs);
            
            for itEle = 1:length(elements)
               elementalStiffnessMatrix = elements(itEle).computeLocalStiffnessMatrix;
               elementalDofIds = elements(itEle).getDofList().getId;
               stiffnessMatrix(elementalDofIds, elementalDofIds) = ...
                   stiffnessMatrix(elementalDofIds, elementalDofIds) + elementalStiffnessMatrix;
            end
            [~, fixedDofs] = femModel.getDofConstraints;
            if ~ isempty(fixedDofs)
                fixedDofIds = fixedDofs.getId();
                reducedStiffnessMatrix = applyMatrixBoundaryConditions(stiffnessMatrix, fixedDofIds);                
            end
        end 
        
        function [forceVector, reducedForceVector] = applyExternalForces(femModel)
            dofs = femModel.getDofArray;
            nDofs = length(dofs);
            forceVector = zeros(1,nDofs);
            [~, fixedDofs] = femModel.getDofConstraints();
            
            % get the external load on the dofs
            for itDof = 1:nDofs
                if ( ~ dofs(itDof).isFixed)
                    forceVector(itDof) = dofs(itDof).getDofLoad;
                end
            end
            
             % get the external forces from every element
%             elements = femModel.getAllElements;
%             for itEle = 1:length(elements)                
%                 elementalForceVector = elements(itEle).computeLocalForceVector;
%                 elementalDofs = elements(itEle).getDofList;
%                 forceVector(elementalDofs.getId) = forceVector(elementalDofs.getId) + elementalForceVector;
%             end
            
            reducedForceVector = forceVector;
            reducedForceVector(fixedDofs.getId()) = [];
            
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
        
        function appendFirstDerivativeValuesToDofs(femModel, values)
            [freeDofs, fixedDofs] = femModel.getDofConstraints();
            if length(freeDofs) ~= length(values)
                error('the arrays of dofs and values are not of the same size')
            end
            
            freeDofs.appendFirstDerivativeValue(values);
            fixedDofs.appendFirstDerivativeValue(zeros(length(fixedDofs),1));
        end
        
        function setSecondDerivativeValuesToDofs(femModel, values)
            [freeDofs, fixedDofs] = femModel.getDofConstraints();
            if length(freeDofs) ~= length(values)
                error('the arrays of dofs and values are not of the same size')
            end
            
            freeDofs.setSecondDerivativeValue(values);
            fixedDofs.setSecondDerivativeValue(zeros(length(fixedDofs),1));
        end
        
        function appendSecondDerivativeValuesToDofs(femModel, values)
            [freeDofs, fixedDofs] = femModel.getDofConstraints();
            if length(freeDofs) ~= length(values)
                error('the arrays of dofs and values are not of the same size')
            end
            
            freeDofs.appendSecondDerivativeValue(values);
            fixedDofs.appendSecondDerivativeValue(zeros(length(fixedDofs),1));
        end
        
        function appendValuesToNodes(femModel, valueName, values)
            [freeDofs, ~] = femModel.getDofConstraints();
            if length(freeDofs) ~= length(values)
                error('the arrays of dofs and values are not of the same size')
            end
            
            %expand the given values to a full vector with length of all
            %system dofs
            valuesFull = zeros(length(femModel.getDofArray),1);
            valuesFull(freeDofs.getId()) = values;
            
            nodes = femModel.getAllNodes();
            for itNode = 1:length(nodes)
                node = nodes(itNode);
                nodeDofIds = node.getDofArray().getId();
                nodeValues = valuesFull(nodeDofIds);
                node.appendStepValue(valueName, nodeValues);
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
        
        function [massMatrix, reducedMassMatrix] = assembleGlobalMassMatrix(femModel)
            elements = femModel.getAllElements;
            ndofs = length(femModel.getDofArray);
            massMatrix = zeros(ndofs);
            
            for itEle = 1:length(elements)
               elementalMassMatrix = elements(itEle).computeLocalMassMatrix;
               elementalDofIds = elements(itEle).getDofList().getId;
               massMatrix(elementalDofIds, elementalDofIds) = ...
                   massMatrix(elementalDofIds, elementalDofIds) + elementalMassMatrix;
            end
            [~, fixedDofs] = femModel.getDofConstraints;
            if ~ isempty(fixedDofs)
                fixedDofIds = fixedDofs.getId();
                reducedMassMatrix = applyMatrixBoundaryConditions(massMatrix, fixedDofIds);                
            end
        end
        
        function dampingMatrix = assembleGlobalDampingMatrix(femModel)
            elements = femModel.getAllElements;
            ndofs = length(femModel.getDofArray);
            dampingMatrix = zeros(ndofs);
            
            for itEle = 1:length(elements)
               elementalDampingMatrix = elements(itEle).computeLocalDampingMatrix;
               elementalDofIds = elements(itEle).getDofList().getId;
               dampingMatrix(elementalDofIds, elementalDofIds) = ...
                   dampingMatrix(elementalDofIds, elementalDofIds) + elementalDampingMatrix;
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
            nodes = femModel.getAllNodes;
            for itNode = 1:length(nodes)
                node = nodes(itNode);
                nodeDofIds = node.getDofArray().getId();
                val = node.getValue(valueName);
                val = val(end,:);
                vals(nodeDofIds) = val;
            end
        end
        
        function vals = assembleValuesVector(femModel, step)
            %ASSEMBLEVALUESVECTOR returns a vector with all dof values of
            %the elements w.r.t. global dof ids
            dofs = femModel.getDofArray;
            nDofs = length(dofs);
            vals = zeros(1,nDofs);
            for itDof = 1:nDofs
                dof = dofs(itDof);
                vals(dof.getId()) = dof.getValue(step);
            end
        end
        
        function vals = assembleValuesVector2(femModel, valueName, step)
            %ASSEMBLEVALUESVECTOR returns a vector with all dof values of
            %the elements w.r.t. global dof ids
            dofs = femModel.getDofArray;
            vals = zeros(1,length(dofs));
            
            for itDof = 1:length(dofs)
               dof = dofs(itDof);
               if strcmp(valueName, dof.getValueType)
                   vals(dof.getId()) = dof.getValue(step);
               end
            end
        end
        
        function vals = assembleFirstDerivativesVector2(femModel, valueName, step)
            %ASSEMBLEFIRSTDERIVATIVESVECTOR returns a vector with all values
            % corresponding to the first derivative of the elements w.r.t.
            % global dof ids
            dofs = femModel.getDofArray;
            vals = zeros(1,length(dofs));
            for itDof = 1:length(dofs)
               dof = dofs(itDof);
               if strcmp(valueName, dof.getValueType)
                   vals(dof.getId()) = dof.getFirstDerivativeValue(step);
               end
            end
        end
        
        function vals = assembleSecondDerivativesVector2(femModel, valueName, step)
            %ASSEMBLEFIRSTDERIVATIVESVECTOR returns a vector with all values
            % corresponding to the first derivative of the elements w.r.t.
            % global dof ids
            dofs = femModel.getDofArray;
            vals = zeros(1,length(dofs));
            for itDof = 1:length(dofs)
               dof = dofs(itDof);
               if strcmp(valueName, dof.getValueType)
                   vals(dof.getId()) = dof.getSecondDerivativeValue(step);
               end
            end
        end
        
        function vals = assembleFirstDerivativesVector(femModel, step)
            %ASSEMBLEFIRSTDERIVATIVESVECTOR returns a vector with all values
            % corresponding to the first derivative of the elements w.r.t.
            % global dof ids
            dofs = femModel.getDofArray;
            vals = zeros(1,length(dofs));
            elements = femModel.getAllElements;
            for itEle = 1:length(elements)
                element = elements(itEle);
                dofIdList = element.getDofList().getId();
                val = element.getFirstDerivativesVector(step);
                vals(dofIdList) = val;
            end
        end
        
        function vals = assembleSecondDerivativesVector(femModel, step)
            %ASSEMBLESECONDDERIVATIVESVECTOR returns a vector with all values
            % corresponding to the second derivative of the elements w.r.t.
            % global dof ids
            dofs = femModel.getDofArray;
            vals = zeros(1,length(dofs));
            elements = femModel.getAllElements;
            for itEle = 1:length(elements)
                element = elements(itEle);
                dofIdList = element.getDofList().getId();
                val = element.getSecondDerivativesVector(step);
                vals(dofIdList) = val;
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
