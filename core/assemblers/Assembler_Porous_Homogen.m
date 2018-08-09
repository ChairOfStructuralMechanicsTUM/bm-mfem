classdef Assembler_Porous_Homogen < Assembler
    %Assembler_Porous_Homogen works as SIMPLEASSEMBLER, but is specified
    %for the application on mixed porous and homogen models
    
    properties
    end
    
    methods %test
        
        function assembling = Assembler_Porous_Homogen(femModel)
            if (nargin < 1)
                error('input model is missing');
            end
        end
    end
    
     methods (Static)
         
        function [stiffnessmatrix,reduced_stiffnessmatrix] = assembleGlobalStiffnessMatrix(femModel)
            if nargin == 1
                elements = femModel.getAllElements;
            elseif nargin == 2
                elements = femModel.getModelPart(modelPartName).getElements;
                if isempty(elements)
                    msg = ['SimpleAssembler: Model part ', modelPartName, ...
                        ' has no elements in it.'];
                    e = MException('MATLAB:bm_mfem:emptyModelPart',msg);
                    throw(e);
                end
            end
            
            ndofs = length(femModel.getDofArray);
            K = zeros(ndofs);
            
             for itEle = 1:length(elements)
               elementalStiffnessMatrix = elements(itEle).computeLocalStiffnessMatrix;
               elementalDofIds = elements(itEle).getDofList().getId;
               K(elementalDofIds, elementalDofIds) = ...
               K(elementalDofIds, elementalDofIds) + elementalStiffnessMatrix;
             end
             
            % Determining intersection (common nodes)
            node_array = femModel.getAllNodes;
            mp1 = femModel.getModelPart('Homogen');
            mp2 = femModel.getModelPart('Porous');
            ids1 = mp1.getNodes().getId();
            ids2 = mp2.getNodes().getId();
            coupling = intersect(ids1,ids2);

            % Applying coupling conditions
            delete = zeros(length(coupling),1);
            
            for iNode = 1:length(coupling)
                node_id = coupling(iNode);
                id_homogen_x = node_array(node_id).getDof('DISPLACEMENT_X').getId;
                id_homogen_y = node_array(node_id).getDof('DISPLACEMENT_Y').getId;
                id_frame_x = node_array(node_id).getDof('FRAME_DISPLACEMENT_X').getId;
                id_frame_y = node_array(node_id).getDof('FRAME_DISPLACEMENT_Y').getId;
                id_fluid_y = node_array(node_id).getDof('FLUID_DISPLACEMENT_Y').getId;

                K(id_homogen_x,:)=K(id_homogen_x,:)+K(id_frame_x,:);
                K(:,id_homogen_x)=K(:,id_homogen_x)+K(:,id_frame_x);
                K(id_homogen_y,:)=K(id_homogen_y,:)+K(id_frame_y,:);
                K(:,id_homogen_y)=K(:,id_homogen_y)+K(:,id_frame_y);
                K(id_homogen_y,:)=K(id_homogen_y,:)+K(id_fluid_y,:);
                K(:,id_homogen_y)=K(:,id_homogen_y)+K(:,id_fluid_y);

                delete([3*iNode-2,3*iNode-1,3*iNode]) = [id_frame_x,id_frame_y,id_fluid_y];
            end
            
            % apply fixed dofs
            [~, fixedDofs] = femModel.getDofConstraints;
            fixedDofIds = fixedDofs.getId();

            stiffnessmatrix=K;
            stiffnessmatrix(delete,:) = [];
            stiffnessmatrix(:,delete) = [];
            
            reduced_stiffnessmatrix=K;
            reduced_stiffnessmatrix([delete' fixedDofIds],:) = [];
            reduced_stiffnessmatrix(:,[delete' fixedDofIds]) = [];
            
        end
        
        function [massmatrix,reduced_massmatrix] = assembleGlobalMassMatrix(femModel)
            if nargin == 1
                elements = femModel.getAllElements;
            elseif nargin == 2
                elements = femModel.getModelPart(modelPartName).getElements;
                if isempty(elements)
                    msg = ['SimpleAssembler: Model part ', modelPartName, ...
                        ' has no elements in it.'];
                    e = MException('MATLAB:bm_mfem:emptyModelPart',msg);
                    throw(e);
                end
            end
            
            ndofs = length(femModel.getDofArray);
            M = zeros(ndofs);
            
             for itEle = 1:length(elements)
               elementalMassmatrix = elements(itEle).computeLocalMassMatrix;
               elementalDofIds = elements(itEle).getDofList().getId;
               M(elementalDofIds, elementalDofIds) = ...
               M(elementalDofIds, elementalDofIds) + elementalMassmatrix;
             end
             
            % Determining intersection (common nodes)
            node_array = femModel.getAllNodes;
            mp1 = femModel.getModelPart('Homogen');
            mp2 = femModel.getModelPart('Porous');
            ids1 = mp1.getNodes().getId();
            ids2 = mp2.getNodes().getId();
            coupling = intersect(ids1,ids2);

            % Applying coupling conditions
            delete = zeros(length(coupling),1);
            
            for iNode = 1:length(coupling)
                node_id = coupling(iNode);
                id_homogen_x = node_array(node_id).getDof('DISPLACEMENT_X').getId;
                id_homogen_y = node_array(node_id).getDof('DISPLACEMENT_Y').getId;
                id_frame_x = node_array(node_id).getDof('FRAME_DISPLACEMENT_X').getId;
                id_frame_y = node_array(node_id).getDof('FRAME_DISPLACEMENT_Y').getId;
                id_fluid_y = node_array(node_id).getDof('FLUID_DISPLACEMENT_Y').getId;

                M(id_homogen_x,:)=M(id_homogen_x,:)+M(id_frame_x,:);
                M(:,id_homogen_x)=M(:,id_homogen_x)+M(:,id_frame_x);
                M(id_homogen_y,:)=M(id_homogen_y,:)+M(id_frame_y,:);
                M(:,id_homogen_y)=M(:,id_homogen_y)+M(:,id_frame_y);
                M(id_homogen_y,:)=M(id_homogen_y,:)+M(id_fluid_y,:);
                M(:,id_homogen_y)=M(:,id_homogen_y)+M(:,id_fluid_y);

                delete([3*iNode-2,3*iNode-1,3*iNode]) = [id_frame_x,id_frame_y,id_fluid_y];
            end
            
            % apply fixed dofs
            [~, fixedDofs] = femModel.getDofConstraints;
            fixedDofIds = fixedDofs.getId();

            massmatrix=M;
            massmatrix(delete,:) = [];
            massmatrix(:,delete) = [];
            
            reduced_massmatrix=M;
            reduced_massmatrix([delete' fixedDofIds],:) = [];
            reduced_massmatrix(:,[delete' fixedDofIds]) = [];
            
        end
        
        function [forceVector, reducedForceVector, delete] = applyExternalForces(femModel)
            dofs = femModel.getDofArray;
            nDofs = length(dofs);
            f = zeros(1,nDofs);
            [~, fixedDofs] = femModel.getDofConstraints();
            
            % get the external load on the dofs
            for itDof = 1:nDofs
                if ( ~ dofs(itDof).isFixed)
                    f(itDof) = dofs(itDof).getDofLoad;
                end
            end
            
            % Determining intersection (common nodes)
            node_array = femModel.getAllNodes;
            mp1 = femModel.getModelPart('Homogen');
            mp2 = femModel.getModelPart('Porous');
            ids1 = mp1.getNodes().getId();
            ids2 = mp2.getNodes().getId();
            coupling = intersect(ids1,ids2);

            % Applying coupling conditions
            delete = zeros(length(coupling),1);
            
            for iNode = 1:length(coupling)
                node_id = coupling(iNode);
                id_homogen_x = node_array(node_id).getDof('DISPLACEMENT_X').getId;
                id_homogen_y = node_array(node_id).getDof('DISPLACEMENT_Y').getId;
                id_frame_x = node_array(node_id).getDof('FRAME_DISPLACEMENT_X').getId;
                id_frame_y = node_array(node_id).getDof('FRAME_DISPLACEMENT_Y').getId;
                id_fluid_y = node_array(node_id).getDof('FLUID_DISPLACEMENT_Y').getId;

                f(id_homogen_x)=f(id_homogen_x)+f(id_frame_x);
                f(id_homogen_y)=f(id_homogen_y)+f(id_frame_y);
                f(id_homogen_y)=f(id_homogen_y)+f(id_fluid_y);

                delete([3*iNode-2,3*iNode-1,3*iNode]) = [id_frame_x,id_frame_y,id_fluid_y];
            end
            
            forceVector = f;
            forceVector(delete) = [];
            
            reducedForceVector = f;
            reducedForceVector([delete' fixedDofs.getId()]) = [];
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
end
end           