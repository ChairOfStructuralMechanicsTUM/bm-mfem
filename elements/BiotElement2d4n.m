classdef BiotElement2d4n < QuadrilateralElement
        
    properties (Access = private)

    end



    
    methods
         
     % Constructor
     function biotElement2d4n = BiotElement2d4n(id,nodeArray)
            
        requiredPropertyNames = cellstr(["OMEGA","NUMBER_GAUSS_POINT", "LAMBDA", ...
                "MU","ETA_S","DENSITY_S","DENSITY_F","ETA_F","PRESSURE_0","HEAT_CAPACITY_RATIO",...
                 "PRANDL_NUMBER","POROSITY","TURTUOSITY","STATIC_FLOW_RESISTIVITY",...
                "VISCOUS_CHARACT_LENGTH","THERMAL_CHARACT_LENGTH"]);
            
        % define the arguments for the super class constructor call
        if nargin == 0
            super_args = {};
        elseif nargin == 2
            if ~(length(nodeArray) == 4 && isa(nodeArray,'Node'))
               error('problem with the nodes in element %d', id);
            end
            super_args = {id, nodeArray, requiredPropertyNames};
        end
            
            % call the super class contructor
            biotElement2d4n@QuadrilateralElement(super_args{:});
            biotElement2d4n.dofNames = cellstr(["DISPLACEMENT_X", ...
                "DISPLACEMENT_Y"]);
        end
    
    
     %Initialization
        function initialize(biotElement2d4n)
            biotElement2d4n.lengthX = computeLength(biotElement2d4n.nodeArray(1).getCoords, ...
                biotElement2d4n.nodeArray(2).getCoords);
            
            biotElement2d4n.lengthY = computeLength(biotElement2d4n.nodeArray(1).getCoords, ...
                biotElement2d4n.nodeArray(4).getCoords);
            
            checkConvexity(biotElement2d4n);
        end
        
       
        function responseDoF = getResponseDofArray(element, step)
            
            responseDoF = zeros(12,1);
            for itNodes = 1:1:4
                nodalDof = element.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.'
                
                
                for itDof = 3:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
        end
    
        
        
            % Shape Function and Derivatives     
        function [N_mat, N, B_b, B_s, J] = computeShapeFunction(biotElement2d4n,xi,eta)

            N = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];
            
            N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4
                -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];
            
            N_mat = sparse(2,8);
            N_mat(1,1:2:end) = N(:);
            N_mat(2,2:2:end) = N(:);
            
            % Coordinates of the nodes forming one element
            ele_coords = zeros(4,2);
            for i=1:4
                ele_coords(i,1) = biotElement2d4n.nodeArray(i).getX;
                ele_coords(i,2) = biotElement2d4n.nodeArray(i).getY;
            end
            
             
            % Jacobian
            J = N_Diff_Par * ele_coords;
            
            % Calculation of B-Matrix
            B=J\N_Diff_Par;%/Jdet;
            Bx=B(1,1:4);
            By=B(2,1:4);
            
            % Calculation of B_bending Matrix
            %B_b = sparse(3,8);
            B_b=[Bx(1),0,Bx(2),0,Bx(3),0,Bx(4),0;
                0,By(1),0,By(2),0,By(3),0,By(4);
                By(1),Bx(1),By(2),Bx(2),By(3),Bx(3),By(4),Bx(4)];
            
            % Calculation of B_shear Matrix
            %B_s = sparse(2,8);
            %If 3 DOF Element: B_s=[Bx(1),N(1),0,Bx(2),N(2),0,Bx(3),N(3),0,Bx(4),N(4),0; By(1),N(1),0,By(2),N(2),0,By(3),N(3),0,By(4),N(4),0];
            B_s=[N(1),0,N(2),0,N(3),0,N(4),0; 
                N(1),0,N(2),0,N(3),0,N(4),0];
        end
        
        
        
%%__________________________________________________________
        
        function K_f = FluidBulkModulus(OMEGA)
        
            HEAT_CAPACITY_RATIO= biotElement2d4n.getPropertyValue('HEAT_CAPACITY_RATIO');    
            PRESSURE_0 = biotElement2d4n.getPropertyValue('PRESSURE_0');
            ETA_F = biotElement2d4n.getPropertyValue('ETA_F');
            PRANDL_NUMBER = biotElement2d4n.getPropertyValue('PRANDL_NUMBER');
            THERMAL_CHARACT_LENGTH = biotElement2d4n.getPropertyValue('THERMAL_CHARACT_LENGTH');
            DENSITY_F = biotElement2d4n.getPropertyValue('DENSITY_F');
         
         
            K_f = (HEAT_CAPACITY_RATIO*PRESSURE_0)/(HEAT_CAPACITY_RATIO-(HEAT_CAPACITY_RATIO-1)*...
               (1+(8*ETA_F)/(i*OMEGA*PRANDL_NUMBER*THERMAL_CHARACT_LENGTH^2*DENSITY_F)*...
               (1+(i*OMEGA*PRANDL_NUMBER*THERMAL_CHARACT_LENGTH^2*DENSITY_F)/(16*ETA_F))^(0.5))^(-1))
        
        end
    



        function [stiffnessMatrix_s, stiffnessMatrix_f] = computeLocalStiffnessMatrix(biotElement2d4n)
        

            ETA_S = biotElement2d4n.getPropertyValue('ETA_S');
            LAMBDA = biotElement2d4n.getPropertyValue('LAMBDA');
            MU = biotElement2d4n.getPropertyValue('MU');
            OMEGA = biotElement2d4n.getPropertyValue('OMEGA');
            POROSITY = biotElement2d4n.getPropertyValue('POROSITY');
            
            LameCoeff = (1+i*ETA_S)*LAMBDA
            ShearModulus = (1+i*ETA_S)*MU
            
            K_f = FluidBulkModulus(OMEGA)
            
            R = POROSITY*K_f
            Q = (1-POROSITY)*K_f
            A = LameCoeff+(1-POROSITY)^2/POROSITY*K_f
            
            
            D_s = [A+2*MU,A,0;A,A+2*MU,0;0,0,MU]
            D_sf = [Q,Q,0;Q,Q,0;0,0,0]
            D_f = [R,R,0;R,R,0;0,0,0]
            

%%__________________________________________________________            
            
 
            %%% get gauss integration factors
            [w,g]=returnGaussPoint(p);
           
            
            %%% empty Stiffnessmatrix
            
            stiffnessMatrix_s=sparse(8,8);
            stiffnessMatrix_f=sparse(8,8);
           
            %%%% WERDEN BENDING- UND SHEARMATRIX BENUTZT? WARUM?
            
            if (biotElement2d4n.getProperties.hasValue('FULL_INTEGRATION')) && (biotElement2d4n.getPropertyValue('FULL_INTEGRATION'))
                for xi = 1 : p
                    for eta = 1 : p
                        [~, ~,B_b, B_s, J] = computeShapeFunction(biotElement2d4n,g(xi),g(eta));

                        stiffnessMatrix_s = stiffnessMatrix_s +             ...
                            B_b' * D_s * B_b * det(J) * w(xi) * w(eta) + ...
                            B_s' * D_sf * B_s * det(J) * w(xi) * w(eta);
                        
                        stiffnessMatrix_f = stiffnessMatrix_f +             ...
                            B_b' * D_f * B_b * det(J) * w(xi) * w(eta) + ...
                            B_s' * D_sf' * B_s * det(J) * w(xi) * w(eta);
                        
                        
                    end
                end
                
            
            else
                for xi = 1 : p
                    for eta = 1 : p
                        [~, ~,B_b, ~, J] = computeShapeFunction(biotElement2d4n,g(xi),g(eta));

                        stiffnessMatrix_s = stiffnessMatrix_s +             ...
                            B_b' * D_s * B_b * det(J) * w(xi) * w(eta) + ...
                            B_b' * D_sf * B_b * det(J) * w(xi) * w(eta);
                        
                        stiffnessMatrix_f = stiffnessMatrix_f +             ...
                            B_b' * D_f * B_b * det(J) * w(xi) * w(eta) + ...
                            B_b' * D_sf' * B_b * det(J) * w(xi) * w(eta);
                 
                       
                    end
                end
                
                [w,g] = returnGaussPoint(1);
              
                for xi = 1 : 1
                    for eta = 1 : 1
                        [~, ~,~, B_s, J] = computeShapeFunction(biotElement2d4n,g(xi),g(eta));
                        
                        stiffnessMatrix_s = stiffnessMatrix_s +             ...
                            B_s' * D_s * B_s * det(J) * w(xi) * w(eta) + ...
                            B_s' * D_sf * B_s * det(J) * w(xi) * w(eta);
                        stiffnessMatrix_f = stiffnessMatrix_f +             ...
                            B_s' * D_f * B_s * det(J) * w(xi) * w(eta) + ...
                            B_s' * D_sf' * B_s * det(J) * w(xi) * w(eta); 
                        
                    end
                end
            end
        end
        
            
        
            
        %%%% WAS IST MIT DER MASSENMATRIX?    
        %%%% WIE KANN SPÄTER NACH SOLID UND FLUID GELOEST WERDEN?
        
        
        function [massMatrix_s, massMatrix_f] = computeLocalMassMatrix(quadrilateralElement2d4n)

            
            massMatrix_s=zeros(8,8);
            massMatrix_f=zeros(8,8);
            
            roh_s = quadrilateralElement2d4n.getPropertyValue('DENSITY_S');
            roh_f = quadrilateralElement2d4n.getPropertyValue('DENSITY_F');
            p = quadrilateralElement2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    
                    [N_mat, ~, ~, ~, J] = computeShapeFunction(biotElement2d4n,xi,eta);
                    
                    massMatrix_s=massMatrix_s + (w(i)*w(j)*roh_s*transpose(N_mat)*N_mat*det(J));
                    
                end
            end
            
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    
                    [N_mat, ~, ~, ~, J] = computeShapeFunction(biotElement2d4n,xi,eta);
                    
                    massMatrix_f=massMatrix_f + (w(i)*w(j)*roh_f*transpose(N_mat)*N_mat*det(J));
                    
                end
            end
            
        end
        
        
        
        
        %%%% FREIHEITSGRADE - WIE DEFINIERT?
        %%%% MÜSSEN BC ZWISCHEN SOLID UND FLUID BZW PORO DEFINIERT WERDEN?
        
        function dofs = getDofList(element)
            dofs([1 3 5 7]) = element.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 4 6 8]) = element.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,8);
            
            vals([1 3 5 7]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 4 6 8]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
        end
             
             
        
        
        %%%% WO MUSS ICH DIE DISPLACEMENT BASED FORMULATION EINBAUEN? 
        %%%% UND WOZU?
    
    
    end
end