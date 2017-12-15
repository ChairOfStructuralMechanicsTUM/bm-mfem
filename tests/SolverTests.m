classdef SolverTests <  matlab.unittest.TestCase
    %SOLVERTESTS Class for solver tests
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Test)
        
        function eigensolverTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            %parameters
            mass = 1.0;
            stiffness = 10.0;
            
            %model
            model = FemModel();
            
            n01 = model.addNewNode(1,1,0,0);
            n02 = model.addNewNode(2,2,0,0);
            base = model.addNewNode(3,0,0,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            s01 = model.addNewElement('SpringDamperElement3d2n', 1, [3 1]);
            s01.setPropertyValue('ELEMENTAL_STIFFNESS', 2*stiffness);
            m01 = model.addNewElement('ConcentratedMassElement3d1n', 2, 1);
            m01.setPropertyValue('ELEMENTAL_MASS', 2*mass);
            s02 = model.addNewElement('SpringDamperElement3d2n', 3, [1 2]);
            s02.setPropertyValue('ELEMENTAL_STIFFNESS', stiffness);
            m02 = model.addNewElement('ConcentratedMassElement3d1n', 4, 2);
            m02.setPropertyValue('ELEMENTAL_MASS', mass);
            
            model.getAllNodes.fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            base.fixDof('DISPLACEMENT_X');
            
            solver = EigensolverStrategy(model);
            solver.solve(2);
            solver.assignModeShapes();
            
            %eigenfrequencies in Hz
            expectedEigenfrequencies = (1 / (2*pi)) .* [sqrt(stiffness/(2*mass)) sqrt(2*stiffness/mass)]';
            actualEigenfrequencies = sort(solver.getEigenfrequencies);
            
            testCase.assertThat(actualEigenfrequencies, IsEqualTo(expectedEigenfrequencies, ...
                    'Within', AbsoluteTolerance(1e-7)))
        end
        
        function newmarkTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            %parameters
            mass = 1.0;
            stiffness = 10.0;
            damping = 2.0;
            uInit = 0.1;
            vInit = 0.0;
            
            %analytical solution
            omega = sqrt(stiffness/mass);
            D = damping / (2 * mass * omega);
            omega_D = omega * sqrt(1 - power(D, 2));
            delta = damping / (2 * mass);
            theta = atan(- (vInit + uInit * delta) / (omega_D * uInit));
            A = sqrt(power(uInit, 2) + power((vInit + uInit * delta) / omega_D, 2));
            
            %model
            model = FemModel();
            
            n01 = model.addNewNode(1,0,1,0);
            n02 = model.addNewNode(2,0,0,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            springEle = model.addNewElement('SpringDamperElement3d2n', 1, [1 2]);
            springEle.setPropertyValue('ELEMENTAL_STIFFNESS', stiffness);
            springEle.setPropertyValue('ELEMENTAL_DAMPING', damping);
            massEle = model.addNewElement('ConcentratedMassElement3d1n', 2, 1);
            massEle.setPropertyValue('ELEMENTAL_MASS', mass);
            
            model.getAllNodes.fixDof('DISPLACEMENT_X');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            n02.fixDof('DISPLACEMENT_Y');
            
            n01.setDofValue('DISPLACEMENT_Y', uInit);
            
            dt = 0.1;
            time = 0;
            endTime = 4;
            solver = NewmarkSolvingStrategy(model, dt);
            while time < endTime
                solver.solve();
                time = time + dt;
                
                actualDisplacementY = n01.getDofValue('DISPLACEMENT_Y','end');
                expectedDisplacementY = A * cos(omega_D * time + theta) * exp(-delta * time);
                testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
                    'Within', AbsoluteTolerance(1e-3)))
                
            end
        end
    end
    
end

