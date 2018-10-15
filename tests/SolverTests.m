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
            
            %eigenfrequencies in rad/s
            expectedEigenfrequencies = sparse([sqrt(stiffness/(2*mass)) sqrt(2*stiffness/mass)]');
            actualEigenfrequencies = solver.getEigenfrequencies;
            
            testCase.assertThat(actualEigenfrequencies, IsEqualTo(expectedEigenfrequencies, ...
                    'Within', AbsoluteTolerance(1e-7)))
        end
        
        function dampedHarmonicTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            a = .01;
            m = createBeam(.4,10,'IY',a^4/12,'IZ',a^4/12,'IT',a^4/12,'CROSS_SECTION',a^2);
            m.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z','ROTATION_X', 'ROTATION_Y', 'ROTATION_Z'});
            
            m.getNode(1).fixAllDofs;
            m.getNode(11).setDofLoad('DISPLACEMENT_Z',100);
            m.getAllElements.addProperty('RAYLEIGH_ALPHA',4.64);
            m.getAllElements.addProperty('RAYLEIGH_BETA',9.1e-5);
            
            s = EigensolverStrategy(m);
            s.harmonicAnalysis(2*pi*logspace(1,4,100),10);
            
            disp = abs(m.getNode(11).getDofValue('DISPLACEMENT_Z','all'));
            load('tests/test_data.mat','damped_harmonic_disp');
            
            testCase.assertThat(disp, IsEqualTo(damped_harmonic_disp, ...
                'Within', AbsoluteTolerance(1e-5)))
            
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
        
        function PODtest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            %solve POD
            POD_model = SolverTests.create12DofModel;
            PODsolver = MORStrategy(POD_model,'POD','nPODmodes',8);
            coarse_sampling = linspace(1e-4,3.188,41);
            fine_sampling = linspace(0,1,120);
            PODsolver.initialize(coarse_sampling);
            PODsolver.solve(fine_sampling);
            disp = abs(POD_model.getNode(12).getDofValue('DISPLACEMENT_X','all'));
            load('tests/test_data.mat','POD_disp');
            
            testCase.assertThat(disp, IsEqualTo(POD_disp, ...
                    'Within', AbsoluteTolerance(1e-5)))
        end
    end
    
    methods (Access = private, Static)
        function model = create12DofModel
            model = FemModel();
            
            %nodes
            n01 = model.addNewNode(1,0,0,0);
            n=12;
            for id = 1:n
                model.addNewNode(1+id,10*id/n,0,0);
            end
            endnode = model.getNode(n+1);
            
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            %elements
            for id = 1:n
                model.addNewElement('SpringDamperElement3d2n', id, [id id+1]);
            end
            model.getAllElements.setPropertyValue('ELEMENTAL_STIFFNESS',1);
            model.getAllElements.setPropertyValue('ELEMENTAL_DAMPING',0.6);
            
            %masses
            mass_id = length(model.getAllElements()) + 1;
            for id = 1:n
                model.addNewElement('ConcentratedMassElement3d1n',id + n, id+1);
            end
            
            model.getElements(mass_id:length(model.getAllElements)).setPropertyValue('ELEMENTAL_MASS',1);
            
            %boundary conditions
            addPointLoad(endnode,10,[1 0 0]);
            model.getAllNodes.fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            n01.fixAllDofs();
        end
    end
    
end

