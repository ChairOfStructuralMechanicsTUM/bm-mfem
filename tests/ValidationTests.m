classdef ValidationTests <  matlab.unittest.TestCase
    %VALIDATIONTESTS Class for larger tests
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Test)
        
        % bridge test taken from IFEM ch. 21 by Felippa
        function bridgeTest(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance

            io = ModelIO('validation_bridge_input.msh');
            model = io.readModel;
            
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            model.getModelPart('fixed_support').fixDof('DISPLACEMENT_X');
            model.getModelPart('fixed_support').fixDof('DISPLACEMENT_Y');
            model.getModelPart('roller_support').fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            
            addPointLoad(model.getNodes([3 5 9 11]),10,[0 -1 0]);
            addPointLoad(model.getNode(7),16,[0 -1 0]);
            
            solver = SimpleSolvingStrategy(model);
            solver.solve();
            
            actualDisplacementX = model.getAllNodes.getDofValue('DISPLACEMENT_X');
            actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');
            
            expectedDisplacementX = [0 0.809536 0.28 0.899001 0.56 0.8475 ...
                0.8475 0.795999 1.135 0.885464 1.415 1.695]';
            expectedDisplacementY = [0 -1.775600 -1.792260 -2.291930 -2.316600 ...
                -2.385940 -2.421940 -2.291930 -2.316600 -1.775600 -1.792260 0]';
            
            testCase.assertThat(actualDisplacementX, IsEqualTo(expectedDisplacementX, ...
                'Within', RelativeTolerance(1e-5)))
            testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
                'Within', RelativeTolerance(1e-5)))
            
            actualElementStress = model.getAllElements.computeElementStress(1);
            expectedElementStress = [28 28 28.75 28.75 28 28 -6.2610 -6.0030 ...
                -6.0300 -6.0300 -6.0030 -6.2610 3.3330 3.0830 4.0000 3.0830 ...
                3.3330 1.6770 3.2020 3.2020 1.6770]';
            
            testCase.assertThat(actualElementStress, IsEqualTo(expectedElementStress, ...
                'Within', RelativeTolerance(1e-3)))
        end
        
        function sdofWithHarmonicExcitation(testCase)
            %SDOFWITHHARMONICEXCITATION sdof spring mass system with
            %harmonic excitiation
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            %parameters
            mass = 80.0;
            stiffness = 200.0;
            F0 = 100.0;
            uInit = 0.0;
            vInit = 0.0;
            
            %analytical solution
            eigenfrequency = sqrt(stiffness/mass);
            excitationFrequency = 0.5 * eigenfrequency;
            beta = excitationFrequency / eigenfrequency;
            
            model = FemModel();
            n01 = model.addNewNode(1,0,0,0);
            n02 = model.addNewNode(2,0,1,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            s = model.addNewElement('SpringDamperElement3d2n',1,[1 2]);
            s.setPropertyValue('ELEMENTAL_STIFFNESS', stiffness);
            m = model.addNewElement('ConcentratedMassElement3d1n',2,2);
            m.setPropertyValue('ELEMENTAL_MASS',mass);
            
            model.getAllNodes.fixDof('DISPLACEMENT_X');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            n01.fixDof('DISPLACEMENT_Y');
            
            dt = 0.05;
            time = 0;
            endTime = 20;
            solver = NewmarkSolvingStrategy(model, dt);
            while time < endTime
                applyHarmonicSineLoad(n02, F0, [0, -1, 0], excitationFrequency, time);
                solver.solve();
                
                actualDisplacementY = n02.getDofValue('DISPLACEMENT_Y','end');
                transientY = uInit * cos(eigenfrequency * time) + vInit / eigenfrequency * ...
                    sin(eigenfrequency * time) - F0 / stiffness * beta / (1 - power(beta,2)) * ...
                    sin(eigenfrequency * time);
                steadyY = F0 / stiffness / (1 - power(beta, 2)) * sin(excitationFrequency * time);
                expectedDisplacementY = - (transientY + steadyY);
                testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
                    'Within', AbsoluteTolerance(1e-2)))
                
                time = time + dt;
            end
        end
        
        function mdofWithDamping(testCase)
            %MDOFWITHDAMPING 2-dof-system with springs, masses, and dampers;
            %taken from HUMAR: Dynamics for Structures (p. 560)
            
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.RelativeTolerance
            tolerance = AbsoluteTolerance(0.07) | RelativeTolerance(0.07);
                        
            model = FemModel();
            n01 = model.addNewNode(1,1,0,0);
            n02 = model.addNewNode(2,2,0,0);
            base = model.addNewNode(3,0,0,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            s01 = model.addNewElement('SpringDamperElement3d2n',1,[3 1]);
            s01.setPropertyValue('ELEMENTAL_STIFFNESS', 2.0);
            s01.setPropertyValue('ELEMENTAL_DAMPING', 0.4);
            m01 = model.addNewElement('ConcentratedMassElement3d1n',2,1);
            m01.setPropertyValue('ELEMENTAL_MASS', 2.0);
            s02 = model.addNewElement('SpringDamperElement3d2n',3,[1 2]);
            s02.setPropertyValue('ELEMENTAL_STIFFNESS', 1.0);
            s02.setPropertyValue('ELEMENTAL_DAMPING', 0.05);
            m02 = model.addNewElement('ConcentratedMassElement3d1n',4,2);
            m02.setPropertyValue('ELEMENTAL_MASS', 1.0);
            s03 = model.addNewElement('SpringDamperElement3d2n',5,[3 2]);
            s03.setPropertyValue('ELEMENTAL_DAMPING', 0.15);
            
            model.getAllNodes.fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            base.fixDof('DISPLACEMENT_X');
            
            n01.setDofValue('DISPLACEMENT_X', 1.0);
            n02.setDofValue('DISPLACEMENT_X', 2.0);
            
            dt = 0.05;
            time = 0;
            endTime = 20;
            solver = NewmarkSolvingStrategy(model, dt);
            while time < endTime
                solver.solve();
                
                actualDisplacementX01 = n01.getDofValue('DISPLACEMENT_X','end');
                actualDisplacementX02 = n02.getDofValue('DISPLACEMENT_X','end');
                
                expectedDisplacementX01 = 1.0064 * exp(-0.08334 * time) * sin(0.70221 * time + atan(9.031));
                expectedDisplacementX02 = 2.0148 * exp(-0.08334 * time) * sin(0.70221 * time + atan(8.153));
                
                testCase.assertThat(actualDisplacementX01, IsEqualTo(expectedDisplacementX01, ...
                    'Within', tolerance))
                testCase.assertThat(actualDisplacementX02, IsEqualTo(expectedDisplacementX02, ...
                    'Within', tolerance))
                
                time = time + dt;
            end
        end
        
        function externalScriptsTest(testCase)
           bridge;
           Bridge_with_inputFile;
        end
        
    end
    
end

