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

            io = GmshInput('validation_bridge_input.msh');
            model = io.readModel;
            
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            model.getModelPart('fixed_support').getNodes.fixAllDofs;
            model.getModelPart('roller_support').getNodes.fixDof('DISPLACEMENT_Y');
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
        
        function twoDofHarmonicAnalysis(testCase)
            %TWODOFHARMONICANALYSIS 2-dof-system with springs and masses;
            %taken from HUMAR: Dynamics for Structures (p. 675)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            model = FemModel();
            support = model.addNewNode(3,0,0,0);
            n01 = model.addNewNode(1,10,0,0);
            model.addNewNode(2,20,0,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            stiffness = 20.0;
            mass = 0.2;
            
            s01 = model.addNewElement('SpringDamperElement3d2n',1,[1 3]);
            s01.setPropertyValue('ELEMENTAL_STIFFNESS', stiffness);
            s02 = model.addNewElement('SpringDamperElement3d2n',2,[1 2]);
            s02.setPropertyValue('ELEMENTAL_STIFFNESS', stiffness/2);
            m01 = model.addNewElement('ConcentratedMassElement3d1n',3,1);
            m01.setPropertyValue('ELEMENTAL_MASS', mass);
            m02 = model.addNewElement('ConcentratedMassElement3d1n',4,2);
            m02.setPropertyValue('ELEMENTAL_MASS', mass/2);
            
            model.getAllNodes.fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            support.fixDof('DISPLACEMENT_X');
            
            addPointLoad(n01,1,[1 0 0]);
            
            exfreq = linspace(.1*sqrt(5),10*sqrt(5),1000);
            
            solver = EigensolverStrategy(model);
            solver.harmonicAnalysis(exfreq,2);
            
            Rd11_actual = model.getNode(1).getDofValue('DISPLACEMENT_X','all');
            Rd12_actual = model.getNode(2).getDofValue('DISPLACEMENT_X','all');
            
            Rd11_expected = (1/3 * stiffness) ./ (0.5 - (exfreq ./ sqrt(stiffness/mass)).^2) ...
                + (1/1.5 * stiffness) ./ (2 - (exfreq ./ sqrt(stiffness/mass)).^2);
            Rd12_expected = (2/3 * stiffness) ./ (0.5 - (exfreq ./ sqrt(stiffness/mass)).^2) ...
                - (2/3 * stiffness) ./ (2 - (exfreq ./ sqrt(stiffness/mass)).^2);
            
            testCase.assertThat(Rd11_actual, IsEqualTo(Rd11_expected / stiffness^2, ...
                'Within', AbsoluteTolerance(1e-7)))
            testCase.assertThat(Rd12_actual, IsEqualTo(Rd12_expected / stiffness^2, ...
                'Within', AbsoluteTolerance(1e-7)))
            
        end
        
        function shellElement3d4nLargeTest(testCase)
            %SHELLELEMENT3D4NLARGETEST plane shell consisting out of 100
            %   elements under static loading at the middle node
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 10, 10, 'elementType', 'ShellElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]);
            
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('THICKNESS', 0.005);
            model.getAllElements.setPropertyValue('DENSITY',7860);
            
            support = [x0 xl y0 yl];
            support.fixAllDofs();
            
            model.getNode(61).setDofLoad('DISPLACEMENT_Z',2500);
            
            solver = SimpleSolvingStrategy(model);
            solver.solve();
            
            actualDisplacement = model.getNode(61).getDofValue('DISPLACEMENT_Z');
            expectedDisplacement = 0.006040637455055775;
            
            testCase.assertThat(actualDisplacement, IsEqualTo(expectedDisplacement, ...
                'Within', RelativeTolerance(1e-7)))
        end
        
        function externalScriptsTest(testCase)
            bridge;
            Bridge_with_inputFile;
        end
        
    end
    
end

