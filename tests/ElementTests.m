classdef ElementTests < matlab.unittest.TestCase
    %ELEMENTTESTS Tests for all elements
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Test)
        
        function testBarElement2d2n(testCase)
            node01 = Node(1,0,0);
            node02 = Node(2,30,40);
            ele01 = BarElement2d2n(1,[node01 node02]);
            ele01.setPropertyValue('CROSS_SECTION',5);
            ele01.setPropertyValue('YOUNGS_MODULUS',1000);
            
            actualSolution = ele01.computeLocalStiffnessMatrix;
            expectedSolution = [36 48 -36 -48; ...
                48 64 -48 -64; ...
                -36 -48 36 48; ...
                -48 -64 48 64];
            testCase.verifyEqual(actualSolution, expectedSolution)
        end
        
        function testBarElement3d2n(testCase)
            node01 = Node(1,0,0,0);
            node02 = Node(2,2,3,6);
            ele01 = BarElement3d2n(1,[node01 node02]);
            ele01.setPropertyValue('CROSS_SECTION',10);
            ele01.setPropertyValue('YOUNGS_MODULUS',343);
            
            actualSolution = ele01.computeLocalStiffnessMatrix;
            expectedSolution = [40 60 120 -40 -60 -120; ...
                60 90 180 -60 -90 -180; ...
                120 180 360 -120 -180 -360; ...
                -40 -60 -120 40 60 120; ...
                -60 -90 -180 60 90 180; ...
                -120 -180 -360 120 180 360];
            testCase.verifyEqual(actualSolution, expectedSolution)
        end
        
        function testSpringDamperElement3d2n(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            model = FemModel();
            model.addNewNode(1,0,0,0);
            model.addNewNode(2,1,1,0);
            model.addNewNode(3,2,0,0);
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            model.addNewElement('SpringDamperElement3d2n',1,[1 2]);
            model.addNewElement('SpringDamperElement3d2n',2,[2 3]);
            model.getAllElements().setPropertyValue('ELEMENTAL_STIFFNESS',100);
            model.getAllElements().setPropertyValue('ELEMENTAL_DAMPING',[1,20,1]);
            
            addPointLoad(model.getNode(2),10,[0 -1 0]);
            
            model.getNode(1).fixDof('DISPLACEMENT_X');
            model.getNode(1).fixDof('DISPLACEMENT_Y');
            model.getNode(3).fixDof('DISPLACEMENT_X');
            model.getNode(3).fixDof('DISPLACEMENT_Y');
            model.getAllNodes.fixDof('DISPLACEMENT_Z');
            
            solver = SimpleSolvingStrategy(model);
            solver.solve();
            
            actualDisplacementX = model.getAllNodes.getDofValue('DISPLACEMENT_X');
            actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');
            actualDisplacementZ = model.getAllNodes.getDofValue('DISPLACEMENT_Z');
            expectedDisplacementX = [0 0 0]';
            expectedDisplacementY = [0 -0.1 0]';
            expectedDisplacementZ = [0 0 0]';
            
            testCase.assertThat(actualDisplacementX, IsEqualTo(expectedDisplacementX, ...
                'Within', RelativeTolerance(1e-7)))
            testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
                'Within', RelativeTolerance(1e-7)))
            testCase.assertThat(actualDisplacementZ, IsEqualTo(expectedDisplacementZ, ...
                'Within', RelativeTolerance(1e-7)))
        end
        
        function testBeamElement3d2nStatic(testCase)
            warning('off','all')
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            model = FemModel();
            
            %nodes
            n01 = model.addNewNode(1,0,0,0);
            for id = 1:4
                model.addNewNode(1+id,id,0,0);
            end
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z', ...
                'ROTATION_X', 'ROTATION_Y', 'ROTATION_Z'});
            
            %elements
            for id = 1:4
                model.addNewElement('BeamElement3d2n', id, [id id+1]);
            end
            model.getAllElements.setPropertyValue('IY',1000);
            model.getAllElements.setPropertyValue('IZ',100);
            model.getAllElements.setPropertyValue('IT',100);
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS',3000);
            model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
            model.getAllElements.setPropertyValue('CROSS_SECTION',3000);
            
            %boundary conditions
            n01.fixAllDofs();
            addPointLoad(model.getNode(5),10,[0 0 -1]);
            
            solver = SimpleSolvingStrategy(model);
            solver.solve();
            
            actualDisplacementX = model.getAllNodes.getDofValue('DISPLACEMENT_X');
            actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');
            actualDisplacementZ = model.getAllNodes.getDofValue('DISPLACEMENT_Z');
            actualRotationX = model.getAllNodes.getDofValue('ROTATION_X');
            actualRotationY = model.getAllNodes.getDofValue('ROTATION_Y');
            actualRotationZ = model.getAllNodes.getDofValue('ROTATION_Z');
            
            expectedDisplacementX = zeros(5,1);
            expectedDisplacementY = zeros(5,1);
            expectedDisplacementZ = [0;-6.11111111111115e-06;-2.22222222222223e-05;-4.50000000000002e-05;-7.11111111111115e-05];
            expectedRotationX = zeros(5,1);
            expectedRotationY = [0;1.16666666666667e-05;2.00000000000001e-05;2.50000000000001e-05;2.66666666666668e-05];
            expectedRotationZ = zeros(5,1);
            
            testCase.assertThat(actualDisplacementX, IsEqualTo(expectedDisplacementX, ...
                'Within', RelativeTolerance(1e-7)))
            testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
                'Within', RelativeTolerance(1e-7)))
            testCase.assertThat(actualDisplacementZ, IsEqualTo(expectedDisplacementZ, ...
                'Within', RelativeTolerance(1e-7)))
            testCase.assertThat(actualRotationX, IsEqualTo(expectedRotationX, ...
                'Within', RelativeTolerance(1e-7)))
            testCase.assertThat(actualRotationY, IsEqualTo(expectedRotationY, ...
                'Within', RelativeTolerance(1e-7)))
            testCase.assertThat(actualRotationZ, IsEqualTo(expectedRotationZ, ...
                'Within', RelativeTolerance(1e-7)))
            
            warning('on','all')
        end
        
%         function testBeamElement3d2nDynamic(testCase)
%             import matlab.unittest.constraints.IsEqualTo
%             import matlab.unittest.constraints.RelativeTolerance
%             model = FemModel();
%             
%             %nodes
%             n01 = model.addNewNode(1,0,0,0);
%             for id = 1:4
%                 model.addNewNode(1+id,id,0,0);
%             end
%             model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z', ...
%                 'ROTATION_X', 'ROTATION_Y', 'ROTATION_Z'});
%             
%             %elements
%             for id = 1:4
%                 model.addNewElement('BeamElement3d2n', id, [id id+1]);
%             end
%             model.getAllElements.setPropertyValue('IY',1000);
%             model.getAllElements.setPropertyValue('IZ',100);
%             model.getAllElements.setPropertyValue('IT',100);
%             model.getAllElements.setPropertyValue('YOUNGS_MODULUS',3000);
%             model.getAllElements.setPropertyValue('SHEAR_MODULUS',3000);
%             model.getAllElements.setPropertyValue('CROSS_SECTION',3000);
%             
%             %boundary conditions
%             n01.fixAllDofs();
%             addPointLoad(model.getNode(5),10,[0 0 -1]);
%             
%             dt = .05;
%             time = 0;
%             endTime = 1.5;
%             ls = linspace(time,endTime,endTime/dt+1);
%             % disp = zeros(1,endTime/dt+1);
%             solver = NewmarkSolvingStrategy(model, dt);
%             % disp = 0;
%             while time < endTime
%                 solver.solve();
%                 time = time + dt;
%                 %     v.plotDeformed();
%                 %     pause(0.1)
%                 %     plot(model.getAllNodes.getDofValue('DISPLACEMENT_Z','end'))
%                 
%                 %     disp(model.getProperties.getValue('STEP')) = n05.getDofValue('DISPLACEMENT_Z','end');
%                 
%             end
%             plot(ls,endnode.getDofValue('DISPLACEMENT_Y','all'))
%             solver = SimpleSolvingStrategy(model);
%             solver.solve();
%             
%             actualDisplacementX = model.getAllNodes.getDofValue('DISPLACEMENT_X');
%             actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');
%             actualDisplacementZ = model.getAllNodes.getDofValue('DISPLACEMENT_Z');
%             
%             expectedDisplacementX = zeros(5,1);
%             expectedDisplacementY = zeros(5,1);
%             expectedDisplacementZ = [0;-6.11111111111115e-06;-2.22222222222223e-05;-4.50000000000002e-05;-7.11111111111115e-05];
%             
%             testCase.assertThat(actualDisplacementX, IsEqualTo(expectedDisplacementX, ...
%                 'Within', RelativeTolerance(1e-7)))
%             testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
%                 'Within', RelativeTolerance(1e-7)))
%             testCase.assertThat(actualDisplacementZ, IsEqualTo(expectedDisplacementZ, ...
%                 'Within', RelativeTolerance(1e-7)))
%         end
        
        function testBeamElement3d2nEigen(testCase)
            warning('off','all')
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            model = FemModel();
            
            %nodes
            n01 = model.addNewNode(1,0,0,0);
            n=10;
            for id = 1:n
                model.addNewNode(1+id,10*id/n,0,0);
            end
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z', ...
                'ROTATION_X', 'ROTATION_Y', 'ROTATION_Z'});
            
            %elements
            for id = 1:n
                model.addNewElement('BeamElement3d2n', id, [id id+1]);
            end
            model.getAllElements.setPropertyValue('IY',0.00001);
            model.getAllElements.setPropertyValue('IZ',0.00001);
            model.getAllElements.setPropertyValue('IT',0.00001);
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS',2.1e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
            model.getAllElements.setPropertyValue('CROSS_SECTION',.01);
            model.getAllElements.setPropertyValue('DENSITY',7850);
            
            %boundary conditions - cantilever
            model.getAllNodes().fixAllDofs();
            model.getAllNodes().unfixDof('DISPLACEMENT_Y');
            model.getAllNodes().unfixDof('ROTATION_Z');
            n01.fixAllDofs();
            
            %solver
            solver = EigensolverStrategy(model);
            solver.solve(5);
            
            %eigenfrequencies in Hz
            actualEigenfrequenciesCantilever = solver.getEigenfrequencies ./ (2*pi);
            expectedEigenfrequenciesCantilever = [0.915242085421142;5.73511288165053;16.0584347115302;31.4797553897702;52.0975276397135];

            testCase.assertThat(actualEigenfrequenciesCantilever, IsEqualTo(expectedEigenfrequenciesCantilever, ...
                'Within', AbsoluteTolerance(1e-5)))
            
            warning('on','all')
        end
        
        function testReissnerMindlinElement3d4nStatic (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 4, 4, ...
                'elementType', 'ReissnerMindlinElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_Z", "ROTATION_X", "ROTATION_Y"]);
            model.getNode(13).setDofLoad('DISPLACEMENT_Z', -0.25);
            support = [x0 xl y0 yl];
            support.fixAllDofs();

            % Properties
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 10920);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('THICKNESS', 0.1);
            model.getAllElements.setPropertyValue('DENSITY', 1);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT', 4);
            model.getAllElements.setPropertyValue('SHEAR_CORRECTION_FACTOR', 5/6);
            
            % Solver
            solver = SimpleSolvingStrategy(model);
            solver.solve();

            actualDisplacementZ = model.getAllNodes.getDofValue('DISPLACEMENT_Z');
            actualRotationX = model.getAllNodes.getDofValue('ROTATION_X');
            actualRotationY = model.getAllNodes.getDofValue('ROTATION_Y');
            
            expectedDisplacementZ = [0;0;0;0;0;0;-0.000422362795098927;-0.000613870097726914;...
                        -0.000422362795098927;0;0;-0.000613870097726913;-0.00168945118039571;...
                        -0.000613870097726913;0;0;-0.000422362795098926;-0.000613870097726913;...
                        -0.000422362795098927;0;0;0;0;0;0];
            expectedRotationX = [0;0;0;0;0;0;0.00254235733839925;4.40124127619258e-18; ...
                        -0.00254235733839923;0;0;0.00473720688683214;1.18261027682316e-17;...
                        -0.00473720688683216;0;0;0.00254235733839923;1.73756544264676e-18;...
                        -0.00254235733839923;0;0;0;0;0;0];
            expectedRotationY = [0;0;0;0;0;0;0.00254235733839925;0.00473720688683215;...
                        0.00254235733839921;0;0;-1.09535396625964e-17;1.06625639192173e-18;...
                        9.26364511108020e-19;0;0;-0.00254235733839922;-0.00473720688683216;...
                        -0.00254235733839923;0;0;0;0;0;0];
            
            testCase.assertThat(actualDisplacementZ, IsEqualTo(expectedDisplacementZ, ...
                'Within', AbsoluteTolerance(1e-12)))
            testCase.assertThat(actualRotationX, IsEqualTo(expectedRotationX, ...
                'Within', AbsoluteTolerance(1e-12)))
            testCase.assertThat(actualRotationY, IsEqualTo(expectedRotationY, ...
                'Within', AbsoluteTolerance(1e-12)))
        end
        
        function testReissnerMindlinElement3d4nDynamic (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 4, 4, ...
                'elementType', 'ReissnerMindlinElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_Z", "ROTATION_X", "ROTATION_Y"]);
            support = [x0 xl y0 yl];
            support.fixDof('DISPLACEMENT_Z');
            
            % Properties
            model.getAllElements.setPropertyValue('THICKNESS', .00125);
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.07e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT', 4);
            model.getAllElements.setPropertyValue('DENSITY', 7850);
            model.getAllElements.setPropertyValue('SHEAR_CORRECTION_FACTOR', 5/6);
            
            % Solver
            model.getNode(13).setDofLoad('DISPLACEMENT_Z', 1);
            
            dt = .01;
            time = 0;
            endTime = .2;
            solver = NewmarkSolvingStrategy(model, dt);
            
            while time < endTime
                solver.solve();
                time = time + dt;                
            end
            
            actualDisplacementZ = model.getNode(13).getDofValue('DISPLACEMENT_Z','all');
            
            expectedDisplacementZ = [0 5.63499025462025e-05 0.000110356168474448 ...
                0.000195291413790720 0.000348022151124683 0.000436442533214557 ...
                0.000549911355006482 0.000638983074213174 0.000633900627028416 ...
                0.000646184028790807 0.000576985852794556 0.000460890315653250 ...
                0.000377492675821615 0.000228575313367927 0.000127279605314402 ...
                7.17349829089917e-05 3.30152776432233e-06 4.13423283962989e-05 ...
                9.48520255923515e-05 0.000166360041372910 0.000315266268863059];
            
            testCase.assertThat(actualDisplacementZ, IsEqualTo(expectedDisplacementZ, ...
                'Within', RelativeTolerance(1e-7)))
        end
        
        function testReissnerMindlinElement3d4nEigen (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
        
            
            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 4, 4, ...
                'elementType', 'ReissnerMindlinElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_Z", "ROTATION_X", "ROTATION_Y"]);
            support = [x0 xl y0 yl];
            support.fixDof('DISPLACEMENT_Z');
            
            % Properties
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 10920);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('THICKNESS', 0.01);
            model.getAllElements.setPropertyValue('DENSITY', 1);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT', 4);
            model.getAllElements.setPropertyValue('SHEAR_CORRECTION_FACTOR', 5/6);

            
            % Solver            
            solver = EigensolverStrategy(model);
            solver.solve(5);
            
            %eigenfrequencies in Hz
            actualEigenfrequencies = solver.getEigenfrequencies('Hz');
            expectedEigenfrequencies = [1.01402840311520;3.22165476249304;3.22165476249320;...
                5.32248978160243;10.5674305646586];
            
            testCase.assertThat(actualEigenfrequencies, IsEqualTo(expectedEigenfrequencies, ...
                'Within', RelativeTolerance(1e-7)))
            
        end
        
        function testQuadrilateralElement2d4n (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            node01 = Node(1,0,0);
            node02 = Node(2,0.5,0);
            node03 = Node(3,0.5,0.25);
            node04 = Node(4,0,0.25);
            
            nodeArray = [node01 node02 node03 node04];
            
            nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
            
            ele01 = QuadrilateralElement2d4n(1,[node01 node02 node03 node04]);
            
            elementArray = [ele01];
            
            elementArray.setPropertyValue('YOUNGS_MODULUS',96);
            elementArray.setPropertyValue('POISSON_RATIO',1/3);
            elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);
            elementArray.setPropertyValue('DENSITY',7860);
            
            actualSolution = ele01.computeLocalStiffnessMatrix;
            
            expectedSolution = [42 18 -6 0 -21 -18 -15 0; ...
                18 78 0 30 -18 -39 0 -69; ...
                -6 0 42 -18 -15 0 -21 18; ...
                0 30 -18 78 0 -69 18 -39; ...
                -21 -18 -15 0 42 18 -6 0; ...
                -18 -39 0 -69 18 78 0 30; ...
                -15 0 -21 18 -6 0 42 -18; ...
                0 -69 18 -39 0 30 -18 78];
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within', AbsoluteTolerance(1e-7)))
            
        end
        
        function testHexahedronElement3d8n (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            node01 = Node(1,-1,-1,-1);
            node02 = Node(2,1,-1,-1);
            node03 = Node(3,1,1,-1);
            node04 = Node(4,-1,1,-1);
            node05 = Node(5,-1,-1,1);
            node06 = Node(6,1,-1,1);
            node07 = Node(7,1,1,1);
            node08 = Node(8,-1,1,1);
            
            nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08];
            
            nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            ele01 = HexahedronElement3d8n(1,[node01 node02 node03 node04 node05 node06 node07 node08]);
            
            ele01.setPropertyValue('YOUNGS_MODULUS',32);
            ele01.setPropertyValue('POISSON_RATIO',1/3);
            ele01.setPropertyValue('NUMBER_GAUSS_POINT',3);
            ele01.setPropertyValue('DENSITY',7860);

            actualSolution = ele01.computeLocalStiffnessMatrix;
            
            expectedSolution = [16,6,6,-8,2,2,-6,-6,1,4,-2,3,4,3,-2,-6,1,-6,-4,-3,-3,0,-1,-1;
                6,16,6,-2,4,3,-6,-6,1,2,-8,2,3,4,-2,-1,0,-1,-3,-4,-3,1,-6,-6;
                6,6,16,-2,3,4,-1,-1,0,3,-2,4,2,2,-8,-6,1,-6,-3,-3,-4,1,-6,-6;
                -8,-2,-2,16,-6,-6,4,2,-3,-6,6,-1,-6,-1,6,4,-3,2,0,1,1,-4,3,3;
                2,4,3,-6,16,6,-2,-8,2,6,-6,1,1,0,-1,-3,4,-2,-1,-6,-6,3,-4,-3;
                2,3,4,-6,6,16,-3,-2,4,1,-1,0,6,1,-6,-2,2,-8,-1,-6,-6,3,-3,-4;
                -6,-6,-1,4,-2,-3,16,6,-6,-8,2,-2,-4,-3,3,0,-1,1,4,3,2,-6,1,6;
                -6,-6,-1,2,-8,-2,6,16,-6,-2,4,-3,-3,-4,3,1,-6,6,3,4,2,-1,0,1;
                1,1,0,-3,2,4,-6,-6,16,2,-3,4,3,3,-4,-1,6,-6,-2,-2,-8,6,-1,-6;
                4,2,3,-6,6,1,-8,-2,2,16,-6,6,0,1,-1,-4,3,-3,-6,-1,-6,4,-3,-2;
                -2,-8,-2,6,-6,-1,2,4,-3,-6,16,-6,-1,-6,6,3,-4,3,1,0,1,-3,4,2;
                3,2,4,-1,1,0,-2,-3,4,6,-6,16,1,6,-6,-3,3,-4,-6,-1,-6,2,-2,-8;
                4,3,2,-6,1,6,-4,-3,3,0,-1,1,16,6,-6,-8,2,-2,-6,-6,-1,4,-2,-3;
                3,4,2,-1,0,1,-3,-4,3,1,-6,6,6,16,-6,-2,4,-3,-6,-6,-1,2,-8,-2;
                -2,-2,-8,6,-1,-6,3,3,-4,-1,6,-6,-6,-6,16,2,-3,4,1,1,0,-3,2,4;
                -6,-1,-6,4,-3,-2,0,1,-1,-4,3,-3,-8,-2,2,16,-6,6,4,2,3,-6,6,1;
                1,0,1,-3,4,2,-1,-6,6,3,-4,3,2,4,-3,-6,16,-6,-2,-8,-2,6,-6,-1;
                -6,-1,-6,2,-2,-8,1,6,-6,-3,3,-4,-2,-3,4,6,-6,16,3,2,4,-1,1,0;
                -4,-3,-3,0,-1,-1,4,3,-2,-6,1,-6,-6,-6,1,4,-2,3,16,6,6,-8,2,2;
                -3,-4,-3,1,-6,-6,3,4,-2,-1,0,-1,-6,-6,1,2,-8,2,6,16,6,-2,4,3;
                -3,-3,-4,1,-6,-6,2,2,-8,-6,1,-6,-1,-1,0,3,-2,4,6,6,16,-2,3,4;
                0,1,1,-4,3,3,-6,-1,6,4,-3,2,4,2,-3,-6,6,-1,-8,-2,-2,16,-6,-6;
                -1,-6,-6,3,-4,-3,1,0,-1,-3,4,-2,-2,-8,2,6,-6,1,2,4,3,-6,16,6;
                -1,-6,-6,3,-3,-4,6,1,-6,-2,2,-8,-3,-2,4,1,-1,0,2,3,4,-6,6,16];
            
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within',AbsoluteTolerance(1e-7)))
        end
        
        
        function testTetrahedronElement3d8n (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            node01 = Node(1,0,0,0);
            node02 = Node(2,0,1,0);
            node03 = Node(3,0,1,1);
            node04 = Node(4,1,1,0);
            
            nodeArray = [node01 node02 node03 node04];
            
            nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
            
            ele01 = TetrahedronElement3d4n(1,[node01 node02 node03 node04]);
            
            ele01.setPropertyValue('YOUNGS_MODULUS',336);
            ele01.setPropertyValue('POISSON_RATIO',1/3);
            ele01.setPropertyValue('DENSITY',2699);

            actualSolution = ele01.computeLocalStiffnessMatrix;
            
            expectedSolution = [ 21,  0,  0,-21, 21,  0,  0,  0,  0,  0,-21,  0;
                                  0, 84,  0, 42,-84, 42,  0,  0,-42,-42,  0,  0;
                                  0,  0, 21,  0, 21,-21,  0,-21,  0,  0,  0,  0;
                                -21, 42,  0,126,-63, 63,-21,  0,-42,-84, 21,-21;
                                 21,-84, 21,-63,126,-63,  0,-21, 42, 42,-21,  0;
                                  0, 42,-21, 63,-63,126,-21, 21,-84,-42,  0,-21;
                                  0,  0,  0,-21,  0,-21, 21,  0,  0,  0,  0, 21;
                                  0,  0,-21,  0,-21, 21,  0, 21,  0,  0,  0,  0;
                                  0,-42,  0,-42, 42,-84,  0,  0, 84, 42,  0,  0;
                                  0,-42,  0,-84, 42,-42,  0,  0, 42, 84,  0,  0;
                                -21,  0,  0, 21,-21,  0,  0,  0,  0,  0, 21,  0;
                                  0,  0,  0,-21,  0,-21, 21,  0,  0,  0,  0, 21];
            
            testCase.assertThat(actualSolution, IsEqualTo(expectedSolution, ...
                'Within',AbsoluteTolerance(1e-7)))
        end
        
        
        function testShellElement3d4nStatic (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            % Model
            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 4, 4, ...
                'elementType', 'ShellElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]);
            addPointLoad(model.getNode(13), 2500, [1 1 2]);
            support = [x0 xl y0 yl];
            support.fixAllDofs();

            % Properties
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('THICKNESS', 0.005);
            model.getAllElements.setPropertyValue('DENSITY',7860);
            
            % Solver
            solver = SimpleSolvingStrategy(model);
            solver.solve();

            % Assertion
            actualDisplacementX = model.getAllNodes.getDofValue('DISPLACEMENT_X');
            actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');
            actualDisplacementZ = model.getAllNodes.getDofValue('DISPLACEMENT_Z');
            actualRotationX = model.getAllNodes.getDofValue('ROTATION_X');
            actualRotationY = model.getAllNodes.getDofValue('ROTATION_Y');
            actualRotationZ = model.getAllNodes.getDofValue('ROTATION_Z');
            
            expectedDisplacementX = [0;0;0;0;0;0;1.33489988272789e-07;...
                9.92693537645756e-08;3.51417111734869e-08;0;0;2.22863312517712e-07;...
                6.21520913991292e-07;2.22863312517712e-07;0;0;3.51417111734869e-08;...
                9.92693537645757e-08;1.33489988272789e-07;0;0;0;0;0;0];
            expectedDisplacementY = [0;0;0;0;0;0;1.33489988272789e-07;...
                2.22863312517712e-07;3.51417111734869e-08;0;0;9.92693537645757e-08;...
                6.21520913991292e-07;9.92693537645757e-08;0;0;3.51417111734870e-08;...
                2.22863312517712e-07;1.33489988272789e-07;0;0;0;0;0;0];
            expectedDisplacementZ = [0;0;0;0;0;0;0.00122309403452578;...
                0.00237730317015391;0.00122309403452578;0;0;0.00237730317015392;...
                0.00544345532928431;0.00237730317015391;0;0;0.00122309403452578;...
                0.00237730317015391;0.00122309403452578;0;0;0;0;0;0];
            expectedRotationX = [0;0;0;0;0;0;0.00720847527292095;...
                0.0160976327860242;0.00720847527292093;0;0;3.94752267012529e-18;...
                6.26418749210311e-18;-1.17566815306664e-17;0;0;-0.00720847527292092;...
                -0.0160976327860242;-0.00720847527292091;0;0;0;0;0;0];
            expectedRotationY = [0;0;0;0;0;0;-0.00720847527292091;...
                -1.73204645887090e-17;0.00720847527292092;0;0;...
                -0.0160976327860242;4.24684987557193e-17;0.0160976327860242;...
                0;0;-0.00720847527292093;1.09450137463053e-18;0.00720847527292091;...
                0;0;0;0;0;0];
            expectedRotationZ = [0;0;0;0;0;0;1.56112458654790e-22;...
                -1.34759401859883e-07;2.53688757107511e-09;0;0;1.34759401859884e-07;...
                1.45250097913107e-22;-1.34759401859883e-07;0;0;-2.53688757107498e-09;...
                1.34759401859884e-07;9.12139527139402e-23;0;0;0;0;0;0];
            
            testCase.assertThat(actualDisplacementX, IsEqualTo(expectedDisplacementX, ...
                'Within', AbsoluteTolerance(1e-12)))
            testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
                'Within', AbsoluteTolerance(1e-12)))
            testCase.assertThat(actualDisplacementZ, IsEqualTo(expectedDisplacementZ, ...
                'Within', AbsoluteTolerance(1e-12)))
            testCase.assertThat(actualRotationX, IsEqualTo(expectedRotationX, ...
                'Within', AbsoluteTolerance(1e-12)))
            testCase.assertThat(actualRotationY, IsEqualTo(expectedRotationY, ...
                'Within', AbsoluteTolerance(1e-12)))
            testCase.assertThat(actualRotationZ, IsEqualTo(expectedRotationZ, ...
                'Within', AbsoluteTolerance(1e-12)))
        end
        
        function testShellElement3d4nDynamic(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            % Model
            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 4, 4, ...
                'elementType', 'ShellElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]);
            model.getNode(13).setDofLoad('DISPLACEMENT_Z', 1);
            support = [x0 xl y0 yl];
            support.fixAllDofs();
            
            % Properties
            model.getAllElements.setPropertyValue('USE_CONSISTENT_MASS_MATRIX',1);
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('THICKNESS', 0.005);
            model.getAllElements.setPropertyValue('DENSITY',7860);
            model.getAllElements.addProperty('RAYLEIGH_ALPHA',3);
            model.getAllElements.addProperty('RAYLEIGH_BETA',4e-4);
            
            % Solver
            dt = .001;
            time = 0;
            endTime = .01;
            solver = NewmarkSolvingStrategy(model, dt);
            
            while time < endTime
                solver.solve();
                time = time + dt;
            end
            
            % Assertion
            actualDisplacementZ = model.getNode(13).getDofValue('DISPLACEMENT_Z','all');
            expectedDisplacementZ = [0,3.21524973685890e-07,...
                9.13911217812143e-07,1.37016296186340e-06,1.75612069115088e-06,...
                2.21377040123328e-06,2.77107997383805e-06,3.36779398166601e-06,...
                3.89553037963957e-06,4.26061622283949e-06,4.44264794892829e-06];
            testCase.assertThat(actualDisplacementZ, IsEqualTo(expectedDisplacementZ, ...
                'Within', RelativeTolerance(1e-7)))
            
        end
        
        function testShellElement3d4nEigen(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            % Model
            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 4, 4, ...
                'elementType', 'ShellElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]);
            support = [x0 xl y0 yl];
            support.fixAllDofs();
            
            % Properties
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('THICKNESS', 0.005);
            model.getAllElements.setPropertyValue('DENSITY',7860);
            
            % Solver
            solver = EigensolverStrategy(model);
            solver.solve(5);
            
            % Assertion
            actualEigenfrequencies = solver.getEigenfrequencies('Hz');
            expectedEigenfrequencies = [41.3911366690480;80.2659687320295;...
                80.2659687320296;103.127323365447;122.970895447928];
            
            testCase.assertThat(actualEigenfrequencies, IsEqualTo(expectedEigenfrequencies, ...
                'Within', RelativeTolerance(1e-7)))
        end
        
        function testDKQElement3d4nStatic (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance

            % Model
            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 4, 4, ...
                'elementType', 'DiscreteKirchhoffElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y"]);
            model.getNode(13).setDofLoad('DISPLACEMENT_Z', 2500);
            support = [x0 xl y0 yl];
            support.fixAllDofs();

            % Properties
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('THICKNESS', 0.005);
            model.getAllElements.setPropertyValue('DENSITY',7860);
            
            % Solver
            solver = SimpleSolvingStrategy(model);
            solver.solve();

            % Assertion
            actualDisplacementZ = model.getAllNodes.getDofValue('DISPLACEMENT_Z');
            actualRotationX = model.getAllNodes.getDofValue('ROTATION_X');
            actualRotationY = model.getAllNodes.getDofValue('ROTATION_Y');
            
            expectedDisplacementZ = [0;0;0;0;0;0;0.00149797814601510;...
                0.00291158986538897;0.00149797814601510;0;0;0.00291158986538897;...
                0.00666684399719017;0.00291158986538897;0;0;0.00149797814601510;...
                0.00291158986538897;0.00149797814601510;0;0;0;0;0;0];
            expectedRotationX = [0;0;0;0;0;0;0.00882854312106295;...
                0.0197154931962282;0.00882854312106307;0;0;4.27500171253438e-18;...
                -8.45148379594905e-18;-2.64711902013052e-17;0;0;-0.00882854312106298;...
                -0.0197154931962282;-0.00882854312106298;0;0;0;0;0;0];
            expectedRotationY = [0;0;0;0;0;0;-0.00882854312106294;...
                1.31064049521836e-17;0.00882854312106299;0;0;-0.0197154931962283;...
                5.76855816486450e-17;0.0197154931962282;0;0;-0.00882854312106298;...
                1.10572958824939e-17;0.00882854312106297;0;0;0;0;0;0];
            
            testCase.assertThat(actualDisplacementZ, IsEqualTo(expectedDisplacementZ, ...
                'Within', AbsoluteTolerance(1e-12)))
            testCase.assertThat(actualRotationX, IsEqualTo(expectedRotationX, ...
                'Within', AbsoluteTolerance(1e-12)))
            testCase.assertThat(actualRotationY, IsEqualTo(expectedRotationY, ...
                'Within', AbsoluteTolerance(1e-12)))
        end
        
        function testDKQElement3d4nDynamic (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance

            % Model
            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 4, 4, ...
                'elementType', 'DiscreteKirchhoffElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y"]);
            model.getNode(13).setDofLoad('DISPLACEMENT_Z', 2500);
            support = [x0 xl y0 yl];
            support.fixAllDofs();

            % Properties
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('THICKNESS', 0.005);
            model.getAllElements.setPropertyValue('DENSITY',7860);
            model.getAllElements.addProperty('RAYLEIGH_ALPHA',3);
            model.getAllElements.addProperty('RAYLEIGH_BETA',4e-4);
            
            % Solver
            dt = .001;
            time = 0;
            endTime = .005;
            solver = NewmarkSolvingStrategy(model, dt);
            
            while time < endTime
                solver.solve();
                time = time + dt;
            end
            
            % Assertion
            actualDisplacementZ = model.getNode(13).getDofValue('DISPLACEMENT_Z','all');
            expectedDisplacementZ = [0 0.000803812434214726 0.00228477804453036...
                0.00342540740465850 0.00439030172787720 0.00553442600308321];
            testCase.assertThat(actualDisplacementZ, IsEqualTo(expectedDisplacementZ, ...
                'Within', RelativeTolerance(1e-7)))
        end
        
        function testDKQElement3d4nEigen(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            
            % Model
            [model, x0, xl, y0, yl] = createRectangularPlate(1, 1, 4, 4, ...
                'elementType', 'DiscreteKirchhoffElement3d4n');
            model.getAllNodes.addDof(["DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y"]);
            support = [x0 xl y0 yl];
            support.fixAllDofs();
            
            % Properties
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.1e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('THICKNESS', 0.005);
            model.getAllElements.setPropertyValue('DENSITY',7860);
            
            % Solver
            solver = EigensolverStrategy(model);
            solver.solve(5);
            
            % Assertion
            actualEigenfrequencies = solver.getEigenfrequencies('Hz');
            expectedEigenfrequencies = [46.1882586911365;103.709269817444;...
                103.709269817444;154.684851632199;197.188070901913];
            
            testCase.assertThat(actualEigenfrequencies, IsEqualTo(expectedEigenfrequencies, ...
                'Within', RelativeTolerance(1e-7)))
        end
    end
end

