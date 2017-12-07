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
%             model.getAllElements.setPropertyValue('SHEAR_MODULUS',2.1e11/2/1.3);
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
            
            actualEigenfrequenciesCantilever = sort(solver.getEigenfrequencies);
            expectedEigenfrequenciesCantilever = [0.915242085421142;5.73511288165053;16.0584347115302;31.4797553897702;52.0975276397135];

            testCase.assertThat(actualEigenfrequenciesCantilever, IsEqualTo(expectedEigenfrequenciesCantilever, ...
                'Within', AbsoluteTolerance(1e-5)))
        end
        
        
        
    end
    
end

