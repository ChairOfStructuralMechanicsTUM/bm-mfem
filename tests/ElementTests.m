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
            
            actualEigenfrequenciesCantilever = sort(solver.getEigenfrequencies);
            expectedEigenfrequenciesCantilever = [0.915242085421142;5.73511288165053;16.0584347115302;31.4797553897702;52.0975276397135];

            testCase.assertThat(actualEigenfrequenciesCantilever, IsEqualTo(expectedEigenfrequenciesCantilever, ...
                'Within', AbsoluteTolerance(1e-5)))
            
            warning('on','all')
        end
        
        function testReissnerMindlinElement3d4nStatic (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
        
            a=linspace(0,1,5);
            % Nodes
            ii = 1;
            for i=1:5
                for j=1:5
                    node(ii) = Node(ii,a(j),a(i),0);
                    ii = ii + 1; 
                end
            end
            nodeArray = node(:)';
            nodeArray.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y'});
            
            % Elements
            ii = 1; 
            for i = 1 : 5 : 16
                for j= i:(i+3)
                    ele(ii) = ReissnerMindlinElement3d4n(ii, [node(j) node(j+1) node(j+6) node(j+5)]);
                    ii = ii + 1;
                end
            end
            elementArray = ele(:)';
            
            % Boundary Conditions 
            ii = 1; 
            for i=1:length(node)
                if node(i).getX == 0 || node(i).getY == 0 || node(i).getX == 1 || node(i).getY == 1
                    boundary(ii) = node(i);
                    ii= ii+1;
                end
            end
            boundary.fixDof('DISPLACEMENT_Z');
            boundary.fixDof('ROTATION_X');
            boundary.fixDof('ROTATION_Y');

            % Properties
            elementArray.setPropertyValue('THICKNESS', .1);
            elementArray.setPropertyValue('YOUNGS_MODULUS', 10920);
            elementArray.setPropertyValue('POISSON_RATIO', 0.3);
            elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
            elementArray.setPropertyValue('DENSITY', 1);
            elementArray.setPropertyValue('SHEAR_CORRECTION_FACTOR', 5/6);
            
            % Solver
            model = FemModel(nodeArray,elementArray);
            model.getNode(13).setDofLoad('DISPLACEMENT_Z', -0.25);
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
                'Within', AbsoluteTolerance(1e-7)))
            testCase.assertThat(actualRotationX, IsEqualTo(expectedRotationX, ...
                'Within', AbsoluteTolerance(1e-7)))
            testCase.assertThat(actualRotationY, IsEqualTo(expectedRotationY, ...
                'Within', AbsoluteTolerance(1e-7)))

        end 
        
        function testReissnerMindlinElement3d4nDynamic (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.AbsoluteTolerance
            
            model = FemModel();
            nelex = 5;
            mid = 13;
            a=linspace(0,1,nelex);
            
            % Nodes
            ii = 1;
            for i=1:nelex
                for j=1:nelex
                    model.addNewNode(ii,a(j),a(i),0);
                    ii = ii + 1;
                end
            end
            model.getAllNodes.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y'});
            nodes = model.getAllNodes();
            
            % Elements
            ii = 1;
            for i = 1 : nelex : (nelex^2-nelex)
                for j= i:(i+nelex-2)
                    model.addNewElement('ReissnerMindlinElement3d4n', ii, [nodes(j) nodes(j+1) nodes(j+nelex+1) nodes(j+nelex)]);
                    ii = ii + 1;
                end
            end
            
            % Boundary Conditions
            for i=1:length(nodes)
                if nodes(i).getX == 0 || nodes(i).getY == 0 || nodes(i).getX == 1 || nodes(i).getY == 1
                    nodes(i).fixDof('DISPLACEMENT_Z');
                end
            end
            
            % Properties
            model.getAllElements.setPropertyValue('THICKNESS', .00125);
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.07e11);
            model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
            model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT', 4);
            model.getAllElements.setPropertyValue('DENSITY', 7850);
            model.getAllElements.setPropertyValue('SHEAR_CORRECTION_FACTOR', 5/6);
            
            % Solver
            model.getNode(mid).setDofLoad('DISPLACEMENT_Z', 1);
            
            dt = .01;
            time = 0;
            endTime = .2;
            solver = NewmarkSolvingStrategy(model, dt);
            
            while time < endTime
                solver.solve();
                time = time + dt;                
            end
            
            actualDisplacementZ = model.getNode(mid).getDofValue('DISPLACEMENT_Z','all');
            
            expectedDisplacementZ = [0 5.63499025462025e-05 0.000110356168474448 ...
                0.000195291413790720 0.000348022151124683 0.000436442533214557 ...
                0.000549911355006482 0.000638983074213174 0.000633900627028416 ...
                0.000646184028790807 0.000576985852794556 0.000460890315653250 ...
                0.000377492675821615 0.000228575313367927 0.000127279605314402 ...
                7.17349829089917e-05 3.30152776432233e-06 4.13423283962989e-05 ...
                9.48520255923515e-05 0.000166360041372910 0.000315266268863059];
            
            testCase.assertThat(actualDisplacementZ, IsEqualTo(expectedDisplacementZ, ...
                'Within', AbsoluteTolerance(1e-7)))
            
        end
        
        function testReissnerMindlinElement3d4nEigen (testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
        
            a=linspace(0,1,5);
            % Nodes
            ii = 1;
            for i=1:5
                for j=1:5
                    node(ii) = Node(ii,a(j),a(i),0);
                    ii = ii + 1; 
                end
            end
            nodeArray = node(:)';
            nodeArray.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y'});
            
            % Elements
            ii = 1; 
            for i = 1 : 5 : 16
                for j= i:(i+3)
                    ele(ii) = ReissnerMindlinElement3d4n(ii, [node(j) node(j+1) node(j+6) node(j+5)]);
                    ii = ii + 1;
                end
            end
            elementArray = ele(:)';
            
            % Boundary Conditions 
            ii = 1; 
            for i=1:length(node)
                if node(i).getX == 0 || node(i).getY == 0 || node(i).getX == 1 || node(i).getY == 1
                    boundary(ii) = node(i);
                    ii= ii+1;
                end
            end
            boundary.fixDof('DISPLACEMENT_Z');

            % Properties
            elementArray.setPropertyValue('THICKNESS', 0.01);
            elementArray.setPropertyValue('YOUNGS_MODULUS', 10920);
            elementArray.setPropertyValue('POISSON_RATIO', 0.3);
            elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
            elementArray.setPropertyValue('DENSITY', 1);
            elementArray.setPropertyValue('SHEAR_CORRECTION_FACTOR', 5/6);
            
            % Solver
            model = FemModel(nodeArray,elementArray);
            
            solver = EigensolverStrategy(model);
            solver.solve(5);
            
            actualEigenfrequencies = sort(solver.getEigenfrequencies);
            expectedEigenfrequencies = [1.01402840311520;3.22165476249304;3.22165476249320;...
                                            5.32248978160243;10.5674305646586];
            
            testCase.assertThat(actualEigenfrequencies, IsEqualTo(expectedEigenfrequencies, ...
                'Within', RelativeTolerance(1e-7)))

        end 
    end
end

