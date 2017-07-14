classdef ElementTests < matlab.unittest.TestCase
    %ELEMENTTESTS Tests for all elements
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Test)
        
        function testBarElement2d2n(testCase)
            node01 = Node(1,0,0);
            node02 = Node(2,30,40);
            mat = Material('test');
            mat.setValue('YOUNGS_MODULUS', 1000);
            ele01 = BarElement2d2n(1,[node01 node02], mat, 5);
            
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
            mat = Material('test');
            mat.setValue('YOUNGS_MODULUS', 343);
            ele01 = BarElement3d2n(1,[node01 node02], mat, 10);
            
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
            
            props = PropertyContainer();
            props.setValue('ELEMENTAL_STIFFNESS',100);
            props.setValue('ELEMENTAL_DAMPING',[1,20,1]);
            model.addNewElement('SpringDamperElement3d2n',1,[1 2],props);
            model.addNewElement('SpringDamperElement3d2n',2,[2 3],props);
            
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
        
        
        
    end
    
end

