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
                0.8475 0.795999 1.135 0.885464 1.415 1.695];
            expectedDisplacementY = [0 -1.775600 -1.792260 -2.291930 -2.316600 ...
                -2.385940 -2.421940 -2.291930 -2.316600 -1.775600 -1.792260 0];
            
            testCase.assertThat(actualDisplacementX, IsEqualTo(expectedDisplacementX, ...
                'Within', RelativeTolerance(1e-5)))
            testCase.assertThat(actualDisplacementY, IsEqualTo(expectedDisplacementY, ...
                'Within', RelativeTolerance(1e-5)))
        end
        
        function externalScriptsTest(testCase)
           bridge;
           Bridge_with_inputFile;
        end
        
    end
    
end

