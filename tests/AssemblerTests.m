classdef AssemblerTests <  matlab.unittest.TestCase
    %ASSEMBLERTESTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Test)
        function testSimpleAssembler(testCase)
            node01 = Node(1,-4,3);
            node02 = Node(2,0,3);
            node03 = Node(3,4,3);
            node04 = Node(4,0,0);
            nodeArray = [node01 node02 node03 node04];
            
            ele01 = BarElement2d2n(1,[node01 node02]);
            ele02 = BarElement2d2n(2,[node02 node03]);
            ele03 = BarElement2d2n(3,[node01 node04]);
            ele04 = BarElement2d2n(4,[node04 node02]);
            ele05 = BarElement2d2n(5,[node03 node04]);
            elementArray = [ele01 ele02 ele03 ele04 ele05];
            
            model = FemModel(nodeArray, elementArray);
            
            model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
            model.getAllElements.setPropertyValue('CROSS_SECTION',2);
            model.getAllElements.setPropertyValue('YOUNGS_MODULUS',3000);
            
            model.initialize;
            
            actualSolution = SimpleAssembler.assembleGlobalStiffnessMatrix(model);
            
            expectedSolution = sparse([2268 -576 -1500 0 0 0 -768 576; ...
                -576 432 0 0 0 0 576 -432; ...
                -1500 0 3000 0 -1500 0 0 0; ...
                0 0 0 2000 0 0 0 -2000; ...
                0 0 -1500 0 2268 576 -768 -576; ...
                0 0 0 0 576 432 -576 -432; ...
                -768 576 0 0 -768 -576 1536 0; ...
                576 -432 0 -2000 -576 -432 0 2864]);
            
            testCase.verifyEqual(actualSolution, expectedSolution)
        end
        
    end
    
end
