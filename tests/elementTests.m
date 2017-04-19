classdef elementTests < matlab.unittest.TestCase
    %ELEMENTTESTS Tests for all elements
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Test)
        
        function testBarElement2d2n(testCase)
            node01 = Node(1,0,0);
            node02 = Node(2,30,40);
            mat = Material('test');
            mat.addParameter('YOUNGS_MODULUS', 1000);
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
            mat.addParameter('YOUNGS_MODULUS', 343);
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
        
    end
    
end

