classdef IOTests <  matlab.unittest.TestCase
    %IOTESTS Tests for input / output of model data
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Test)
        
        function mdpaInputTest(testCase)
        %MDPAINPUTTEST tests the input from a mdpa file using properties,
        %   line comments, different element types, and model parts
            import matlab.unittest.constraints.IsTrue
            
            io = MdpaInput('tests/mdpatest.mdpa');
            model = io.readModel;
            
            testCase.verifyEqual(model.getNode(14).getX,2.213);
            testCase.verifyEqual(model.getElement(9).getPropertyValue('DENSITY'),7850);
            testCase.verifyEqual(model.getElement(11).getPropertyValue('POISSON_RATIO'),0.3);
            testCase.verifyThat(isa(model.getElement(11),'ReissnerMindlinElement3d4n'),IsTrue);
            testCase.verifyEqual(model.getModelPart('PointLoad').getNodes.getId(),11);
        end
        
        function ansysInputTest(testCase)
        %ANSYSINPUTTEST tests input from an ANSYS model in a modal analysis
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
        
            io = AnsysInput('examples/model_qa.txt', getPaths('ansys'));
            model = io.readModel;
            
            solver = EigensolverStrategy(model);
            solver.solve(10);
            
            actualEigenfrequencies = solver.getEigenfrequencies('Hz');
            expectedEigenfrequencies = [83.17041059064;519.8081004268;1217.833700837;...
                1463.709928893;2903.635403021;3695.100847286;4870.092778333;...
                6297.883310161;7337.460973345;9109.774679695];
            
            testCase.assertThat(actualEigenfrequencies, IsEqualTo(expectedEigenfrequencies, ...
                'Within', RelativeTolerance(1e-7)))
            
        end
    end
    
end

