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
            
            model = io.readModel;
            model.getNode(10).setDofLoad('DISPLACEMENT_Y', -1);
            solver = SimpleSolvingStrategy(model);
            solver.solve();
            
            actualDisplacement = model.getAllNodes().getDofValue('DISPLACEMENT_Y');
            expectedDisplacement = [0.0000E+00;-1.9071E-06;-4.3293E-08;-1.6443E-07;...
                -3.5260E-07;-5.9656E-07;-8.8517E-07;-1.2072E-06;-1.5517E-06;-1.9077E-06;...
                -1.9072E-06;0.0000E+00;-1.5515E-06;-1.2073E-06;-8.8516E-07;-5.9656E-07;...
                -3.5260E-07;-1.6443E-07;-4.3294E-08;1.9033E-09;-4.1353E-08;-1.6284E-07;...
                -3.5126E-07;-5.9549E-07;-8.8436E-07;-1.2067E-06;-1.5514E-06];
            
            testCase.assertThat(actualDisplacement, IsEqualTo(expectedDisplacement, ...
                'Within', RelativeTolerance(1e-4)))
            
        end
    end
    
end

