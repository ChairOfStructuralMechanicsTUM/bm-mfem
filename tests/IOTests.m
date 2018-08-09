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
            
            testCase.assumeEqual(exist('getPaths','file'),2, ...
                'Script getPaths not found. Please create it to use ANSYS import.');
            testCase.assumeEqual(exist(getPaths('ansys'),'file'),2, ...
                'ANSYS executable not found. Please check your getPaths script.');

            %eigenvalue analysis
            io = AnsysInput('tests/ansys_model.inp', getPaths('ansys'));
            model = io.readModel;

            solver = EigensolverStrategy(model);
            solver.solve(10);

            actualEigenfrequencies = solver.getEigenfrequencies('Hz');
            expectedEigenfrequencies = [83.17041059064;519.8081004268;1217.833700837;...
                1463.709928893;2903.635403021;3695.100847286;4870.092778333;...
                6297.883310161;7337.460973345;9109.774679695];

            testCase.assertThat(actualEigenfrequencies, IsEqualTo(expectedEigenfrequencies, ...
                'Within', RelativeTolerance(1e-7)))

            % static analysis
            model = io.readModel;
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

            % dynamic analysis
            model = io.readModel;

            dt = .00005;
            time = 0;
            endTime = .005;
            solver = NewmarkSolvingStrategy(model, dt);

            while time < endTime
                solver.solve();
                time = time + dt;
            end

            actualDisplacement = model.getNode(10).getDofValue('DISPLACEMENT_Y','all');
            load('tests/test_data.mat','ansys_input_dynamic');

            testCase.assertThat(actualDisplacement, IsEqualTo(ansys_input_dynamic, ...
                'Within', RelativeTolerance(1e-4)))
            
        end
    end
    
end

