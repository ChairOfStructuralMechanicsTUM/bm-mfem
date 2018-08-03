classdef VariousTests <  matlab.unittest.TestCase
    %VARIOUSTESTS Class for various small tests
    %   e.g. IO, properties, etc.
    
    properties
    end
    
    methods (Test)
        
        function elementPropertyTest(testCase)
        %ELEMENTPROPERTYTEST tests the correct assignment and copy behavior
        %   of properties assigned to an element
            model = FemModel();
            model.addNewNode(1,0,0,0);
            model.addNewNode(2,1,0,0);
            model.addNewNode(3,2,0,0);
            ele01 = model.addNewElement('BarElement3d2n',1,[1 2]);
            ele02 = model.addNewElement('BarElement3d2n',2,[3 2]);
            
            prop = PropertyContainer();
            prop.addValue('CROSS_SECTION',1);
            prop.addValue('YOUNGS_MODULUS',1);
            
            ele01.setProperties(prop);
            ele02.setProperties(prop);
            
            ele02.setPropertyValue('YOUNGS_MODULUS',2);
            
            testCase.verifyEqual(ele01.getPropertyValue('CROSS_SECTION'), 1)
            testCase.verifyEqual(ele01.getPropertyValue('YOUNGS_MODULUS'), 1)
            testCase.verifyEqual(ele02.getPropertyValue('CROSS_SECTION'), 1)
            testCase.verifyEqual(ele02.getPropertyValue('YOUNGS_MODULUS'), 2)
        end
        
        function checkConvexityTestQuads(testCase)
        %CHECKCONVEXITYTEST tests the convexity check for quads
            nodeArray = [Node(1,0,0) Node(3,2,1) Node(2,2,0) Node(4,0,1)];
            nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
            
            testCase.verifyError(@() QuadrilateralElement2d4n(1,nodeArray), ...
                'MATLAB:bm_mfem:elementNotConvex');
        end
        
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
        
        function modelPartTest(testCase)
        %MODELPARTTEST tests the behavior of modelparts
            import matlab.unittest.constraints.IsTrue
            
            model = FemModel();
            model.addNewNode(1,0,0,0);
            model.addNewNode(2,1,0,0);
            model.addNewNode(3,2,0,0);
            model.addNewNode(4,3,0,0);
            model.addNewElement('BarElement3d2n',1,[1 2]);
            model.addNewElement('BarElement3d2n',2,[3 2]);
            model.addNewElement('BarElement3d2n',3,[3 4]);
            model.getAllNodes.addDof(["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]);
            model.getAllElements().setPropertyValue('DENSITY',100);
            
            model.addNewModelPart('mp1',[1 2],1);
            model.addNewModelPart('mp2',[2 3 4],[2 3]);
            model.getModelPart('mp2').getElements().setPropertyValue('DENSITY',200);
            model.getModelPart('mp1').getNodes().fixAllDofs();
            
            testCase.verifyEqual(model.getElement(1).getPropertyValue('DENSITY'),100);
            testCase.verifyEqual(model.getElement(2).getPropertyValue('DENSITY'),200);
            testCase.verifyEqual(model.getNode(1).getDofArray().isFixed(),true(1,6));
            testCase.verifyEqual(model.getNode(3).getDofArray().isFixed(),false(1,6));
        end
    end
    
end

