classdef VariousTests <  matlab.unittest.TestCase
    %VARIOUSTESTS Class for various small tests
    %   e.g. IO, properties, etc.
    
    properties
    end
    
    methods (Test)
        
        function elementPropertyTest(testCase)
            %ELEMENTPROPERTYTEST tests the correct assignment and copy
            %behavior of properties assigned to an element
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
        
        function checkConvexityTest(testCase)
            nodeArray = [Node(1,0,0) Node(3,2,1) Node(2,2,0) Node(4,0,1)];
            nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
            
            testCase.verifyError(@() QuadrilateralElement2d4n(1,nodeArray), ...
                'MATLAB:bm_mfem:elementNotConvex');
        end
    end
    
end

