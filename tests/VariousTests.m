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
            prop.setValue('test',1);
            prop.setValue('test_copy',1);
            
            ele01.setProperties(prop);
            ele02.setProperties(prop);
            
            ele02.setPropertyValue('test_copy',2);
            
            testCase.verifyEqual(ele01.getPropertyValue('test'), 1)
            testCase.verifyEqual(ele01.getPropertyValue('test_copy'), 1)
            testCase.verifyEqual(ele02.getPropertyValue('test'), 1)
            testCase.verifyEqual(ele02.getPropertyValue('test_copy'), 2)
        end
    end
    
end

