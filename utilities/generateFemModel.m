function femModel = generateFemModel(nodeArray, elementType, elementArray, materialData, propertyData)
%GENERATEFEMMODEL Generates a FEM model from node and element data
%   femModel = GENERATEFEMMODEL(NODEARRAY, 'ELEMENTTYPE', ELEMENTARRAY, MATERIALDATA, PROPERTYDATA)
%       NODEARRAY : nodes in the form id, x, y(, z)
%       'ELEMENTTYPE' : name of the desired element; available elements:
%           BarElement2d2n, BarElement3d2n
%       ELEMENTPROPERTIES : map of all neccessary properties
%       ELEMENTARRAY : elements in the form id, node...

femNodes(1, length(nodeArray)) = Node;

for itNode = 1:length(nodeArray)
    cNode = nodeArray(itNode,:);
    node = Node(cNode(1), cNode(2), cNode(3), cNode(4));
    femNodes(itNode) = node;
end

switch elementType
    case 'BarElement2d2n'
        femElements(1, length(elementArray)) = BarElement2d2n;
        for itEle = 1:length(elementArray)
            cEle = elementArray(itEle,:);
            nProperty = cEle(2);
            cProperties = propertyData(nProperty);
            nMaterial = cProperties.getValue('material');
            
            element = BarElement2d2n(cEle(1), ...
                [femNodes(cEle(3)) femNodes(cEle(4))], ...
                materialData(nMaterial), ...
                cProperties.getValue('crossSectionArea'));
            femElements(itEle) = element;
        end
        
    case 'BarElement3d2n'
        femElements(1, length(elementArray)) = BarElement3d2n;
        for itEle = 1:length(elementArray)
            cEle = elementArray(itEle,:);
            nProperty = cEle(2);
            cProperties = propertyData(nProperty);
            nMaterial = cProperties.getValue('material');
            
            element = BarElement3d2n(cEle(1), ...
                [femNodes(cEle(3)) femNodes(cEle(4))], ...
                materialData(nMaterial), ...
                cProperties.getValue('crossSectionArea'));
            femElements(itEle) = element;
        end
    otherwise
        disp('unknown element type')
end

femModel = FemModel(femNodes, femElements);

end

