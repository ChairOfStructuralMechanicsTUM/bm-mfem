io=MdpaInput('test.mdpa');
model=io.readMdpa
No appropriate method, property, or field 'readMdpa' for class 'MdpaInput'.
 
model=io.readModel();
Array indices must be positive integers or logical values.

Error in MdpaInput.readProperties (line 77)
            props(nProp) = property;

Error in MdpaInput/readModel (line 37)
                    [fid,props] = obj.readProperties(fid,props,nProp);
 
model=io.readModel();
Array indices must be positive integers or logical values.

Error in MdpaInput.readProperties (line 77)
            props(nProp) = property;

Error in MdpaInput/readModel (line 37)
                    [fid,props] = obj.readProperties(fid,props,nProp);
 
io=MdpaInput('mdpatest.mdpa');
model=io.readModel()

model = 

  FemModel with no properties.

model.getNodeArray
No appropriate method, property, or field 'getNodeArray' for class 'FemModel'.
 
model.getAllNodes

ans = 

  1×14 Node array with no properties.

model.getSubModelParts
No appropriate method, property, or field 'getSubModelParts' for class 'FemModel'.
 
model.getModelParts
No appropriate method, property, or field 'getModelParts' for class 'FemModel'.
 
Did you mean:
clear
io=MdpaInput('mdpatest.mdpa');
model=io.readModel()

model = 

  FemModel with no properties.

v=Visualization(model);
v.plotUndeformed
Index exceeds array bounds.

Error in Visualization/plotUndeformed (line 45)
                    text(c(1), c(2), c(3), elemStr,'Interpreter','latex', ...
 
model.getModelParts
No appropriate method, property, or field 'getModelParts' for class 'FemModel'.
 
Did you mean:
model.getAllModelParts

ans =

  1×6 cell array

    {'DISPLACEMENT'}    {'DISPLACEMENT_al…'}    {'Parts_structure'}    {'PointLoad'}    {'ROTATION'}    {'empty'}

a=model.getModelPart('Parts_structure')

a = 

  FemModelPart with properties:

              name: 'Parts_structure'
         nodeArray: [1×11 Node]
      elementArray: [1×10 BeamElement3d2n]
          dofArray: []
         fixedDofs: []
          freeDofs: []
    parentFemModel: [1×1 FemModel]

