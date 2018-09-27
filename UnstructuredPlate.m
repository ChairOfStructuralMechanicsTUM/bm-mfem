            io = MdpaInput('test.mdpa');
            model = io.readModel;
            model.getAllNodes.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y','DISPLACEMENT_X', 'DISPLACEMENT_Y', 'ROTATION_Z'})
           
            v = Visualization(model);
            v.plotUndeformed();