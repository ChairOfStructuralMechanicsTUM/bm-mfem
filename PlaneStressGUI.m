function varargout = PlaneStressGUI(varargin)
% PLANESTRESSGUI MATLAB code for PlaneStressGUI.fig
%      PLANESTRESSGUI, by itself, creates a new PLANESTRESSGUI or raises the existing
%      singleton*.
%
%      H = PLANESTRESSGUI returns the handle to a new PLANESTRESSGUI or the handle to
%      the existing singleton*.
%
%      PLANESTRESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLANESTRESSGUI.M with the given input arguments.
%
%      PLANESTRESSGUI('Property','Value',...) creates a new PLANESTRESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlaneStressGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlaneStressGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlaneStressGUI

% Last Modified by GUIDE v2.5 03-Jun-2018 16:46:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlaneStressGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PlaneStressGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PlaneStressGUI is made visible.
function PlaneStressGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlaneStressGUI (see VARARGIN)

% Choose default command line output for PlaneStressGUI
handles.output = hObject;

bm_fem;
guidata(hObject, handles);
set(findall(handles.uipanel_visibility, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uipanel_properties, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uipanel_bC, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uipanel_solver, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uibuttongroup_vis, '-property', 'enable'), 'enable', 'off');
set(handles.edit1,'string','Quad4_4x40.msh');
set(handles.axes1,'Visible', 'off');

cla(handles.axes1);
colorbar off

fieldType = {'Select Field' 'displacement_x' 'displacement_y' 'displacement_absolute'...
    'sigma_xx' 'sigma_yy' 'sigma_xy' 'prin_I' 'prin_II' 'vm_stress'};
set(handles.popupmenu2,'String',fieldType);

% UIWAIT makes PlaneStressGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PlaneStressGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Browse.
function pushbutton_Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    file = uigetfile('tests/*.msh');
    set(handles.edit1,'string', file);

% --- Executes on button press in pushbutton_ImportModel.
function pushbutton_ImportModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ImportModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    file = get(handles.edit1,'string');
    io = ModelIO(file);
    model = io.readModel;
    model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
    handles.model = model;
    vis = VisualizationGUI(model);
    handles.vis = vis;
    handles.loaded = 'false';
    handles.constrained = false;
    
    %populate popupmenu with modelParts
    modelParts = keys(getAllModelParts(model));
    modelPartsDefault = {'Select Model Part'};
    modelPartsList = [modelPartsDefault, modelParts];
    set(handles.popupmenu3,'string',modelPartsList);
    set(handles.popupmenu4,'string',modelPartsList);
    
    axes1 = handles.axes1;
    cla(axes1);
    reset(axes1);
    axes1.Color = [0.94,0.94,0.94];
    axis off
    colorbar off
    
    vis = handles.vis;
    vis.plotUndeformed();

    set(handles.popupmenu2,'Value',1);
    set(handles.checkbox3,'Value',1);
    set(handles.checkbox4,'Value',0);
    set(findall(handles.uipanel_bC, '-property', 'enable'), 'enable', 'off');
    set(findall(handles.uipanel_solver, '-property', 'enable'), 'enable', 'off');
    set(findall(handles.uibuttongroup_vis, '-property', 'enable'), 'enable', 'off');
    set(findall(handles.uipanel_visibility, '-property', 'enable'), 'enable', 'on');
    set(findall(handles.uipanel_visibility, '-property', 'enable'), 'Value', 0);
    set(findall(handles.uipanel_properties, '-property', 'enable'), 'enable', 'on');

    guidata(hObject,handles);
% --- Executes on button press in checkbox_nodeID.
function checkbox_nodeID_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_nodeID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_nodeID
tic;
    if(get(hObject,'Value') == get(hObject,'Max'))
        handles.vis.plotNumbering('nodes','deformed');                   
    else
        handles.vis.clearNumbering;
%         set(findobj(gcf,'tag','NodeNum'),'Visible','off');
    end
toc

% --- Executes on button press in checkbox_elemID.
function checkbox_elemID_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_elemID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_elemID

    if(get(hObject,'Value') == get(hObject,'Max'))
        set(findobj(gcf,'tag','ElemNum'),'Visible','on');                     
    else
        set(findobj(gcf,'tag','ElemNum'),'Visible','off');
    end



function edit_youngsModulus_Callback(hObject, eventdata, handles)
% hObject    handle to edit_youngsModulus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_youngsModulus as text
%        str2double(get(hObject,'String')) returns contents of edit_youngsModulus as a double


% --- Executes during object creation, after setting all properties.
function edit_youngsModulus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_youngsModulus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_thickness_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thickness as text
%        str2double(get(hObject,'String')) returns contents of edit_thickness as a double


% --- Executes during object creation, after setting all properties.
function edit_thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_prxy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_prxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_prxy as text
%        str2double(get(hObject,'String')) returns contents of edit_prxy as a double


% --- Executes during object creation, after setting all properties.
function edit_prxy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_prxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_numberGaussPoint_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numberGaussPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numberGaussPoint as text
%        str2double(get(hObject,'String')) returns contents of edit_numberGaussPoint as a double


% --- Executes during object creation, after setting all properties.
function edit_numberGaussPoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numberGaussPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_submit.
function pushbutton_submit_Callback(hObject, eventdata, handles)
Properties();

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LoadBC();


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ConstrainBC();


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_solve.
function pushbutton_solve_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_solve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(findall(handles.uibuttongroup_vis, '-property', 'enable'), 'enable', 'off');

model = handles.model;
nodes = model.getAllNodes();
loaded = handles.loaded;

if strcmp(loaded,'static')
    solver = SimpleSolvingStrategy(model);
    solver.solve();
    set(findall(handles.uipanel11, '-property', 'enable'), 'enable', 'off');
    
else strcmp(loaded,'dynamic')
    endTime = handles.endTime;
    time = 0;
    dt = handles.dt;
    load = handles.load;
    direction = handles.direction;
    excitationFrequency = handles.excitationFrequency;
    node = nodes(handles.nodeId);
    
    solver = NewmarkSolvingStrategy(model, dt);
    while time < endTime
        applyHarmonicSineLoad(node, load, direction, excitationFrequency, time);
        solver.solve();
        time = time + dt;
    end
    totalSteps = endTime/dt;
    set(handles.text9,'String',totalSteps);
    set(handles.edit9,'String',1);
    set(findall(handles.uipanel11, '-property', 'enable'), 'enable', 'off');
end

% Scale Deflection to 5 % of maximum Model Length 
maxModelLength = max(abs([max(nodes.getX)-min(nodes.getX),max(nodes.getY)-min(nodes.getY)]));
maxDeformation = max(max(abs(nodes.getDofValue('DISPLACEMENT_X'))),max(abs(nodes.getDofValue('DISPLACEMENT_Y'))));
scaling = maxModelLength * 0.05 / maxDeformation;
handles.vis.setScaling(scaling);
set(handles.edit7,'String',scaling);

handles.vis.plotField('Select Field',1);
% handles.vis.plotDeformed(10);
handles.vis.plotConstrain('deformed');
handles.vis.plotLoad('deformed');

set(findobj(gcf,'Tag','deformed'),'Visible','on');
set(handles.popupmenu2,'Value',1);
set(handles.checkbox3,'Value',1);
set(handles.checkbox4,'Value',1);
set(handles.checkbox6,'Value',1);
set(handles.checkbox7,'Value',1);
set(findall(handles.uibuttongroup_vis, '-property', 'enable'), 'enable', 'on');
set(findall(handles.uipanel_properties, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uipanel_bC, '-property', 'enable'), 'enable', 'off');
set(handles.pushbutton11,'Enable','on');
set(hObject,'Enable','off');
guidata(hObject, handles);
% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
    if(get(hObject,'Value') == get(hObject,'Max'))
        set(findobj(gcf,'Tag','undeformed'),'Visible','on');                     
    else
        set(findobj(gcf,'Tag','undeformed'),'Visible','off');
    end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
contents = cellstr(get(handles.popupmenu2,'String'));
fieldType = contents{get(handles.popupmenu2,'Value')};

if(get(hObject,'Value') == get(hObject,'Max'))
    
    set(handles.text6,'Enable','on');
    set(handles.text7,'Enable','on');
    set(handles.edit7,'Enable','on');
    set(handles.popupmenu2,'Enable','on');
    set(findobj(gcf,'Tag',fieldType),'Visible','on');

    
else
    set(findobj(gcf,'Tag',fieldType),'Visible','off');
    set(handles.text6,'Enable','off');
    set(handles.text7,'Enable','off');
    set(handles.edit7,'Enable','off');
    set(handles.popupmenu2,'Enable','off');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
% axes1 = handles.axes1;
% children = axes1.Children;
% set(children,'Visible','off');

vis = handles.vis;
contents = cellstr(get(hObject,'String'));
fieldType = contents{get(hObject,'Value')};
vis.plotField(fieldType);
uistack(findobj(gcf,'Tag','load'),'top');
uistack(findobj(gcf,'Tag','constrain'),'top');
uistack(findobj(gcf,'Tag','undeformed'),'top');
uistack(findobj(gcf,'Tag','NodeNum'),'top');
uistack(findobj(gcf,'Tag','ElemNum'),'top');

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

scaling = str2num(get(handles.edit7,'String'));

if scaling ~= handles.deformedScaling
    handles.vis.setScaling(scaling);
    handles.vis.plotDeformed;
end

if(get(hObject,'Value') == get(hObject,'Max'))
    
    set(handles.checkbox3,'Enable','on')
    set(handles.checkbox4,'Enable','on')
    set(handles.popupmenu2,'Enable','off')
    axes1 = handles.axes1;
    children = axes1.Children;
    set(children,'Visible','off');
    colorbar off

    if(get(handles.checkbox3,'Value') == get(handles.checkbox3,'Max'))
        set(findobj(gcf,'Tag','undeformed'),'Visible','on');                     
    else
        set(findobj(gcf,'Tag','undeformed'),'Visible','off');
    end
    
    if(get(handles.checkbox4,'Value') == get(handles.checkbox4,'Max'))
        set(findobj(gcf,'Tag','deformed'),'Visible','on');                     
    else
        set(findobj(gcf,'Tag','deformed'),'Visible','off');
    end
    
end
    

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
scaling = str2num(get(handles.edit7,'String'));

if(get(hObject,'Value') == get(hObject,'Max'))
    
    axes1 = handles.axes1;
    children = axes1.Children;
    set(children,'Visible','off');
    
    set(handles.popupmenu2,'Enable','on')
    set(handles.checkbox3,'Enable','off')
    set(handles.checkbox4,'Enable','off')
    
    contents = cellstr(get(handles.popupmenu2,'String'));
    fieldType = contents{get(handles.popupmenu2,'Value')};
    
    if ~strcmp(fieldType,'Select Field')
        if scaling == handles.fieldScaling
            set(findobj(gcf,'Tag',fieldType),'Visible','on');
            colorbar
        else
            handles.vis.plotField(fieldType);
        end
    end
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton_solve, 'enable', 'off');

model = handles.model;
nodes = model.getAllNodes();
nodes.setDofLoad('DISPLACEMENT_X',0);
nodes.setDofLoad('DISPLACEMENT_Y',0);
handles.loaded = 'false';
guidata(hObject,handles);

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton_solve, 'enable', 'off');

model = handles.model;
nodes = model.getAllNodes();
nodes.unfixDof('DISPLACEMENT_X');
nodes.unfixDof('DISPLACEMENT_Y');
handles.constrained = false;
guidata(hObject,handles);


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5

if(get(hObject,'Value') == get(hObject,'Max'))
    axis on               
else
    axis off
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

scaling = str2num(get(hObject,'String'));
contents = cellstr(get(handles.popupmenu2,'String'));
fieldType = contents{get(handles.popupmenu2,'Value')};

if scaling ~= handles.vis.getScaling
    handles.vis.setScaling(scaling);
    handles.vis.plotField(fieldType);
    if (get(handles.checkbox6,'Value') == get(hObject,'Max'))
        handles.vis.plotLoad('deformed');
    end
end
    
guidata(hObject,handles);
% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if(get(hObject,'Value') == get(hObject,'Max'))
    set(findobj(gcf,'Tag','load'),'Visible','on');
else
    set(findobj(gcf,'Tag','load'),'Visible','off');
end

% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7
if(get(hObject,'Value') == get(hObject,'Max'))
    set(findobj(gcf,'Tag','constrain'),'Visible','on');
else
    set(findobj(gcf,'Tag','constrain'),'Visible','off');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1);
colorbar off
handles.vis.plotUndeformed;
handles.vis.plotLoad('undeformed');
handles.vis.plotConstrain('undeformed');
set(handles.pushbutton_solve,'Enable','on');
set(handles.checkbox_elemID,'Enable','on');
set(handles.checkbox_nodeID,'Enable','on');
set(findall(handles.uibuttongroup_vis, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uipanel_properties, '-property', 'enable'), 'enable', 'on');
set(findall(handles.uipanel_bC, '-property', 'enable'), 'enable', 'on');


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
contents = cellstr(get(handles.popupmenu2,'String'));
fieldType = contents{get(handles.popupmenu2,'Value')};

step = get(hObject,'Value');
handles.vis.plotField(fieldType,step);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_density_Callback(hObject, eventdata, handles)
% hObject    handle to edit_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_density as text
%        str2double(get(hObject,'String')) returns contents of edit_density as a double


% --- Executes during object creation, after setting all properties.
function edit_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton13,'Enable','on');
step = str2num(get(handles.edit9,'String'));
step = step + 1;
set(handles.edit9,'string',step);
totalSteps = str2num(get(handles.text9,'String'));
if step == totalSteps
    set(hObject,'Enable','off');
end
contents = cellstr(get(handles.popupmenu2,'String'));
fieldType = contents{get(handles.popupmenu2,'Value')};

handles.vis.plotField(fieldType,step+1);

function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
contents = cellstr(get(handles.popupmenu2,'String'));
fieldType = contents{get(handles.popupmenu2,'Value')};

totalSteps = str2num(get(handles.text9,'String'));
step = str2num(get(hObject,'String'));

if step == totalSteps
    set(handles.pushbutton12,'Enable','off');
elseif step == 1
    set(handles.pushbutton12,'Enable','on');
else
    set(handles.pushbutton12,'Enable','on');
end
    
handles.vis.plotField(fieldType,step+1);
% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton12,'Enable','on');
step = str2num(get(handles.edit9,'String'));
step = step - 1;
if step == 1
    set(hObject,'Enable','off');
end
set(handles.edit9,'string',step);
contents = cellstr(get(handles.popupmenu2,'String'));
fieldType = contents{get(handles.popupmenu2,'Value')};

handles.vis.plotField(fieldType,step+1);



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.popupmenu2,'String'));
fieldType = contents{get(handles.popupmenu2,'Value')};

step = str2num(get(handles.edit9,'String'));

model = handles.model;
nodes = model.getAllNodes;
if ~isempty(get(handles.edit10,'String'))
    startNodeId = str2double(get(handles.edit10,'String'));
else
    startContents = cellstr(get(handles.popupmenu3,'String'));
    modelPart = startContents{get(handles.popupmenu3,'Value')};
    startNodeId = model.getModelPart(modelPart).getId();
end

if ~isempty(get(handles.edit11,'String'))
    endNodeId = str2double(get(handles.edit11,'String'));
else
    endContents = cellstr(get(handles.popupmenu4,'String'));
    modelPart = endContents{get(handles.popupmenu4,'Value')};
    endNodeId = model.getModelPart(modelPart).getId();
end

selectedNodeIds = getNodeIdsAlongLine(startNodeId,endNodeId,nodes);

handles.vis.plotLineData(selectedNodeIds,fieldType,step);
% if strcmp(dir,'horizontal')
%     x = selectedNodes.getX';
%     value = selectedNodes.getDofValue('DISPLACEMENT_X');
%     line(x,value);
% else
%     y = selectedNodes.getY';
%     value = [1;2;3;4;5];
%     [y,value] = sortrows([y value]);
%     figure
%     line(value,y);
% end

    


% --- Executes on selection change in popupmenu13.
function popupmenu13_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu13 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu13


% --- Executes during object creation, after setting all properties.
function popupmenu13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double
contents = cellstr(get(handles.popupmenu2,'String'));
fieldType = contents{get(handles.popupmenu2,'Value')};

step = str2num(get(handles.edit9,'String'));

model = handles.model;
if ~isempty(get(handles.edit20,'String'))
    nodeId = str2double(get(handles.edit20,'String'));
else
    startContents = cellstr(get(handles.popupmenu13,'String'));
    modelPart = startContents{get(handles.popupmenu13,'Value')};
    nodeId = model.getModelPart(modelPart).getId();
end
handles.vis.plotLineData(nodeId,fieldType,step); 

% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.popupmenu2,'String'));
fieldType = contents{get(handles.popupmenu2,'Value')};

step = str2num(get(handles.edit9,'String'));

model = handles.model;
if ~isempty(get(handles.edit20,'String'))
    nodeId = str2double(get(handles.edit20,'String'));
else
    startContents = cellstr(get(handles.popupmenu13,'String'));
    modelPart = startContents{get(handles.popupmenu13,'Value')};
    nodeId = model.getModelPart(modelPart).getId();
end
handles.vis.plotLineData(nodeId,fieldType,step);    