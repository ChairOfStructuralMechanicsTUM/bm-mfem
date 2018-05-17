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

% Last Modified by GUIDE v2.5 17-May-2018 16:29:18

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

% Update handles structure
guidata(hObject, handles);
set(findall(handles.uipanel_visibility, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uipanel_properties, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uipanel_bC, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uipanel_solver, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uibuttongroup_vis, '-property', 'enable'), 'enable', 'off');
set(handles.edit1,'string','Quad4_4x40.msh');
set(handles.edit_youngsModulus,'string','2e+05');
set(handles.edit_thickness,'string','0.5');
set(handles.edit_prxy,'string','0.3');
set(handles.edit_numberGaussPoint,'string','3');
set(handles.axes1,'Visible', 'off');

cla(handles.axes1);
colorbar off

fieldType = {'Select Field' 'sigma_xx' 'sigma_yy' 'sigma_xy' 'prin_I' 'prin_II' 'vm_stress'};
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
    handles.loaded = false;
    handles.constrained = false;
    guidata(hObject,handles);
    
    axes1 = handles.axes1;
    cla(axes1);
    reset(axes1);
    axes1.Color = [0.94,0.94,0.94];
    axis off
    colorbar off
    
    vis = VisualizationGUI(model);
    vis.plotUndeformed();

    set(handles.popupmenu2,'Value',1);
    set(handles.checkbox3,'Value',1);
    set(handles.checkbox4,'Value',0);
    set(handles.radiobutton1,'Value',1);
    set(findall(handles.uipanel_bC, '-property', 'enable'), 'enable', 'off');
    set(findall(handles.uipanel_solver, '-property', 'enable'), 'enable', 'off');
    set(findall(handles.uibuttongroup_vis, '-property', 'enable'), 'enable', 'off');
    set(findall(handles.uipanel_visibility, '-property', 'enable'), 'enable', 'on');
    set(findall(handles.uipanel_visibility, '-property', 'enable'), 'Value', 0);
    set(findall(handles.uipanel_properties, '-property', 'enable'), 'enable', 'on');


% --- Executes on button press in checkbox_NodeNum.
function checkbox_NodeNum_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_NodeNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_NodeNum

    if(get(hObject,'Value') == get(hObject,'Max'))
        set(findobj(gcf,'tag','NodeNum'),'Visible','on');                     
    else
        set(findobj(gcf,'tag','NodeNum'),'Visible','off');
    end


% --- Executes on button press in checkbox_ElemNum.
function checkbox_ElemNum_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ElemNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ElemNum

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
% hObject    handle to pushbutton_submit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model = handles.model;
elementArray = model.getAllElements();
youngsModulus = str2num(get(handles.edit_youngsModulus,'string'));
thickness = str2num(get(handles.edit_thickness,'string'));
prxy = str2num(get(handles.edit_prxy,'string'));
nrGaussPoint = str2num(get(handles.edit_numberGaussPoint,'string'));

elementArray.setPropertyValue('THICKNESS', thickness);
elementArray.setPropertyValue('YOUNGS_MODULUS', youngsModulus);
elementArray.setPropertyValue('POISSON_RATIO', prxy);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', nrGaussPoint);

set(findall(handles.uipanel_bC, '-property', 'enable'), 'enable', 'on');

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
nodeArray = model.getAllNodes();
solver = SimpleSolvingStrategy(model);
solver.solve();
% evaluation_point = model.getModelPart('LMPoint').getId();
% dy = nodeArray(evaluation_point).getDofValue('DISPLACEMENT_Y');
% fprintf("Displacement at Node %i: %f\n", evaluation_point, dy);

% set(findall(handles.uipanel7, '-property', 'enable'), 'enable', 'on');

cla(handles.axes1);
colorbar off
vis = VisualizationGUI(model);
vis.plotUndeformed();
vis.plotDeformed();

set(handles.popupmenu2,'Value',1);
set(handles.checkbox3,'Value',1);
set(handles.checkbox4,'Value',0);
set(handles.radiobutton1,'Value',1);
set(findall(handles.uibuttongroup_vis, '-property', 'enable'), 'enable', 'on');
set(handles.popupmenu2,'enable','off');
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
    if(get(hObject,'Value') == get(hObject,'Max'))
        set(findobj(gcf,'Tag','deformed'),'Visible','on');                     
    else
        set(findobj(gcf,'Tag','deformed'),'Visible','off');
    end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
axes1 = handles.axes1;
children = axes1.Children;
set(children,'Visible','off');

model = handles.model;
vis = VisualizationGUI(model);
contents = cellstr(get(hObject,'String'));
fieldType = contents{get(hObject,'Value')};

h = findobj(gcf,'Tag',fieldType);
if ~isempty(h)
    set(h,'Visible','on');
    colorbar
else   
    if ~strcmp(fieldType,'Select Field')
        vis.plotField(fieldType);
    else
        colorbar off
    end
end

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
        set(findobj(gcf,'Tag',fieldType),'Visible','on');
        colorbar
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
handles.loaded = false;
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
