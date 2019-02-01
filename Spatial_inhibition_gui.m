function varargout = Spatial_inhibition_gui(varargin)
% SPATIAL_INHIBITION_GUI MATLAB code for Spatial_inhibition_gui.fig
%      SPATIAL_INHIBITION_GUI, by itself, creates a new SPATIAL_INHIBITION_GUI or raises the existing
%      singleton*.
%
%      H = SPATIAL_INHIBITION_GUI returns the handle to a new SPATIAL_INHIBITION_GUI or the handle to
%      the existing singleton*.
%
%      SPATIAL_INHIBITION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPATIAL_INHIBITION_GUI.M with the given input arguments.
%
%      SPATIAL_INHIBITION_GUI('Property','Value',...) creates a new SPATIAL_INHIBITION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Spatial_inhibition_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Spatial_inhibition_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Spatial_inhibition_gui

% Last Modified by GUIDE v2.5 01-Feb-2019 05:42:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Spatial_inhibition_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @Spatial_inhibition_gui_OutputFcn, ...
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

global data;


% --- Executes just before Spatial_inhibition_gui is made visible.
function Spatial_inhibition_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Spatial_inhibition_gui (see VARARGIN)

% Choose default command line output for Spatial_inhibition_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Spatial_inhibition_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Spatial_inhibition_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global data;
update(data,handles);
% spine = data.spineCmpt(get(handles.listbox1,'Value'));
% excLevel = interp1(data.stimRange,get(handles.slider1,'Value'),'nearest');
% shaft = get(handles.listbox2,'Value');
% if length(data.inhRange)>1
%     inhLevel = interp1(data.inhRange,get(handles.slider2,'Value'),'nearest');
% else
%     inhLevel = data.inhRange;
% end
% 
% set(handles.text2,'String',['Excitation intensity = ' num2str(excLevel)]);
% set(handles.text3,'String',['Inhibition intensity = ' num2str(inhLevel)]);
% 
% 
% spineIndex = find(data.spineCmpt==spine);
% excIndex = find(data.stimRange==excLevel);
% inhIndex = find(data.inhRange==inhLevel);
% 
% cla(handles.axes1);
% map = jet(1000);
% high = max(data.inhDiff(spineIndex,:,excIndex,inhIndex));
% low = min(data.inhDiff(spineIndex,:,excIndex,inhIndex));
% h = plot(handles.axes1,data.g);
% for c = 1:length(data.shaftCmpt)
%     inhCOI = data.shaftCmpt(c);
%     difference = data.inhDiff(spineIndex,inhCOI,excIndex,inhIndex);
%     if (high-low)>0
%         relDiff = round(999*(difference-low)/(high-low)+1);
%     else
%         relDiff = 1;
%     end
%     cmptColor = map(relDiff,:);
%     highlight(h,inhCOI);
%     highlight(h,inhCOI,'NodeColor',cmptColor);
%     title(['Excitation on ' num2str(spine) ', intensity = ' num2str(excLevel) ', inhLevel = ' num2str(inhLevel)]);
% end
% highlight(h,spine,'Marker','s');
% highlight(h,spine);
% highlight(h,spine);
% 
% cla(handles.axes2);
% plot(handles.axes2,squeeze(data.noInhibition(spineIndex,excIndex,:)));
% hold(handles.axes2,'on');
% for i = 1:length(shaft)
%     plot(handles.axes2,squeeze(data.yesInhibition(spineIndex,shaft(i),excIndex,inhIndex,:)));
% end
% list = get(handles.listbox2,'String');
% legend(handles.axes2,'No inh',list{shaft});


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
global data;
update(data,handles);
% spine = data.spineCmpt(get(handles.listbox1,'Value'));
% excLevel = interp1(data.stimRange,get(handles.slider1,'Value'),'nearest');
% shaft = get(handles.listbox2,'Value');
% if length(data.inhRange)>1
%     inhLevel = interp1(data.inhRange,get(handles.slider2,'Value'),'nearest');
% else
%     inhLevel = data.inhRange;
% end
% 
% set(handles.text2,'String',['Excitation intensity = ' num2str(excLevel)]);
% set(handles.text3,'String',['Inhibition intensity = ' num2str(inhLevel)]);
% 
% 
% spineIndex = find(data.spineCmpt==spine)
% excIndex = find(data.stimRange==excLevel)
% inhIndex = find(data.inhRange==inhLevel);
% 
% cla(handles.axes1);
% map = jet(1000);
% high = max(data.inhDiff(spineIndex,:,excIndex,inhIndex))
% low = min(data.inhDiff(spineIndex,:,excIndex,inhIndex))
% h = plot(handles.axes1,data.g);
% for c = 1:length(data.shaftCmpt)
%     inhCOI = data.shaftCmpt(c);
%     difference = data.inhDiff(spineIndex,inhCOI,excIndex,inhIndex);
%     if (high-low)>0
%         relDiff = round(999*(difference-low)/(high-low)+1);
%     else
%         relDiff = 1;
%     end
%     cmptColor = map(relDiff,:)
%     highlight(h,inhCOI);
%     highlight(h,inhCOI,'NodeColor',cmptColor);
%     title(['Excitation on ' num2str(spine) ', intensity = ' num2str(excLevel) ', inhLevel = ' num2str(inhLevel)]);
% end
% highlight(h,spine,'Marker','s');
% highlight(h,spine);
% highlight(h,spine);
% 
% cla(handles.axes2);
% plot(handles.axes2,squeeze(data.noInhibition(spineIndex,excIndex,:)));
% hold(handles.axes2,'on');
% for i = 1:length(shaft)
%     plot(handles.axes2,squeeze(data.yesInhibition(spineIndex,shaft(i),excIndex,inhIndex,:)));
% end
% list = get(handles.listbox2,'String');
% legend(handles.axes2,'No inh',list{shaft});


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
global data;
update(data,handles);
% spine = data.spineCmpt(get(handles.listbox1,'Value'));
% excLevel = interp1(data.stimRange,get(handles.slider1,'Value'),'nearest');
% shaft = get(handles.listbox2,'Value');
% if length(data.inhRange)>1
%     inhLevel = interp1(data.inhRange,get(handles.slider2,'Value'),'nearest');
% else
%     inhLevel = data.inhRange;
% end
% 
% set(handles.text2,'String',['Excitation intensity = ' num2str(excLevel)]);
% set(handles.text3,'String',['Inhibition intensity = ' num2str(inhLevel)]);
% 
% 
% spineIndex = find(data.spineCmpt==spine);
% excIndex = find(data.stimRange==excLevel);
% inhIndex = find(data.inhRange==inhLevel);
% 
% cla(handles.axes1);
% map = jet(1000);
% high = max(data.inhDiff(spineIndex,:,excIndex,inhIndex));
% low = min(data.inhDiff(spineIndex,:,excIndex,inhIndex));
% h = plot(handles.axes1,data.g);
% for c = 1:length(data.shaftCmpt)
%     inhCOI = data.shaftCmpt(c);
%     difference = data.inhDiff(spineIndex,inhCOI,excIndex,inhIndex);
%     if (high-low)>0
%         relDiff = round(999*(difference-low)/(high-low)+1);
%     else
%         relDiff = 1;
%     end
%     cmptColor = map(relDiff,:);
%     highlight(h,inhCOI);
%     highlight(h,inhCOI,'NodeColor',cmptColor);
%     title(['Excitation on ' num2str(spine) ', intensity = ' num2str(excLevel) ', inhLevel = ' num2str(inhLevel)]);
% end
% highlight(h,spine,'Marker','s');
% highlight(h,spine);
% highlight(h,spine);
% 
% cla(handles.axes2);
% plot(handles.axes2,squeeze(data.noInhibition(spineIndex,excIndex,:)));
% hold(handles.axes2,'on');
% for i = 1:length(shaft)
%     plot(handles.axes2,squeeze(data.yesInhibition(spineIndex,shaft(i),excIndex,inhIndex,:)));
% end
% list = get(handles.listbox2,'String');
% legend(handles.axes2,'No inh',list{shaft});

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global data;
update(data,handles);
% spine = data.spineCmpt(get(handles.listbox1,'Value'));
% excLevel = interp1(data.stimRange,get(handles.slider1,'Value'),'nearest');
% shaft = get(handles.listbox2,'Value');
% if length(data.inhRange)>1
%     inhLevel = interp1(data.inhRange,get(handles.slider2,'Value'),'nearest');
% else
%     inhLevel = data.inhRange;
% end
% 
% set(handles.text2,'String',['Excitation intensity = ' num2str(excLevel)]);
% set(handles.text3,'String',['Inhibition intensity = ' num2str(inhLevel)]);
% 
% 
% spineIndex = find(data.spineCmpt==spine);
% excIndex = find(data.stimRange==excLevel);
% inhIndex = find(data.inhRange==inhLevel);
% 
% cla(handles.axes1);
% map = jet(1000);
% high = max(data.inhDiff(spineIndex,:,excIndex,inhIndex));
% low = min(data.inhDiff(spineIndex,:,excIndex,inhIndex));
% h = plot(handles.axes1,data.g);
% for c = 1:length(data.shaftCmpt)
%     inhCOI = data.shaftCmpt(c);
%     difference = data.inhDiff(spineIndex,inhCOI,excIndex,inhIndex);
%     if (high-low)>0
%         relDiff = round(999*(difference-low)/(high-low)+1);
%     else
%         relDiff = 1;
%     end
%     cmptColor = map(relDiff,:);
%     highlight(h,inhCOI);
%     highlight(h,inhCOI,'NodeColor',cmptColor);
%     title(['Excitation on ' num2str(spine) ', intensity = ' num2str(excLevel) ', inhLevel = ' num2str(inhLevel)]);
% end
% highlight(h,spine,'Marker','s');
% highlight(h,spine);
% highlight(h,spine);
% 
% cla(handles.axes2);
% plot(handles.axes2,squeeze(data.noInhibition(spineIndex,excIndex,:)));
% hold(handles.axes2,'on');
% for i = 1:length(shaft)
%     plot(handles.axes2,squeeze(data.yesInhibition(spineIndex,shaft(i),excIndex,inhIndex,:)));
% end
% list = get(handles.listbox2,'String');
% legend(handles.axes2,'No inh',list{shaft});

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data;

[file folder] = uigetfile('/Users/Aaron/Documents');
data = load([folder '/' file]);
% 

set(handles.slider1,'value',data.stimRange(1));
set(handles.slider1,'min',data.stimRange(1));
set(handles.slider1,'max',data.stimRange(end));
set(handles.text2,'String',['Excitation intensity = ' num2str(data.stimRange(1))]);

set(handles.slider2,'value',data.inhRange(1));
set(handles.slider2,'min',data.inhRange(1));
if length(data.inhRange)==1
    set(handles.slider2,'max',data.inhRange(end)+1);
else
    set(handles.slider2,'max',data.inhRange(end));
end

set(handles.text3,'String',['Inhibition intensity = ' num2str(data.inhRange(1))]);

for i = 1:size(data.spineCmpt,2)
    spineList{i} = num2str(data.spineCmpt(i));
end
set(handles.listbox1,'String',spineList);

for i = 1:size(data.shaftCmpt,2)
    shaftList{i} = num2str(data.shaftCmpt(i));
    if data.shaftCmpt(i) == 1
        shaftList{i} = 'Soma';
    end
end
set(handles.listbox2,'String',shaftList);
set(handles.listbox2,'max',1000,'min',0);


function update(dat, hand)

spine = dat.spineCmpt(get(hand.listbox1,'Value'));
excLevel = interp1(dat.stimRange,dat.stimRange,get(hand.slider1,'Value'),'nearest');
shaft = get(hand.listbox2,'Value');
if length(dat.inhRange)>1
    inhLevel = interp1(dat.inhRange,dat.inhRange,get(hand.slider2,'Value'),'nearest');
else
    inhLevel = dat.inhRange;
end

set(hand.text2,'String',['Excitation intensity = ' num2str(excLevel)]);
set(hand.text3,'String',['Inhibition intensity = ' num2str(inhLevel)]);


spineIndex = find(dat.spineCmpt==spine);
excIndex = find(dat.stimRange==excLevel);
inhIndex = find(dat.inhRange==inhLevel);

cla(hand.axes1);
map = jet(1000);
high = max(dat.inhDiff(spineIndex,:,excIndex,inhIndex));
low = min(dat.inhDiff(spineIndex,:,excIndex,inhIndex));
h = plot(hand.axes1,dat.g);
for c = 1:length(dat.shaftCmpt)
    inhCOI = dat.shaftCmpt(c);
    difference = dat.inhDiff(spineIndex,inhCOI,excIndex,inhIndex);
    if (high-low)>0
        relDiff = round(999*(difference-low)/(high-low)+1);
    else
        relDiff = 1;
    end
    cmptColor = map(relDiff,:);
    highlight(h,inhCOI);
    highlight(h,inhCOI,'NodeColor',cmptColor);
%     title(['Excitation on ' num2str(spine) ', intensity = ' num2str(excLevel) ', inhLevel = ' num2str(inhLevel)]);
end
highlight(h,spine,'Marker','s');
highlight(h,spine);
highlight(h,spine);

cla(hand.axes2);
plot(hand.axes2,squeeze(dat.noInhibition(spineIndex,excIndex,:)));
hold(hand.axes2,'on');
for i = 1:length(shaft)
    plot(hand.axes2,squeeze(dat.yesInhibition(spineIndex,shaft(i),excIndex,inhIndex,:)));
end
list = get(hand.listbox2,'String');
legend(hand.axes2,'No inh',list{shaft});
