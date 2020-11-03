function varargout = SelectBehavioralStateGUI_eLife2020(varargin)
% SELECTBEHAVIORALSTATEGUI_IOS_eLife2020 MATLAB code for SelectBehavioralStateGUI_IOS_eLife2020.fig
%      SELECTBEHAVIORALSTATEGUI_IOS_eLife2020, by itself, creates a new SELECTBEHAVIORALSTATEGUI_IOS_eLife2020 or raises the existing
%      singleton*.
%
%      H = SELECTBEHAVIORALSTATEGUI_IOS_eLife2020 returns the handle to a new SELECTBEHAVIORALSTATEGUI_IOS_eLife2020 or the handle to
%      the existing singleton*.
%
%      SELECTBEHAVIORALSTATEGUI_IOS_eLife2020('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTBEHAVIORALSTATEGUI_IOS_eLife2020.M with the given input arguments.
%
%      SELECTBEHAVIORALSTATEGUI_IOS_eLife2020('Property','Value',...) creates a new SELECTBEHAVIORALSTATEGUI_IOS_eLife2020 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SelectBehavioralStateGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SelectBehavioralStateGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SelectBehavioralStateGUI_IOS_eLife2020

% Last Modified by GUIDE v2.5 28-Feb-2020 21:03:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SelectBehavioralStateGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SelectBehavioralStateGUI_OutputFcn, ...
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


% --- Executes just before SelectBehavioralStateGUI is made visible.
function SelectBehavioralStateGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectBehavioralStateGUI_IOS_eLife2020 (see VARARGIN)

% Choose default command line output for SelectBehavioralStateGUI_IOS_eLife2020
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SelectBehavioralStateGUI_IOS_eLife2020 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SelectBehavioralStateGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ButtonSelect_IOS_eLife2020
% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ButtonSelect_IOS_eLife2020
% Hint: get(hObject,'Value') returns toggle state of togglebutton2


% --- Executes on button press in togglebutton3.
function togglebutton3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ButtonSelect_IOS_eLife2020
% Hint: get(hObject,'Value') returns toggle state of togglebutton3
