function varargout = BehaviorStateButtonGUI_SVM(varargin)
% BEHAVIORSTATEBUTTONGUI_SVM MATLAB code for BehaviorStateButtonGUI_SVM.fig
%      BEHAVIORSTATEBUTTONGUI_SVM, by itself, creates a new BEHAVIORSTATEBUTTONGUI_SVM or raises the existing
%      singleton*.
%
%      H = BEHAVIORSTATEBUTTONGUI_SVM returns the handle to a new BEHAVIORSTATEBUTTONGUI_SVM or the handle to
%      the existing singleton*.
%
%      BEHAVIORSTATEBUTTONGUI_SVM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEHAVIORSTATEBUTTONGUI_SVM.M with the given input arguments.
%
%      BEHAVIORSTATEBUTTONGUI_SVM('Property','Value',...) creates a new BEHAVIORSTATEBUTTONGUI_SVM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BehaviorStateButtonGUI_SVM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BehaviorStateButtonGUI_SVM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BehaviorStateButtonGUI_SVM

% Last Modified by GUIDE v2.5 29-Jul-2019 12:44:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BehaviorStateButtonGUI_SVM_OpeningFcn, ...
                   'gui_OutputFcn',  @BehaviorStateButtonGUI_SVM_OutputFcn, ...
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

% --- Executes just before BehaviorStateButtonGUI_SVM is made visible.
function BehaviorStateButtonGUI_SVM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BehaviorStateButtonGUI_SVM (see VARARGIN)

% Choose default command line output for BehaviorStateButtonGUI_SVM

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BehaviorStateButtonGUI_SVM wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = BehaviorStateButtonGUI_SVM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ButtonSelect_SVM

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ButtonSelect_SVM

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ButtonSelect_SVM

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ButtonSelect_SVM


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1
