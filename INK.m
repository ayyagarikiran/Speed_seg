function varargout = INK(varargin)
% INK MATLAB code for INK.fig
%      INK, by itself, creates a new INK or raises the existing
%      singleton*.
%
%      H = INK returns the handle to a new INK or the handle to
%      the existing singleton*.
%
%      INK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INK.M with the given input arguments.
%
%      INK('Property','Value',...) creates a new INK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before INK_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to INK_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help INK

% Last Modified by GUIDE v2.5 28-Jan-2016 12:01:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @INK_OpeningFcn, ...
                   'gui_OutputFcn',  @INK_OutputFcn, ...
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


% --- Executes just before INK is made visible.
function INK_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to INK (see VARARGIN)

% Choose default command line output for INK
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes INK wait for user response (see UIRESUME)
% uiwait(handles.figure1);


global drawing;
drawing =0;
set(gcf,'WindowButtonDownFcn',@mouseDown)
set(gcf,'WindowButtonMotionFcn',@mouseMove)
set(gcf,'WindowButtonUpFcn',@mouseUp)

global pnt
global Npnt
pnt = zeros(1000,3);
Npnt = 0;
tic

% --- Outputs from this function are returned to the command line.
function varargout = INK_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ClearButton.
function ClearButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla
global pnt
global Npnt
pnt = zeros(1000,3);
Npnt = 0;

% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pnt
global Npnt
global speed1
global speed2
global disti
global A
global arclen
global B
global segpoints
if Npnt<1000 
pnt(Npnt+1:end,:) =[];
end
dlmwrite('speed.txt',speed1)
avgspd = disti(end)/(pnt(end,3)-pnt(1,3));
threshold=avgspd/3;
curvthr = 50;
[pk,loc] = findpeaks(-speed2);
seg = find(pk>-threshold);
locs = loc(seg);
% figure()
% plot(speed2)
B = smooth(abs(diff(A)./diff(arclen)),'moving');
[pks,loccur] = findpeaks(B);
seg = find(pks>curvthr);
locur =loccur(seg)+5;
% figure()
% plot(B)
segments=unique(sort([locs;locur]));
for i=1:numel(segments)-1
    if abs(segments(i)-segments(i+1))<5
        Pnts(i) = 0;
    else
        Pnts(i) = segments(i);
    end
    
    if i == numel(segments)-1
        if abs(Pnts(i)-segments(i+1))<5
            break
        else
            Pnts(i+1) = segments(i+1);
        end
    end
end

segpoints=nonzeros(Pnts);
segpoints=[1;segpoints;length(pnt)];
for i=1:length(segpoints)-1
            N=segpoints(i+1)-segpoints(i);
            part=pnt(segpoints(i):segpoints(i+1),:);
            xs =sum(part(:,1));
            ys = sum(part(:,2));
            xys =sum(part(:,1).*part(:,2));
            xsq = sumsqr(part(:,1));
           slope = ([xsq xs;xs 5])\[xys;ys];
           if abs(slope(1))>90
                err = (1/5)*sum(abs(((part(:,2)-slope(2))./slope(1))-part(:,1)));
            else
                err = (1/5)*sum(abs((slope(1)*(part(:,1)+slope(2))-part(:,2))));
            end
                xs = sum(part(:,1));
                ys = sum(part(:,2));
                xys = sum(part(:,1).*part(:,2));
                xsq = sumsqr(part(:,1));
                ysq = sumsqr(part(:,2));
                r1 = sum(-((part(:,1).^2+part(:,2).^2).*part(:,1)));
                r2 = sum(-(part(:,1).^2+part(:,2).^2).*part(:,2));
                r3 = sum(-((part(:,1).^2+part(:,2).^2)));
                cirfit = ([2*xsq 2*xys xs;2*xys 2*ysq ys;2*xs 2*ys 5])\[r1;r2;r3];
                radius = sqrt(cirfit(1)^2+cirfit(2)^2-cirfit(3));
                err1 = (1/N)*(sum(sqrt((part(:,1)+cirfit(1)).^2+(part(:,2)+cirfit(2)).^2)-radius));
                if abs(err)<abs(err1)
                    p_fit=polyfit(part(:,1),part(:,2),1);
                    xi=linspace(part(1,1),part(end,1));
                    w=p_fit(1).*xi+p_fit(2)
                    figure()
                   
                    %plot(xi,w)
                    plot(xi,w,'r','linewidth',4)
                    hold on
                    
                   
                else    
                    p_fit=polyfit(part(:,1),part(:,2),2);
                    xi=linspace(part(1,1),part(end,1));
                    w=p_fit(1).*xi.^2++p_fit(2).*xi+p_fit(3);
                    hold on
                    %plot(xi,w)
                    plot(xi,w,'r','linewidth',4)
                    
                end
                
end

           


function mouseDown(hObject, eventdata, handles) 
global drawing
drawing = 1;

function mouseUp(hObject, eventdata, handles) 
global drawing
drawing = 0;

function mouseMove(hObject, eventdata, handles) 
global drawing
global Npnt
global pnt
global speed
global speed1
global speed2
global dist
global time
global timei
global disti
global A
global arclen

if drawing
    C = get(gca,'CurrentPoint');
    if C(1,1)<1 && C(1,1)>0 && C(1,2)<1 && C(1,2)>0
        Npnt = Npnt+1;
        pnt(Npnt,1) = C(1,1);
        pnt(Npnt,2) = C(1,2);
        pnt(Npnt,3) = toc;
        plot(C(1,1),C(1,2),'k','marker','o','MarkerFaceColor','r')
        hold on
        xlim([0 1]); ylim([0 1]);
        set(gca,'XTick',[],'YTick',[])
        box on
        dist(1)=0;
        time(1)=0;
        if Npnt>1
            dist(Npnt)=distance(pnt(Npnt,1:2),pnt((Npnt-1),1:2));
            time(Npnt)=timediff(pnt(Npnt,3),pnt((Npnt-1),3));
        end 
        disti=cumsum(dist);
        timei=cumsum(time);
        if Npnt>4
            xs = sum(pnt(Npnt-4:Npnt,1));
            ys = sum(pnt(Npnt-4:Npnt,2));
            xys =sum(pnt(Npnt-4:Npnt,1).*pnt(Npnt-4:Npnt,2));
            xsq = sumsqr(pnt(Npnt-4:Npnt,1));
           slope = ([xsq xs;xs 5])\[xys;ys];

            if abs(slope(1))>90
                err = (1/5)*sum(abs(((pnt(Npnt-4:Npnt,2)-slope(2))./slope(1))-pnt(Npnt-4:Npnt,1)));
            else
                err = (1/5)*sum(abs((slope(1)*pnt(Npnt-4:Npnt,1)+slope(2))-pnt(Npnt-4:Npnt,2)));
            end
            arclen(Npnt-4) = disti(Npnt);
            if err <= 0.1*arclen(Npnt-4)
                A(Npnt-4) = (atan(slope(1)));
            else
                xs = sum(pnt(Npnt-4:Npnt,1));
                ys = sum(pnt(Npnt-4:Npnt,2));
                xys = sum(pnt(Npnt-4:Npnt,1).*pnt(Npnt-4:Npnt,2));
                xsq = sumsqr(pnt(Npnt-4:Npnt,1));
                ysq = sumsqr(pnt(Npnt-4:Npnt,2));
                r1 = sum(-((pnt(Npnt-4:Npnt,1).^2+pnt(Npnt-4:Npnt,2).^2).*pnt(Npnt-4:Npnt,1)));
                r2 = sum(-((pnt(Npnt-4:Npnt,1).^2+pnt(Npnt-4:Npnt,2).^2).*pnt(Npnt-4:Npnt,2)));
                r3 = sum(-((pnt(Npnt-4:Npnt,1).^2+pnt(Npnt-4:Npnt,2).^2)));
                cirfit = ([2*xsq 2*xys xs;2*xys 2*ysq ys;2*xs 2*ys 5])\[r1;r2;r3];
                x = pnt(Npnt-3,1);
                y = (-cirfit(2)+sqrt((cirfit(2)^2)-(cirfit(3)-(x^2)-2*cirfit(1)*(x))))/(2*cirfit(1));
                A(Npnt-4) = (atan(-(x+cirfit(1))/(y+cirfit(2))));
            end
        end
    end
    
end

speed = (disti(3:end)-disti(1:end-2))./(timei(3:end)-timei(1:end-2));
speed1(1)=0;
if isempty(disti)==0
    speed1(length(disti)) = 0;
end
speed1(2:length(disti)-1)=speed;
speed2 = smooth(speed1,'moving');



function d=distance(pnt1,pnt2)
d=norm(pnt1(1,1:2)-pnt2(1,1:2));

function t=timediff(pnt1,pnt2)
t=pnt1-pnt2;



