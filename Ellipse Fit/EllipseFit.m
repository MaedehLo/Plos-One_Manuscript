clear
clc
%close all

% A should be a vertical vector as input to have a corroct answers

%A=xlsread('Corner Edge.xlsx');

%A=xlsread('Middle Edge.xlsx');

%A=xlsread('Middle Center.xlsx');

AA=xlsread('test5.xlsx');
A=AA';
%A(150,:)=[];
%%
Bin_Num=180/5;
figure(1)
h = histogram(A,Bin_Num,'Normalization','probability');
hold on

counts = h.Values;
th = 0:pi/Bin_Num:pi*(Bin_Num-1)/(Bin_Num);
r=[counts];
x = 1*r.*cos(th);
y= 1*r.*sin(th);

%%
figure(2)
for i=1:Bin_Num
  
 %   x(i)=cos(A(i,1)*pi/180);
 %   y(i)=sin(A(i,1)*pi/180);
    
   XY(i,1)=x(i);
   XY(i,2)=y(i);
   XY(i+Bin_Num,1)=-x(i);
   XY(i+Bin_Num,2)=-y(i);
    
    %plot([-x(i) x(i)],[-y(i) y(i)],'LineWidth',2)
   % plot([0 x(i)],[0 y(i)],'LineWidth',2)
    %hold on
   %plot(x(i),y(i),'*')
    hold on
    %plot(-x(i),-y(i),'*')
    hold on
end

%%
% Coefficiant in quadratic equation of Ellipse
ellipse=EllipseCoeff(XY);  

%%
% Plot Ellipse 1 (direct method)
[a,b,c,d,e,f]=deal(ellipse(1),ellipse(2),ellipse(3),ellipse(4),ellipse(5),ellipse(6));

delta=b*b-4*a*c;

Fun=@(X,Y)a*(X^2) + b*X*Y + c*(Y^2) +d*X + e*Y + f;
fimplicit(Fun,'black')
hold on
%ezplot( a*(X1^2) + b*X1*Y1 + c*(Y1^2) +d*X1 + e*Y1 + f)
%hold on

%%
% Ellipse major(M) and moinor(N) axis, Eccentricity(E) and orientation(p)
Axis=EllipseAxis(ellipse);
[m,n,E,phi]=deal(Axis(1),Axis(2),Axis(3),Axis(4));

 %%
% Plot Ellipse 2 (based on axis)
% Parameterize the equation
t = linspace(0, 360,1000);
phaseShift = 0; 
xAmplitude =  m;  
yAmplitude =  n ;  

G= [m n];
MM= max(G); % major axis
NN= min(G); % minor axis
s=MM/NN; % the ratio of major axis to minor axis

% Plot the rotated ellipse.
rotationAngle = phi*180/pi;    % input orientation angle

%% NOTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% angle of major axis (angle of mean axis in 0-180) for not symmetry fibers
if rotationAngle>0
alpha=180-rotationAngle;    
else 
alpha=-rotationAngle;
end

%%%% OR
% alpha is the average of the fibers for symmetry distribution of fibers
% alpha=mean(A);

%%
% Index Orientation S
S=Index_Orientation(A,alpha)

transformMatrix = [cosd(rotationAngle), sind(rotationAngle);...
  -sind(rotationAngle), cosd(rotationAngle)]; 

xAligned = xAmplitude * cosd(t);
yAligned = yAmplitude * sind(t);
xyAligned = [xAligned; yAligned];

xyRotated = transformMatrix * xyAligned;
xyRotated=xyRotated';
xRotated = xyRotated(:,1);
yRotated = xyRotated(:,2);

plot(x,y,'o','MarkerSize',7,'MarkerFaceColor','b');
hold on
plot(-x,-y,'o','MarkerSize',7,'MarkerFaceColor','b');
hold on
plot(xRotated, yRotated,'--r', 'LineWidth', 1);
hold on

%axis([-1 1 -1 1])
axis equal
grid on


