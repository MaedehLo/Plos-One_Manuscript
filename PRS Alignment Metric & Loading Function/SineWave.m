close all
clear
clc

Loc=22;   % number of locations we use from Raman data
A_Loc=11;  % number of locations with good Amplitude Alignment Factor
xc=75;   % phase angle 
alpha=[0,30,60,90,120,150,180];
S=xlsread('SineParameters.xlsx')
B=xlsread('fin_PC1.xlsx')
A=S(:,1);
y0=S(:,2);
x=0:1:180;
figure
for i=1:1:A_Loc
    y=y0(i)+A(i)*sin(pi*(x-xc)/90);
    plot(x,y)
    hold on
end
  for j=1:1:Loc
    plot(alpha,B(:,j),'*black')
    hold on
  end
axis([-10 190 -3 3])
xlabel('Polarized Angle')
ylabel('Adjused PC1 Score')