close all
clear all
clc
%%
%Preprocess Raman spectra (11th sub) and PCA analysis
%Read Raman text files respectively from 0 to 180 degree(7 polarized angle)
%Each polarized angle consists of ROI1(11 locations)&ROI2(11 locations)
%Putting (0 degree) data from loc 1 to loc 22 into Matrix Y
for ii=0:30:180                                                             
filename1 = sprintf('M%d_1.txt',ii);       % %d should replace with diffrernt angle 0,30, ... in the file name
filename2 = sprintf('M%d_2.txt',ii);
filename3 = sprintf('M%d_pbs.txt',ii);
M1 = dlmread(filename1,'',1,0);       %read raw ROI1 into M1
M2 = dlmread(filename2,'',1,0);       %read raw ROI2 into M2
PBS = dlmread(filename3,'',1,0);      %read raw PBS into matrix PBS

X= M1(751:1153,3);             %read wavenumber from 1400cm-1 to 1800 cm-1

%%
%Consider for each location and change based on each Raman data direction
% it can  be 0 or -12 based on direction
%adjust spcetral data order for each polarized angle:
%"0" means keeping original order while "-12" means keep opposite order
%adj1=[0 0 0 0 -12 -12 -12];    %adj1:adjust ROI1 for 7 polarized angles
%adj2=[0 0 0 0 -12 -12 -12];    %adj2:adjust ROI2 for 7 polarized angles

adj1=[0 0 0 0 0 0 0];
adj2=[0 0 0 0 0 0 0];

%%
% read ROI1, ROI2 and PBS from 1400-1800cm-1 range
% Y is for all 22 location (RoI 1 & RoI 2)

% RoI 1
for i =1:1:11 
    i0 = abs(i+adj1(ii/30+1)); 
    Y(:,i) = M1((750+1600*(i0-1)):(1152+1600*(i0-1)),4); %ROI1
end

% RoI 2
for j=1:1:11 
    j0= abs(j+adj2(ii/30+1));               
    Y(:,11+j)= M2((750+1600*(j0-1)):(1152+1600*(j0-1)),4); %ROI2
end

%PBS (k: number of iteration for measurments in Raman)
for k = 1:5
    P(:,k) = PBS((750+1600*(k-1)):(1152+1600*(k-1)),3);  % read pbs Intensity (see it is column 2 or 3 !)
end
P_ave = mean(P,2);             % Average of PBS background for each row, 2 is related to mean
YS=Y-P_ave;                    % Subtraction of PBS background

%%
% Smoothing by Savitzky-Golay filtering with 3th order, window size=11
order = 3;
framelen = 11;
sgfYS = sgolayfilt(YS,order,framelen);
plot(X,YS,':')                 % Plot 11th sub spectra from Raman
hold on
plot(X,sgfYS,'.-')             % Plot smoothed spectra
hold on
legend('signal','sgolay')
title('Unnormalized')
drawnow

%%
% SNV normalization for 2 ROI*11 locations= 22 spectra
hFig=figure;
Mean=mean(sgfYS);
Std=std(sgfYS);
SNV_Y= (sgfYS-Mean)./Std;
plot(X,SNV_Y,'-');
legend show
SNVY(1:403,1+22*(ii/30):22*(ii/30+1))=SNV_Y;
title('SNV Normalized') 
end

xlswrite('SNVY.xlsx',SNVY)      % Save preprocessed 22 spectra in excel
close(hFig)

% plot preprocessed Raman spectra of 22 locations for each angle
for t= 1:7
    figure(8+t)
    plot(X,SNVY(:,1+22*(t-1):22*t))
    t1=30*(t-1);
    title(sprintf('polarized angle%d',t1))
    legend show
end

%%
% remove outlier
% m1=144;m2=11;
% n=22;
% if m1>n
%     o1=rem(m1,n);
% else
%     o1=m1
% end
% if m2>n
%     o2=rem(m2,n);
% else
%     o2=m2
% end
% if o1<o2
%     SNVY(:,o1:o1:end)=[];
% end

%%
%cut off outliers e.g.v, w and ...(v<w) in sequence
% v1 & v2 & ... are the numbers of bad location
e=22;   % total location of RoI1 & RoI2
%v1=1;
%v2=2;
%SNVY(:,[v1 v2 v1+e v2+e v1+2*e v2+2*e v1+3*e v2+3*e v1+4*e v2+4*e v1+5*e v2+5*e v1+6*e v2+6*e])=[];
%v3=3;
%v4=4;
%SNVY(:,[v3 v4 v3+e v4+e v3+2*e v4+2*e v3+3*e v4+3*e v3+4*e v4+4*e v3+5*e v4+5*e v3+6*e v4+6*e])=[];

%for j=1:11
%    i=12-j;
%    SNVY(:,[i i+e i+2*e i+3*e i+4*e i+5*e i+6*e])=[];
%    e=e-1;
%end
%%
%plot mean spectra for 7 polarized angles
%Loc=22-off_Loc;
Loc=22;
figure
for l=1:7
    s=mean(SNVY(:,1+Loc*(l-1):Loc*l),2);
    plot(X,s);
    hold on
end
legend('0','30','60','90','120','150','180')
xlabel('Wavenumber(cm^{-1})')
ylabel('Normalized Raman Intensity (a.u.)')
%%
% PCA analysis
% Get PC1 Score by subtracting mu from SNV_Y data and multiplying by coeff
coeff = xlsread('Loading coeff.xlsx');
mu = mean(SNVY,2);
SNVYS = SNVY-mu;
coeffPC = coeff.';
PC1_Score = coeffPC*SNVYS;
u = mean(PC1_Score); % standarized PC1 Score
SD = std(PC1_Score);
PC1_std = (PC1_Score-u)./SD;
PC1 = zeros(7,Loc);
for p = 1:1:Loc
   PC1(:,p) = PC1_std(1,Loc*((0:6))+p);
end 
angle = [0;30;60;90;120;150;180];
figure
plot(angle, PC1, '+b');
hold on
axis([-20 220 -3.5  3.5]);
xlabel('Polarized Angle (\circ)');                                            
ylabel('PC1 Score (a.u.)');
ave_PC1 =mean(PC1,2);
fin_PC1= [PC1 ave_PC1];
xlswrite('fin_PC1.xlsx',fin_PC1);%-->Save PC1 in excel, LAST COLUMN!!! is averaged PC1 Score per angle

%% Save PC1 of different locations in different excels for sine fitting
% in originlab
xlRange = 'A1:B7';
sheet = 1;
Angle =[0; 30 ;60 ;90; 120; 150; 180];
for q= 1:1:Loc  %--> 22 locations if w.o/ excluding outlines
   PC=[Angle PC1(1+7*(q-1):7*q).'];
   filename = sprintf('PC1_Score%3d.xls', q);
   xlswrite(filename,PC,sheet,xlRange);
end
%%
%Find best angle Xc to fit sine fit 
diff_PC1=[1/2*(abs(ave_PC1(1)-ave_PC1(4))+abs(ave_PC1(4)-ave_PC1(7)));
abs(ave_PC1(2)-ave_PC1(5));
abs(ave_PC1(3)-ave_PC1(6))];
[Max,Index] = max(diff_PC1);
if (Index == 1)&&(ave_PC1(1)>ave_PC1(4))&&(abs(ave_PC1(1)-ave_PC1(4))>abs(ave_PC1(4)-ave_PC1(7)))
    disp('phase angle should be set: xc=-45 degree')
elseif (Index == 1)&&(ave_PC1(1)<ave_PC1(4))&&(abs(ave_PC1(1)-ave_PC1(4))>abs(ave_PC1(4)-ave_PC1(7)))
    disp('phase angle should be set: xc=45 degree')
elseif (Index == 1)&&(ave_PC1(4)<ave_PC1(7))&&(abs(ave_PC1(1)-ave_PC1(4))<abs(ave_PC1(4)-ave_PC1(7)))
    disp('phase angle should be set: xc=-45 degree')
elseif (Index == 1)&&(ave_PC1(4)>ave_PC1(7))&&(abs(ave_PC1(1)-ave_PC1(4))<abs(ave_PC1(4)-ave_PC1(7)))
     disp('phase angle should be set: xc=45 degree')
elseif (Index == 2)&&(ave_PC1(2)>ave_PC1(5))
    disp('phase angle should be set: xc=-15 degree')
elseif (Index == 2)&&(ave_PC1(2)<ave_PC1(5))
    disp('phase angle should be set: xc=75 degree')
elseif (Index == 3)&&(ave_PC1(3)>ave_PC1(6))
    disp('phase angle should be set: xc=15 degree')
else (Index == 3)&&(ave_PC1(3)<ave_PC1(6))
     disp('phase angle should be set: xc=-75 degree')
end



