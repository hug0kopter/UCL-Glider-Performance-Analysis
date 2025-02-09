clc
clear all
close all

rho_SW=1027;
WingArea=1.112;%change base on your wing
folderPath='Put_The_Data_Here/CFD/';
FileList=dir(folderPath);
FileList=FileList(4:end);
SpeedList=[0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
NumberOfSpeed=length(SpeedList);
TemplateAOA=transpose(linspace(-7,7,15));
% length(TemplateAOA)

%BlankSpeedList=zeros(1,NumberOfFiles);
TableSpaceRef=zeros(length(TemplateAOA),NumberOfSpeed);
% AoAdebuglist=zeros(length(TemplateAOA),1);
% CDatSandAOA=zeros(length(TemplateAOA),1);
% CLatSandAOA=zeros(length(TemplateAOA),1);

FileName1=append(folderPath,FileList(1).name);
FileName2=append(folderPath,FileList(2).name);
FileName3=append(folderPath,FileList(3).name);


AoAdebuglist=importdata(FileName1);
CDatSandAOA=importdata(FileName2);
CLatSandAOA=importdata(FileName3);
    


TempDatSandAOA=1/2.*CDatSandAOA.*rho_SW*WingArea;
TempLatSandAOA=1/2.*CLatSandAOA.*rho_SW*WingArea;
DatSandAOA=[];
LatSandAOA=[];
for i=1:length(SpeedList)
DatSandAOA=[DatSandAOA,TempDatSandAOA(:,i).*SpeedList(i).^2];
LatSandAOA=[LatSandAOA,TempLatSandAOA(:,i).*SpeedList(i).^2];
end

NBatSandAOA=(DatSandAOA.^2+LatSandAOA.^2).^(1/2);
%% 
NB_Speed_Relationship=repmat(struct(),[length(TemplateAOA),1]);
UsefullAOARange=[];
for i=1:length(TemplateAOA)
NBForAOAi=NBatSandAOA(i,:);
UsefullAOARange=[UsefullAOARange,i];
y=transpose(SpeedList);
y2=transpose(LatSandAOA(i,:));
x=transpose(NBatSandAOA(i,:));
if i==7
xu=[x(1);x(2);x(3);x(5)];
y2u=[y2(1);y2(2);y2(3);y2(5)];
end 
f = fittype('a1*exp(b1*x) + a2*exp(b2*x) ');
start_point = [-0.4, 0.5, 0, 0]; 
fit_options = fitoptions('Method','NonlinearLeastSquares','StartPoint', start_point,'Algorithm','Levenberg-Marquardt'); 
[fit_result, gof] = fit(x, y, f, fit_options); 
disp(fit_result);
plot(fit_result,x,y);
NB_Speed_Relationship(i).FitFunction=fit_result;
disp('请检查拟合结果，Please check the fit result')

f2 = fittype('poly1');
if i==7
[fit_result2, gof2] = fit(xu, y2u, f2); 
else
[fit_result2, gof2] = fit(x, y2, f2); 
end
disp(fit_result2);
plot(fit_result2,x,y2);
NB_Speed_Relationship(i).FitFunction2=fit_result2;
disp('请检查拟合结果2，Please check the fit result2')
end

close all 
save('DAW_MST1CFD_RAM')
clear
%%
clear
load('DAW_MST1CFD_RAM.mat')
if true % to run the entire plot section
close all
Total_NB=14;
NB=[3,4,5,6,7,8,9,10,11];% here is the NB value you want to have a look can change and play
NBdown=Total_NB-NB;
Depth1=200;
UsefullAOA=TemplateAOA;%(UsefullAOARange);
UsefullNBSR=NB_Speed_Relationship;%(UsefullAOARange);
SpeedatNBandAOA=zeros(length(UsefullAOA),length(NB));
LatNBandAOA=zeros(length(UsefullAOA),length(NB));
BetaatNBandAOA=zeros(length(UsefullAOA),length(NB));
thetaatNBandAOA=zeros(length(UsefullAOA),length(NB));
UatNBandAOA=zeros(length(UsefullAOA),length(NB));
WatNBandAOA=zeros(length(UsefullAOA),length(NB));
movingDistance=zeros(length(UsefullAOA),length(NB));
movingDperL=zeros(length(UsefullAOA),length(NB));

for i=1:length(NB)
for j=1:length(UsefullAOA)
    SpeedatNBandAOA(j,i)=feval(UsefullNBSR(j).FitFunction, NB(i));
    LatNBandAOA(j,i)=feval(UsefullNBSR(j).FitFunction2, NB(i));
    BetaatNBandAOA(j,i)=acosd(LatNBandAOA(j,i)./NB(i));
    UatNBandAOA(j,i)=cosd(BetaatNBandAOA(j,i))*SpeedatNBandAOA(j,i);
    WatNBandAOA(j,i)=sind(BetaatNBandAOA(j,i))*SpeedatNBandAOA(j,i);
    thetaatNBandAOA(j,i)=BetaatNBandAOA(j,i)-UsefullAOA(j);
    movingDistance(j,i)=Depth1/tand(BetaatNBandAOA(j,i));
    movingDperL(j,i)=movingDistance(j,i)/(NB(i)*0.1*2/rho_SW*1000);%in L

end
end

% %put the useless AOA range here
% UsefullAOA=round(UsefullAOA,3);
uselessAOARange=[-1,0];
LB=find(UsefullAOA==uselessAOARange(1));
UB=find(UsefullAOA==uselessAOARange(2));

UsefullAOA=[UsefullAOA(1:LB-1);UsefullAOA(UB:end)];
SpeedatNBandAOA=[SpeedatNBandAOA(1:LB-1,:);SpeedatNBandAOA(UB:end,:)];
LatNBandAOA=[LatNBandAOA(1:LB-1,:);LatNBandAOA(UB:end,:)];
BetaatNBandAOA=[BetaatNBandAOA(1:LB-1,:);BetaatNBandAOA(UB:end,:)];
thetaatNBandAOA=[thetaatNBandAOA(1:LB-1,:);thetaatNBandAOA(UB:end,:)];
UatNBandAOA=[UatNBandAOA(1:LB-1,:);UatNBandAOA(UB:end,:)];
WatNBandAOA=[WatNBandAOA(1:LB-1,:);WatNBandAOA(UB:end,:)];
movingDistance=[movingDistance(1:LB-1,:);movingDistance(UB:end,:)];
movingDperL=[movingDperL(1:LB-1,:);movingDperL(UB:end,:)];

checkAOA0=find(UsefullAOA==0);
%general change
generalYlim=[0,0.2];
generalXlim=[0,1];
generalXlim2=[-1,0];
F1P2start=0;
%positive part

figure(1) 
hold on
grid on
%change if need negtive part
%xlim([-1,1])
xlim(generalXlim)
ylim(generalYlim)
Leg=cell(1,1);
PNB=[];

for i=1:length(NB)
plot(UatNBandAOA(LB:end,i),WatNBandAOA(LB:end,i))
STRNB2=append('NetB=',num2str(NB(i)),'N');
Leg{i}=STRNB2;
end

%change if need negtive part
istart=LB+F1P2start;
for i=istart:istart+5
plot([0,UatNBandAOA(i,:)],[0,WatNBandAOA(i,:)],'--')
STRAOA2=append('AOA=',num2str(UsefullAOA(i)),'deg');
Leg{i-istart+1+length(NB)}=STRAOA2;
end
%legend(Leg,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg,'NumColumns',1,'Location','northwest','Box','off')
title('Glide polar')
ylabel('Vertical Velocity w m/s')
xlabel('Horizontal Velocity u m/s')

figure(2)
hold on
grid on
Leg2=cell(1,1);
for i=1:length(NB)
plot(UsefullAOA(LB:end),thetaatNBandAOA(LB:end,i))
STR2=append('NetB=',num2str(NB(i)),'N');
Leg2{i}=STR2;
end
yline(0,'b--')
Leg2{length(NB)+1}='zero line';
%legend(Leg2,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg2,'NumColumns',1,'Box','off')
title('Pitch angle plot')
ylabel('Pitch Angle Theta deg')
xlabel('AOA alpha deg')
%% 

figure(3)
hold on
Leg3=cell(1,1);
contourf(NB,UsefullAOA(LB:end),movingDistance(LB:end,:))
colorbar
title('The distance per dive m')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(4)
hold on
Leg4=cell(1,1);
contourf(NB,UsefullAOA(LB:end),movingDperL(LB:end,:))
colorbar
title('The distance per dive every Liter change m/L')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(5)
hold on
grid on
Leg5=cell(1,1);
AOA1=UsefullAOA(LB:end);
MDL1=movingDperL(LB:end,:)./1000;
SpeedNBAOA1=SpeedatNBandAOA(LB:end,:);
plot3(MDL1,SpeedNBAOA1,AOA1)
ylim([0,1])
xlabel('The distance per dive every Liter change km')
ylabel('General Speed m/s')
zlabel('AOA alpha deg')
Leg2=Leg2(1:length(Leg2)-1);
%legend(Leg2,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg2,'NumColumns',1,'Box','off')
%%
figure(11)
hold on
Leg2=cell(1,1);
for i=1:length(NB)
plot(UsefullAOA(LB:end),BetaatNBandAOA(LB:end,i))
STR2=append('NetB=',num2str(NB(i)),'N');
Leg2{i}=STR2;
end
yline(0,'b--')
Leg2{length(NB)+1}='zero line';
%legend(Leg2,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg2,'NumColumns',1,'Box','off')
title('Glide angle plot')
ylabel('Glide Angle Beta deg')
xlabel('AOA alpha deg')

%%
figure(13)
hold on
grid on
Leg13=cell(1,1);
AOA1Start=LB;
AOA1end=LB+7;
AOA1=UsefullAOA(AOA1Start:AOA1end);
for i13=1:length(AOA1)
STRAOA13=append('AOA=',num2str(AOA1(i13)),'deg');
Leg13{i13}=STRAOA13;
end
L1=LatNBandAOA(AOA1Start:AOA1end,:);
NB1=NB;
SpeedNBAOA1=SpeedatNBandAOA(AOA1Start:AOA1end,:);
%[Xk,Yk]=meshgrid(NB1,AOA1);
plot3(L1,SpeedNBAOA1,NB1)
%ylim([0,1])
xlabel('Lift')
ylabel('General Speed m/s')
%zlabel('AOA alpha deg')
zlabel('NB')
title('Positive Part')
legend(Leg13,'NumColumns',1,'Location','northwest','Box','off')

%% 

%negtive part
%general change
generalYlim=[0,0.2];
generalXlim=[0,1];
generalXlim2=[-1,0];
F1P2start=0;
%positive part
figure(6) 
hold on
grid on
%change if need negtive part
%xlim([-1,1])
xlim(generalXlim)
ylim(generalYlim)
Leg=cell(1,1);
PNB=[];

for i=1:length(NB)
plot(-UatNBandAOA(1:LB-1,i),WatNBandAOA(1:LB-1,i))
STRNB2=append('NetB=',num2str(NB(i)),'N');
Leg{i}=STRNB2;
end

%change if need negtive part
istart=1+F1P2start;
for i=LB-istart-5:LB-istart
plot([0,-UatNBandAOA(i,:)],[0,WatNBandAOA(i,:)],'--')
STRAOA2=append('AOA=',num2str(UsefullAOA(i)),'deg');
Leg{i-(LB-istart-5)+1+length(NB)}=STRAOA2;
end
%legend(Leg,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg,'NumColumns',1,'Location','northwest','Box','off')
title('Glide polar')
ylabel('Vertical Velocity w m/s')
xlabel('Horizontal Velocity u m/s')










figure(7)
hold on
grid on
Leg2=cell(1,1);
for i=1:length(NB)
plot(UsefullAOA(1:LB-1),180-thetaatNBandAOA(1:LB-1,i))
STR2=append('NetB=',num2str(NB(i)),'N');
Leg2{i}=STR2;
end
yline(0,'b--')
Leg2{length(NB)+1}='zero line';
%legend(Leg2,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg2,'NumColumns',1,'Box','off')
title('Pitch angle plot')
ylabel('Pitch Angle Theta deg')
xlabel('AOA alpha deg')
%% 

figure(8)
hold on
Leg3=cell(1,1);
contourf(NB,UsefullAOA(1:LB-1),-movingDistance(1:LB-1,:))
colorbar
title('The distance per rise m')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(9)
hold on
Leg4=cell(1,1);
contourf(NB,UsefullAOA(1:LB-1),-movingDperL(1:LB-1,:))
colorbar
title('The distance per rise every Liter change m/L')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(10)
hold on
grid on
Leg5=cell(1,1);
AOA1=UsefullAOA(1:LB-1);
MDL1=-movingDperL(1:LB-1,:)./1000;
SpeedNBAOA1=SpeedatNBandAOA(1:LB-1,:);
plot3(MDL1,SpeedNBAOA1,AOA1)
ylim([0,1])
xlabel('The distance per rise every Liter change km')
ylabel('General Speed m/s')
zlabel('AOA alpha deg')
Leg2=Leg2(1:length(Leg2)-1);
%legend(Leg2,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg2,'NumColumns',1,'Box','off')
%%
figure(12)
hold on
Leg2=cell(1,1);
for i=1:length(NB)
plot(UsefullAOA(1:LB-1),180-BetaatNBandAOA(1:LB-1,i))
STR2=append('NetB=',num2str(NB(i)),'N');
Leg2{i}=STR2;
end
yline(0,'b--')
Leg2{length(NB)+1}='zero line';
%legend(Leg2,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg2,'NumColumns',1,'Box','off')
title('Glide angle plot')
ylabel('Glide Angle Beta deg')
xlabel('AOA alpha deg')

%%
figure(14)
hold on
grid on
Leg13=cell(1,1);
AOA1Start=LB-6;
AOA1end=LB-1;
AOA1=UsefullAOA(AOA1Start:AOA1end);
for i13=1:length(AOA1)
STRAOA13=append('AOA=',num2str(AOA1(i13)),'deg');
Leg13{i13}=STRAOA13;
end
L1=LatNBandAOA(AOA1Start:AOA1end,:);
NB1=NB;
SpeedNBAOA1=SpeedatNBandAOA(AOA1Start:AOA1end,:);
%[Xk,Yk]=meshgrid(NB1,AOA1);
plot3(L1,SpeedNBAOA1,NB1)
%ylim([0,1])
xlabel('Lift')
ylabel('General Speed m/s')
%zlabel('AOA alpha deg')
zlabel('NB')
title('Positive Part')
legend(Leg13,'NumColumns',1,'Location','northeast','Box','off')

%% Energy Calc
P1=1.2;
TPCHalf=Depth1./WatNBandAOA;
DPCHalf=abs(movingDistance);
DPLHalf=abs(movingDperL)./2;
AX_EHalf=P1.*TPCHalf;%Ws
PressureSurface=1e5;
PressureBottom=20e5;
BEVolume=(Total_NB.*0.1.*2./rho_SW);%m3
EtaBE=0.35;
BE_E=(PressureSurface+PressureBottom).*BEVolume./EtaBE;%Nm

AOAP=UsefullAOA(LB:end);
AOAN=UsefullAOA(1:LB-1);

AXEHalf_AOAP_NB=AX_EHalf(LB:end,:);
AXEHalf_AOAN_NB=AX_EHalf(1:LB-1,:);
DPCHalf_AOAP_NB=DPCHalf(LB:end,:);
DPCHalf_AOAN_NB=DPCHalf(1:LB-1,:);
DPLHalf_AOAP_NB=DPLHalf(LB:end,:);
DPLHalf_AOAN_NB=DPLHalf(1:LB-1,:);
TPCHalf_AOAP_NB=TPCHalf(LB:end,:);
TPCHalf_AOAN_NB=TPCHalf(1:LB-1,:);

for iAXE=1:length(AOAN)
%debugAXE1=AXEHalf_AOAN_NB(iAXE,:);
tempAXE=AXEHalf_AOAP_NB+fliplr(AXEHalf_AOAN_NB(iAXE,:));
tempDPC=DPCHalf_AOAP_NB+fliplr(DPCHalf_AOAN_NB(iAXE,:));
tempDPL=DPLHalf_AOAP_NB+fliplr(DPLHalf_AOAN_NB(iAXE,:));
tempTPC=TPCHalf_AOAP_NB+fliplr(TPCHalf_AOAN_NB(iAXE,:));
AXE_AOAP_AOAN_NB3DMatrix(:,iAXE,:)=tempAXE;
TotalE_AOAP_AOAN_NB3DMatrix(:,iAXE,:)=tempAXE+BE_E;
DPC_AOAP_AOAN_NB3DMatrix(:,iAXE,:)=tempDPC;
DPL_AOAP_AOAN_NB3DMatrix(:,iAXE,:)=tempDPL;
TPC_AOAP_AOAN_NB3DMatrix(:,iAXE,:)=tempTPC;
end


BatteryE=2000*3600;%Ws
NumOfC_AOAP_AOAN_NB3DMatrix=BatteryE./TotalE_AOAP_AOAN_NB3DMatrix;
TotalD_AOAP_AOAN_NB3DMatrix=DPC_AOAP_AOAN_NB3DMatrix.*NumOfC_AOAP_AOAN_NB3DMatrix./1000;
[XAOAN,YAOAP,ZNB]=meshgrid(AOAN,AOAP,NB);


TotalTime_AOAP_AOAN_NB3DMatrix=TPC_AOAP_AOAN_NB3DMatrix.*NumOfC_AOAP_AOAN_NB3DMatrix./3600./168;
Power1_AOAP_AOAN_NB3DMatrix=TotalE_AOAP_AOAN_NB3DMatrix./TPC_AOAP_AOAN_NB3DMatrix;


figure(15)
hold on
contourslice(XAOAN,YAOAP,ZNB,TotalD_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c15=colorbar;
c15.Label.String='Total Distant per battery charge km';
title('Total Distant 4D contourplot')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change N')
figure(16)
hold on
contourslice(XAOAN,YAOAP,ZNB,AXE_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c16=colorbar;
c16.Label.String='AUX Energy Ws';
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change N')
title('AUX Energy 4D contourplot')

figure(18)
hold on
contourslice(XAOAN,YAOAP,ZNB,TotalTime_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c18=colorbar;
c18.Label.String='Total Time per battery charge week';
title('Total Time 4D contourplot')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change N')
figure(19)
hold on
contourslice(XAOAN,YAOAP,ZNB,Power1_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c19=colorbar;
c19.Label.String='Power W';
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change N')
title('Power 4D contourplot')
figure(20)
hold on
contourslice(XAOAN,YAOAP,ZNB,NumOfC_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c20=colorbar;
c20.Label.String='Number Of Cycles per battery charge';
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change N')
title('Num of cycles 4D contourplot')
%%
figure(21)
hold on
contourf(AOAN,AOAP,TotalD_AOAP_AOAN_NB3DMatrix(:,:,6))
c21=colorbar;
c21.Label.String='Total Distant per battery charge km';
title('Total Distant 4D contourplot section on NB=8N')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
figure(22)
hold on
contourf(AOAN,AOAP,TotalTime_AOAP_AOAN_NB3DMatrix(:,:,6))
c22=colorbar;
c22.Label.String='Total Time per battery charge week';
title('Total Time 4D contourplot section on NB=8N')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
figure(23)
hold on
contourf(AOAN,AOAP,NumOfC_AOAP_AOAN_NB3DMatrix(:,:,5))
c23=colorbar;
c23.Label.String='Number of Cycles per battery charge';
title('NumOfCycles 4D contourplot section on NB=10N')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
figure(24)
hold on
contourf(AOAN,AOAP,Power1_AOAP_AOAN_NB3DMatrix(:,:,5))
c24=colorbar;
c24.Label.String='Average Power of glider W';
title('Power 4D contourplot section on NB=10N')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
figure(25)
hold on
contourslice(XAOAN,YAOAP,ZNB,DPC_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c25=colorbar;
c25.Label.String='Distance per Cycles m';
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change N')
title('DPC 4D contourplot')
figure(26)
hold on
contourf(AOAN,AOAP,DPC_AOAP_AOAN_NB3DMatrix(:,:,end))
c26=colorbar;
c26.Label.String='Distance per Cycles m';
title('DPC 4D contourplot section on NB=30N')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
figure(27)
hold on
contourf(AOAN,AOAP,DPC_AOAP_AOAN_NB3DMatrix(:,:,5))
c27=colorbar;
c27.Label.String='Distance per Cycles m';
title('DPC 4D contourplot section on NB=10N')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
figure(28)
hold on
contourf(AOAN,AOAP,DPL_AOAP_AOAN_NB3DMatrix(:,:,1))
c28=colorbar;
c28.Label.String='Distance per Cycles every Liter change m/L';
title('DPL 4D contourplot section on NB=6N')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
figure(29)
hold on
contourf(AOAN,AOAP,DPL_AOAP_AOAN_NB3DMatrix(:,:,5))
c29=colorbar;
c29.Label.String='Distance per Cycles every Liter change m/L';
title('DPL 4D contourplot section on NB=10N')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
figure(30)
hold on
contourf(AOAN,AOAP,TotalTime_AOAP_AOAN_NB3DMatrix(:,:,1))
c30=colorbar;
c30.Label.String='Total Time per battery charge week';
title('Total Time 4D contourplot section on NB=6N')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
% %%
% figure(17);
% cmap = parula(256);
% SliceD17=[0 0 1];
% %SF17=[300 130 250];%Easy in looking
% SF17=[1 1 1];%More accurate
% 
% sliceViewer(TotalD_AOAP_AOAN_NB3DMatrix,"Colormap",cmap,"Parent",figure(17),"SliceDirection",SliceD17,'ScaleFactors',SF17)
% %Change the following according to SliceDirection
% 
% if SliceD17==[0 1 0]
% title('AOAN for rise alphaN deg')
% ylabel('Net Buoyance change N')
% elseif SliceD17==[1 0 0]
% title('Net Buoyance change N')
% ylabel('AOAP for dive alphaP deg')
% elseif SliceD17==[0 0 1]
% title('AOAN for rise alphaN deg')
% ylabel('AOAP for dive alphaP deg')
% end
% set(gcf, 'Position', [50, 50, 800, 600])
% %set(gca,'position',[0.1,0.1,0.9,100] )
% c17=colorbar;
% c17.Label.String='Total Distant per battery charge km';
% %volumeViewer(TotalD_AOAP_AOAN_NB3DMatrix)


end