clc
clear all
close all

rho_SW=1025;
WingArea=0.848;%change base on your wing
folderPath='Put_The_Data_Here/OPT1/';
FileList=dir(folderPath);
FileList=FileList(4:end);
NumberOfFiles=length(FileList);
TemplateAOA=transpose(linspace(-10,10,101));
length(TemplateAOA)
SpeedList=[];
BlankSpeedList=zeros(1,NumberOfFiles);
TableSpaceRef=zeros(length(TemplateAOA),NumberOfFiles);
AoAdebuglist=zeros(length(TemplateAOA),1);
CDatSandAOA=zeros(length(TemplateAOA),1);
CLatSandAOA=zeros(length(TemplateAOA),1);

for i=1:NumberOfFiles
   
    FileName=append(folderPath,FileList(i).name);
    F1=importdata(FileName);
    Temp1=split(FileList(i).name,'-');
    Temp2=split(Temp1(2),' ');
    Temp3=split(Temp2(1),'_');
    speed=str2num(cell2mat(append(Temp3(1),'.',Temp3(2))));
    
    TempAOA=F1.data(:,1);
    TempCD=F1.data(:,6);
    TempCL=F1.data(:,3);
    for j=1:101
        R1=round(TemplateAOA(j),3);
        R2=round(TempAOA(1),3);
        if R1==R2
            j1=j;
            break
        end    
    end
  
    SpeedList=[SpeedList,speed];
    AoAdebuglist=[AoAdebuglist,[AoAdebuglist(1:j1-1,i);TempAOA;AoAdebuglist(j1+length(TempCD):length(TemplateAOA),i)]];
    CDatSandAOA=[CDatSandAOA,[CDatSandAOA(1:j1-1,i);TempCD;CDatSandAOA(j1+length(TempCD):length(TemplateAOA),i)]];
    CLatSandAOA=[CLatSandAOA,[CLatSandAOA(1:j1-1,i);TempCL;CLatSandAOA(j1+length(TempCD):length(TemplateAOA),i)]];
    
end
AoAdebuglist(:,1)=[];
CDatSandAOA(:,1)=[];
CLatSandAOA(:,1)=[];
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
if NBForAOAi ~= BlankSpeedList
kN0=find(NBForAOAi);
if length(kN0)>=1/2*length(SpeedList)
UsefullAOARange=[UsefullAOARange,i];
y=transpose(SpeedList);
y2=transpose(LatSandAOA(i,:));
x=transpose(NBatSandAOA(i,:));
f = fittype('a1*exp(b1*x) + a2*exp(b2*x) ');
start_point = [0.895, 0.1672, -0.1, -1]; 
fit_options = fitoptions('Method','NonlinearLeastSquares','StartPoint', start_point,'Algorithm','Levenberg-Marquardt'); 
[fit_result, gof] = fit(x, y, f, fit_options); 
disp(fit_result);
plot(fit_result,x,y);
NB_Speed_Relationship(i).FitFunction=fit_result;
disp('请检查拟合结果，Please check the fit result')

f2 = fittype('poly1');
[fit_result2, gof2] = fit(x, y2, f2); 
disp(fit_result2);
plot(fit_result2,x,y2);
NB_Speed_Relationship(i).FitFunction2=fit_result2;
disp('请检查拟合结果2，Please check the fit result2')

end
end
end
close all %we dont need the graph after checking
save("DAW_OPT1_RAM.mat")
clear

%% now is to plotting, you can run this section seperate from the main after you done the fit
clear
load('DAW_OPT1_RAM.mat')
close all
NB=[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,30];% here is the NB value you want to have a look can change and play
UsefullAOA=TemplateAOA(UsefullAOARange);
UsefullNBSR=NB_Speed_Relationship(UsefullAOARange);
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
    movingDistance(j,i)=200/tand(BetaatNBandAOA(j,i));
    movingDperL(j,i)=movingDistance(j,i)/(NB(i)*0.1*2/rho_SW*1000);%in L

end
end

%put the useless AOA range here
uselessAOARange=[0,0.8];
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
F1P2start=0;
if checkAOA0~=[]


figure(1) 

hold on
grid on
%change if need negtive part
%xlim([-1,1])
xlim(generalXlim)
ylim(generalYlim)
Leg=cell(1,1);
PNB=[];
for i=1:length(UsefullAOA)
    if round(UsefullAOA(i),3)==0
        AOA0=i;
    end
end
for i=1:length(NB)
plot(UatNBandAOA(AOA0:end,i),WatNBandAOA(AOA0:end,i))
STRNB2=append('NetB=',num2str(NB(i)),'N');
Leg{i}=STRNB2;
end

%change if need negtive part
istart=AOA0+F1P2start;
for i=istart:istart+6
plot([0,UatNBandAOA(i,:)],[0,WatNBandAOA(i,:)],'--')
STRAOA2=append('AOA=',num2str(UsefullAOA(i)),'deg');
Leg{i-istart+1+length(NB)}=STRAOA2;
end
%legend(Leg,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg,'NumColumns',1,'Location','northwest', 'Box','off')
title('Glide polar')
ylabel('Vertical Velocity w m/s')
xlabel('Horizontal Velocity u m/s')

figure(2)
hold on
grid on
Leg2=cell(1,1);
for i=1:length(NB)
plot(UsefullAOA(AOA0:end),thetaatNBandAOA(AOA0:end,i))
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
contourf(NB,UsefullAOA(AOA0:end),movingDistance(AOA0:end,:))
colorbar
title('The distance per dive or rise m')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(4)
hold on
Leg4=cell(1,1);
contourf(NB,UsefullAOA(AOA0:end),movingDperL(AOA0:end,:))
colorbar
title('The distance per dive or rise every Liter change m/L')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(5)
hold on
grid on
Leg5=cell(1,1);
AOA1=UsefullAOA(AOA0:end);
MDL1=movingDperL(AOA0:end,:)./1000;
SpeedNBAOA1=SpeedatNBandAOA(AOA0:end,:);
plot3(MDL1,SpeedNBAOA1,AOA1)
ylim([0,1])
xlabel('The distance per dive or rise every Liter change km')
ylabel('General Speed m/s')
zlabel('AOA alpha deg')
Leg2=Leg2(1:length(Leg2)-1);
%legend(Leg2,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg2,'NumColumns',1,'Box','off')


figure(6)
hold on
grid on
Leg2=cell(1,1);
for i=1:length(NB)
plot(UsefullAOA(AOA0:end),BetaatNBandAOA(AOA0:end,i))
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
else
figure(1) 
hold on
%change if need negtive part
%xlim([-1,1])
xlim(generalXlim)
ylim(generalYlim)
Leg=cell(1,1);
PNB=[];
% for i=1:length(UsefullAOA)
%     if round(UsefullAOA(i),3)==0
%         AOA0=i;
%     end
% end
for i=1:length(NB)
plot(UatNBandAOA(LB:end,i),WatNBandAOA(LB:end,i))
STRNB2=append('NetB=',num2str(NB(i)),'N');
Leg{i}=STRNB2;
end

%change if need negtive part
istart=LB+F1P2start;
for i=istart:istart+6
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
title('The distance per dive or rise m')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(4)
hold on
Leg4=cell(1,1);
contourf(NB,UsefullAOA(LB:end),movingDperL(LB:end,:))
colorbar
title('The distance per dive or rise every Liter change m/L')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(5)
hold on
Leg5=cell(1,1);
AOA1=UsefullAOA(LB:end);
MDL1=movingDperL(LB:end,:)./1000;
SpeedNBAOA1=SpeedatNBandAOA(LB:end,:);
plot3(MDL1,SpeedNBAOA1,AOA1)
ylim([0,1])
xlabel('The distance per dive or rise every Liter change km')
ylabel('General Speed m/s')
zlabel('AOA alpha deg')
Leg2=Leg2(1:length(Leg2)-1);
%legend(Leg2,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg2,'NumColumns',1,'Box','off')

figure(6)
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
figure(7)
hold on
grid on
Leg5=cell(1,1);
AOA1Start=LB+11;
AOA1end=LB+26;
AOA1=UsefullAOA(AOA1Start:AOA1end);
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
Leg2=Leg2(1:length(Leg2));
%legend(Leg2,'NumColumns',1,'Location','bestoutside','Box','off')
%legend(Leg2,'NumColumns',1,'Box','off')
%% 
P1=1.2;
TPC=200./WatNBandAOA(LB:end,:).*2;
DPC=movingDistance(LB:end,:).*2;
AXE=P1.*TPC;%Ws
PS=1e5;
PB=20e5;
BEVolume=(NB.*0.1.*2./rho_SW);%m3
EtaBE=0.5;
BEE=(PS+PB).*BEVolume./EtaBE;%Nm

TotalE=AXE;
for i=1:length(AXE)
TotalE(i,:)=AXE(i,:)+BEE;
end

BatteryE=2000*3600;
NumOfC=BatteryE./TotalE;
TotalD=DPC.*NumOfC./1000;
TotalTime=TPC.*NumOfC./3600./168;
Power1=TotalE./TPC;

figure(8)
hold on
contourf(NB,UsefullAOA(LB:end),TotalD)
colorbar
title('The Total Distant per battery charge km')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
figure(9)
hold on
contourf(NB,UsefullAOA(LB:end),AXE)
colorbar
title('The AX Energy Ws')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
figure(10)
hold on
contourf(NB,UsefullAOA(LB:end),TotalTime)
colorbar
title('Total Time per battery charge weeks')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
figure(11)
hold on
contourf(NB,UsefullAOA(LB:end),Power1)
colorbar
title('Power of Glider W')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
figure(12)
hold on
contourf(NB,UsefullAOA(LB:end),NumOfC)
colorbar
title('Number of cycles per battery charge')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
end