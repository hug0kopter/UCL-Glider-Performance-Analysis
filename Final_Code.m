clc
clear all
close all
%% Input Constents:(Can run this section seperate to update some constents)

%general
rho_SW=1027;
WingArea=1.112;%change base on your wing

%Change NB
%After this line the constents are changeable without rerun the fitting
Total_NB=14;% Set total capacity of BE
NB=[3,4,5,6,7,8,9,10,11];% here is the NB value you want to have a look can change and play
Depth1=200;% Set target depth

%Data Process
uselessAOARange=[-2.8,-1.2];% Some aoa around the transit between Pos and Neg will be useless

%Energy calculation
P1=1.2;%the aux power in w
PressureSurface=1e5;%the surf pressure in pa
PressureBottom=21e5;%the bottom pressure in pa
EtaBE=0.35;%BE efficiency - estimate
BatteryE=2000*3600;%total battery cap in Ws - estimates

%4D slicer
SliceControler=6;%want slice base on which NB index

%% Read Data
folderPath='Put_The_Data_Here/MECHSPACE_XFLR5/';
FileList=dir(folderPath);
FileList=FileList(4:end);
NumberOfFiles=length(FileList);
%Declear Spaces
TemplateAOA=transpose(linspace(-10,10,101));
SpeedList=[];
BlankSpeedList=zeros(1,NumberOfFiles);
TableSpaceRef=zeros(length(TemplateAOA),NumberOfFiles);
AoAdebuglist=zeros(length(TemplateAOA),1);
CDatSandAOA=zeros(length(TemplateAOA),1);
CLatSandAOA=zeros(length(TemplateAOA),1);
%Read Each File
for i=1:NumberOfFiles
    FileName=append(folderPath,FileList(i).name);
    F1=importdata(FileName);
    %Reconstruct the Speed from Filename
    Temp1=split(FileList(i).name,'-');
    Temp2=split(Temp1(2),' ');
    Temp3=split(Temp2(1),'_');
    speed=str2num(cell2mat(append(Temp3(1),'.',Temp3(2))));
    
    TempAOA=F1.data(:,1);
    TempCD=F1.data(:,6);
    TempCL=F1.data(:,3);
    %Find the starting point of the AOA data
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
%Remove the first term that is 0
AoAdebuglist(:,1)=[];
CDatSandAOA(:,1)=[];
CLatSandAOA(:,1)=[];

%Calculate L & D
TempDatSandAOA=1/2.*CDatSandAOA.*rho_SW*WingArea;
TempLatSandAOA=1/2.*CLatSandAOA.*rho_SW*WingArea;
DatSandAOA=[];
LatSandAOA=[];
for i=1:length(SpeedList)
DatSandAOA=[DatSandAOA,TempDatSandAOA(:,i).*SpeedList(i).^2];
LatSandAOA=[LatSandAOA,TempLatSandAOA(:,i).*SpeedList(i).^2];
end
%Find the NetBuoyancy at steady state
NBatSandAOA=(DatSandAOA.^2+LatSandAOA.^2).^(1/2);
%% Reverse relationship calc
%Set the struct to save relationships
NB_Speed_Relationship=repmat(struct(),[length(TemplateAOA),1]);
%Find the valid aoa range
UsefullAOARange=[];
for i=1:length(TemplateAOA)
    NBForAOAi=NBatSandAOA(i,:);
    %Non blank check
    if NBForAOAi ~= BlankSpeedList
        %Find Zero term index
        kN0=find(NBForAOAi);
        %If more than half is zero
        if length(kN0)>=1/2*length(SpeedList)
            %Extra limitation for quantitiy normally is +7~-7
            if i>16 && i<86
                %Save the valid aoa range
                UsefullAOARange=[UsefullAOARange,i];
                %Set up the curve fitting environment
                y=transpose(SpeedList);
                y2=transpose(LatSandAOA(i,:));
                x=transpose(NBatSandAOA(i,:));
                %EXP fit for Speed and NB relationship at different AOA
                f = fittype('a1*exp(b1*x) + a2*exp(b2*x) ');

                start_point = [0.895, 0.1672, -0.1, -1]; %Here the staring point need to use the toolbox to find for different dataset
                
                fit_options = fitoptions('Method','NonlinearLeastSquares','StartPoint', start_point,'Algorithm','Levenberg-Marquardt');
                [fit_result, gof] = fit(x, y, f, fit_options);
                disp(fit_result);
                plot(fit_result,x,y);
                %Save the fit function
                NB_Speed_Relationship(i).FitFunction=fit_result;
                disp('请检查拟合结果，Please check the fit result')
                %Linear fit to get Speed and Lift relationship at different AOA
                f2 = fittype('poly1');
                [fit_result2, gof2] = fit(x, y2, f2);
                disp(fit_result2);
                plot(fit_result2,x,y2);
                %Save the fit function
                NB_Speed_Relationship(i).FitFunction2=fit_result2;
                disp('请检查拟合结果2，Please check the fit result2')
            end
        end
    end
end
%Save everything into RAM file to save time on doing fitting
close all 
save('Final_Code_RAM')
clear

%% Run this part seperatly
clear
load('Final_Code_RAM')
if true %force to run the entire plot section
    if true %force to run before the process part
%% Process data
close all

UsefullAOARange=UsefullAOARange(1:end);%debug if imag
UsefullAOA=TemplateAOA(UsefullAOARange);
UsefullNBSR=NB_Speed_Relationship(UsefullAOARange);
%Declear the space
SpeedatNBandAOA=zeros(length(UsefullAOA),length(NB));
LatNBandAOA=zeros(length(UsefullAOA),length(NB));
BetaatNBandAOA=zeros(length(UsefullAOA),length(NB));
thetaatNBandAOA=zeros(length(UsefullAOA),length(NB));
UatNBandAOA=zeros(length(UsefullAOA),length(NB));
WatNBandAOA=zeros(length(UsefullAOA),length(NB));
movingDistance=zeros(length(UsefullAOA),length(NB));
movingDperL=zeros(length(UsefullAOA),length(NB));
NBinLforDPLCalc=zeros(length(UsefullAOA),length(NB));
%calculate for each term 
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
    NBinLforDPLCalc(j,i)=(NB(i)*0.1*2/rho_SW*1000);%in L
end
end

%Find the seperation between Pos and Neg
UsefullAOA=round(UsefullAOA,3);
LB=find(UsefullAOA==uselessAOARange(1));
UB=find(UsefullAOA==uselessAOARange(2));

%Split and reform the data
UsefullAOA=[UsefullAOA(1:LB-1);UsefullAOA(UB:end)];
SpeedatNBandAOA=[SpeedatNBandAOA(1:LB-1,:);SpeedatNBandAOA(UB:end,:)];
LatNBandAOA=[LatNBandAOA(1:LB-1,:);LatNBandAOA(UB:end,:)];
BetaatNBandAOA=[BetaatNBandAOA(1:LB-1,:);BetaatNBandAOA(UB:end,:)];
thetaatNBandAOA=[thetaatNBandAOA(1:LB-1,:);thetaatNBandAOA(UB:end,:)];
UatNBandAOA=[UatNBandAOA(1:LB-1,:);UatNBandAOA(UB:end,:)];
WatNBandAOA=[WatNBandAOA(1:LB-1,:);WatNBandAOA(UB:end,:)];
movingDistance=[movingDistance(1:LB-1,:);movingDistance(UB:end,:)];
movingDperL=[movingDperL(1:LB-1,:);movingDperL(UB:end,:)];
NBinLforDPLCalc=[NBinLforDPLCalc(1:LB-1,:);NBinLforDPLCalc(UB:end,:)];
checkAOA0=find(UsefullAOA==0);
%% Graph stage
%general change
generalYlim=[0,0.2];
generalXlim=[0,1];
generalXlim2=[-1,0];
F1P2start=0;%Figure1 series part 2 that draw the aoa range starting point adjustment
    end % here to run the process stage
%% Positive part which is for dive
%%
figure(1) 
hold on
grid on
% change if need negtive part
% xlim([-1,1])
xlim(generalXlim)
ylim(generalYlim)
Leg=cell(1,1);%Set for legend use
PNB=[];

for i=1:length(NB)
plot(UatNBandAOA(LB:end,i),WatNBandAOA(LB:end,i))
STRNB2=append('NetB=',num2str(NB(i)),'N');
Leg{i}=STRNB2;
end

% change if need negtive part
istart=LB+F1P2start;
for i=istart:istart+5% Here can chane for the amount of AOA lines needed
plot([0,UatNBandAOA(i,:)],[0,WatNBandAOA(i,:)],'--')
STRAOA2=append('AOA=',num2str(UsefullAOA(i)),'deg');
Leg{i-istart+1+length(NB)}=STRAOA2;% Continue the Legend after the part1
end
%legend(Leg,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg,'NumColumns',1,'Location','northwest','Box','off')
title('Glide polar for dive')
ylabel('Vertical Velocity w m/s')
xlabel('Horizontal Velocity u m/s')
%%
figure(2)
hold on
grid on
Leg2=cell(1,1);%Set for legend use
for i=1:length(NB)
plot(UsefullAOA(LB:end),thetaatNBandAOA(LB:end,i))
STR2=append('NetB=',num2str(NB(i)),'N');
Leg2{i}=STR2;
end
yline(0,'b--')% Mark the zero pitch
Leg2{length(NB)+1}='zero line';
%legend(Leg2,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg2,'NumColumns',1,'Box','off')
title('Pitch angle plot for dive')
ylabel('Pitch Angle Theta deg')
xlabel('AOA alpha deg')
%%
figure(3)
hold on
contourf(NB,UsefullAOA(LB:end),movingDistance(LB:end,:))
colorbar
title('The distance per dive m')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(4)
hold on
contourf(NB,UsefullAOA(LB:end),movingDperL(LB:end,:))
colorbar
title('The distance per dive every Liter change m/L')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(5)
hold on
grid on

AOA1=UsefullAOA(LB:end);
MDL1=movingDperL(LB:end,:)./1000;
SpeedNBAOA1=SpeedatNBandAOA(LB:end,:);
plot3(MDL1,SpeedNBAOA1,AOA1)
ylim([0,1])
xlabel('The distance per dive every Liter change km')
ylabel('General Speed m/s')
zlabel('AOA alpha deg')
Leg5=Leg2(1:length(Leg2)-1);%need to run graph2 first
%legend(Leg5,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg5,'NumColumns',1,'Box','off')
%%
figure(6)
hold on
Leg6=cell(1,1);
for i=1:length(NB)
plot(UsefullAOA(LB:end),BetaatNBandAOA(LB:end,i))
STR2=append('NetB=',num2str(NB(i)),'N');
Leg6{i}=STR2;
end
yline(0,'b--')
Leg6{length(NB)+1}='zero line';
%legend(Leg6,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg6,'NumColumns',1,'Box','off')
title('Glide angle plot for dive')
ylabel('Glide Angle Beta deg')
xlabel('AOA alpha deg')
%%
figure(7)
hold on
grid on
Leg7=cell(1,1);
AOA1Start=LB+10;
AOA1end=LB+17;
AOA1=UsefullAOA(AOA1Start:AOA1end);
for i7=1:length(AOA1)
STRAOA7=append('AOA=',num2str(AOA1(i7)),'deg');
Leg7{i7}=STRAOA7;
end
L1=LatNBandAOA(AOA1Start:AOA1end,:);
NB1=NB;
SpeedNBAOA1=SpeedatNBandAOA(AOA1Start:AOA1end,:);
[Xk,Yk]=meshgrid(NB1,AOA1);
plot3(L1,SpeedNBAOA1,NB1)
xlabel('Lift')
ylabel('General Speed m/s')
zlabel('NB')
title('Positive Part')
legend(Leg7,'NumColumns',1,'Location','northwest','Box','off')

%% Negtive part which is rising
%%
figure(8) 
hold on
grid on
% xlim([-1,1])
xlim(generalXlim)
ylim(generalYlim)
Leg8=cell(1,1);
PNB=[];

for i=1:length(NB)
plot(-UatNBandAOA(1:LB-1,i),WatNBandAOA(1:LB-1,i))
STRNB2=append('NetB=',num2str(NB(i)),'N');
Leg8{i}=STRNB2;
end

istart=1+F1P2start;
for i=LB-istart-5:LB-istart
plot([0,-UatNBandAOA(i,:)],[0,WatNBandAOA(i,:)],'--')
STRAOA2=append('AOA=',num2str(UsefullAOA(i)),'deg');
Leg8{i-(LB-istart-5)+1+length(NB)}=STRAOA2;
end
%legend(Leg,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg8,'NumColumns',1,'Location','northwest','Box','off')
title('Glide polar for rise')
ylabel('Vertical Velocity w m/s')
xlabel('Horizontal Velocity u m/s')
%%
figure(9)
hold on
grid on
Leg9=cell(1,1);
for i=1:length(NB)
plot(UsefullAOA(1:LB-1),180-thetaatNBandAOA(1:LB-1,i))
STR2=append('NetB=',num2str(NB(i)),'N');
Leg9{i}=STR2;
end
yline(0,'b--')
Leg9{length(NB)+1}='zero line';
%legend(Leg9,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg9,'NumColumns',1,'Box','off')
title('Pitch angle plot for rise')
ylabel('Pitch Angle Theta deg')
xlabel('AOA alpha deg')
%% 
figure(10)
hold on
contourf(NB,UsefullAOA(1:LB-1),-movingDistance(1:LB-1,:))
colorbar
title('The distance per rise m')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(11)
hold on
contourf(NB,UsefullAOA(1:LB-1),-movingDperL(1:LB-1,:))
colorbar
title('The distance per rise every Liter change m/L')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(12)
hold on
grid on
AOA1=UsefullAOA(1:LB-1);
MDL1=-movingDperL(1:LB-1,:)./1000;
SpeedNBAOA1=SpeedatNBandAOA(1:LB-1,:);
plot3(MDL1,SpeedNBAOA1,AOA1)
ylim([0,1])
xlabel('The distance per rise every Liter change km')
ylabel('General Speed m/s')
zlabel('AOA alpha deg')
Leg12=Leg9(1:length(Leg9)-1);
%legend(Leg12,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg12,'NumColumns',1,'Box','off')
%%
figure(13)
hold on
Leg13=cell(1,1);
for i=1:length(NB)
    plot(UsefullAOA(1:LB-1),180-BetaatNBandAOA(1:LB-1,i))
    STR2=append('NetB=',num2str(NB(i)),'N');
    Leg13{i}=STR2;
end
yline(0,'b--')
Leg13{length(NB)+1}='zero line';
%legend(Leg13,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg13,'NumColumns',1,'Box','off')
title('Glide angle plot')
ylabel('Glide Angle Beta deg')
xlabel('AOA alpha deg')
%%
figure(14)
hold on
grid on
Leg14=cell(1,1);
AOA1Start=LB-17;
AOA1end=LB-10;
AOA1=UsefullAOA(AOA1Start:AOA1end);
for i14=1:length(AOA1)
STRAOA14=append('AOA=',num2str(AOA1(i14)),'deg');
Leg14{i14}=STRAOA14;
end
L1=LatNBandAOA(AOA1Start:AOA1end,:);
NB1=NB;
SpeedNBAOA1=SpeedatNBandAOA(AOA1Start:AOA1end,:);
[Xk,Yk]=meshgrid(NB1,AOA1);
plot3(L1,SpeedNBAOA1,NB1)
%ylim([0,1])
xlabel('Lift')
ylabel('General Speed m/s')
zlabel('NB')
title('Positive Part')
legend(Leg14,'NumColumns',1,'Location','northeast','Box','off')
%% Energy Calc
TPCHalf=Depth1./WatNBandAOA;
DPCHalf=abs(movingDistance);
% DPLHalf=abs(movingDperL)./2;%This need to change as the NB now is coupled
NBinLHalf=NBinLforDPLCalc;%in L
AX_EHalf=P1.*TPCHalf;%Ws
%general calculating the terms
BEVolume=(Total_NB.*0.1.*2./rho_SW);%m3

BE_E=(PressureSurface+PressureBottom).*BEVolume./EtaBE;%Nm
%Energy calculation
%Setting up the 4D dataset
AOAP=UsefullAOA(LB:end);
AOAN=UsefullAOA(1:LB-1);

AXEHalf_AOAP_NB=AX_EHalf(LB:end,:);
AXEHalf_AOAN_NB=AX_EHalf(1:LB-1,:);
DPCHalf_AOAP_NB=DPCHalf(LB:end,:);
DPCHalf_AOAN_NB=DPCHalf(1:LB-1,:);
% DPLHalf_AOAP_NB=DPLHalf(LB:end,:);
% DPLHalf_AOAN_NB=DPLHalf(1:LB-1,:);
TPCHalf_AOAP_NB=TPCHalf(LB:end,:);
TPCHalf_AOAN_NB=TPCHalf(1:LB-1,:);

for iAXE=1:length(AOAN)
%debugAXE1=AXEHalf_AOAN_NB(iAXE,:);
tempAXE=AXEHalf_AOAP_NB+fliplr(AXEHalf_AOAN_NB(iAXE,:));
tempDPC=DPCHalf_AOAP_NB+fliplr(DPCHalf_AOAN_NB(iAXE,:));
% tempDPL=DPLHalf_AOAP_NB+fliplr(DPLHalf_AOAN_NB(iAXE,:));
tempDPL=tempDPC./Total_NB;
tempTPC=TPCHalf_AOAP_NB+fliplr(TPCHalf_AOAN_NB(iAXE,:));
AXE_AOAP_AOAN_NB3DMatrix(:,iAXE,:)=tempAXE;
TotalE_AOAP_AOAN_NB3DMatrix(:,iAXE,:)=tempAXE+BE_E;
DPC_AOAP_AOAN_NB3DMatrix(:,iAXE,:)=tempDPC;
DPL_AOAP_AOAN_NB3DMatrix(:,iAXE,:)=tempDPL;
TPC_AOAP_AOAN_NB3DMatrix(:,iAXE,:)=tempTPC;
end
%calculating other further terms

NumOfC_AOAP_AOAN_NB3DMatrix=BatteryE./TotalE_AOAP_AOAN_NB3DMatrix;
TotalD_AOAP_AOAN_NB3DMatrix=DPC_AOAP_AOAN_NB3DMatrix.*NumOfC_AOAP_AOAN_NB3DMatrix./1000;
[XAOAN,YAOAP,ZNB]=meshgrid(AOAN,AOAP,NB);

TotalTime_AOAP_AOAN_NB3DMatrix=TPC_AOAP_AOAN_NB3DMatrix.*NumOfC_AOAP_AOAN_NB3DMatrix./3600./168;
Power1_AOAP_AOAN_NB3DMatrix=TotalE_AOAP_AOAN_NB3DMatrix./TPC_AOAP_AOAN_NB3DMatrix;
%%Graphic part2
%%
%THIS IS A 4D GRAPH
figure(15)
hold on
contourslice(XAOAN,YAOAP,ZNB,TotalD_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c15=colorbar;
c15.Label.String='Total Distant per battery charge km';
title('Total Distant 4D contourplot')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change for down N')
%%
%THIS IS A 4D GRAPH
figure(16)
hold on
contourslice(XAOAN,YAOAP,ZNB,AXE_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)
c16=colorbar;
c16.Label.String='AUX Energy Ws';
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change for down N')
title('AUX Energy 4D contourplot')
%%
%Figure 17 is a CT scanning software that is not functioning well in this code,
%if need it to determin the operation region, release it in the end
%Hide it if not use, it will cost alot of time
%%
%THIS IS A 4D GRAPH
figure(18)
hold on
contourslice(XAOAN,YAOAP,ZNB,TotalTime_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c18=colorbar;
c18.Label.String='Total Time per battery charge week';
title('Total Time 4D contourplot')
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change for down N')
%%
%THIS IS A 4D GRAPH
figure(19)
hold on
contourslice(XAOAN,YAOAP,ZNB,Power1_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c19=colorbar;
c19.Label.String='Power W';
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change for down N')
title('Power 4D contourplot')
%%
%THIS IS A 4D GRAPH
figure(20)
hold on
contourslice(XAOAN,YAOAP,ZNB,NumOfC_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c20=colorbar;
c20.Label.String='Number Of Cycles per battery charge';
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change for down N')
title('Num of cycles 4D contourplot')
%%
%THIS IS A 4D GRAPH
figure(21)
hold on
contourslice(XAOAN,YAOAP,ZNB,DPC_AOAP_AOAN_NB3DMatrix,AOAN,AOAP,NB,5)

c21=colorbar;
c21.Label.String='Distance per Cycles m';
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
zlabel('Net Buoyance change N')
title('DPC 4D contourplot')
%%
%This is a key section view graph that need to notice
figure(22)
hold on
contourf(AOAN,AOAP,TotalD_AOAP_AOAN_NB3DMatrix(:,:,SliceControler))
c22=colorbar;
c22.Label.String='Total Distant per battery charge km';
title(append('Total Distant 4D contourplot section on NB=',num2str(NB(SliceControler)),'N for dive'))
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
%%
%This is a key section view graph that need to notice
figure(23)
hold on
contourf(AOAN,AOAP,TotalTime_AOAP_AOAN_NB3DMatrix(:,:,SliceControler))
c23=colorbar;
c23.Label.String='Total Time per battery charge week';
title(append('Total Time 4D contourplot section on NB=',num2str(NB(SliceControler)),'N for dive'))
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
%%
figure(24)
hold on
contourf(AOAN,AOAP,NumOfC_AOAP_AOAN_NB3DMatrix(:,:,SliceControler))
c24=colorbar;
c24.Label.String='Number of Cycles per battery charge';
title(append('NumOfCycles 4D contourplot section on NB=',num2str(NB(SliceControler)),'N for dive'))
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
%%
figure(25)
hold on
contourf(AOAN,AOAP,Power1_AOAP_AOAN_NB3DMatrix(:,:,SliceControler))
c25=colorbar;
c25.Label.String='Average Power of glider W';
title(append('Power 4D contourplot section on NB=',num2str(NB(SliceControler)),'N for dive'))
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')

%%
figure(26)
hold on
contourf(AOAN,AOAP,DPC_AOAP_AOAN_NB3DMatrix(:,:,SliceControler))
c26=colorbar;
c26.Label.String='Distance per Cycles m';
title(append('DPC 4D contourplot section on NB=',num2str(NB(SliceControler)),'N for dive'))
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
%%
figure(27)
hold on
contourf(AOAN,AOAP,DPL_AOAP_AOAN_NB3DMatrix(:,:,SliceControler))
c27=colorbar;
c27.Label.String='Distance per Cycles every Liter change m/L';
title(append('DPL 4D contourplot section on NB=',num2str(NB(SliceControler)),'N for dive'))
xlabel('AOAN for rise alphaN deg')
ylabel('AOAP for dive alphaP deg')
%%
% figure(17);
% cmap = parula(256);
% SliceD17=[0 0 1];%change to choose the direction of slice
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