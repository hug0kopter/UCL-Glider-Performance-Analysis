clc
clear all
close all

rho_SW=1025;% sea water density
WingArea=0.848;%change base on your wing
folderPath='Put_The_Data_Here/UOPT1/';%Put your folder path here just the folder!
FileList=dir(folderPath);%checking what is in the folder
FileList=FileList(4:end);%first 3 is . .. and Matlab drive so dont want them
NumberOfFiles=length(FileList);%creat a index
TemplateAOA=transpose(linspace(-10,10,101));%set up a AOA range to set the matrix
length(TemplateAOA)% i was try to creat a index but then i just copy every time so
SpeedList=[];%set a empty space for speed
BlankSpeedList=zeros(1,NumberOfFiles);% setup a empty sample for later compare
TableSpaceRef=zeros(length(TemplateAOA),NumberOfFiles);%these are all declear for empty space
AoAdebuglist=zeros(length(TemplateAOA),1);
CDatSandAOA=zeros(length(TemplateAOA),1);
CLatSandAOA=zeros(length(TemplateAOA),1);
%then is to take data from csv
for i=1:NumberOfFiles
   
    FileName=append(folderPath,FileList(i).name);%generate file name
    F1=importdata(FileName);%data import
    Temp1=split(FileList(i).name,'-');%this is to deconstruct the speed number from its name
    Temp2=split(Temp1(2),' ');
    Temp3=split(Temp2(1),'_');
    speed=str2num(cell2mat(append(Temp3(1),'.',Temp3(2))));%construct speed data
    
    TempAOA=F1.data(:,1);% temp saving the AOA CD CL data
    TempCD=F1.data(:,6);
    TempCL=F1.data(:,3);
    for j=1:101 % find the starting point of data as it is not full AOA range
        R1=round(TemplateAOA(j),3);%No idea there is some minor noise in the data so just round it
        R2=round(TempAOA(1),3);
        if R1==R2
            j1=j;% set up a pointer
            break
        end    
    end
  % Save every data in list for every file.
    SpeedList=[SpeedList,speed];
    AoAdebuglist=[AoAdebuglist,[AoAdebuglist(1:j1-1,i);TempAOA;AoAdebuglist(j1+length(TempCD):length(TemplateAOA),i)]];
    CDatSandAOA=[CDatSandAOA,[CDatSandAOA(1:j1-1,i);TempCD;CDatSandAOA(j1+length(TempCD):length(TemplateAOA),i)]];
    CLatSandAOA=[CLatSandAOA,[CLatSandAOA(1:j1-1,i);TempCL;CLatSandAOA(j1+length(TempCD):length(TemplateAOA),i)]];
    
end
AoAdebuglist(:,1)=[];% clean the first col that is full of zero
CDatSandAOA(:,1)=[];% clean the first col that is full of zero
CLatSandAOA(:,1)=[];% clean the first col that is full of zero
TempDatSandAOA=1/2.*CDatSandAOA.*rho_SW*WingArea;
TempLatSandAOA=1/2.*CLatSandAOA.*rho_SW*WingArea;
DatSandAOA=[];
LatSandAOA=[];
for i=1:length(SpeedList)
DatSandAOA=[DatSandAOA,TempDatSandAOA(:,i).*SpeedList(i).^2];
LatSandAOA=[LatSandAOA,TempLatSandAOA(:,i).*SpeedList(i).^2];
end
% clculate L and D
NBatSandAOA=(DatSandAOA.^2+LatSandAOA.^2).^(1/2);
% clculate NetBouyancy
%% Here is begine the Fitting process to reverse the relationship between Speed and NB
NB_Speed_Relationship=repmat(struct(),[length(TemplateAOA),1]);% creat a space to save fitted functions
UsefullAOARange=[];% declear for AOA the have atleast 1/2 of the speed have data.
for i=1:length(TemplateAOA)
NBForAOAi=NBatSandAOA(i,:); % read NB
if NBForAOAi ~= BlankSpeedList % check if empty
kN0=find(NBForAOAi); % find none zero term
if length(kN0)>=1/2*length(SpeedList) %check if > 1/2 of the length that means have enough data to fit
UsefullAOARange=[UsefullAOARange,i]; % set the index of AOA data
y=transpose(SpeedList); % fiting toolbox need to trans all the matrix to vertical
y2=transpose(LatSandAOA(i,:)); % fiting toolbox need to trans all the matrix to vertical
x=transpose(NBatSandAOA(i,:)); % fiting toolbox need to trans all the matrix to vertical
f = fittype('a1*exp(b1*x) + a2*exp(b2*x) '); % setting fittype to be 2 order exp
start_point = [0.895, 0.1672, -0.1, -1]; % here depends on your matlab version and your pc may need to be change and play around.
fit_options = fitoptions('Method','NonlinearLeastSquares','StartPoint', start_point,'Algorithm','Levenberg-Marquardt'); %set the method and algorithm used
[fit_result, gof] = fit(x, y, f, fit_options); %just~ do it
disp(fit_result);% show the parameter
plot(fit_result,x,y);% plot out for check if the curve is nice, if not nice please abort and change the starting point
NB_Speed_Relationship(i).FitFunction=fit_result; % save the fit function
disp('请检查拟合结果，Please check the fit result') % 请检查拟合结果

f2 = fittype('poly1');% doing linear fit to check the L and NB relationship
[fit_result2, gof2] = fit(x, y2, f2); %just~ do it
disp(fit_result2);% show the parameter
plot(fit_result2,x,y2);% plot out for check if the curve is nice, if not nice just pray it be nice.
NB_Speed_Relationship(i).FitFunction2=fit_result2;% save the fit function
disp('请检查拟合结果2，Please check the fit result2')% 请检查拟合结果

end
end
end
close all %we dont need the graph after checking
save("DAW_UOPT1_RAM.mat")
clear

%% now is to plotting, you can run this section seperate from the main after you done the fit
clear
load('DAW_UOPT1_RAM.mat')
close all
NB=[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,30];% here is the NB value you want to have a look can change and play
UsefullAOA=TemplateAOA(UsefullAOARange);%use index to get data
UsefullNBSR=NB_Speed_Relationship(UsefullAOARange);%use index to get function
SpeedatNBandAOA=zeros(length(UsefullAOA),length(NB));%declear
LatNBandAOA=zeros(length(UsefullAOA),length(NB));%declear
BetaatNBandAOA=zeros(length(UsefullAOA),length(NB));%declear
thetaatNBandAOA=zeros(length(UsefullAOA),length(NB));%declear
UatNBandAOA=zeros(length(UsefullAOA),length(NB));%declear
WatNBandAOA=zeros(length(UsefullAOA),length(NB));%declear
movingDistance=zeros(length(UsefullAOA),length(NB));%declear
movingDperL=zeros(length(UsefullAOA),length(NB));%declear
% save data from raw data
for i=1:length(NB)
for j=1:length(UsefullAOA)
    SpeedatNBandAOA(j,i)=feval(UsefullNBSR(j).FitFunction, NB(i));%use fit function to clc
    LatNBandAOA(j,i)=feval(UsefullNBSR(j).FitFunction2, NB(i));%use fit function to clc
    BetaatNBandAOA(j,i)=acosd(LatNBandAOA(j,i)./NB(i));%clc
    UatNBandAOA(j,i)=cosd(BetaatNBandAOA(j,i))*SpeedatNBandAOA(j,i);%clc
    WatNBandAOA(j,i)=sind(BetaatNBandAOA(j,i))*SpeedatNBandAOA(j,i);%clc
    thetaatNBandAOA(j,i)=BetaatNBandAOA(j,i)-UsefullAOA(j);%clc
    movingDistance(j,i)=200/tand(BetaatNBandAOA(j,i));%clc
    movingDperL(j,i)=movingDistance(j,i)/(NB(i)*0.1*2/rho_SW*1000);%in m/L

end
end

%put the useless AOA range here
uselessAOARange=[0,0.8]; % that is you change the AOA range that you do not want, normally the middle part around zero
LB=find(UsefullAOA==uselessAOARange(1));% find the index of two boundary of the useless region
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
%cut the data
checkAOA0=find(UsefullAOA==0);% check if you cut out zero incase you want zero to be in though is no reason for that
%general change
generalYlim=[0,0.2];
generalXlim=[0,1];
F1P2start=0;% the starting point of your AOA line in F1, 0 is from the start, 1 is the 1st one and so on.
if checkAOA0~=[]% if you somehow want to have zero use AOA0 as index to seperate positive and negative


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
%%
figure(7)
hold on
grid on
Leg5=cell(1,1);
AOA1Start=AOA0+11;
AOA1end=AOA0+26;
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
TPC=200./WatNBandAOA(AOA0:end,:).*2;
DPC=movingDistance(AOA0:end,:).*2;
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


figure(8)
hold on
contourf(NB,UsefullAOA(AOA0:end),TotalD)
colorbar
title('The Total Distant per battery charge km')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
figure(9)
hold on
contourf(NB,UsefullAOA(AOA0:end),AXE)
colorbar
title('The AX Energy Ws')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
else % here is the normal way that use the LB index as the seperation
figure(1) 
hold on
grid on
%change if need negtive part
%xlim([-1,1])
xlim(generalXlim)
ylim(generalYlim)
Leg=cell(1,1);% this is for legend in massive way
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

%change if need negtive part but idont think you need anymore, just use MS
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
grid on
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
grid on
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