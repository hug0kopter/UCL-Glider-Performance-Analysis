clc
clear all
close all
%% Input Constants:(Run this section if any of the constants need updating)

%general
rho_SW=1027; %salt water density
WingArea=1.112; % Planar wing area looking from the top (reference area)

% Inputting the Net Buoyancy values for the code to interpolate for
%If any of the constants need to be changed, just rerun this section. Its
%not neccesary to rerun the fitting part when the constants here are
%changed.

Total_NB = 14;% Total volume of the VBD piston
NB =[3,4,5,6,7,8,9,10,11];% Net buoyancy values for which the lines will be plotted
Depth1 = 200;% Set target depth

%Data Process
uselessAOARange=[-1,0];% Removing the AoA between -2.8 and -1.2 
% because the lift generated is too low causing our glider to fall vertically

%Energy calculation
P1 = 1.2; %ambient power draw from electronics (not including actuators) [W]
PandRPC1=71928; %Pitch and Roll Energy per dive and rise [Ws]
PressureSurface = 1e5; %Atmospheric pressure at water surface
PressureBottom = 21e5; %Water pressure at maximum depth
EtaBE = 0.35; %buoyancy engine efficiency - our estimate
BatteryE = 1512*3600 ; %total battery capacity in joules - our estimate:
% 1512 [Wh] * 3600 [s/h]

%4D slicer
SliceControler = 6; %Index of the Net Buoyancy we want to use when descending
% Net Buoyancy when ascending is calculated by subtracting this from the
% value of Total_NB, declared above

%% Read Data
folderPath='Put_The_Data_Here/CFD/';
FileList=dir(folderPath);
FileList=FileList(4:end);
SpeedList=[0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];%list of all velocities for which data was collected

NumberOfSpeed=length(SpeedList);
TemplateAOA=transpose(linspace(-7,7,15));%AoA range -7~+7

TableSpaceRef=zeros(length(TemplateAOA),NumberOfSpeed);
% Assemble the file name for read
FileName1=append(folderPath,FileList(1).name);
FileName2=append(folderPath,FileList(2).name);
FileName3=append(folderPath,FileList(3).name);
%Read the prepared tables
AoAdebuglist=importdata(FileName1);
CD_matrix=importdata(FileName2);
CL_matrix=importdata(FileName3);

%Calculating temporary matricies to be multiplied by velocity at a later
%stage
Temporary_Drag_Matrix = 1/2.*CD_matrix.*rho_SW*WingArea; %(1/2)*CD*density*Wing Area
Temporary_Lift_Matrix = 1/2.*CL_matrix.*rho_SW*WingArea; %(1/2)*CL*density*Wing Area
drag_matrix=[]; %matrix of drag where the dimensions are velocity and AoA
lift_matrix=[]; %matrix of lift where the dimensions are velocity and AoA

%multiplying by velocity to calculate actual lift and drag
for i=1:length(SpeedList)
drag_matrix=[drag_matrix,Temporary_Drag_Matrix(:,i).*SpeedList(i).^2];
lift_matrix=[lift_matrix,Temporary_Lift_Matrix(:,i).*SpeedList(i).^2];
end

% Calculating net buoyancy for terminal velocity at each of the velocities
% and AoAs in the dataset
net_buoyancy_matrix=(drag_matrix.^2+lift_matrix.^2).^(1/2);

%% Reverse relationship calc

%Set the structure to save relationships
NB_Speed_Relationship = repmat(struct(),[length(TemplateAOA),1]);

% Processing Data
UsefullAOARange=[];
for i=1:length(TemplateAOA)
    NB_line = net_buoyancy_matrix(i,:); % read each line of the net buoyancy matrix
    %Saving the valid aoa range
    UsefullAOARange=[UsefullAOARange,i];
    %Set up the curve fitting environment
    y = transpose(SpeedList);
    y2 = transpose(lift_matrix(i,:));
    x = transpose(net_buoyancy_matrix(i,:));
    if i==7 %Manual overwrite correction for fit this case
        xu=[x(1);x(2);x(3);x(5)];
        y2u=[y2(1);y2(2);y2(3);y2(5)];
    end
    %Exponential fit for a relationship between Speed and net
    %buoyancy, for the current angle of attack at current
    %iteration
    f = fittype('a1*exp(b1*x) + a2*exp(b2*x) ');

    %Defining the starting points of the fit for our dataset,
    %this was determined using datafitting toolbox
    start_point = [-0.4, 0.5, 0, 0]; 

    fit_options = fitoptions('Method','NonlinearLeastSquares','StartPoint', start_point,'Algorithm','Levenberg-Marquardt');
    [fit_result, gof] = fit(x, y, f, fit_options); %Fitting the curve
    disp(fit_result); %displaying the equation of the currently fitted curve
    plot(fit_result,x,y); %showing the fitting plot for visual inspection

    %Save the fit function
    NB_Speed_Relationship(i).FitFunction = fit_result;
    disp('Please check the fit result') % printing a line after each curve fit for debugging

    %Linear fit to find the relationship between lift and
    %velocity at the angle of attack of the current loop
    %iteration
    f2 = fittype('poly1');
    if i==7
        [fit_result2, gof2] = fit(xu, y2u, f2); 
    else
        [fit_result2, gof2] = fit(x, y2, f2); 
    end
    disp(fit_result2); %displaying the equation of the second fit
    plot(fit_result2,x,y2); %showing the fitting plot for visual inspection

    %Save the fit function
    NB_Speed_Relationship(i).FitFunction2 = fit_result2;
    disp('Please check the fit result2') % printing a line after each curve fit for debugging
end
%Save everything into RAM file to save time on doing fitting
close all 
save('CFD_Final_Code_RAM')
clear

%%
% loading the fitted curves into the code so that it can be run without 
% reruning the curve fit every time
clear
load('CFD_Final_Code_RAM.mat')
if true %run this section to reload the data and run the entire plotting section

    if true %run this section only to recalculate all data neccesary for the plots without plotting
%% Process data
close all

UsefullAOARange = UsefullAOARange(1:end);% for debugging
UsefullAOA = TemplateAOA;%(UsefullAOARange); % Range of AoA we interested in (-7 to 7) and excluding the empty results
UsefullNBSR = NB_Speed_Relationship;%(UsefullAOARange); % Taking only the fitted curves which correspond with AoA's of interest

%Declare empty matricies for data to be appended later
velocity_matrix = zeros(length(UsefullAOA),length(NB)); % matrix of velocities with dimensions of net buoyancy and Aoa
L_NBandAOA = zeros(length(UsefullAOA),length(NB)); % lift matrix
B_NBandAOA = zeros(length(UsefullAOA),length(NB)); % gliding angle (Beta) matrix
theta_NBandAOA = zeros(length(UsefullAOA),length(NB)); % pitch angle matrix
U_NBandAOA = zeros(length(UsefullAOA),length(NB)); % horizontal velocity matrix
W_NBandAOA = zeros(length(UsefullAOA),length(NB)); % vertical velocity matrixok
horizontal_distance_travelled = zeros(length(UsefullAOA),length(NB)); %distance per dive/rise matrix
movingDperL = zeros(length(UsefullAOA),length(NB)); % distance per litre for each AoA and NB
NB_in_litres = zeros(length(UsefullAOA),length(NB)); % NB in litres for energy calculations

% Calculating all performance characteristics using the trigonometric
% relationships from the free body diagram
for i=1:length(NB) % for each net buoyancy value
    for j=1:length(UsefullAOA)% for each angle of attack of interest
        velocity_matrix(j,i) = feval(UsefullNBSR(j).FitFunction, NB(i)); % calculating velocity for each NB from fitted curves
        L_NBandAOA(j,i) = feval(UsefullNBSR(j).FitFunction2, NB(i)); % calculating lift for each NB from fitted curves 
        B_NBandAOA(j,i) = acosd(L_NBandAOA(j,i)./NB(i)); % calculating lift for each NB from fitted curves = acos(Lift/Net Buoyancy) 
        U_NBandAOA(j,i) = cosd(B_NBandAOA(j,i))*velocity_matrix(j,i); % calculating horizontal velocity for each NB from fitted curves = cos(Beta)* Velocity) 
        W_NBandAOA(j,i) = sind(B_NBandAOA(j,i))*velocity_matrix(j,i); % calculating vertical velocity for each NB from fitted curves = sin(Beta)* Velocity)
        theta_NBandAOA(j,i) = B_NBandAOA(j,i)-UsefullAOA(j); % calculating theta (pitch angle) for each NB from fitted curves = glide angle - AoA
        horizontal_distance_travelled(j,i) = Depth1/tand(B_NBandAOA(j,i)); % calculating distance travelled horizontally = depth/tan(Beta)
        movingDperL(j,i) = horizontal_distance_travelled(j,i)/(NB(i)*0.1*2/rho_SW*1000);%  % calculating distance travelled horizontally per litre
        NB_in_litres(j,i) = (NB(i)*0.1*2/rho_SW*1000);% converting NB to litres
    end
end

%Find the seperation between Pos and Neg
%UsefullAOA = round(UsefullAOA,3);
LB = find(UsefullAOA==uselessAOARange(1)); % Lower boundary of the AoA range (Separation point right before the lowest AoA when diving)
UB = find(UsefullAOA==uselessAOARange(2)); % Upper boundary of the AoA range

% Splitting the dataset, removing the AoAs when the glider drops vertically
% and then combining them again:
UsefullAOA=[UsefullAOA(1:LB-1);UsefullAOA(UB:end)];

velocity_matrix = [velocity_matrix(1:LB-1,:);velocity_matrix(UB:end,:)];

L_NBandAOA = [L_NBandAOA(1:LB-1,:);L_NBandAOA(UB:end,:)];

B_NBandAOA = [B_NBandAOA(1:LB-1,:);B_NBandAOA(UB:end,:)];

theta_NBandAOA = [theta_NBandAOA(1:LB-1,:);theta_NBandAOA(UB:end,:)];

U_NBandAOA = [U_NBandAOA(1:LB-1,:);U_NBandAOA(UB:end,:)];

W_NBandAOA = [W_NBandAOA(1:LB-1,:);W_NBandAOA(UB:end,:)];

horizontal_distance_travelled = [horizontal_distance_travelled(1:LB-1,:);horizontal_distance_travelled(UB:end,:)];

movingDperL = [movingDperL(1:LB-1,:);movingDperL(UB:end,:)];

NB_in_litres = [NB_in_litres(1:LB-1,:);NB_in_litres(UB:end,:)];

checkAOA0 = find(UsefullAOA==0);
%% Graph stage
% Specifying the domain and range for the plots (x and y values)
generalYlim=[0,0.2];
generalXlim=[0,1];
generalXlim2=[-1,0];

F1P2start = 0; %To adjust the starting point for the 5 polar curves on figure 1 and 8,
% alter this number, the curves plotted will be in 0.2 degrees increment
% if u want to have a look at a certain range of angles just keep
% increasing this variable until you see it.

    end % here to run the process stage
%% The following section and all its plots are for positive AoA's which means the glider is diving
%
%% Figure1: glider polar in dive
figure(1) 
hold on
grid on
% xlim([-1,1])
xlim(generalXlim)
ylim(generalYlim)
Leg=cell(1,1);%Set for legend use
PNB=[];

for i=1:length(NB)
plot(U_NBandAOA(LB:end,i),W_NBandAOA(LB:end,i))
STRNB2=append('NetB=',num2str(NB(i)),'N');
Leg{i}=STRNB2;
end

istart=LB+F1P2start;
for i=istart:istart+5% Here can chane for the amount of AOA lines needed
plot([0,U_NBandAOA(i,:)],[0,W_NBandAOA(i,:)],'--')
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
plot(UsefullAOA(LB:end),theta_NBandAOA(LB:end,i))
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
contourf(NB,UsefullAOA(LB:end),horizontal_distance_travelled(LB:end,:))
colorbar
title('The distance per dive m')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(4)
hold on
contourf(NB,UsefullAOA(LB:end),movingDperL(LB:end,:))
colorbar
title('The distance per dive per VBD volume change [m/L]')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(5)
hold on
grid on

AOA1=UsefullAOA(LB:end);
MDL1=movingDperL(LB:end,:)./1000;
SpeedNBAOA1=velocity_matrix(LB:end,:);
plot3(MDL1,SpeedNBAOA1,AOA1)
ylim([0,1])
xlabel('The distance per dive per VBD volume change [km/L]')
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
plot(UsefullAOA(LB:end),B_NBandAOA(LB:end,i))
STR2=append('NetB=',num2str(NB(i)),'N');
Leg6{i}=STR2;
end
yline(0,'b--')
Leg6{length(NB)+1}='zero line';
%legend(Leg6,'NumColumns',1,'Location','bestoutside','Box','off')
legend(Leg6,'NumColumns',1,'Box','off')
title('Glide angle vs. AoA plot for dive')
ylabel('Glide Angle Beta deg')
xlabel('AOA alpha deg')
%%
figure(7)
hold on
grid on
Leg7=cell(1,1);
AOA1Start=LB;
AOA1end=LB+7;
AOA1=UsefullAOA(AOA1Start:AOA1end);
for i7=1:length(AOA1)
STRAOA7=append('AOA=',num2str(AOA1(i7)),'deg');
Leg7{i7}=STRAOA7;
end
L1=L_NBandAOA(AOA1Start:AOA1end,:);
NB1=NB;
SpeedNBAOA1=velocity_matrix(AOA1Start:AOA1end,:);
%[Xk,Yk]=meshgrid(NB1,AOA1);
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
plot(-U_NBandAOA(1:LB-1,i),W_NBandAOA(1:LB-1,i))
STRNB2=append('NetB=',num2str(NB(i)),'N');
Leg8{i}=STRNB2;
end

istart=1+F1P2start;
for i=LB-istart-5:LB-istart
plot([0,-U_NBandAOA(i,:)],[0,W_NBandAOA(i,:)],'--')
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
plot(UsefullAOA(1:LB-1),180-theta_NBandAOA(1:LB-1,i))
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
contourf(NB,UsefullAOA(1:LB-1),-horizontal_distance_travelled(1:LB-1,:))
colorbar
title('The distance per rise m')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(11)
hold on
contourf(NB,UsefullAOA(1:LB-1),-movingDperL(1:LB-1,:))
colorbar
title('The distance per rise per VBD volume change [m/L]')
ylabel('AOA alpha deg')
xlabel('Net Buoyance change N')
%%
figure(12)
hold on
grid on
AOA1=UsefullAOA(1:LB-1);
MDL1=-movingDperL(1:LB-1,:)./1000;
SpeedNBAOA1=velocity_matrix(1:LB-1,:);
plot3(MDL1,SpeedNBAOA1,AOA1)
ylim([0,1])
xlabel('The distance per rise per VBD volume change [km/L]')
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
    plot(UsefullAOA(1:LB-1),180-B_NBandAOA(1:LB-1,i))
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
AOA1Start=LB-6;
AOA1end=LB-1;
AOA1=UsefullAOA(AOA1Start:AOA1end);
for i14=1:length(AOA1)
STRAOA14=append('AOA=',num2str(AOA1(i14)),'deg');
Leg14{i14}=STRAOA14;
end
L1=L_NBandAOA(AOA1Start:AOA1end,:);
NB1=NB;
SpeedNBAOA1=velocity_matrix(AOA1Start:AOA1end,:);
%[Xk,Yk]=meshgrid(NB1,AOA1);
plot3(L1,SpeedNBAOA1,NB1)
%ylim([0,1])
xlabel('Lift')
ylabel('General Speed m/s')
zlabel('NB')
title('Positive Part')
legend(Leg14,'NumColumns',1,'Location','northeast','Box','off')
%% Energy Calculations

TPC_matrix = Depth1./W_NBandAOA; %Time per dive/rise
DPC_matric = abs(horizontal_distance_travelled); %Distance per dive/rise
AX_E = P1.*TPC_matrix;% Power consumed due to ambient current draw (no actuators) in Joules

%Flow Work calculations
BEVolume = (Total_NB.*0.1.*2./rho_SW);% buoyancy engine volume in m3
BE_E = (PressureSurface+PressureBottom).*BEVolume./EtaBE;% Energy consumed by VBD per cycle (dive+rise) in Joules

%Energy calculation
%Setting up the 4D dataset
AOAP = UsefullAOA(LB:end); %positive angles of attack (meaning in dive)
AOAN = UsefullAOA(1:LB-1);%negative angles of attack (meaning in rise)
AXE_AOAP_NB = AX_E(LB:end,:); % ambient energy drawn when going down
AXE_AOAN_NB = AX_E(1:LB-1,:); % ambient energy drawn when going up
DPC_AOAP_NB = DPC_matric(LB:end,:);% distance per cycle going down (diving)
DPC_AOAN_NB = DPC_matric(1:LB-1,:);% distance per cycle going up (rising)
TPC_AOAP_NB = TPC_matrix(LB:end,:);% time per dive
TPC_AOAN_NB = TPC_matrix(1:LB-1,:);% time per rise

% This for loop adds each of the matricies with the same matrix but flipped
% around NB axis, this means that we can calculate total energy consumption
% per cycle, making sure that all combination of NB in dive and NB in rise
% add to the total VBD volume.
for iAXE=1:length(AOAN)
    %debugAXE1=AXEHalf_AOAN_NB(iAXE,:);
    tempAXE = AXE_AOAP_NB + fliplr(AXE_AOAN_NB(iAXE,:)); % summing 
    tempDPC = DPC_AOAP_NB + fliplr(DPC_AOAN_NB(iAXE,:));
    tempDPL = tempDPC./Total_NB;
    tempTPC = TPC_AOAP_NB + fliplr(TPC_AOAN_NB(iAXE,:));
    AXE_AOAP_AOAN_NB3DMatrix(:,iAXE,:) = tempAXE;
    TotalE_AOAP_AOAN_NB3DMatrix(:,iAXE,:) = tempAXE + BE_E + PandRPC1;
    DPC_AOAP_AOAN_NB3DMatrix(:,iAXE,:) = tempDPC;
    DPL_AOAP_AOAN_NB3DMatrix(:,iAXE,:) = tempDPL;
    TPC_AOAP_AOAN_NB3DMatrix(:,iAXE,:) = tempTPC;
end
%calculating other further terms

NumOfC_AOAP_AOAN_NB3DMatrix = BatteryE./TotalE_AOAP_AOAN_NB3DMatrix; %number of cycles 3D matrix
TotalD_AOAP_AOAN_NB3DMatrix = DPC_AOAP_AOAN_NB3DMatrix.*NumOfC_AOAP_AOAN_NB3DMatrix./1000;% total distance 3D matrix

TotalTime_AOAP_AOAN_NB3DMatrix = TPC_AOAP_AOAN_NB3DMatrix.*NumOfC_AOAP_AOAN_NB3DMatrix./3600./168;% total mission time 3D matrix
Power1_AOAP_AOAN_NB3DMatrix = TotalE_AOAP_AOAN_NB3DMatrix./TPC_AOAP_AOAN_NB3DMatrix; % Average power draw matrix for each configuration [W]
[XAOAN,YAOAP,ZNB] = meshgrid(AOAN,AOAP,NB);
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