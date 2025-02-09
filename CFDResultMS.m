clc
clear
close all
folderPath='Put_The_Data_Here/';
FileName1=append(folderPath,'CFDSimulationData.CSV');
DataField1=importdata(FileName1);
Velocity1=DataField1.data(:,1);
AOA1=DataField1.data(:,2);
NB1=DataField1.data(:,12);
Time1=DataField1.data(:,13);
Distance1=DataField1.data(:,14);
% Table1=readtable(FileName1);
AOA2=[AOA1(1:6);AOA1(8:15)];
NB2=[NB1(1:6);NB1(8:15)];
Time2=[Time1(1:6);Time1(8:15)];
Distance2=[Distance1(1:6);Distance1(8:15)];
AOA3=[AOA1(16:21);AOA1(23:30)];
NB3=[NB1(16:21);NB1(23:30)];
Time3=[Time1(16:21);Time1(23:30)];
Distance3=[Distance1(16:21);Distance1(23:30)];
figure(1)
scatter3(AOA1,Time1,NB1,5,Velocity1)
title('No idea')
ylabel('NB N')
xlabel('AOA deg')
zlabel('Time per div or rise (s)')
colorbar
figure(2)
scatter3(AOA1,Distance1,NB1,5,Velocity1)
title('No idea')
ylabel('NB N')
xlabel('AOA deg')
zlabel('Distance per div or rise (m)')
figure(3)
hold on
scatter(AOA2,NB2,30,Distance2,'x')
scatter(AOA3,NB3,5,Distance3)

ylabel('NB N')
xlabel('AOA deg')
legend('0.3','0.4')
c3=colorbar;
c3.Label.String='Distance per div or rise (m)';