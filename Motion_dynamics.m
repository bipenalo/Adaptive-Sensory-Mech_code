% Contrast sensitivity vs exposure duration


%% input
close all
% 3.5% 77 ms
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics35_77.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics35_77_times_.mat');
Mdynamic35_77 = MDynamics35_77';
Mdynamic35_77_time = MDynamics35_77_times_';
% 100 ms
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics243_100.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics243_100_times_.mat');
Mdynamic_100 = MDynamics243_100';
Mdynamic_100_time = MDynamics243_100_times_';
% 149 ms
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics18_149_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics18_149_times_.mat');
Mdynamic_149_time = MDynamics18_149_times_';

% 1.5% 203 ms
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics15_203.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics15_203_times_.mat');
Mdynamic_203 = MDynamics15_203';
Mdynamic_203_time = MDynamics15_203_times_';

% 390 ms
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics13_390.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/Motion Dynamics/MDynamics13_390_times_.mat');
Mdynamic_390 = MDynamics13_390';
Mdynamic_390_time = MDynamics13_390_times_';
% % 1.07% 390 ms
% load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDynamic107_390.mat');
% load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDynamic107_390_times_.mat');
% Mdynamic107_390 = MDynamic107_390';
% Mdynamic107_390_time = MDynamic107_390_times_';
% 788 ms
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDynamic1_788.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDynamic1_788_times_.mat');
Mdynamic_788 = MDynamic_788';
Mdynamic_788_time = MDynamic_788_times_';
%%
close all;
th = 31.26; % Fixed threshold required to reach the decision criterion.

% 77 ms 
Mdynamic35_77(Mdynamic35_77 < 0) =0 ;
y77 = Mdynamic35_77(:, 1800);
t77 = Mdynamic35_77_time;
q77 = cumtrapz(t77, y77); % cumulative integral
figure(2);
hold on
title('3.5%');
xlabel('time arb units');
ylabel('activity arb units');
plot(t77, y77);
plot(t77, q77);
refline([0 th]);
hold off
% 100 ms
Mdynamic_100(Mdynamic_100 < 0) =0 ;
y100 = Mdynamic_100(:, 1800);
t100 = Mdynamic_100_time;
q100 = cumtrapz(t100, y100); % cumulative integral
figure(3);
hold on
title('2.43%');
xlabel('time arb units');
ylabel('activity arb units');
plot(t100, y100);
plot(t100, q100);
refline([0 th]);
hold off
% 149 ms
Mdynamic_149(Mdynamic_149 < 0) =0 ;
y149 = Mdynamic_149(:, 1800);
t149 = Mdynamic_149_time;
q149 = cumtrapz(t149, y149); % cumulative integral
figure(4);
hold on
title('1.9%');
xlabel('time arb units');
ylabel('activity arb units');
plot(t149, y149);
plot(t149, q149);
refline([0 th]);
hold off

% 203ms
Mdynamic_203(Mdynamic_203 < 0) =0 ;
y203 = Mdynamic_203(:, 1800);
t203 = Mdynamic_203_time;
q203 = cumtrapz(t203, y203); % cumulative integral
figure(6);
hold on
title('1.5%');
xlabel('time arb units');
ylabel('activity arb units');
plot(t203, y203);
plot(t203, q203);
refline([0 th]);
hold off

% 390ms
Mdynamic_390(Mdynamic_390 < 0) =0 ;
y390 = Mdynamic_390(:, 1800);
t390 = Mdynamic_390_time;
q390 = cumtrapz(t390, y390); % cumulative integral
figure(8);
hold on
title('1.4%');
xlabel('time arb units');
ylabel('activity arb units');
plot(t390, y390);
plot(t390, q390);
refline([0 th]);
hold off
% 788ms
Mdynamic_788(Mdynamic_788 < 0) =0 ;
y788 = Mdynamic_788(:, 1800);
t788 = Mdynamic_788_time;
q788 = cumtrapz(t788, y788); % cumulative integral
figure(10);
hold on
title('1%');
xlabel('time arb units');
ylabel('activity arb units');
plot(t788, y788);
plot(t788, q788);
refline([0 th]);
hold off



%% Final plot
% importing data
Burr_data1 = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Burr_data.csv');
Burr_data2 = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Burr_data2.csv');
Burr_avg = (Burr_data1(:,2) + Burr_data2(:,2))/2;
T_Burr= [Burr_data1(:,2) Burr_data2(:,2) ]';
Burr_std = std(T_Burr)/sqrt(2);
sensitivity = [28.2 38.46 52.63 66.67 72 72 72 72 72 72];
time = Burr_data1(:,1);
figure();
set(gca,'FontSize',14,'LineWidth',1.5,'XMinorTick','on','XScale','log',...
    'XTickLabel',{'100','1000', '10000'},'YMinorTick','on','YScale','log','YTickLabel',...
    {'10','100', '1000'});
hold on
shadedErrorBar(time, Burr_avg', Burr_std,'lineprops','-r','patchSaturation',0.33);
% Create loglog
mark = 'none';
% loglog(Burr_data(:,1), Burr_data(:,2),'MarkerFaceColor',[0 0 0],'MarkerSize',10,'Marker', mark,...
%     'LineWidth',3,...
%     'Color',[250 128 114]/255);
loglog(time, sensitivity,'MarkerFaceColor',[0 0 0],'MarkerSize',10,'Marker', mark,...
    'LineWidth',3,...
    'Color',[0 0 0]);

% Create ylabel
ylabel({'Contrast Sensitivity'},'FontWeight','bold','FontSize',16);
ylim([10 1000]);
% Create xlabel
xlabel({'Duration (ms)'},'FontWeight','bold','FontSize',16);
xlim([50 10000]);
% loglog(Burr_data(:,1), Burr_data(:,2), 'r');
% loglog(time, sensitivity);
% ylim([10 1000]);
% xlim([50 10000]);
