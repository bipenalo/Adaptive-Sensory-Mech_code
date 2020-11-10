% File containing the code for the Motion Discrimination task

% import data from Tadin's paper
close all;
Tadin_MDtask = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/paper2_data/motion discrimination tadin.csv');
Tadin_MDtask_std1 = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/paper2_data/motion discrimination tadin std1.csv');
Tadin_MDtask_std2 = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/paper2_data/motion discrimination tadin std2.csv');

Tadin_subject = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/paper2_data/motion discrimination tadin individual sub.csv');
Tadin_subject_std1 = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/paper2_data/motion discrimination tadin individual sub std1.csv');
Tadin_subject_std2 = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/paper2_data/motion discrimination tadin individual sub std2.csv');

% 9.09 contrast
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask9_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask9_times.mat');
Yi_T_MDtask_9 = MDtask9_';
t_T_MDtask_9 = MDtask9_times_';

% 11.13 contrast
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask11_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask11_times.mat');
Yi_T_MDtask_11 = MDtask11_';
t_T_MDtask_11 = MDtask11_times_';

% 17.65 contrast
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask15_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask15_times.mat');
Yi_T_MDtask_15 = MDtask15_';
t_T_MDtask_15 = MDtask15_times_';

% 24 contrast
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask22_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask22_times.mat');
Yi_T_MDtask_22 = MDtask22_';
t_T_MDtask_22 = MDtask22_times_';

% 46 contrast
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask46_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask46_times.mat');
Yi_T_MDtask_46 = MDtask46_';
t_T_MDtask_46 = MDtask46_times_';

% load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask92_.mat'); %MDtask92_test MDtask92_
% load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask/MDtask92_times.mat'); %MDtask92_test_times_ MDtask92_times
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask92_test.mat'); %MDtask92_test MDtask92_
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MDtask92_test_times_.mat'); %MDtask92_test_times_ MDtask92_times
% Yi_T_MDtask_92 = MDtask92_';  %MDtask92_
% t_T_MDtask_92 = MDtask92_times_';  %MDtask92_times_
Yi_T_MDtask_92 = MDtask_test';  %MDtask92_
t_T_MDtask_92 = MDtask_test_times_';

contrast = [9.09 11.13 16.65 24 46 92];
th = 3; % 3.0
close all;
y9 = Yi_T_MDtask_9(:, 1800) + 0;
y9(y9 < 0.0) = 0.0;
t9 = t_T_MDtask_9;
q9 = cumtrapz(t9, y9); % cumulative integral
[q9_th, index] = unique(q9);
th9 = interp1(q9_th, t9(index), th);
figure();
hold on
title('9.09%');
xlabel('time arb units');
ylabel('activity arb units');
plot(t9, y9);
plot(t9, q9);
hold off

y11 = Yi_T_MDtask_11(:, 1800);
y11(y11 < 0.0) = 0.0;
t11 = t_T_MDtask_11;
q11 = cumtrapz(t11, y11);
[q11_th, index] = unique(q11);
th11 = interp1(q11_th, t11(index), th);
figure();
hold on
title('11.13 %');
xlabel('time arb units');
ylabel('activity arb units');
plot(t11, y11);
plot(t11, q11);
hold off

y15 = Yi_T_MDtask_15(:, 1800);
y15(y15 < 0.0) = 0.0;
t15 = t_T_MDtask_15;
q15 = cumtrapz(t15, y15);
[q15_th, index] = unique(q15);
th15 = interp1(q15_th, t15(index), th);
figure();
hold on
title('17.65 %');
xlabel('time arb units');
ylabel('activity arb units');
plot(t15, y15);
plot(t15, q15);
hold off

y22 = Yi_T_MDtask_22(:, 1800) ;
y22(y22 < 0.0) = 0.0;
t22 = t_T_MDtask_22;
q22 = cumtrapz(t22, y22);
[q22_th, index] = unique(q22);
th22 = interp1(q22_th, t22(index), th);
figure();
hold on
title('24 %');
xlabel('time arb units');
ylabel('activity arb units');
plot(t22, y22);
plot(t22, q22);
hold off

y46 = Yi_T_MDtask_46(:, 1800) ;
y46(y46 < 0.0) = 0.0;
t46 = t_T_MDtask_46;
q46 = cumtrapz(t46, y46);
[q46_th, index] = unique(q46);
th46 = interp1(q46_th, t46(index), th);
figure();
hold on
title('46 %');
xlabel('time arb units');
ylabel('activity arb units');
plot(t46, y46);
plot(t46, q46);
hold off

y92 = Yi_T_MDtask_92(:, 1800);
y92(y92 < 0.0) = 0.0;
t92 = t_T_MDtask_92;
q92 = cumtrapz(t92, y92);
[q92_th, index] = unique(q92);
th92 = interp1(q92_th, t92(index), th);
figure();
hold on
title('92 %');
xlabel('time arb units');
ylabel('activity arb units');
plot(t92, y92);
plot(t92, q92);
hold off

figure();
plot(Yi_T_MDtask_92(150,:))
xlabel('space arb units');
ylabel('activity arb units');

scale = 23.28/0.7; % 23.28/0.7;
ths = scale * [th9, th11, th15, th22, th46, th92]; 
fig = figure();
hold on
ypos1= Tadin_MDtask_std1(:,2) - Tadin_MDtask(:,2);
yneg1= Tadin_MDtask(:,2) - Tadin_MDtask_std2(:,2);
p1 = shadedErrorBar(contrast, Tadin_MDtask(:,2), ypos1, 'lineprops','-r','patchSaturation',0.33);
p2 = plot(contrast, Tadin_MDtask(:,2), '-r', 'MarkerFaceColor','r','MarkerSize',10,'Marker','o'); 
set(gca, 'XScale','log', 'YScale','log');
p3= createfig_MDtask(contrast, ths, 'none');
legend([p2 p3],"Tadin's data", 'Model');
legend('boxon');


% Chi-square Goodness of fit test
stdv = ypos1*2;
var = stdv.^2;
chi2 = sum((ths' - Tadin_MDtask(:,2)).^2 ./ var);
 
