% Segmentation task
close all;

% import the data from Tadin's paper
Tadin = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/paper2_data/Motion Segregation Tadin.csv');
tadin_std1 = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/paper2_data/Motion Segregation Tadin std1.csv');
tadin_std2 = csvread('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Adative_CS/Data/paper2_data/Motion Segregation Tadin std2.csv');

%%
% 9.09% contrast
% channel 1
close all;

load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/MStask9_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/MStask9_times.mat');
t_T_ch1_9 = MStask9_times_';
Yi_T_ch1_9 = MStask9_';
space = 1:3600;
figure(5);
hold on
h = surface(space, t_T_ch1_9, Yi_T_ch1_9);
view(0,0);
colormap summer;
title('Adaptive C-S output at 9.09%');
xlabel('Space');
ylabel('Time');
zlabel('Activity');
set(h,'LineStyle','none') % removes the grid lines 
% channel 2
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask9_ch2.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask9_times_ch2.mat');
t_T_ch2_9 = MStask9_times_ch2';
Yi_T_ch2_9 = MStask9_ch2';
space = 1:3600;
%figure(6);
h = surface(space, t_T_ch2_9, Yi_T_ch2_9);
view(0,0);
%colormap cool;
title('Adaptive C-S output at 9.09%');
xlabel('Space');
ylabel('Time');
zlabel('Activity');
set(h,'LineStyle','none') % removes the grid lines 
hold off

t9 = size(Yi_T_ch1_9, 1);
ch1_9x = zeros(uint8(t9/1), 1);
sep = 16.5;
for i = 1:uint8(t9/1)
    
    ch1_9x(i,1) = EstimateCell(Yi_T_ch1_9, i, sep, [1 1800], 1);
    
end
t9_ch2 = size(Yi_T_ch2_9, 1);
ch2_9x = zeros(uint8(t9_ch2/1), 1);
for i = 1:uint8(t9_ch2/1)
    
    ch2_9x(i,1) = EstimateCell(Yi_T_ch2_9, i, sep, [1300 1800], 2);
    
end



figure()
hold on
plot(t_T_ch2_9, ch2_9x)
plot(t_T_ch1_9, ch1_9x)
hold off
d9 = ch1_9x(100) - ch2_9x(100);
SS9 = zeros(size(ch1_9x, 1), 1);

SS9(ch1_9x > 0) = d9;
th = 10;
y9 = SS9;
q9 = cumtrapz(t_T_ch1_9, y9); % cumulative integral
[q9_th, index] = unique(q9);
th9 = interp1(q9_th, t_T_ch1_9(index), th);
figure();
hold on 
title('9%');
plot(t_T_ch1_9, y9)
plot(t_T_ch1_9, q9)
hold off
%%
% 11.13% contrast
% channel 1
close all;
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask11_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask11_times.mat');
t_T_ch1_11 = MStask11_times_';
Yi_T_ch1_11 = MStask11_';

% channel 2
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask11_ch2.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask11_times_ch2.mat');
t_T_ch2_11 = MStask11_times_ch2';
Yi_T_ch2_11 = MStask11_ch2';
t11 = size(Yi_T_ch1_11, 1);
ch1_11x = zeros(t11, 1);

% 3D plot
space = 1:3600;
figure(51);
hold on
h = surface(space, t_T_ch1_11, Yi_T_ch1_11);
view(0,0);
colormap cool;
title('Adaptive C-S output at 11%');
xlabel('Space');
ylabel('Time');
zlabel('Activity');
set(h,'LineStyle','none') % removes the grid lines 

h = surface(space, t_T_ch2_11, Yi_T_ch2_11);
view(0,0);
%colormap cool;
title('Adaptive C-S output at 11%');
xlabel('Space');
ylabel('Time');
zlabel('Activity');
set(h,'LineStyle','none') % removes the grid lines 
hold off
sep = 16.5;
for i = 1:t11
    
    ch1_11x(i,1) = EstimateCell(Yi_T_ch1_11, i, sep, [1 1800], 1);
    
end
t11_ch2 = size(Yi_T_ch2_11, 1);
ch2_11x = zeros(t11_ch2, 1);
for i = 1:t11_ch2
    
    ch2_11x(i,1) = EstimateCell(Yi_T_ch2_11, i, sep, [1300 1800], 2);
    
end

%dif = uint8(t11_ch2/2) - uint8(t11/2) + 1;
% SS11 = ch1_11x(2:end) - ch2_11x;% Only here the # of elements in ch1 is larger than in ch2. 
figure()
hold on
plot(t_T_ch2_11, ch2_11x)
plot(t_T_ch1_11, ch1_11x)
hold off
d11 = ch1_11x(100) - ch2_11x(100);

SS11 = zeros(size(ch1_11x, 1), 1);
SS11(ch1_11x > 0) = d11;
th = 10;
y11 = SS11;
q11 = cumtrapz(t_T_ch1_11, y11); % cumulative integral
[q11_th, index] = unique(q11);
th11 = interp1(q11_th, t_T_ch1_11(index), th);
figure();
hold on 
title('11.13%');
plot(t_T_ch1_11, y11)
plot(t_T_ch1_11, q11)
hold off
%%
% 15.65% contrast
% channel 1
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/MStask15_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/MStask15_times.mat');
t_T_ch1_15 = MStask15_times_';
Yi_T_ch1_15 = MStask15_';

% channel 2
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask15_ch2.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask15_times_ch2.mat');
t_T_ch2_15 = MStask15_times_ch2';
Yi_T_ch2_15 = MStask15_ch2';

t15_ch1 = size(Yi_T_ch1_15, 1);
ch1_15x = zeros(t15_ch1, 1);
sep = 16.5;

for i = 1:t15_ch1
    
    ch1_15x(i,1) = EstimateCell(Yi_T_ch1_15, i, sep, [1 1800], 1);
    
end
t15_ch2 = size(Yi_T_ch2_15, 1);
ch2_15x = zeros(t15_ch2, 1);
for i = 1:t15_ch2
    
    ch2_15x(i,1) = EstimateCell(Yi_T_ch2_15, i, sep, [1300 1800], 2);
    
end

% dif = uint8(t15_ch2/2) - uint8(t15_ch1/2) + 1; 
% SS15 = ch1_15x - ch2_15x(dif:end);% Spatial 
% plot separation as a function of time
figure()
hold on
plot(t_T_ch2_15, ch2_15x)
plot(t_T_ch1_15, ch1_15x)
hold off
d15 = ch1_15x(100) - ch2_15x(100);

SS15 = zeros(size(ch1_15x, 1), 1);
SS15(ch1_15x > 0) = d15;
th = 10;
y15 = SS15;
q15 = cumtrapz(t_T_ch1_15, y15); % cumulative integral
[q15_th, index] = unique(q15);
th15 = interp1(q15_th, t_T_ch1_15(index), th);
figure();
hold on 
title('15.65%');
plot(t_T_ch1_15, y15)
plot(t_T_ch1_15, q15)
hold off
%%
% 22% contrast
% channel 1
close all;
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask22_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask22_times.mat');
t_T_ch1_22 = MStask22_times_';
Yi_T_ch1_22 = MStask22_';

% channel 2
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask22_ch2.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask22_times_ch2.mat');
t_T_ch2_22 = MStask22_times_ch2';
Yi_T_ch2_22 = MStask22_ch2';

figure()
hold on
h = surface(space, t_T_ch1_22, Yi_T_ch1_22);
view(0,0);
colormap hot;
title('Adaptive C-S output at 22%');
xlabel('Space');
ylabel('Time');
zlabel('Activity');
set(h,'LineStyle','none') % removes the grid lines 
h = surface(space, t_T_ch2_22, Yi_T_ch2_22);
view(0,0);
colormap hot;
title('Adaptive C-S output at 22%');
xlabel('Space');
ylabel('Time');
zlabel('Activity');
set(h,'LineStyle','none') % removes the grid lines 
hold off

sep = 14;
t22_ch1 = size(Yi_T_ch1_22, 1);
ch1_22x = zeros(t22_ch1, 1);
for i = 1:t22_ch1
    
    ch1_22x(i,1) = EstimateCell(Yi_T_ch1_22, i, sep, [1 1800], 1);
    
end
t22_ch2 = size(Yi_T_ch2_22, 1);
ch2_22x = zeros(t22_ch2, 1);
for i = 1:t22_ch2
    
    ch2_22x(i,1) = EstimateCell(Yi_T_ch2_22, i, sep, [1300 1800], 2);
    
end

% plot separation as a function of time
figure()
hold on
plot(t_T_ch2_22, ch2_22x)
plot(t_T_ch1_22, ch1_22x)
hold off
d22 = ch1_22x(100) - ch2_22x(100);

SS22 = zeros(size(ch1_22x, 1), 1);
SS22(ch1_22x > 0) = d22;
th = 10;
y22 = SS22;
q22 = cumtrapz(t_T_ch1_22, y22); % cumulative integral
[q22_th, index] = unique(q22);
th22 = interp1(q22_th, t_T_ch1_22(index), th);
figure();
hold on 
title('24%');
plot(t_T_ch1_22, y22)
plot(t_T_ch1_22, q22)
hold off

%%
% 46% contrast
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/MStask92_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/MStask92_times.mat');
t_T_ch1_46 = MStask92_times_';
Yi_T_ch1_46 = MStask92_';
% channel 2
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask46_ch2.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask46_times_ch2.mat');
t_T_ch2_46 = MStask46_times_ch2';
Yi_T_ch2_46 = MStask46_ch2';

t46_ch1 = size(Yi_T_ch1_46, 1);
ch1_46x = zeros(t46_ch1, 1);
sep = 11;

for i = 1:t46_ch1
    
    ch1_46x(i,1) = EstimateCell(Yi_T_ch1_46, i, sep, [1 1800], 1);
    
end

t46_ch2 = size(Yi_T_ch2_46, 1);
ch2_46x = zeros(t46_ch2, 1);
for i = 1:t46_ch2
    
    ch2_46x(i,1) = EstimateCell(Yi_T_ch2_46, i, sep, [1300 1800], 2);
    
end

figure();
hold on
plot(t_T_ch1_46, ch1_46x);
plot(t_T_ch2_46, ch2_46x);
hold off;
th = 10;
SS46 = zeros(size(ch1_46x, 1), 1);
d46 = ch1_46x(100) - ch2_46x(100);

SS46(ch1_46x > 0) = d46;

y46 = SS46;
q46 = cumtrapz(t_T_ch1_46, y46); % cumulative integral
[q46_th, index] = unique(q46);
th46 = interp1(q46_th, t_T_ch1_46(index), th);
figure();
hold on 
title('46%');
plot(t_T_ch1_46, y46)
plot(t_T_ch1_46, q46)
hold off

%%
% 92% contrast
close all;
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/MStask92_.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/MStask92_times.mat');
t_T_ch1_92 = MStask92_times_';
Yi_T_ch1_92 = MStask92_';
space = 1:3600;
figure(6);
hold on
h = surface(space, t_T_ch1_92, Yi_T_ch1_92);
view(0,0);
colormap hot;
title('Adaptive C-S output at 92%');
xlabel('Space');
ylabel('Time');
zlabel('Activity');
set(h,'LineStyle','none') % removes the grid lines 

% channel 2
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask92_ch2.mat');
load('/Users/Boris/Library/Mobile Documents/com~apple~CloudDocs/Boris_Documents/Second Paper/Modeling paper/Data/MStask/Trial2/MStask92_times_ch2.mat');
t_T_ch2_92 = MStask92_times_ch2';
Yi_T_ch2_92 = MStask92_ch2';
space = 1:3600;
%figure(6);
h = surface(space, t_T_ch2_92, Yi_T_ch2_92);
view(0,0);
colormap hot;
title('Adaptive C-S output at 92%');
xlabel('Space');
ylabel('Time');
zlabel('Activity');
set(h,'LineStyle','none') % removes the grid lines 
hold off

t92_ch1 = size(Yi_T_ch1_92, 1);
ch1_92x = zeros(t92_ch1, 1);
sep = 11;

for i = 1:t92_ch1
    
    ch1_92x(i,1) = EstimateCell(Yi_T_ch1_92, i, sep, [1 1800], 1);
    
end
t92_ch2 = size(Yi_T_ch2_92, 1);
ch2_92x = zeros(t92_ch2, 1);
for i = 1:t92_ch2
    
    ch2_92x(i,1) = EstimateCell(Yi_T_ch2_92, i, sep, [1300 1800], 2);
    
end

figure();
hold on
plot(t_T_ch1_92, ch1_92x);
plot(t_T_ch2_92, ch2_92x);
hold off;
d92 = ch1_92x(100) - ch2_92x(100);

th = 10;
SS92 = zeros(size(ch1_92x, 1), 1);

SS92(ch1_92x > 0) = d92;

y92 = SS92;
q92 = cumtrapz(t_T_ch1_92, y92); % cumulative integral
[q92_th, index] = unique(q92);
th92 = interp1(q92_th, t_T_ch1_92(index), th);
figure();
hold on 
title('92%');
plot(t_T_ch1_92, y92)
plot(t_T_ch1_92, q92)
hold off

%% final plot

close all;
contrast = [9.09 11.13 16.65 24.31 46.3908 92];
scale = 23.28/0.7;
ths = scale * [th9, th11, th15, th22, th92, th92];
figure()
hold on
ypos= tadin_std1(:,2) - Tadin(:,2);
yneg= Tadin(:,2) - tadin_std2(:,2);

p1=shadedErrorBar(contrast, Tadin(:,2), ypos, 'lineprops','-r','patchSaturation',0.33);
p2=plot(contrast, Tadin(:,2), '-r', 'MarkerFaceColor','r','MarkerSize',10,'Marker','o'); 
set(gca, 'XScale','log', 'YScale','log');
p3=createfig_MDtask(contrast, ths, 'none');
% plot(contrast, scaled_th, 'LineWidth',2,'Color',[0 0 0]);
legend([p2 p3], "Tadin's data", 'Model'); 
legend('boxon');
hold off

% Chi-square Goodness of fit test
stdv = ypos*2;
var = stdv.^2;
chi2 = sum((ths' - Tadin(:,2)).^2 ./ var);

