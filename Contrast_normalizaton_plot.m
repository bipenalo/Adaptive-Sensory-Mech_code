% Background normalization plot
close all

% values obtained with the NormalizationStage.jl file 
L0_1 = [0.2047139, 0.2126243, 0.2231004, 0.2654416, 0.3592341, 0.6306032];
L1 = [0.2142464, 0.2325554, 0.2426315, 0.2873892, 0.3820787, 0.6381485];
L10 = [0.216066, 0.234756, 0.2523096, 0.2897853, 0.384524, 0.638913];
L100 = [0.2166249, 0.2349783, 0.2525395, 0.2900271, 0.3847703, 0.6389895];
L1000 = [0.2166468,0.2350006, 0.2525625, 0.2900513, 0.3847949, 0.6389972];
L10000 = [0.216649, 0.2350028, 0.2525648,  0.2900537, 0.3847974, 0.638998];
L100000 = [ 0.2166492, 0.235003, 0.252565, 0.2900739, 0.3848976, 0.639998];

contrast = [2.8, 5.5, 11, 22, 46, 92];
hold on
% plot(contrast, L0_1);
% plot(contrast, L10);
% plot(contrast, L100);
plot(contrast, L1000);
% plot(contrast, L10000);
% plot(contrast, L100000);
% legend('0.1', '1', '10', '100', '1000', '10000', '1000000');
lgd1 = legend('10000');
title(lgd1, 'Luminance (cd m^{-2})', 'Interpreter','tex');
xlabel('Contrast (%)');
ylabel('Normalized output');
hold off

figure(2);
hold on
val2_8 = [0.20471 0.21424 0.21662 0.21664 0.21665; ... 
          0.2126243 0.23255 0.23497 0.23500 0.23501; ...
          0.22310 0.24263 0.25253 0.25256 0.25257; ...
          0.26544 0.28739 0.29002 0.29005 0.29007; ...
          0.3592 0.3820 0.3847 0.38479 0.38489; ...
          0.63060 0.63814 0.63898 0.63899 0.63999]; % luminance 0.1, 1, 100, 1000, 1000000
X = categorical({'0.1','1','100','1000', '1000000'});
X = reordercats(X,{'0.1','1','100','1000', '1000000'}); % reorder them
bar(X, val2_8', 'stacked');
xlabel('Luminance (cd m^{-2})', 'Interpreter','tex');
ylabel('Cumulative normalized output');
lgd = legend({'2.8', '5.5', '11', '22', '46', '92'}, 'Orientation','horizontal');
title(lgd, 'Contrast (%)');
hold off

%% single plot
subplot(1,2,1);
hold on
val2_8 = [0.20471 0.21424 0.21662 0.21664 0.21665; ... 
          0.2126243 0.23255 0.23497 0.23500 0.23501; ...
          0.22310 0.24263 0.25253 0.25256 0.25257; ...
          0.26544 0.28739 0.29002 0.29005 0.29007; ...
          0.3592 0.3820 0.3847 0.38479 0.38489; ...
          0.63060 0.63814 0.63898 0.63899 0.63999]; % luminance 0.1, 1, 100, 1000, 1000000
X = categorical({'0.1','1','100','1000', '1000000'});
X = reordercats(X,{'0.1','1','100','1000', '1000000'}); % reorder them
bar(X, val2_8', 'stacked');
xlabel('Luminance (cd m^{-2})', 'Interpreter','tex');
ylabel('Cumulative normalized output');
lgd = legend({'2.8', '5.5', '11', '22', '46', '92'}, 'Orientation','horizontal');
title(lgd, 'Contrast (%)');
hold off


subplot(1,2,2);
hold on
plot(contrast, L1000);
lgd1 = legend('10000');
title(lgd1, 'Luminance (cd m^{-2})', 'Interpreter','tex');
xlabel('Contrast (%)');
ylabel('Normalized output');
hold off
