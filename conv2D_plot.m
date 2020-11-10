% code to plot figure 2b
close all
% load('workspa_23Jun');

% 2D Gaussian --> Low-pass filter 
x = linspace(0, 499, 499); % 499
y = linspace(0, 499, 499);
[X, Y] = meshgrid(x, y);
sigxS= 250; sigyS= 250;
Rfs = Gaussian2D(0.5, X, Y, sigxS, sigyS, 250, 'S');
figure();
h1 = surf(X,Y,Rfs);
set(h1,'LineStyle','none') % removes the grid lines 

sigxC= 2; sigyC= 2;
Rfc = Gaussian2D(2, X, Y, sigxC, sigyC, 250, 'C');
figure();
h2 = surf(X,Y,Rfc);
set(h2,'LineStyle','none') % removes the grid lines 

DoG2d = Rfc - Rfs;
figure();
inicio = 215; fin = 285;
h3 = surf(X(inicio:fin, inicio:fin),Y(inicio:fin,inicio:fin),DoG2d(inicio:fin,inicio:fin));
axis tight; view(90,90);
set(h3,'LineStyle','none') % removes the grid lines 


% load image of Ruben Blades 2
ruben2 = imread('Ruben2.jpeg');
ruben2 = imresize(ruben2, [512, 512]);
ruben2 = rgb2gray(ruben2);
Ruben2 = im2double(ruben2);

% convert pixel images to percentage.
RubenPerc = (Ruben2./max(max(Ruben2)));
Bckg = [0.1, 1, 10, 100, 1000, 10000, 1000000]; % background intensities
n = 7;
Ruben_norm = cell(n, 1);
for k = 1:n
    
    Ruben_Bckg = 255*RubenPerc * Bckg(k) + Bckg(k);

    % perform the convolution operation
    Ruben_Ex = conv2((Ruben_Bckg), Rfc, 'same'); % high-pass filter
    Ruben_Inh = conv2((Ruben_Bckg), Rfs, 'same');

    % calculate steady-state response
    A = 0.5; B = 1; D = 0; scale = 255;
    Ruben_norm{k} = scale * (B * Ruben_Ex - D * Ruben_Inh) ./(A + Ruben_Ex + Ruben_Inh);
end

figure(4)
subplot(1,2,1)
hold on
imshow(ruben2);
title('Original scene', 'FontSize', 16);
hold off
subplot(1,2,2)
hold on
imshow(uint8(255 * Ruben_norm{2}));
title('Low-passed output', 'FontSize', 16);
hold off

figure()
subplot(2,1,1);
plot(Ruben_norm{5}(257,:), 'Linewidth', 2);
% plot(g * Ruben_Xij_B100_(257,:), 'Linewidth', 2);
legend({'1000'}, 'FontSize',10);
xlabel('Pixels', 'FontSize', 17);
ylabel('Normalized Activity', 'FontSize', 17);
lgd = legend({'1000'});
title(lgd, 'Luminance (cd m^{-2})', 'Interpreter','tex');
axis tight;
hold off

means = [mean(mean(Ruben_norm{1})), mean(mean(Ruben_norm{2})), mean(mean(Ruben_norm{4})), ...
    mean(mean(Ruben_norm{5})), mean(mean(Ruben_norm{7}))];
subplot(2,1,2);
X = categorical({'1','10','100','1000', '1000000'});
bar(X, means);
xlabel('Luminance (cd m^{-2})', 'Interpreter','tex');
ylabel('Normalized average pixel output');


%%
function res = Gaussian2D(A, x, y, sigx, sigy, c, type)
        if type == 'S'        
            res = A * exp(-((x - c).^2/(2*sigx^2) + (y - c).^2/(2*sigy^2))) + 0*0.009;
        else
            res = A * exp(-((x - c).^2/(2*sigx^2) + (y - c).^2/(2*sigy^2)));
        end
end
