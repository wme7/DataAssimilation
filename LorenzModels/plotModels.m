% Make preaty plots

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultTextFontName','Times',...
'DefaultTextFontSize',20,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',20,...
'DefaultLineLineWidth',1.5,...
'DefaultAxesBox','on',...
'defaultAxesLineWidth',1.5,...
'DefaultFigureColor','w',...
'DefaultLineMarkerSize',7.75)

% Example 1
figure(1);
[x,y,z] = Lorenz63(); 
plot3(x,y,z); 
view(-138,33);
title('Lorenz-63');
xlabel('x');
ylabel('y');
zlabel('z');
grid on
fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 6 6];
print('-dpng' ,'Lorenz63.png');
fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 6 6];
print('-depsc' ,'Lorenz63.eps');

% Example 2
figure(2); 
[x,y,z] = Lorenz84(); 
plot3(x,y,z); 
view(-59,37);
title('Lorenz-84');
xlabel('x');
ylabel('y');
zlabel('z');
grid on
fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 6 6];
print('-dpng' ,'Lorenz84.png');
fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 6 6];
print('-depsc' ,'Lorenz84.eps');

% Example 2
figure(3); 
X = Lorenz96();
plot3(X(:,1),X(:,2),X(:,3)); 
view(-47,25);
title('Lorenz-96')
xlabel('x');
ylabel('y');
zlabel('z');
grid on
fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 6 6];
print('-dpng' ,'Lorenz96.png');
fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 6 6];
print('-depsc' ,'Lorenz96.eps');