% Run Lorenz Models

% Example 1
figure(1);
[x,y,z] = Lorenz63(); 
plot3(x,y,z);
grid on

% Example 2
figure(2); 
[x,y,z] = Lorenz84(); 
plot3(x,y,z);
grid on

% Example 2
figure(3); 
X = Lorenz96();
plot3(X(:,1),X(:,2),X(:,3));
grid on