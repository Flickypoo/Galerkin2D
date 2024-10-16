function [u,u_x,u_y] = computeExactSol(pts,solData)

% Data
alpha = solData.alpha;
beta = solData.beta;

% Coordinates
x = pts(:,1);
y = pts(:,2); 

u = sin(alpha*pi*x).*sin(beta*pi*y); 
u_x = alpha*pi*cos(alpha*pi*x).*sin(beta*pi*y);
u_y = beta*pi*sin(alpha*pi*x).*cos(beta*pi*y);