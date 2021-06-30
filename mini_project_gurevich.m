 %%% Research Skills, Mini project
%%% Michelle Gurevich

% mini_project_gurevich.m
% Description - Chaos and fractals
% This script provides three examples of chaotic dynamical systems and
% plots them with respect to the appropriate phase spaces. These are the
% Lorenz attractor, the Julia set, and a special case of the latter: the
% Mandelbrot set. All required background information is given in the
% accompanying document.

% Author: Michelle Gurevich
% Date: 15 Mar 2021

%% Part I: Lorenz attractor

% instantiate system parameters to those Lorenz used, though these can be
% adjusted (but should remain positive)
sigma = 10;
beta = 8/3;
rho = 28;

% set up the ODEs which describe fluid motion of the system (here v is the
% vector [x,y,z])
f = @(t,v) [sigma * (v(2) - v(1)); v(1) * (rho - v(3)) - v(2); ... 
    v(1) * v(2) - beta * v(3)];

% use ODE solver to solve system with intitial condition vector v = [0,1,0]
% (i.e. x = z = 0, y = 1) and over the range of 0 to 120 in t; second 
% vector will be used to demonstrate sensitivty to initial conditions
t = [0 120];
[t,v] = ode45(f,t,[0 1 0]);
[t,v2] = ode45(f,t,[0 .99999 0]);

% plot results using comet animation (note the z and y axes are switched
% for clarified animation)
[xmat,ymat,zmat] = peaks(100); 
xvec = v(:,1);
yvec = v(:,2);
zvec = v(:,3);
comet3(xvec,zvec,yvec)
title('Animated Lorenz attractor with \sigma=10, \beta=8/3, \rho=28')

% display plot with result normal to YZ plane (view is from directly
% overhead, i.e. azimuth = 0 and elevation = 90) thereby generating the 
% familiar 'butterfly' image
figure()
plot3(v(:,1),v(:,3),v(:,2))
hold on

% demonstrate sensitivity to initial conditions by plotting result of small
% initial perturbation (on the order of 10^-5) alongside original plot
plot3(v2(:,1),v2(:,3),v2(:,2))
view(2)
hold off
legend('original y condition','perturbed y')
title('Lorenz attractor with \sigma=10, \beta=8/3, \rho=28 seen from above')

%% Part II: Julia set

% call plot_julia_set() and pass through the c value corresponding to the
% Douady rabbit fractal 
c_douady = -0.123 + 0.745i;
plot_julia_set(c_douady)

%% Part III: Mandelbrot zoom

% call plot_mandelbrot_set() to generate plot of the Mandelbrot fractal
plot_mandelbrot_set()

%% Function definitions

% define julia_set function to take as argument a value c and use it to
% solve the recursive relation for the corresponding Julia set
function [z,f] = julia_set(c)
    % define x and y vectors, z
    [x,y] = meshgrid(-1.5:0.001:1.5,-1.5:0.001:1.5);
    z = x + 1i * y;
    % define recursive relation for Julia set
    f = @(s) s.^2 + c;
end

% define plot_julia_set function to plot the results of julia_set()
function plot_julia_set(c)
    [z,f] = julia_set(c);
    % iterate the recursive formula 20 times and assign output to z
    for k = 1:20
        z = f(z);
    end   
    % remove unbounded values so all values converge
    r = max([c,3]);
    z(abs(z) >= r) = 0;
    
    % plot results
    contour(abs(z))
    xlabel('Re(z)')
    ylabel('Im(z)')
    % get and use real and imaginary values of c for title
    cReIm = split(string(c),'+');
    title(sprintf('Julia set for c = %s + %s', cReIm(1), cReIm(2)));
end

% define plot_mandelbrot_set function to take no arguments but to generate
% and plot the Mandelbrot fractal
function plot_mandelbrot_set()
	% define x and y vectors, z
    [x,y] = meshgrid(-2:.0009:1,-1.5:.0009:1.5);
    z = x + 1i * y;
    c = x + 1i * y;
    % define recursive relation for Julia set
    f = @(s) s.^2 + c;
    k = zeros(size(c));
    % iterate the recursive formula 40 times and assign output to z
    for j = 1:40
        z = f(z);
        % remove unbounded values so all values converge
        k(abs(z) > 2 & k == 0) = 40 - j;
    end
    % plot results
    imagesc(k)
    colormap parula
    xlabel('Re(z)')
    ylabel('Im(z)')
    title('Mandelbrot set');
end
