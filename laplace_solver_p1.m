%% Adapted from Prof. Mark Lawrence
% Jerry Gao
% Prof. Lawrence
% FL25 ESE 3300/330 Engineering Electromagnets Principles
% Midterm Case Study 1
% 18 October 2025
% Problem 1
%-------------------------------------------------------------------------%
%  This simple program computes the Electric Fields due to 
%  Parallel plate Capacitors using the Finite difference method (FDM)  
%-------------------------------------------------------------------------%
clc;
close all;
clear all;

%-------------------------------------------------------------------------%
%                   SYMBOLS USED IN THIS CODE                             
%-------------------------------------------------------------------------%
% E = Total electric field matrix using Poisson's equation
% V = Potential matrix
% Nx = Number of grid points in X-Direction
% Ny = Number of grid points in Y-Direction
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                         INITIALIZATION                                  
%          Here, all the grid, size, charges, etc. are defined
%-------------------------------------------------------------------------%
% Enter the dimensions
Nx = 201;     % Number of X-grids
Ny = 201;     % Number of Y-grids
mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y
   
Ni = 100;  % Number of iterations for the Poisson solver (original 750)
V = zeros(Nx,Ny);   % Potential (Voltage) matrix
T = 0;            % Top-wall potential
B = 0;            % Bottom-wall potential
L = 0;            % Left-wall potential
R = 0;            % Right-wall potential

%-------------------------------------------------------------------------%
% Initializing edges potentials (DO NOT MODIFY)
%-------------------------------------------------------------------------%
V(1,:) = L;
V(Nx,:) = R;
V(:,1) = B;
V(:,Ny) = T;

%-------------------------------------------------------------------------%
% Initializing Corner potentials (DO NOT MODIFY)
%-------------------------------------------------------------------------%
V(1,1) = 0.5 * (V(1,2) + V(2,1));
V(Nx,1) = 0.5 * (V(Nx-1,1) + V(Nx,2));
V(1,Ny) = 0.5 * (V(1,Ny-1) + V(2,Ny));
V(Nx,Ny) = 0.5 * (V(Nx,Ny-1) + V(Nx-1,Ny));

%-------------------------------------------------------------------------%
% 201 grids
% 200 grid spacings
% each grid spacing = 10 nm
% total grid side length = 2000 nm
length_plate = 51;  % Length of plate in terms of number of grids (500 nm)
lp = floor(length_plate/2);
% KEY VARIABLE: separation distance s >= 52 (waveguide & nano wires cannot
% physically touch)
s = 52; % 52 pixels = 520 nm
% KEY VARIABLE: nanowire width d
d = 1; % 1 pixel = 10 nm
position_plate = s / 2; % Position of plate on x axis
pp1 = mpx + position_plate;
pp2 = mpx - position_plate;

% Forcing equipotential wires
conductor = false(Nx,Ny);
conductor(pp1:pp1+d,mpy-lp:mpy+lp) = true; % Inside nanowire
conductor(pp2-d:pp2,mpy-lp:mpy+lp) = true; % Inside nanowire

% Permittivity (epsilon) "field"
% Note on the section below: I did not realize that epsilon_r = 1 til I was
% almost done with implementation, i.e. the waveguide's permittivity
% constant is exactly that of free space/vacuum so doing this in this case
% is lowkey overkill, but this neatly sets up the next problem
epsilon_0 = 8.854e-12;
epsilon_r = 1;
epsilon_waveguide = epsilon_0 * epsilon_r;
epsilon_field = ones(Nx,Ny) * epsilon_0;
epsilon_field(76:126,76:126) = epsilon_waveguide;

for z = 1:Ni    % Number of iterations
    for i=2:Nx-1
        for j=2:Ny-1
            % The next two lines are meant to force the matrix to hold the 
            % potential values for all iterations
            V(pp1:pp1+d,mpy-lp:mpy+lp) = 0.5; % ACTUAL BCS (V(x=s) = 1)
            V(pp2-d:pp2,mpy-lp:mpy+lp) = -0.5; % ACTUAL BCS (V(x=0) = 0)
            % BC Problem: original BCs run from 0 to 1, but grounding
            % everywhere else (edges, corners, etc.) to 0V by default
            % results in 0V wire disappearing
            % Solution: offset electric potentials so to range from
            % -0.5 to +0.5
            %V(mpy-lp2:mpy+lp2,mpy-lp2:mpy+lp2) = -1;
            if ~(conductor(i,j))
                epsilon_x_positive = 0.5 * (epsilon_field(i,j) + epsilon_field(i+1,j));
                epsilon_x_negative = 0.5 * (epsilon_field(i,j) + epsilon_field(i-1,j));
                epsilon_y_positive = 0.5 * (epsilon_field(i,j) + epsilon_field(i,j+1));
                epsilon_y_negative = 0.5 * (epsilon_field(i,j) + epsilon_field(i,j-1));
                V(i,j) = (epsilon_x_positive * V(i+1,j) + epsilon_x_negative * V(i-1,j) + ...
                    epsilon_y_positive * V(i,j+1) + epsilon_y_negative * V(i,j-1)) / ( ...
                    epsilon_x_positive + epsilon_x_negative + epsilon_y_positive + epsilon_y_negative);
            end
        end
    end
end

% Take transpose for proper x-y orientation
V = V';
[Ex,Ey]=gradient(V);
Ex = -Ex * 10e8; % rescale from V/pixel back to SI unit V/m
Ey = -Ey * 10e8;

% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2);
x = (1:Nx)-mpx;
y = (1:Ny)-mpy;

%% Contour Display for electric potential
figure(1)
contour_range_V = -101:0.5:101;
contourf(x,y,V,50,'LineColor','none')%,contour_range_V,'linewidth',0.5);%contour(x,y,V,contour_range_V,'linewidth',0.5)
%hold on, quiver(x,y,Ex,Ey,3,"k")
title('Electric Potential Distribution V(x,y) in Volts','fontsize',14,'color','black');
axis([min(x) max(x) min(y) max(y)]);
cb = colorbar('location','eastoutside','fontsize',14,'color','black');
new_voltage_labels = arrayfun(@(x) sprintf('%.1f',x),cb.Ticks+0.5, 'uniformoutput', false);
cb.TickLabels = new_voltage_labels;
xlabel('X-Axis in Nanometers (10^{-9} m)','fontsize',14,'color','black');
ylabel('Y-Axis in Nanometers (10^{-9} m)','fontsize',14,'color','black');
h1=gca;
h1.XTickLabel = h1.XTick * 10;
h1.YTickLabel = h1.YTick * 10;
h1.XAxis.TickLabelColor = 'black';
h1.YAxis.TickLabelColor = 'black';
h1.Color = 'black';
set(h1,'fontsize',14);
fh1 = figure(1);
set(fh1, 'color', 'white')

%% Contour Display for electric field
figure(2)
contour_range_E = -20:0.05:20;
contourf(x,y,E,50,'LineColor', 'None');%,contour_range_E,'linewidth',0.5);%contour(x,y,E,contour_range_E,'linewidth',0.5)
title('Electric Field Distribution E(x,y) in V/m','fontsize',14,'color','black');
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14,'color','black');
xlabel('X-Axis in Nanometers (10^{-9} m)','fontsize',14,'color','black');
ylabel('Y-Axis in Nanometers (10^{-9} m)','fontsize',14,'color','black');
h2=gca;
h2.XTickLabel = h2.XTick * 10;
h2.YTickLabel = h2.YTick * 10;
h2.XAxis.TickLabelColor = 'black';
h2.YAxis.TickLabelColor = 'black';
h2.Color = 'black';
set(h2,'fontsize',14);
fh2 = figure(2);
set(fh2, 'color', 'white')

%% Quiver Display for electric field Lines
figure(3)
contour(x,y,E,'linewidth',0.5);
hold on, quiver(x,y,Ex,Ey,6)%quiver(x,y,Ex,Ey,2)
title('Electric Field Vector Lines E(x,y) in V/m','fontsize',14,'color','black');
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14,'color','black');
xlabel('X-Axis in Nanometers (10^{-9} m)','fontsize',14,'color','black');
ylabel('Y-Axis in Nanometers (10^{-9} m)','fontsize',14,'color','black');
h3=gca;
h3.XTickLabel = h3.XTick * 10;
h3.YTickLabel = h3.YTick * 10;
h3.XAxis.TickLabelColor = 'black';
h3.YAxis.TickLabelColor = 'black';
h3.Color = 'black';
set(h3,'fontsize',14);
fh3 = figure(3);
set(fh3, 'color', 'white')

%% (Approximate) Modulation strength calculation
waveguide = E(76:126,76:126); % 500 nm by 500 nm centered at the origin
average_field_strength = mean(waveguide(:));
disp(average_field_strength);

%% Maximum speed of modulation calculation (based on r & c given by tuned s & d)
conductivity = 6.3e7;
cross_sectional_area = 500e-9 * d * 10e-9;
resistance = 50e-6 / (conductivity * cross_sectional_area);
disp(resistance);

surface_area = (cross_sectional_area + 500e-9 * 50e-6 + d * 10e-9 * 50e-6) * 2;
capacitance = epsilon_waveguide * surface_area / (s * 10e-9);
disp(capacitance);

modulation_speed = 1 / (resistance * capacitance);
disp(modulation_speed);

%-------------------------------------------------------------------------%
% REFERENCE
%           SADIKU, ELEMENTS OF ELECTROMAGNETICS, 4TH EDITION, OXFORD
%-------------------------------------------------------------------------%

%% NEXT STEPS:
% 1. Hand-tune values of s & d to find largest average magnitude of e field
% through the 500x500 waveguide (modulation strength) (done for now; can
% include plot for varying values later)
% 2. Consider RC properties for modulation speed (C found via surface
% area) (TODO: consider bcs for current inside wire