clc;
close all;
clear all;

%% HEADER (Problem 1)
% Jerry Gao
% Prof. Lawrence
% FL25 ESE 3300/330 Engineering Electromagnets Principles
% Midterm Case Study 1
% 18 October 2025

%% Solve for electric field between nanowire electrodes
% CONT: Solve 3D Laplace's Equation in Cartesian Coords to find V then E
% BCs: V(y=0) = 0; V(y=s) = 1
% New BCs Oct 20 in class: inside wire at width center vs one side

%% Determine how electro-optic performance varies w/ respect to nanowire
% width "d" and nanowire separation distance "s"


%% Estimate maximum speed of modulation based on resisitive & capacitive
% properties


%% TODO: UPDATED QUESTIONS for 10/21 office hours
% 1. Is E field direction incorrect in Figure 2? (ANS: yes but does not
% matter in terms of solving the problems)
% 2. Is there a specific quantity for max speed of modulation? Same as
% bandwidth for RC low pass filter? Look up material property values for R
% & C e.g. silver nanowire resistivity/conductivity? ANS: 1/(rc)
% 3. "Account for a nonuniform dielectric constant" but epsilon_r for both
% problems are constant? (ANS: higher epsilon_r for problem 2 only applies
% to orange material; everything else same thus now non uniform)
% 4. Correct extension from 2d to 3d in laplace_solver for FDM? How do we
% handle Y & Z boundary conditions (only x-direction defined for 0V & 1V;
% z-direction -> V independent of wire length) ANS: (no extension; need
% depth 50 micrometers for surface area calculation for capacitance)
% 5. Include all calculations in report including solved by MATLAB? Submit
% code separately? Hard to fit everything into 2-3 pages

% modulation strength -> e field / area; try diff values of d & s