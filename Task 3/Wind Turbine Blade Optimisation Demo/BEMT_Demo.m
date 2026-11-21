function CP = BEMT_Demo(B, radiusLimits, lambda, chordProfile, AoAProfile, noElements, V, kinematicViscosity, convergenceTolerance, iterationLimit, plot)
% Clear Command Window
clc
% Clear figures
close all

%% Notes on model inputs
% Input 1: number of blades

% Input 2: start and end radius for blade in metres e.g. [0.2, 0.5]

% Input 3: tip-speed ratio, e.g. 4

% Input 4: chord lengths at root and tip in metres e.g [0.08, 0.038]

% Input 5: Angle of Attack (AoA) at root and tip in degrees e.g [5, 2]

% Input 6: number of blade elements (at least 10)

% Input 7: velocity of free-stream in m/s (incoming air)

% Input 8: kinematic viscosity of the air e.g. 1.5e-5

% Input 9: the maximum allowable fluctuation in a and aPrime to conclude
%          convergence

% Input 10: The limit on the number of iterations for the iterative solve.

% Input 11: plot toggle - "1" to enable plotting of blade, "0" to suppress
%           will greatly affect speed, so be sure to input "0" for use in optimisation loops!
%           note this will also save the blade section coordinates as CSV
%           files - feel free to modufy/suppress this in plotOutputPoints.m

%% Introductory matter

% Load previously created lookup table for S822 airfoil CL and CD values
load('S822_CL_Lookup.mat');
load('S822_CD_Lookup.mat');
load('S822_Alpha_Values.mat');
load('S822_Re_Values.mat');

%loads the S822 aerofoil coordinates, used for plotting
load('S822_Points.mat');


% Array of radii for each blade element
r = linspace(radiusLimits(1), radiusLimits(2), noElements);

%Omega
omega = V*lambda/radiusLimits(2);


% Chord length - you can enter 2 values as a vector [C1 C2] and it will
% linearly interpolate from C1 at the blade root to C2 at the blade tip
if length(chordProfile) == 2
    c = linspace(chordProfile(1), chordProfile(2), noElements);
end

% Angle of attack - you can enter 2 values as a vector [A1 A2] and it will
% lienarly interpolate from C1 at the blade root to C2 at the blade tip
if length(AoAProfile) == 2
    AoA = deg2rad(linspace(AoAProfile(1), AoAProfile(2), noElements));
end


%% Iterative BEMT Solution
for i = 1:noElements %loop through each local blade section, from root to tip
    
    % guess the axial and angular induction factors - the solver will then
    % iterate from this reasonable starting condition to find the correct
    % values that balance the Blade Element Momentum Theory equations.

    %calculate local blade solidity
    sigma = B * c(i) / (2*pi*r(i));

    %set the initial conditions of the solver
    if i == 1
        a = 0.33;
        aPrime = 0;
    end


    
    %initialise changes (delta) in axial and angular induction factors
    %this iterative solver tracks these changes to check for convergence -
    %i.e. when the solution stops changing appreciably, the while loop
    %below will stop
    deltaA = inf;
    deltaAPrime = inf;
    
    %initialise iteration count - if the while loop fails to converge, a
    %stop condition is reached at iter = iterationLimit
    iter = 0;    
   

    while (abs(deltaA) > convergenceTolerance || abs(deltaAPrime) > convergenceTolerance) && iter < iterationLimit
        
        % Compute the inflow angle
        phi = atan((1-a)/(lambda * (r(i)/r(end)) * (1+aPrime)));
        
        %calculate the required, local blade pitch angle to acheive the
        %angle of attack specified, given the inflow angle
        gamma = phi - AoA(i);
        
        %calulate the relative velocity at the loacl blade section - this
        %will affect the amount of force the blade section produces
        relWind = V * (1-a) / sin(phi);
        
        %caclulate the local Reynolds number at the local blade section
        Re = relWind * c(i) / kinematicViscosity;

        %all blade sections must be within the bounds of the data set
        %provided - the folowing code saturates the Reynolds number value
        %at the lower or upper bound respectively if the local Reynolds
        %number goes outside these bounds.
        %This is just to keep your code running - but this should be
        %avoided, as you are not looking up the correct CL and CD values
        %and hence will not be simulating the right physics
        if Re < min(S822_Re_Values) 
            disp('Warning - Local Reynolds number is outside lower bound of supplied aerofoil data - solver to use data at lower bound');
            Re = min(S822_Re_Values);
        elseif Re > max(S822_Re_Values) 
            disp('Warning - Local Reynolds number is outside upper bound of supplied aerofoil data - solver to use data at upper bound');
            Re = max(S822_Re_Values);
        end
        
        % Look up the CL and CD through linear interpolation of the input data
        CL = griddata(S822_Alpha_Values, S822_Re_Values, S822_CL_Lookup', rad2deg(AoA(i)), Re);
        CD = griddata(S822_Alpha_Values, S822_Re_Values, S822_CD_Lookup', rad2deg(AoA(i)), Re);
        
        %transform CL and CD to normal and tangential forces
        Cn = CL * cos(phi) + CD * sin(phi);
        Ct = CL * sin(phi) - CD * cos(phi);

        %calcutlate axial induction factor
        aNew = 1 / (4 * sin(phi)^2 / (sigma * Cn) + 1);

 
        %calculate angular induction factor
        aPrimeNew = 1 / (4 * sin(phi) * cos(phi) / (sigma * Ct) - 1);
        
        
        %calculate deltas (used for convergence of while loop)
        deltaA = aNew - a;
        deltaAPrime = aPrimeNew - aPrime;
        
        %update induction factors
        a = a + deltaA;
        aPrime = aPrime + deltaAPrime;
        
        %update itation count
        iter = iter + 1;
        
    end

    %store blade values after each blade section is converged
    FtStore(i)      = 0.5 * Ct * 1.225 * c(i) * relWind^2; %force per unit span in the tangential direction (useful direction that pulls the turbine blade round)
    FnStore(i)      = 0.5 * Cn * 1.225 * c(i) * relWind^2; %force per unit span in the normal direction (direction that bends the blades)
    a_store(i)      = a; %axial induction factor
    aPrime_store(i) = aPrime; %axial induction factor
    phiStore(i)     = phi; %inflow angle
    gammaStore(i)   = gamma; %blade twist angle
    ReStore(i)      = Re; %Reynolds number
    CLStore(i)      = CL; %Coefficient of lift
    CDStore(i)      = CD; %Coefficient of drag
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Build the set of points to represent the current blade element %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plot == true %toggle plotting on and off with a "1" or "0" in the "plot" input argument
        plotOutputPoints(c(i), gammaStore(i), r(i), i);
    end
    
end


%% Perform numerical integration across blade elements
%so far we have force per unit span across the blade - to get this into a
%torque value, we need to integrate this force per unit span 

%calculate radial spacing between blade elements (equal spacing)
h = r(2) - r(1);

%integrate the force per unit span * radial position across the blade, then
%multiply by the number of blades
Torque = B * trapz(r, r.*FtStore);

%Power is Torque multiplied by angluar velocity, omega
Power = omega * Torque;

%calculate wind power in Watts
windPower = 0.5 * 1.225 * pi * r(end)^2 * V^3;

%calculate coefficient of power
CP = Power / windPower;

%calculate coefficient of thrust
CT = Torque / (windPower/V);

%The function outputs CP only currently, but this can be updated using
%square bracket notation
%[CP, CT, ...] = BEMT_Demo()






