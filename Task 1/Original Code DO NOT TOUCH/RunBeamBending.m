clc, clear all, close all

N = 100; %the beam solver discretises beam into 100 nodes
t = 1e-3*ones(1,N); %set the skin thickess of the blade at every node
c = 0.15; %define the chord length - KEEP FIXED

I = secondMomentArea(t,c); %caclulate the second moment of area

%submit to beamBending:
% number of nodes for numerical solve, N
% vector of second moment of area values across the span, I 
% (I must have N elements)
[delta, M, V, L, P]  = beamBending(N,I); 
%beamBending function returns:
%beam deflection, delta
%bending moment distribution, M
%shear force distriubition, V
%blade length is returned, L
%distributed aerodynamic loading, P

%define the spanwise position of every node (in meters)
nodeLocations = linspace(0,L,N);

%plot the results
subplot(2,1,1) %tells MATLAB you want 2 subplots in one figure (check the docs for more info) 
plot(nodeLocations,delta,'ko-') %plot the deflection (delta) at each node location across the blade
xlabel('distance, m') %set labels
ylabel('deflection, mm')

subplot(2,1,2) %plot on the second subplot
plot(nodeLocations,(M*(0.06*c)./I)/(1e6),'ko-') %plot the bending stress (M*y / I) at every node location
xlabel('distance, m') %set labels
ylabel('Bending Stress, MPa')


%define the second moment of area of each blade cross section
%see approaches outlined in the theory document on moodle
function I = secondMomentArea(t,c)

    %for the sake of having a working script, I am going to make a simple
    %model of how the second moment of area scales with thickness of
    %material (say skin or box section) and the chord length of the blade

    I = (10e-2)*c.*t.^2;

    %This is a placeholder for you to add in a real second moment of area
    %calculation.

end
