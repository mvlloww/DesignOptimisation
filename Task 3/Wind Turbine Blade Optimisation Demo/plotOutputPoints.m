function plotOutputPoints(c, gamma, r, i)
    %plots blade section and outputs coordinates
    %Inputs:
    % 1) chord length
    % 2) gamma angle (angle of blade chordline to rotor plane
    % 3) radius of blade element
    % 4) iteration count (i.e. blade element number)

    load('S822_Points.mat');
    
    % Scale the points using the current chord length
    profile_Xa = S822.x;
    profile_Ya = S822.y;

    % Mirror y-coordinates
    profile_Xa = -profile_Xa;

    % Scale the points by multiplying the unit chord length by the actual
    % chord length
    profile_Xa = profile_Xa * c;
    profile_Ya = profile_Ya * c;

    profile_Xb = 0.5 * profile_Xa;
    profile_Yb = 0.5 * profile_Ya;

    profile_Xb = profile_Xb - (max(profile_Xb) - min(profile_Xb))/2;

    % Rotate points to get correct gamma angle
    profile_X2a = profile_Xa * cos(-gamma) - profile_Ya * sin(-gamma);
    profile_Y2a = profile_Xa * sin(-gamma) + profile_Ya * cos(-gamma);

    profile_X2b = profile_Xb * cos(-gamma) - profile_Yb * sin(-gamma);
    profile_Y2b = profile_Xb * sin(-gamma) + profile_Yb * cos(-gamma);

    % Plot airfoil section in correct orientation for each blade element
    hold on
    figure(1),
    plot3(profile_X2a, profile_Y2a, r*ones(size(profile_X2a)));
    axis square
    axis equal
    
    if i == 1
        quiver3(0, -0.1, r(1), 0, 0.1, 0)
    end
    xlim([-0.3, 0.1])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view([45,45])
    grid minor
    drawnow
    hold off
    
    %%% Write profile coordinates to one CSV file per blade element
    fileName1 = strcat('Foil_Points_', num2str(i), 'a.csv');
    fileName2 = strcat('Foil_Points_', num2str(i), 'b.csv');

    csvwrite(fileName1,100*[profile_X2a(2:end-1), profile_Y2a(2:end-1), r*ones(length(profile_X2a(2:end-1)), 1)]);
    csvwrite(fileName2,100*[profile_X2b(2:end-1), profile_Y2b(2:end-1), r*ones(length(profile_X2b(2:end-1)), 1)]);
   
end

