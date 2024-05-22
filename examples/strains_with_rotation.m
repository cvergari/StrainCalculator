% This script uses StrainCalculator to compare differnt types of strain definition.
% A fake data series is generated (using :meth:`~StrainCalculator.generateFakeData`) 
% to obtain a uniform horizontal tension field *and a rigid rotation*. 
%
% The Green deformation field is graphed at each step, while the average
% strain over the whole surface is plotted at the end, to compare the
% different computation methods.
%
% The latter plot highlights how engineering strain is sensitive to
% rotation: it increases and then decreases suggesting compression , 
% although the sample still undergoes tension (but it is rotated).
% Engineering strain and Almansi strain also show a "phantom" vertical
% component of strain.
% 
% .. figure:: images/strain_rotation_1.png
%
%    Last strain step, with the superimposition of a horizontal strain and
%    a 60° counterclockwise rotation
%
% .. figure:: images/strain_rotation_2.png
%
%    Comparison of the different strain calculations.
%


clear
close all
clc


% Add path of parent folder to Matlab path
addpath([fileparts(mfilename('fullpath')) '/../']);

% Prepare an empty object
obj = StrainCalculator();

% Generate fake data for Horizotnal tension
strainType = 'HorizontalTension';
obj.generateFakeData('type', strainType, 'maxStrain', 0.5);

% Add an increasing rotation component at each step
max_angle = 60;  % degrees
obj = add_rotation(obj, max_angle);

% Recalculate displacements
obj.calcDisplacements();

Methods = {'Engineering', 'LogStrain', 'True', 'Almansi', 'Green'};

for k = 1 : length(Methods)

    obj.calcStrain(Methods{k});
    
    % Strain is uniform, so take a single sample for each step and each
    % method.
    ex(1,k) = 0;
    for n = 2 : length(obj.strain.X)
        ex(n,k) = mean(obj.strain.X{n}(1));
    end
    ey(1,k) = 0;
    for n = 2 : length(obj.strain.X)
        ey(n,k) = mean(obj.strain.Y{n}(1));
    end
    
    % Principal strains 
    eI(1,k) = 0;
    for n = 2 : length(obj.strain.X)
        eI(n,k) = mean(obj.strain.PrincipalStrainI{n}(1));
    end
    eII(1,k) = 0;
    for n = 2 : length(obj.strain.X)
        eII(n,k) = mean(obj.strain.PrincipalStrainII{n}(1));
    end
        
end

% Plot strains at each step using Green strain
obj.calcStrain('Green');
figure
for k = 1 : obj.nSteps
    cla
    obj.plotGrid(k)
    hold on;
    if k < 2  ;   continue;   end
    
    % Plot strain from second step
    h = obj.plotStrain('step', k, 'direction', 'X');
    set(h, 'FaceAlpha', 0.5)
    axis([-1 1.2 0 1.8])

    if k == 2
        c = colorbar(gca);
        c.Label.String = 'Strain';
        caxis(obj.strain.limits.('X'))
    end
    title(['Step n. ' num2str(k) '/' num2str(obj.nSteps)])
    
    drawnow
    pause(0.3)    
    
end

figure
subplot(2,2,1)
plot(ex)
xlabel('Step n.')
ylabel('Horizontal strain')
legend(Methods, 'Location', 'NorthWest');
title('Horizontal strains')

subplot(2,2,2)
plot(ey)
xlabel('Step n.')
ylabel('Vertical strain')
title('Vertical strains')

subplot(2,2,3)
plot(eI)
xlabel('Step n.')
ylabel('Strain')
title('Principal strain I')

subplot(2,2,4)
plot(eII)
xlabel('Step n.')
ylabel('Strain')
title('Principal strain II')



function obj = add_rotation(obj, max_angle)
    % Add a rotational component to the positions at each step
    step_angles = linspace(0, max_angle, obj.nSteps);
    
    % Add increasing angle to the positions
    for step = 1 : length(step_angles)
         theta = step_angles(step);
         R = [cosd(theta), -sind(theta);
              sind(theta), cosd(theta)];
         % extract coordinates
         xy = [obj.Xpositions{step}(:), ...
               obj.Ypositions{step}(:)];
         % Rotate
         xy = R * xy';
         % Set them back in the object
         obj.Xpositions{step} = reshape(xy(1,:), size(obj.Xpositions{step})); 
         obj.Ypositions{step} = reshape(xy(2,:), size(obj.Ypositions{step})); 
         
    end
    
    
end

