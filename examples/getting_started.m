% StrainCalculator calculates strains from a displacement field. One way to
% use it, is to virtually generate a grid and the position of each grid element at 
% each time step. In this example we will do just that, and then we will 
% let StrainCalculator compute displacements and
% strains, and we will compare the strains to known values.
%
% Let's first create the object, and a simple grid:
%
% .. code:: matlab
%
%    % Create object
%    obj = StrainCalculator();
%
%    % Generate grid
%    [X0, Y0] = meshgrid(0:10, 0:10);
%    % Center grid for convenience
%    X0 = X0 - 5;  
%    Y0 = Y0 - 5;
%    % Plot
%    plot(X0(:), Y0(:),'b.')
%
% The grid spacing is 1 both in the horizontal and vertical direction. This
% will be automatically computed by StrainCalculator, as it assumes an
% homogeneous grid.
%
% Now we have to stretch the grid. Let's simulate a tension test in the 
% horizontal direction , from rest to 10% strain (as engineering strain) 
% in 10 loading steps.
%
% .. code:: matlab
%
%   % Store grid elements position at each step
%   X = cell(11,1); 
%   Y = cell(11,1);
%
%   max_step = 10;
%   % Notice that *step* starts from zero to store the resting step
%   for step = 0 : max_step
%       current_strain = 0.1 * step / max_step;
%       X{step + 1} = X0 + X0 .* current_strain;
%       Y{step + 1} = Y0 - Y0 .* current_strain .* 0.5;  % Assume Poisson's ratio = -0.5;
%    end
%
%    % Plot final grid
%    hold on;
%    plot(X{end}(:), Y{end}(:),'r.')
%
% We computed the position of each grid element at each step, so not we are
% ready to calculate displacements and strains:
%
% .. code:: matlab
%
%   obj.calcStrain('Engineering')
%   
% Done! Now you can check that the final horizontal engineering strain. In
% the final step, strain is 0.1 everywhere, as expected.
%
% .. code:: matlab
% 
%   % Assert equality tolerance for round-off errors
%   assert(all(obj.strain.X{end}(:) - 0.1 < 1e-7))  
%
% The script also includes the plot of all loading steps. At the end, you
% should get something like this:
%
% .. figure:: images/tension_last_step.png
%  
%    Last step of a tension test in the horizontal direction. Note that the 
%    maximum horizontal engineering strain is 0.1, and the vertical one is
%    0.05 (Poisson's ratio = -0.5).
   
clear
close all
clc

% Add path of parent folder to Matlab path
addpath([fileparts(mfilename('fullpath')) '/../']);

obj = StrainCalculator();

% Generate grid
[X0, Y0] = meshgrid(0:10, 0:10);
% Center grid for convenience
X0 = X0 - 5;  
Y0 = Y0 - 5;
% Plot
plot(X0(:), Y0(:),'b.')

% Store grid elements position at each step
X = cell(11,1); 
Y = cell(11,1);

max_step = 10;
% Notice that *step* starts from zero to store the resting step
for step = 0 : max_step
    current_strain = 0.1 * step / max_step;
    X{step + 1} = X0 + X0 .* current_strain;
    Y{step + 1} = Y0 - Y0 .* current_strain .* 0.5;  % Assume Poisson's ratio = -0.5;
end

% Plot final grid
hold on;
plot(X{end}(:), Y{end}(:),'r+')

obj.Xpositions = X;
obj.Ypositions = Y;

obj.calcStrain('Engineering')

% Animate tension test
figure
for k = 1 : obj.nSteps
   subplot(1,2,1)
   cla
   obj.plotGrid(k)
   hold on;

   h = obj.plotStrain('step', k, 'direction', 'X');
   set(h, 'FaceAlpha', 0.5)
   axis([-6 6 -6 6])

   subplot(1,2,2)
   cla
   obj.plotGrid(k)
   hold on;

   h = obj.plotStrain('step', k, 'direction', 'Y');
   set(h, 'FaceAlpha', 0.5)
   axis([-6 6 -6 6])

   
   if k == 1
       subplot(1,2,1)
       title('X direction')
       colorbar(gca)
       subplot(1,2,2)
       title('Y direction')
       colorbar(gca)

   end
   drawnow
   pause(0.3)
   
end

