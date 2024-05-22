% Help of plot strains

clear
close all
clc

% Add path of parent folder to Matlab path
addpath([fileparts(mfilename('fullpath')) '/../']);

obj = StrainCalculator();
% Choose a strain type
strainType = 'HorizontalTension';
% strainType = 'HorizontalCompression';
% strainType = 'VerticalTension';
% strainType = 'VerticalCompression';
% strainType = 'HorizontalShear';
% strainType = 'PureShear';


obj.generateFakeData(strainType);

obj.calcStrain('Green');

figure
for k = 2 : obj.nSteps  % First step is empty, start from 2
    cla
    % Plot grid
    obj.plotGrid(k)
    hold on;
    % Choose interesting direction to plot
    if contains(strainType, 'Shear')
       plotDirection = 'XY';
    elseif contains(strainType, 'Horizontal')
       plotDirection = 'X';
    else
       plotDirection = 'Y';
    end
    
    % Plot strain
	h = obj.plotStrain('step', k, 'direction', plotDirection);
    set(h, 'FaceAlpha', 0.5)
    axis([0 1.2 0 1.2])
    title(['Step n. ' num2str(k) '/' num2str(obj.nSteps)])
    
    % Add colorbar
    if k == 2
        c = colorbar(gca);
        c.Label.String = 'Strain';
        caxis(obj.strain.limits.(plotDirection))
    end
    drawnow
    pause(0.3)
   
end


