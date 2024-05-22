% Help of plot strains

clear
close all
clc

% Add path of parent folder to Matlab path
addpath([fileparts(mfilename('fullpath')) '/../']);

obj = StrainCalculator();

nElements = 5;
XinitialPosition = repmat(linspace(0,(nElements-1)*10,nElements), nElements, 1);
YinitialPosition = flipud(XinitialPosition');

Xdisp = [-2,-1,0,1,2;-2,-1,0,1,2;-2,-1,0,1,2;-2,-1,0,1,2;-2,-1,0,1,2];
Ydisp = [-1,-1,-1,-1,-1;-0.5,-0.5,-0.5,-0.5,-0.50;0,0,0,0,0;0.5,0.5,0.5,0.5,0.5;1,1,1,1,1];
obj.Xpositions = {XinitialPosition, XinitialPosition + Xdisp};
obj.Ypositions = {YinitialPosition, YinitialPosition + Ydisp};

nSteps = 2;

obj.calcStrain('Log');
plotDir = 'X';

figure
for k = 2 : obj.nSteps  % First step is empty, start from 2
    cla
    % Plot grid
    obj.plotGrid(k)
    hold on;

    % Plot strain
	h = obj.plotStrain('step', k, 'direction', plotDir);
    set(h, 'FaceAlpha', 0.5)
%     axis([0 1.2 0 1.2])
    
    % Add colorbar
    if k == 2
        c = colorbar(gca);
        c.Label.String = 'Strain';
        if max(obj.strain.limits.(plotDir)) == min(obj.strain.limits.(plotDir))
            caxis([-1, 1] + max(obj.strain.limits.(plotDir)));
        else
            caxis(obj.strain.limits.(plotDir));
        end
    end
    drawnow
    pause(0.3)
   
end


