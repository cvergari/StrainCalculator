% This script uses StrainCalculator to compare differnt types of strain definition.
% A fake data series is generated (using :meth:`~StrainCalculator.generateFakeData`) 
% to obtain a uniform deformation field. The deformation field is computed using the different
% methods available in :meth:`~StrainCalculator.calcStrain`. Finally, all loading steps and all
% methods are plotted to be compared.
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
obj.generateFakeData('type', strainType, 'maxStrain', 0.3);

Methods = {'Engineering', 'LogStrain', 'True', 'Almansi', 'Green', 'Infinitesimal'};

failed = [];
for k = 1 : length(Methods)

    obj.calcStrain(Methods{k});
    
    % Strain is uniform, so take a single sample for each step and each
    % method.
    e(1,k) = 0;
    for n = 2 : length(obj.strain.X)
        e(n,k) = obj.strain.X{n}(1);
    end
    
    % Principal strains should be identical
    eII(1,k) = 0;
    for n = 2 : length(obj.strain.X)
        eII(n,k) = obj.strain.PrincipalStrainII{n}(1);
    end
    
    
end

figure
subplot(1,2,1)
plot(e)
xlabel('Step n.')
ylabel('Strain')
legend(Methods);
title('Strains')


subplot(1,2,2)
plot(eII)
xlabel('Step n.')
ylabel('Strain')
legend(Methods);
title('Principal strains')
