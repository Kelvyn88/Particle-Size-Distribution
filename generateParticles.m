function [domainSize area areaFraction time] = generateParticles( figureNumber, size, matrix, maxAttempt)
% Plot particles in a box of size 'size*size'.
% Input:  figure  - new figure for plot number
%         size    - domain size
%         matrix  - particles to fit
%         maxAttempt - maximum attempts
% Output: domainSize  - domain size
%         area   - area used
%         areaFraction   - fraction of area populated
%         time  - tic - toc
%% 
matrix(find(matrix>90)) = 90; % Take particles size only between 1-90
matrix(find(matrix<1)) = 1;

sortParticles = sort(matrix, 'descend'); % sort decreasing

domainSize = size;          % size of domain
tic;                        % measure time until toc
while (1)
    fig_f = figure(figureNumber);
    clf(figureNumber);               % clear figure for next figures
    axis([0 domainSize 0 domainSize]);
    dx = zeros([1 length(sortParticles)]); %X coordinates
    dy = zeros([1 length(sortParticles)]); %Y coordinates
    dr = zeros([1 length(sortParticles)]); %radii
    for i = 1:length(sortParticles)
        attempt=0;
        particleMin = sortParticles(i)/2;
        particleSize = domainSize - sortParticles(i);
        while (1) %% TEST OVERLAP SECTION
            x = particleMin + rand()*particleSize;
            y = particleMin + rand()*particleSize;
            overlap=0;
            for j=1:i-1 %for every existing particle: test overlap
                dist = sqrt((dx(j)-x)^2+(dy(j)-y)^2);
                if (dist<(particleMin+dr(j)))
                    overlap=1;
                    break;
                end
            end
            
            % Professor code, print particle
            if (overlap==0) % create particle
                dx(i) = x; % xparticle 
                dy(i) = y; % yparticle
                dr(i) = particleMin; %radius
                rectangle('position', [x-particleMin, y-particleMin, 2*particleMin, 2*particleMin], ...
                'Curvature', [1,1], ...
                'EdgeColor', 'g', 'FaceColor', 'c');
                % pause; close;
                break;
            else
                attempt = attempt+1;
                if (attempt >= maxAttempt) 
                    break;
                end
            end
        end
        
        % If overlap
        if (overlap == 1)
            break;
        end
    end
    if (overlap == 0)
        break;
    end
end

% Title labels
time = toc;
ptitle = sprintf('Number of particles=%d \t  Size=%dx%d domain', ...
    length(sortParticles), domainSize, domainSize);
title(ptitle);
stime = sprintf('elapsed time=%.2f s', toc);
xlabel(stime);

% Area fraction calc (3.b)
area = 1/4*pi * sum(sortParticles .^2); % calculate total particle surface
areaFraction =area/(domainSize^2);    % Area Fraction calculation
