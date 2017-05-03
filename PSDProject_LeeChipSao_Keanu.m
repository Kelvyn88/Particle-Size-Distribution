% Keanu Lee Chip Sao
% Particle Size Distribution
% 4/10/2017

%% (1)
%(a) Convert to Second Moment distribution
FirstMomentPSD= readtable('firstmomentpsd.txt');
MaxSize = FirstMomentPSD.MaxSize;
MinSize = FirstMomentPSD.MinSize;
count= FirstMomentPSD.Count;


width = MaxSize - MinSize %Bin width
d = (MinSize + MaxSize)/2 ; % Median diameter
r = d/2 ; % Radius
n_i=count;
n_all = sum(n_i);
MD_1 = (n_i/n_all)./width;

a = 4*pi*(r.^2); % Surface
a_i = count.*a ; % Partial surface
a_all = sum(a_i); % Total Area
MD_2 = (a_i/a_all)./width; % 2nd moment distribution


%(b) Third Moment Distribution
V = (4/3)*pi*(r.^3);
V_i = count.*V; % Partial Volume
V_all = sum(V_i);
MD_3 = (V_i/V_all)./width;

% (c) Average Diameters
% Count Mean Diameter
n_i= count; % Partial Counts
d_i=d;
N=sum(count);
d_cm=sum(n_i.*d_i)/N


%Area Mean Diameter: d_am
a_i=a_i; % Partial Areas
d_i=d;
A=a_all;
d_am=sum(a_i.*d_i)/A


p1=sum(n_i.*d_i.^3);
p2=sum(n_i.*d_i.^2);
d_am=p1./p2 

%Mass mean diameter: d_mm
p1=sum(n_i.*d_i.^4);
p2=sum(n_i.*d_i.^3);
d_mm=p1./p2 

% % (d) Plot functions
subplot(2,1,1), hold on; 
title('Second Moment Distribution'); xlabel('Bindwidth'); ylabel('Area Mean Diameter');
% %fprintf('CMD %f',CMD);
bar(MD_2,1);
% 
subplot(2,1,2), hold on; 
title('Third Moment Distribution'); xlabel('Bindwidth'); ylabel('Mass Mean Diameter');
bar(MD_3,1);
pause; close;

%% (2)
sigma = sqrt(sum(count .* (d-d_cm).^2)./(N-1))
norm_count=90*MD_1;

s22 = 2* sigma^2;
s2p = sigma * sqrt(2*pi); 
% section b
x = linspace(0,90);
fx = 90 * exp( -(x-d_cm).^2 / s22) / s2p;

fig2 = figure(2);
bar(d,norm_count, 1);
hold on
xlabel(gca, 'Particle size');
ylabel(gca, 'Count');
title('Mean diameter by Count');
plot(x, fx, 'g'); pause; close;
hold off
%Frequncy function
frequency = d/sum(d) % Events/Total
%% (3) 
% Random list of radii
RandomNumbers=randn(1,100)*sigma+d_cm
bar(RandomNumbers);
pause, close;

% Packing of particles (b), (c), (d)
% Initial values
size = 4000; %Domain size
particlesNumber = 100;
retry = 4; % Try 'retry' times until you get it.


fprintf('Particles    Domain Size      Area fraction    Time\n');
for i=1:10 % 100 and double each interval calculation
    attempt = retry; % Maximum overlap allowance
    if (attempt < 100)
        attempt = 100;
    end
    RandomNumbers =randn([1 particlesNumber])*sigma + d_cm; % Use our gaussian distribution
    [domainSize area areaFraction time] = generateParticles(i+2, size, RandomNumbers, attempt);
    pause; close;
    fprintf('  %4.2f\t     %dx%d\t   %4.2f%%\t    %4.2f\n', particlesNumber, domainSize, domainSize, areaFraction*100, time);
    numberp(i) = particlesNumber;
    timep(i) = time;
    areaFractionp(i) = areaFraction;
    retry = 3 * retry;
    particlesNumber = 2*particlesNumber; % Double number of particles
    i = i+1;
    if (domainSize>size) % When bigger than size, end.
        break;
    end
end

%% part 4

figure();
% plot time as a function of number of particles
subplot(2,1,1);
loglog(numberp, timep); set(gca, 'Fontsize', 11); title('T=f(N_p)');
set(gca, 'Fontsize', 10);xlabel('Number of particles');ylabel('Time(s)');

%Time vs Area
subplot(2,1,2);
loglog(areaFractionp, timep);set(gca, 'Fontsize', 11); title('T=f(A)');
set(gca, 'Fontsize', 10); xlabel('Area fraction');ylabel('Time(s)');


