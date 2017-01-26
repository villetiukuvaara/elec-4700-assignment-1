%%Question 1

%% Question 1
% To calculate the thermal energy, it is noted that there are two degrees
% of freedom for the electrons. Using Maxwell's principle of equipartition
% of energy,
%
% $$\overline{KE}=\frac{1}{2}kT = 2(\frac{1}{2}m\overline{v^2}) \Rightarrow \overline{v^2}=\frac{2kT}{m}$$

vth = sqrt(2*k*T/m);

%%
% Or 18.7 um/s. The mean free path, $l$, is simply

l = vth*0.2e-12;

%%
% Thus the mean free path is about 37.4 nm.

% Set up the simulation

m0 = 9.10938356e-31;
m = 0.26*m0;
T = 300;
k = 1.38064852e-23;


%%
% Set up parameters for the simulation:
width = 100e-9;
length = 200e-9;
population_size = 1000;
plot_population = 10;
time_step = width/vth/100;
iterations = 1000;
show_movie = 1;

% Each row corresponds to an electron with the position and velocities
% [x y vx vy]
state = zeros(population_size, 4);
trajectories = zeros(iterations, plot_population*2);
temperature = zeros(iterations,1);

% Generate an initial population
for i = 1:population_size
    angle = rand*2*pi;
    state(i,:) = [length*rand width*rand vth*cos(angle) vth*sin(angle)];
end

for i = 1:iterations
    %Update positions
    state(:,1:2) = state(:,1:2) + time_step.*state(:,3:4);
    
    j = state(:,1) > length;
    state(j,1) = state(j,1) - length;
    
    j = state(:,1) < 0;
    state(j,1) = state(j,1) + length;
    
    j = state(:,2) > width;
    state(j,2) = 2*width - state(j,2);
    state(j,4) = -state(j,4);
    
    j = state(:,2) < 0;
    state(j,2) = -state(j,2);
    state(j,4) = -state(j,4);
    
    temperature(i) = k/m/population_size*(sum(state(:,3).^2) + sum(state(:,4).^2));
    
    for j=1:plot_population
        trajectories(i, (2*j):(2*j+1)) = state(j, 1:2);
    end 
    
    if show_movie && mod(i,10) == 0
        figure(1);
        subplot(2,1,1);
        hold off;
        plot(state(1:plot_population,1), state(1:plot_population,2), 'o');
        axis([0 length 0 width]);
        title(sprintf('Trajectories for %d of %d Electrons with Fixed Velocity',...
        plot_population, population_size));
        xlabel('x (m)');
        ylabel('y (m)');
        
        %figure(2);
        if i > 1
            subplot(2,1,2);
            hold off;
            plot(time_step*(1:i), temperature(1:i));
            axis([0 time_step*iterations 0 max(temperature)*1.2]);
            title('Semiconductor Temperature');
            xlabel('Time (s)');
            ylabel('Temperature (K)');
        end
        pause(0.1);
    end
end

if ~show_movie
    figure(1);
    title(sprintf('Electron Trajectories for %d of %d Electrons with Fixed Velocity',...
        plot_population, population_size));
    xlabel('x (m)');
    ylabel('y (m)');
    hold on;
    for i=1:plot_population
        plot(trajectories(:,i*2), trajectories(:,i*2+1));
    end
end
% Store the positions and velocities of the electrons