%% Reaction-limited stochastic biosensor
% Simulates the ideal case where molecule concentration stays constant with time
% All receptors/traps are approachable via all molecules 

%% System specifications
L = 1e-4;               % cm
W = 1e-4;               % cm 

A = L*W;                % cm^2  
l = L/5;               % Length of smaller domain
w = W/5;  
h = l;                  % Need to compute volume
v = l*w*h;              % Volume of domain, for caluclation of density
a = l*w;                % Area of domain

lsites = L/l;
wsites = W/w;

% Introduce mul_fac here

kfm = 5e-7;              % Forward reaction rate s^-1
krm = 0.01;             % Reverse reaction rate s^-1

kfs = 0.1;
krs = 0.1;

N0 = 200;               % #Surface traps/receptors 
Ps = 20000;              % 

%% Simulation setup
time = [1e-6];
for i = -6:1:2
    time = cat(2, time , 2*10^i:10^i:10*10^i);
end

output = time;
t = time(1);

for j = 1:500
    
    %% Trap and target molecule distribution
    rec_dist = zeros(lsites,wsites,2);   % Stores number of unbound(1) and bound(2) #receptors

    lind = ceil(rand(N0,1)*lsites);      % Generate random array along length
    wind = ceil(rand(N0,1)*wsites);      % Generate random array along width

    ind1 = sub2ind(size(rec_dist),lind,wind,ones(N0,1));    

    for i = 1:N0
        rec_dist(ind1(i)) = rec_dist(ind1(i)) + 1;
    end

    %targ_dist = rec_dist(:,:,1);

    forward_cum = kfm * cumsum( reshape( rec_dist(:,:,1)*Ps, numel(targ_dist), 1 ) ) ;
    reverse_cum = krm * cumsum( reshape( rec_dist(:,:,2), numel(rec_dist(:,:,2)), 1 ) );

    % SSA loop

     while t < time(end) 

        kf_total = forward_cum(end);        % kf*(N0-N)*ps
        kr_total = reverse_cum(end);        % kr*N  

        r_total = kf_total + kr_total;

        temp1 = kf_total/r_total;
        temp2 = temp1 + kr_total/r_total;

        r1 = rand();
        r2 = rand();

        tau = 1/r_total*log(1/r2);

        t = t + tau;

        if  r1 <= temp1 

            ind = find(forward_cum > r1*r_total, 1);
            [sub1 sub2] = ind2sub(size(rec_dist(:,:,1)),ind);

            rec_dist(sub1,sub2,1) = rec_dist(sub1,sub2,1) - 1;
            rec_dist(sub1,sub2,2) = rec_dist(sub1,sub2,2) + 1;

            nf = rec_dist(sub1,sub2,1);
            nc = rec_dist(sub1,sub2,2);

            N(time>=t) = N(time>=t) + 1;

            forward_cum(ind:end) = forward_cum(ind:end) - kfm*Ps;
            reverse_cum(ind:end) = reverse_cum(ind:end) + krm;

        elseif r1 <= temp2

            r1 = r1 - temp1;
            ind = find(reverse_cum > r1*r_total, 1);
            [sub1 sub2] = ind2sub(size(rec_dist(:,:,1)),ind);

            rec_dist(sub1,sub2,1) = rec_dist(sub1,sub2,1) + 1;
            rec_dist(sub1,sub2,2) = rec_dist(sub1,sub2,2) - 1;

            nf = rec_dist(sub1,sub2,1);
            nc = rec_dist(sub1,sub2,2);

            N(time>=t) = N(time>=t) - 1;

            forward_cum(ind:end) = forward_cum(ind:end) + kfm*Ps;
            reverse_cum(ind:end) = reverse_cum(ind:end) - krm;

        end


    %      imagesc(rec_dist(:,:,2),[0,5]);
    %      colormap (flipud(autumn));
    %      colorbar;
    %      axis off;
    %      drawnow;

     end    
     
     output = cat(1,output,N);
     t = time(1);
     N = zeros(size(time));
end

ideal_sol = N0 * 1/(1 + krm/(kfm*Ps)) * (1 - exp(-time*(kfm*Ps + krm)));
output = cat(1,output,ideal_sol);
