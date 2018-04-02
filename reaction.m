%% Reaction limited stochastic biosensor
% Simulates the non-ideal target receptor binding 
% 0K: No surface diffusion of target molecules

%% System specifications
L = 1e-4;               % cm
W = 1e-4;               % cm 
Na = 6.022e23;
A = L*W;                % cm^2  

lsites = 10;
wsites = 10;

l = L/lsites;
w = W/wsites;
a = l*w;                % Area of domain

% Introduce mul_fac here
kf = 4e-7;                        % Reference stochastic reaction rates for entire volume of biosensor
kr = 0.002;

%  Steady state value = kf*Ps/(kf*Ps+kr) assuming Ps stays constant

kfm = kf*(lsites*wsites);          % Stochastic Forward reaction rate s^-1 = kf_deterministic/(Na*V)
                                   % Smaller the volume larger the stochastic rate                     
krm = kr;                          % Stochastic Reverse reaction rate s^-1 = kr_deterministic

N0 = 200;                          % #Surface traps/receptors 
Ps = 20000;                        % 

%% Simulation setup
time = [1e-6];
for i = -6:1:2
    time = cat(2, time , 2*10^i:10^i:10*10^i);
end

output = time;


for j = 1:100                       % Remove this for loop if single simulation is required

    %% Trap and target molecule distribution

    rec_dist = zeros(lsites,wsites,2);        % Stores number of unbound(1) and bound(2) #receptors
    targ_dist = zeros(lsites,wsites);         % Stores # target molecules

    lind = ceil(rand(N0,1)*lsites);           % Generate random array along length
    wind = ceil(rand(N0,1)*wsites);           % Generate random array along width

    ind1 = sub2ind(size(rec_dist),lind,wind,ones(N0,1));

    lind = ceil(rand(Ps,1)*lsites);
    wind = ceil(rand(Ps,1)*wsites);

    ind2 = sub2ind(size(targ_dist),lind,wind);

    for i = 1:N0
        rec_dist(ind1(i)) = rec_dist(ind1(i)) + 1;
    end

    for i = 1:Ps
        targ_dist(ind2(i)) = targ_dist(ind2(i)) + 1;
    end

    N = zeros(size(time));
    t = time(1);

    forward_cum = kfm * cumsum( reshape( rec_dist(:,:,1).*targ_dist, numel(targ_dist), 1 ) ) ;
    reverse_cum = krm * cumsum( reshape( rec_dist(:,:,2), numel(rec_dist(:,:,2)), 1 ) );

    %% SSA loop

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
            targ_dist(sub1,sub2) = targ_dist(sub1,sub2) - 1;

            nf = rec_dist(sub1,sub2,1);
            nc = rec_dist(sub1,sub2,2);
            ps = targ_dist(sub1,sub2);

            N(time>=t) = N(time>=t) + 1;

            forward_cum(ind:end) = forward_cum(ind:end) - kfm*(nf + ps + 1);
            reverse_cum(ind:end) = reverse_cum(ind:end) + krm;

        elseif r1 <= temp2

            r1 = r1 - temp1;
            ind = find(reverse_cum > r1*r_total, 1);
            [sub1 sub2] = ind2sub(size(rec_dist(:,:,1)),ind);

            rec_dist(sub1,sub2,1) = rec_dist(sub1,sub2,1) + 1;
            rec_dist(sub1,sub2,2) = rec_dist(sub1,sub2,2) - 1;
            targ_dist(sub1,sub2) = targ_dist(sub1,sub2) + 1;   

            nf = rec_dist(sub1,sub2,1);
            ps = targ_dist(sub1,sub2);
            nc = rec_dist(sub1,sub2,2);

            N(time>=t) = N(time>=t) - 1;

            forward_cum(ind:end) = forward_cum(ind:end) + kfm*(nf + ps - 1);
            reverse_cum(ind:end) = reverse_cum(ind:end) - krm;

        end

     end
     
     output = cat(1,output,N);
end
 
%      imagesc(rec_dist(:,:,2),[0,10]);
%      colormap (flipud(autumn));
%      colorbar;
%      axis off;
%      drawnow;
%     title('Occupation on biosensor surface')
     
      ideal_sol = N0 * kf*Ps/(kf*Ps+kr) * (1 - exp(-time*(kf*Ps+kr)));
      output = cat(1,output,ideal_sol);
      
%      plot(time, ideal_sol);
%      hold on;
%      stairs(time, N);
%      xlabel('Time');
%      ylabel('#Bound receptors');
     