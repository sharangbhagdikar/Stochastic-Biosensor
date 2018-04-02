%% Diffusion limited stochastic biosensor
% Simulates the non-ideal target receptor binding 
% Uniform molecule density assumed initially

%% System specifications
L = 1e-4;               % cm
W = 1e-4;               % cm 
H = 1e-4;
Na = 6.022e23;
A = L*W;                % cm^2  
l = L/10;               % Length of smaller domain
w = W/10;  
h = H/10;               
v = l*w*h;              % Volume of domain, for calculation of density
a = l*w;                % Area of domain

lsites = L/l;
wsites = W/w;
hsites = H/h;
% Introduce mul_fac here

%kf = 4e-7;                        % Reference stochastic reaction rates for entire volume of biosensor
%kr = 0.002;
% Steady state value = kf*Ps/(kf*Ps+kr) assuming Ps stays constant

D = 1e-9;                             % Diffusion constant
kx = D/l^2;
ky = D/w^2;
kz = D/h^2;

%kfm = kf*(lsites*wsites);          % Stochastic Forward reaction rate s^-1 = kf_deterministic/(Na*V)
                                    % Smaller the volume larger the stochastic rate                     
%krm = kr;                          % Stochastic Reverse reaction rate s^-1 = kr_deterministic

N0 = 200;                           % #Surface traps/receptors 
Ps = 2000;                          % 

%% Simulation setup
time = [1e-6];
for i = -6:1:-2
    time = cat(2, time , 2*10^i:10^i:10*10^i);
end

%output = time;
output = [];

%for j = 1:25           % Remove if only single simulation required
%% Trap and target molecule distribution

    targ_dist = zeros(lsites,wsites,hsites);         % Stores # target molecules
    targ_dist(:,:,end) = 4;                          % Add for fixed Ps at boundary
    %targ_dist(:,:,1) = 0;

    % First array adjacent to receptor is assumed to contain zero particles
    % Remove if random distribution not required
%     lind = ceil(rand(Ps,1)*lsites);
%     wind = ceil(rand(Ps,1)*wsites);
%     hind = 1 + ceil(rand(Ps,1)*(hsites-1));
%     ind2 = sub2ind(size(targ_dist),lind,wind,hind);  
%     for i = 1:Ps
%         targ_dist(ind2(i)) = targ_dist(ind2(i)) + 1;
%     end

    N = zeros(size(time));
    t = time(1);
    
    tempmat = targ_dist;
    tempmat(end,:,:) = 0;
    tempmat(:,:,end) = 0; % Remove
    right_cum = kx * reshape ( cumsum ( reshape( tempmat, numel(tempmat), 1 ) ), size(tempmat) );
    
    tempmat = targ_dist;
    tempmat(1,:,:) = 0;
    tempmat(:,:,end) = 0; % Remove
    left_cum = kx * reshape ( cumsum ( reshape( tempmat, numel(tempmat), 1 ) ), size(tempmat) );
    
    tempmat = targ_dist;
    tempmat(:,:,end) = 0;
    up_cum = kz * reshape ( cumsum( reshape( tempmat, numel(tempmat), 1 ) ), size(tempmat) );
    
    tempmat = targ_dist;
    tempmat(:,:,1) = 0;
    down_cum = kz * reshape ( cumsum( reshape( tempmat, numel(tempmat), 1 ) ), size(tempmat) );
    
    tempmat = targ_dist;
    tempmat(:,end,:) = 0;
    tempmat(:,:,end) = 0; % Remove
    front_cum = ky * reshape ( cumsum( reshape( tempmat, numel(tempmat), 1 ) ), size(tempmat) );
        
    tempmat = targ_dist;
    tempmat(:,1,:) = 0;
    tempmat(:,:,end) = 0; % Remove
    back_cum = ky * reshape ( cumsum( reshape( tempmat, numel(tempmat), 1 ) ), size(tempmat) );
    tempmat = [];
    
%% SSA loop

%     while t < time(end) 
         
        r_arr = [right_cum(end) left_cum(end) up_cum(end) down_cum(end) front_cum(end) back_cum(end)];        
        r_total = sum(r_arr);

        temp = cumsum(r_arr)/r_total;
        
        r1 = rand();
        r2 = rand();

        tau = 1/r_total*log(1/r2);

        t = t + tau;

        
%--------------------------------------------------------------------------------------%
%---------------------------------------RIGHT--------------------------------------=---%
%--------------------------------------------------------------------------------------%

        
        if  r1 <= temp(1)
            
            disp('right')
            ind = find(right_cum > r1*r_total, 1);
            
            [sub1 sub2 sub3] = ind2sub(size(targ_dist),ind);
            
            if sub3 ~= hsites       % Remove for free source diffusion
                
                indt = sub2ind(size(targ_dist), sub1+1, sub2, sub3);
                
                targ_dist(ind) = targ_dist(ind) - 1;
                targ_dist(indt) = targ_dist(indt) + 1;
                
                if sub3 ~= hsites
                    up_cum(ind:indt-1) = up_cum(ind:indt-1) - kz;
                end
                
                down_cum(ind:indt-1) = down_cum(ind:indt-1) - kz;
                
                if sub2 ~= wsites
                    front_cum(ind:indt-1) = front_cum(ind:indt-1) - ky;
                end
                if sub2 ~= 1
                    back_cum(ind:indt-1) = back_cum(ind:indt-1) - ky;
                end
                
                if sub1 == 1
                    
                    right_cum(ind:indt-1) = right_cum(ind:indt-1) - kx;
                    left_cum(indt:end) = left_cum(indt:end) + kx;
                    
                elseif sub1 == (lsites - 1)
                    
                    right_cum(ind:end) = right_cum(ind:end) - kx;
                    left_cum(ind:indt-1) = left_cum(ind:indt-1) - kx;
                    
                else
                    right_cum(ind:indt-1) = right_cum(ind:indt-1) - kx;
                    left_cum(ind:indt-1) = left_cum(ind:indt-1) - kx;
                    
                end
            end
           
%--------------------------------------------------------------------------------------%
%----------------------------------------LEFT------------------------------------------%
%--------------------------------------------------------------------------------------%

        elseif r1 <= temp(2)
            disp('left')
            r1 = r1 - temp(1);
            
            ind = find(left_cum > r1*r_total, 1);
            [sub1 sub2 sub3] = ind2sub(size(targ_dist),ind);
            
            if sub3 ~= hsites  % Remove    
            
                indt = sub2ind(size(targ_dist), sub1-1, sub2, sub3);
                
                targ_dist(ind) = targ_dist(ind) - 1;
                targ_dist(indt) = targ_dist(indt) + 1;
                
                if sub3 ~= hsites
                    up_cum(indt:ind-1) = up_cum(indt:ind-1) + kz;
                end
                
                down_cum(indt:ind-1) = down_cum(indt:ind-1) + kz;
                
                if sub2 ~= wsites
                    front_cum(indt:ind-1) = front_cum(indt:ind-1) + ky;
                end
                
                if sub2 ~= 1
                    back_cum(indt:ind-1) = back_cum(indt:ind-1) + ky;
                end
                
                
                if sub1 == 2
                   
                    right_cum(indt:ind-1) = right_cum(indt:ind-1) + kx;
                    left_cum(ind:end) = left_cum(ind:end) - kx;
                    
                elseif sub1 == lsites
                    
                    right_cum(indt:end) = right_cum(indt:end) + kx;
                    left_cum(indt:ind-1) = left_cum(indt:ind-1) + kx;
                    
                else
                    
                    right_cum(indt:ind-1) = right_cum(indt:ind-1) + kx;
                    left_cum(indt:ind-1) = left_cum(indt:ind-1) + kx;
                    
                end
            end
                
%--------------------------------------------------------------------------------------%
%----------------------------------------UP--------------------------------------------%
%--------------------------------------------------------------------------------------%

        elseif r1 < temp(3)
            disp('up')
            r1 = r1 - temp(2);
            
            ind = find(up_cum > r1*r_total, 1);
            [sub1 sub2 sub3] = ind2sub(size(targ_dist),ind);
            
                indt = sub2ind(size(targ_dist), sub1, sub2, sub3+1);
                
                targ_dist(ind) = targ_dist(ind) - 1;
                targ_dist(indt) = targ_dist(indt) + 1;
                
                if sub1 ~= lsites
                    right_cum(ind:indt-1) = right_cum(ind:indt-1) - kx;
                end
                if sub1 ~= 1
                    left_cum(ind:indt-1) = left_cum(ind:indt-1) - kx;
                end
                if sub2 ~= wsites
                    front_cum(ind:indt-1) = front_cum(ind:indt-1) - ky;
                end
                if sub2 ~= 1
                    back_cum(ind:indt-1) = back_cum(ind:indt-1) - ky;
                end
                
                down_cum(ind:indt-1) = down_cum(ind:indt-1) - kz;
                up_cum(ind:indt-1) = up_cum(ind:indt-1) - kz;
                
                if sub3 == 1
                    disp('Something is not right');
                    
                elseif sub3 == hsites-1
                    
                    up_cum(indt:end) = up_cum(indt:end) - kz;
                    
                    down_cum(indt:end) = down_cum(indt:end) - kz;     %Remove
                
                    if sub1 ~= lsites                                     %Remove  
                        right_cum(indt:end) = right_cum(indt:end) - kx;
                    end
                    if sub1 ~= 1                                          %Remove
                        left_cum(indt:end) = left_cum(indt:end) - kx;
                    end
                    if sub2 ~= wsites                                     %Remove
                        front_cum(indt:end) = front_cum(indt:end) - ky;
                    end
                    if sub2 ~= 1                                          %Remove
                        back_cum(indt:end) = back_cum(indt:end) - ky;
                    end

                    targ_dist(indt) = targ_dist(indt) - 1;            %Remove
                    
                                                           
                end
                
            
            
 
%--------------------------------------------------------------------------------------%
%----------------------------------------DOWN------------------------------------------%
%--------------------------------------------------------------------------------------%
            
            
        elseif r1 <= temp(4)
            disp('down')
            r1 = r1 - temp(3);
            
            ind = find(down_cum > r1*r_total, 1);
            [sub1 sub2 sub3] = ind2sub(size(targ_dist),ind);
                        
            indt = sub2ind(size(targ_dist), sub1, sub2, sub3-1);

            targ_dist(ind) = targ_dist(ind) - 1;

            if sub3 ~= 2
                targ_dist(indt) = targ_dist(indt) + 1;
            end

            up_cum(indt:ind-1) = up_cum(indt:ind-1) + kz;
            down_cum(indt:ind-1) = down_cum(indt:ind-1) + kz;

            if sub1 ~= lsites
                right_cum(indt:ind-1) = right_cum(indt:ind-1) + kx;
            end
            if sub1 ~= 1
                left_cum(indt:ind-1) = left_cum(indt:ind-1) + kx;
            end
            if sub2 ~= wsites
                front_cum(indt:ind-1) = front_cum(indt:ind-1) + ky;
            end
            if sub2 ~= 1
                back_cum(indt:ind-1) = back_cum(indt:ind-1) + ky;
            end

            
            if sub3 == 2
                
                up_cum(indt:end) = up_cum(indt:end) - kz;
                down_cum(indt:end) = down_cum(indt:end) - kz;

                if sub1 ~= lsites  
                    right_cum(indt:end) = right_cum(indt:end) - kx;
                end
                if sub1 ~= 1 
                    left_cum(indt:end) = left_cum(indt:end) - kx;
                end
                if sub2 ~= wsites
                    front_cum(indt:end) = front_cum(indt:end) - ky;
                end
                if sub2 ~= 1
                    back_cum(indt:end) = back_cum(indt:end) - ky;
                end


            elseif sub3 == hsites

                up_cum(ind:end) = up_cum(ind:end) + kz;
                down_cum(ind:end) = down_cum(ind:end) + kz;        %Remove
                
                targ_dist(ind) = targ_dist(ind) + 1;                %Remove
                
                if sub1 ~= lsites  
                    right_cum(ind:end) = right_cum(ind:end) + kx; %Remove
                end
                if sub1 ~= 1 
                    left_cum(ind:end) = left_cum(ind:end) + kx;   %Remove
                end
                if sub2 ~= wsites
                    front_cum(ind:end) = front_cum(ind:end) + ky; %Remove
                end
                if sub2 ~= 1
                    back_cum(ind:end) = back_cum(ind:end) + ky;   %Remove
                end

            end

                
                
            
%--------------------------------------------------------------------------------------%
%----------------------------------------FRONT-----------------------------------------%
%--------------------------------------------------------------------------------------%
            
        elseif r1 < temp(5)
            disp('front')
            r1 = r1 - temp(4);
            
            ind = find(front_cum > r1*r_total, 1);
            [sub1 sub2 sub3] = ind2sub(size(targ_dist),ind);
            
            if sub3 ~= hsites %Remove
                
                indt = sub2ind(size(targ_dist), sub1, sub2+1, sub3);
                
                targ_dist(ind) = targ_dist(ind) - 1;
                targ_dist(indt) = targ_dist(indt) + 1;
                
                if sub1 ~= lsites
                    right_cum(ind:indt-1) = right_cum(ind:indt-1) - kx;
                end
                if sub1 ~= 1
                    left_cum(ind:indt-1) = left_cum(ind:indt-1) - kx;
                end
                if sub3 ~= hsites
                    up_cum(ind:indt-1) = up_cum(ind:indt-1) - kz;
                end
                
                down_cum(ind:indt-1) = down_cum(ind:indt-1) - kz;
                
                if sub2 == 1
                   
                    front_cum(ind:indt-1) = front_cum(ind:indt-1) - ky;
                    back_cum(indt:end) = back_cum(indt:end) + ky;
                    
                elseif sub2 == wsites-1
                    
                    front_cum(ind:end) = front_cum(ind:end) - ky;
                    back_cum(ind:indt-1) = back_cum(ind:indt-1) - ky;
                    
                else
                    
                    front_cum(ind:indt-1) = front_cum(ind:indt-1) - ky;
                    back_cum(ind:indt-1) = back_cum(ind:indt-1) - ky;
                    
                end
                
            end

%--------------------------------------------------------------------------------------%
%----------------------------------------BACK------------------------------------------%
%--------------------------------------------------------------------------------------%
            
        elseif r1 <= temp(6)
            disp('back')
            r1 = r1 - temp(5);
            
            ind = find(back_cum > r1*r_total, 1);
            [sub1 sub2 sub3] = ind2sub(size(targ_dist),ind);
            
            if sub3 ~= hsites %Remove

                indt = sub2ind(size(targ_dist), sub1, sub2-1, sub3);
                
                targ_dist(ind) = targ_dist(ind) - 1;
                targ_dist(indt) = targ_dist(indt) + 1;
                
                if sub3 ~= hsites
                    up_cum(indt:ind-1) = up_cum(indt:ind-1) + kz;
                end
                
                down_cum(indt:ind-1) = down_cum(indt:ind-1) + kz;
                
                if sub1 ~= lsites
                    right_cum(indt:ind-1) = right_cum(indt:ind-1) + kx;
                end
                
                if sub1 ~= 1
                    left_cum(indt:ind-1) = left_cum(indt:ind-1) + kx;
                end
                
                if sub2 == 2
                   
                    front_cum(indt:ind-1) = front_cum(indt:ind-1) + ky;
                    back_cum(ind:end) = back_cum(ind:end) - ky;
                    
                elseif sub2 == wsites
                    
                    front_cum(indt:end) = front_cum(indt:end) + ky;
                    back_cum(indt:ind-1) = back_cum(indt:ind-1) + ky;
                    
                else
                    
                    front_cum(indt:ind-1) = front_cum(indt:ind-1) + ky;
                    back_cum(indt:ind-1) = back_cum(indt:ind-1) + ky;
                    
                end
            end
                
        end
        
%     end
     
     ran = reshape(sum(sum(targ_dist,1),2),hsites,1);   
     output = cat(2,output,ran);
%end

% plot(1:hsites,output);
% hold on
% plot(1:hsites,mean(output,2),'r','linewidth',3);
% secr = cat(2,secr,mean(output,2));



%      imagesc(rec_dist(:,:,2),[0,10]);
%      colormap (flipud(autumn));
%      colorbar;
%      axis off;
%      drawnow;
%     title('Occupation on biosensor surface')
     
%      ideal_sol = N0 * kf*Ps/(kf*Ps+kr) * (1 - exp(-time*(kf*Ps+kr)));
%      output = cat(1,output,ideal_sol);
      
%      plot(time, ideal_sol);
%      hold on;
%      stairs(time, N);
%      xlabel('Time');
%      ylabel('#Bound receptors');
     