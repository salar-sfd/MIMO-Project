function [X, min_distance_history, max_distance_history, mean_distance_history] = sphericalCode(M, N, dt, damping)
    X_old = randn(M, N)+1j*randn(M, N);
    X_old = X_old./vecnorm(X_old, 2, 1);
    V_old = zeros(M, N);
    F_old = zeros(M, N);
    min_distance_history = [];
    max_distance_history = [];
    mean_distance_history = [];
    for steps=1:500
        F_new = zeros(M, N, N);
        min_distance = inf;
        max_distance = 0;
        mean_distance = 0;
        for i=1:N
            for j=(i+1):N
                xi = X_old(:, i);
                xj = X_old(:, j);

                % Calculate inner product
                distance = norm(xi-xj);
                mean_distance = mean_distance + distance/(N*(N-1)/2);
                if distance<min_distance
                    min_distance = distance;
                elseif distance>max_distance
                    max_distance = distance;
                end
                
                force = (xi-xj)/(((xi-xj)'*(xi-xj))+0.01);
                tforcei = force-(xi'*force)*xi;
                tforcej = -force+(xj'*force)*xj;
                F_new(:, i, j) = tforcei;
                F_new(:, j, i) = tforcej;
            end
        end
        min_distance_history = [min_distance_history, min_distance];
        max_distance_history = [max_distance_history, max_distance];
        mean_distance_history = [mean_distance_history, mean_distance];
        F_new = sum(F_new, 3);
    
        for i=1:N
            xi = X_old(:, i);
            velocity =  (1-damping)*V_old(:, i)+F_old(:,i)*dt;
            tvelocity = velocity-(xi'*velocity)*xi;
            V_new(:, i) = tvelocity;
        end

        for i=1:N
            X_new(:, i) = X_old(:, i)+V_old(:, i)*dt;
        end

        F_old = F_new;
        V_old = V_new;
        X_old = X_new./vecnorm(X_new, 2, 1);
        

        % ----------------
        % Plotting section
        % ----------------
        
        if M==3
            figure(1);
            clf; hold on;
            quiver3(zeros(1,N), zeros(1,N), zeros(1,N), ...
                    X_old(1,:), X_old(2,:), X_old(3,:), 0, 'LineWidth',1.5);
            axis equal;
            xlabel('x'); ylabel('y'); zlabel('z');
            xlim([-1,1]); ylim([-1,1]); zlim([-1,1]);
            title('3D Vectors on Sphere');

            % Add a unit sphere
            r = 1;
            n = 20; % Number of faces for the sphere
            phi = linspace(0,2*pi,n);
            theta = linspace(-pi/2,pi/2,n);
            x = r*cos(theta)'*cos(phi);
            y = r*cos(theta)'*sin(phi);
            z = r*sin(theta)'*ones(1,n);
            surf(x,y,z,'FaceAlpha',0.2,'EdgeColor','none'); % Transparent sphere
            
            axis equal;
            xlabel('x'); ylabel('y'); zlabel('z');
            xlim([-1.2,1.2]); ylim([-1.2,1.2]); zlim([-1.2,1.2]);
            title('3D Vectors on Sphere');
            grid on;
            view(30, 45); % Set a good initial viewing angle
        elseif M==4
            figure(1);
            clf; hold on;
            quiver(zeros(1,N), zeros(1,N), ...
                   X_old(1,:), X_old(2,:), 0, 'LineWidth',1.5);
            axis equal;
            xlabel('x'); ylabel('y');
            xlim([-1,1]); ylim([-1,1]);
            title('2D Vectors on Circle');
        else
            figure(1);
            hold on
            plot(min_distance_history, 'b-');
            plot(max_distance_history, 'g-');
            plot(mean_distance_history, 'r-');
            legend("Minimum Distance", "Max Distance", "Mean Distance");
            xlabel('Time Step');
            ylabel('Distance');
            title('Convergence of Distance');
            grid on;
            drawnow;
        end
        drawnow;
        disp(min_distance)
    end
    X = X_old;
end