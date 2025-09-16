clc, clear, close all;

%% Preprocessing

M = 8;

%% Question 1

    %% First solution
    
    lr = 0.001;
    [D, u, losses] = frameDesigner(M, 2, lr); 
    
    figure
    
    quiver(zeros(1, M), zeros(1, M), D(1, :), D(2, :));
    axis equal
    title('Frame')
    
    disp(min(u))


%% Question 2

[D, u, losses] = frameDesigner(M, 3, lr); 

figure

quiver3(zeros(1, M), zeros(1, M), zeros(1, M), D(1, :), D(2, :), D(3, :));
axis equal

disp(min(u))

%%
N = 8;
d = 3;

% Step 1: Generate random matrix and orthogonalize its rows
A = randn(N, d);  % Random N x d matrix
[Q, ~] = qr(A, 0);  % QR decomposition to get orthonormal columns
Q = Q';  % Now Q is d x N

% Step 2: Construct tight frame synthesis operator
% Normalize the columns to make it a Parseval tight frame (A = I)
F = sqrt(d/N) * Q;  % Tight frame of size d x N
frame_vectors = F'; % Each row is a vector in R^d (N x d)

% Step 3: Plot the frame vectors (only for d = 2 or 3)
figure; hold on; axis equal;
if d == 2
    for i = 1:N
        quiver(0, 0, frame_vectors(i,1), frame_vectors(i,2), 0, 'LineWidth', 1.5);
    end
    xlabel('x'); ylabel('y');
    title(['Tight Frame in \mathbb{R}^2 with N = ', num2str(N)]);
elseif d == 3
    for i = 1:N
        quiver3(0, 0, 0, frame_vectors(i,1), frame_vectors(i,2), frame_vectors(i,3), 0, 'LineWidth', 1.5);
    end
    xlabel('x'); ylabel('y'); zlabel('z');
    title(['Tight Frame in \mathbb{R}^3 with N = ', num2str(N)]);
    view(3);
else
    error('Plotting is only supported for d = 2 or 3.');
end

grid on;

%%
clc, clear, close all;

[X, min_distance_history, max_distance_history, mean_distance_history] = sphericalCode(128, 129, 0.02, 0.1);
