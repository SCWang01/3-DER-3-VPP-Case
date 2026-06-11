% This script runs the iteration solving method.
clc;
% Initial guesses.
x_1=0.1;
x_2=0.1;
x_3=0.1;
% Step size and convergence threshold.
epis=0.1;
conver=1e-10;
% Loop control and history.
storage=zeros(1000,3);
afflication_1=zeros(1000,3);
afflication_2=zeros(1000,3);
afflication_3=zeros(1000,3);

% Seed the history to allow the first difference computation.
storage(1,:)=[x_1,x_2,x_3];
for iteration=2:200
    % Solve subproblem 1 with current x_2.
    solution_1=Sub_VPP_1(x_2,x_3, epis);
    x_1=solution_1.x_1;
    afflication_1(iteration,:)=solution_1.y_1_afflication;
    % Solve subproblem 2 with updated x_1.
    solution_2=Sub_VPP_2(x_1,x_3, epis);
    x_2=solution_2.x_2;
    afflication_2(iteration,:)=solution_2.y_2_afflication;
    % Solve subproblem 3 with updated x_2.
    solution_3=Sub_VPP_3(x_1,x_2, epis);
    x_3=solution_3.x_3;
    afflication_3(iteration,:)=solution_3.y_3_afflication;
    % Store iteration results.
    storage(iteration,:)=[x_1,x_2,x_3];
    % L2-norm change between consecutive iterations.
    difference=norm(storage(iteration,:) - storage(iteration-1,:), 2);
    if difference<=conver
        break;
    end
end
