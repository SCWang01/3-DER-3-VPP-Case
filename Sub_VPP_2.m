function solution_1=Sub_VPP_2(x_1,x_3,eps)
    switching_1=0;
    switching_2=30;
    switching_3=60;
    %DER-1+DER-2:100 ,DER-1+DER-3:120,DER-2+DER-3:80 
    solution_1=struct();
    F=[];
    x_2=sdpvar(1,1);
    y_1=sdpvar(3,1);
    y_2=sdpvar(3,1);
    y_3=sdpvar(3,1);
    M=1e4;
    K_1=binvar(3,1);
    K_2=binvar(3,1);
    K_3=binvar(3,1);
    mu_1=sdpvar(3,1);
    mu_2=sdpvar(3,1);
    mu_3=sdpvar(3,1);
    lambda_1=sdpvar(1,1);
    lambda_2=sdpvar(1,1);
    lambda_3=sdpvar(1,1);
    EPR_1=sdpvar(3,1);
    EPR_2=sdpvar(3,1);
    EPR_3=sdpvar(3,1);
    % Polytope constraints
    F=[F 0<=x_2<=1];
    % Transferred KKT Conditions 
    %DER1
    F=[F EPR_1(1)==0];
    F=[F EPR_1(2)==120*(1-x_2)-switching_1];
    F=[F EPR_1(3)==240*(1-x_3)-switching_1];
    for i=1:3
        F=[F y_1(i)<=K_1(i)*M];
        F=[F mu_1(i)<= (1-K_1(i))*M];
    end
    F=[F eps*y_1-EPR_1+lambda_1*ones(3,1)-mu_1==0];
    F=[F sum(y_1)==1];
    F=[F mu_1>=0];
    F=[F y_1 >=0];

    %DER2
    EPR_2=sdpvar(3,1);
    F=[F EPR_2(1)==120*(1-x_1)-switching_2];
    F=[F EPR_2(2)==0];
    F=[F EPR_2(3)==72*(1-x_3)-switching_2];
    for i=1:3
        F=[F y_2(i)<=K_2(i)*M];
        F=[F mu_2(i)<= (1-K_2(i))*M];
    end
    F=[F eps*y_2-EPR_2+lambda_2*ones(3,1)-mu_2==0];
    F=[F sum(y_2)==1];
    F=[F mu_2>=0];
    F=[F y_2 >=0];

    %DER3
    EPR_3=sdpvar(3,1);
    F=[F EPR_3(1)==240*(1-x_1)-switching_3];
    F=[F EPR_3(2)==72*(1-x_2)-switching_3];
    F=[F EPR_3(3)==0];
    for i=1:3
        F=[F y_3(i)<=K_3(i)*M];
        F=[F mu_3(i)<= (1-K_3(i))*M];
    end
    F=[F eps*y_3-EPR_3+lambda_3*ones(3,1)-mu_3==0];
    F=[F sum(y_3)==1];
    F=[F mu_3>=0];
    F=[F y_3 >=0];



    obj=x_2*(240*y_1(2)+144*y_3(2))+switching_2*(1-y_2(2));

    %optimization settings
    ops=sdpsettings('solver','gurobi');
    %obj=-switching_1*y_1;
    optimize(F,-obj,ops);  
    solution_1.x_2=double(x_2);
    solution_1.y_1=double(y_1);
    solution_1.y_2=double(y_2);
    solution_1.y_3=double(y_3);
    solution_1.y_2_afflication=[y_1(2),y_2(2),y_3(2)];
end
