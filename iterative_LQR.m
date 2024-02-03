close all; clear all; format compact; beep off;%#ok<CLALL> 

Ts = 0.1;  %サンプリング周期[sec]
t_sim_end = 20;   %シミュレーション時間
tsim = 0:Ts:t_sim_end;   %シミュレーションの基準となる時刻

%コントローラーのパラメーター
iLQR.Ts = Ts; %制御周期
iLQR.tf = 1.0; %予測ホライズンの長さ[sec]
iLQR.N = 10; %予測区間の分割数
iLQR.iter = 10; %繰り返し回数の上限値
iLQR.dt = iLQR.tf/iLQR.N; %予測ホライズンの分割幅
iLQR.torelance = 1; %評価関数のズレの許容範囲

% 評価関数中の重み
iLQR.Q = 100 * [1 0 0; 0 1 0; 0 0 0];
iLQR.R = 1*eye(2);
iLQR.S = 100 * [1 0 0; 0 1 0; 0 0 0];
iLQR.b = 10; %バリア関数の重み
iLQR.e = 0;

%車のパラメーター
car.R = 0.05; %車輪半径[m]
car.T= 0.2; %車輪と車輪の間の幅(車体の横幅)[m]
car.r = car.R/2; %よく使う値を先に計算
car.rT = car.R/car.T;

% 初期条件
car.x = [0;0;0];
iLQR.u = [0;0];

iLQR.len_x = length(car.x);
iLQR.len_u = length(iLQR.u);

iLQR.U = zeros(iLQR.len_u,iLQR.N);    %コントローラに与える初期操作量

%目標地点
iLQR.x_ob = [10,0,0]';

%目標速度
iLQR.u_ob = [0;0];

%拘束条件
iLQR.umax = [15;15];
iLQR.umin = [-15;-15];

%バリア関数用のパラメータ
iLQR.del_bar = 0.1;
iLQR.bar_C = [eye(iLQR.len_u); -eye(iLQR.len_u)];
iLQR.bar_d = [iLQR.umax; -iLQR.umin];
iLQR.con_dim = length(iLQR.bar_d);

%回避関数用のパラメータ
iLQR.xe = [5;0.4;0]; %回避地点
iLQR.xe2 = [5;-0.4;0];
iLQR.xe3 = [5.4;0;0];
iLQR.E = diag([1 1 0]);
iLQR.radius = 0.2;

% シミュレーション(ode45)の設定
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

%臨時追加
X = zeros(iLQR.len_x,iLQR.N+1);
log.X = zeros(iLQR.len_x,iLQR.N+1,401);

norm_ob = norm(car.x-iLQR.xe);


tic;
%シミュレーションループ
for i = 1:length(tsim)

    log.u(:,i) = reshape(iLQR.U,[iLQR.len_u*iLQR.N,1]);
    log.x(:,i) = car.x;
    log.X(:,:,i) = X;
    car.x
    
    %シミュレーション計算（最後だけは計算しない）
    if i ~= length(tsim)
        [t,xi] = ode45( @(t,xi) two_wheel_car(t,xi,iLQR.u,car),...
            [tsim(i) tsim(i+1)], car.x, opts);
        %xi = func_risan(car.x,iLQR.u,iLQR.dt,car);
    else
        break
    end

    %{
    norm_ob_new = norm(car.x-iLQR.xe);

    if norm_ob_new > norm_ob
        iLQR.b = 0;
    end
    
    norm_ob = norm_ob_new;
    %}

    %iLQR法コントローラーを関数として実装
    [U,X] = iLQR_control(iLQR,car);

    iLQR.u = U(:,1);
    iLQR.U = U;
    
    car.x = xi(end,:)';

end
toc;


%
%グラフ化
f=figure();

subplot(5,1,1)
plot(tsim,log.x(1,:),'LineWidth',1); ylabel('x1'); xlabel('Time[s]');
ylim([0,3])
%yticks(0:1:3)
grid on;
set(gca, 'FontSize', 9.5)
subplot(5,1,2);
plot(tsim,log.x(2,:),'LineWidth',1); ylabel('x2'); xlabel('Time[s]');
ylim([0,2])
grid on;
set(gca, 'FontSize', 9.5)
subplot(5,1,3);
plot(tsim,log.x(3,:),'LineWidth',1); ylabel('Φ'); xlabel('Time[s]');
%yticks(0:0.1:0.7)
grid on;
set(gca, 'FontSize', 9.5)
subplot(5,1,4);
plot(tsim,log.u(1,:),'LineWidth',1); ylabel('u1'); xlabel('Time[s]');
%yticks(0:1:11)
grid on;
set(gca, 'FontSize', 9.5)
subplot(5,1,5);
plot(tsim,log.u(2,:),'LineWidth',1); ylabel('u2'); xlabel('Time[s]');
%yticks(0:1:10)
grid on;
set(gca, 'FontSize', 9.5)
%}

%
figure;
plot(log.x(1,:),log.x(2,:),'LineWidth',1); ylabel('y'); xlabel('x');
xlim([-1,11]); ylim([-2,2]);
viscircles(iLQR.xe(1:2)',iLQR.radius,'Color','b');
viscircles(iLQR.xe2(1:2)',iLQR.radius,'Color','b');
viscircles(iLQR.xe3(1:2)',iLQR.radius,'Color','b');
pbaspect([12 4 1])
%


%iLQRコントローラー
function [U,X] = iLQR_control(iLQR,car)
    U = iLQR.U; %前ステップの入力を初期解として使う
    X = Predict(car.x,U,iLQR,car); %状態変数の将来値を入力初期値から予測
    J = CalcJ(X,U,iLQR); %評価関数の初期値を計算
    
    loop = 0;
    while true
        [K,d,dV1,dV2] = Backward(X,U,iLQR,car); %Backward Passの計算
    
        [X,U,J_new] = Forward(X,U,K,d,dV1,dV2,J,iLQR,car); %Forward Passの計算

        loop = loop + 1;

        if abs(J_new-J) <= iLQR.torelance %評価関数値が収束して来たら
            break
        end

        if loop == iLQR.iter %繰り返し回数が限界に来たら
            break 
        end

        J = J_new;
    end
end


%状態変数の初期予測関数
function X = Predict(x,U,iLQR,car)
    xk = x;
    xk = func_risan(xk,U(:,1),iLQR.dt,car);
    X = zeros(iLQR.len_x,iLQR.N+1);
    X(:,1) = xk;
    for i = 1:1:iLQR.N
        xk = func_risan(xk,U(:,i),iLQR.dt,car);
        X(:,i+1) = xk;
    end
end


%評価関数の計算
function J = CalcJ(X,U,iLQR)
    %終端コストの計算
    phi = (1/2)*(X(:,end) - iLQR.x_ob)' * iLQR.S * (X(:,end) - iLQR.x_ob);
    
    %途中のコストの計算
    L = 0;
    for i = 1:1:iLQR.N
        bar_val = barrier(U(:,i),iLQR);
        eva_val = evasion(X(:,i),iLQR.xe,iLQR) + evasion(X(:,i),iLQR.xe2,iLQR) + evasion(X(:,i),iLQR.xe3,iLQR);
        L = L + (1/2)*(X(:,i) - iLQR.x_ob)' * iLQR.Q * (X(:,i) - iLQR.x_ob)...
            + (1/2)*U(:,i)' * iLQR.R * U(:,i)...
            + iLQR.b * bar_val...
            + iLQR.e * eva_val;
    end
    L = L*iLQR.dt; %最後にまとめてdtをかける
    J = phi + L;
end


%Backward Pass
function [K,d,dV1,dV2] = Backward(X,U,iLQR,car)
    sk = (X(:,end)-iLQR.x_ob)'*iLQR.S; %Vxの初期値
    Sk = iLQR.S; %Vxxの初期値
    Qk = iLQR.Q;
    Rk = iLQR.R;

    K = zeros(iLQR.len_u,iLQR.len_x,iLQR.N);
    d = zeros(iLQR.len_u,iLQR.N);
    dV1 = 0; %dが1次の項を入れる変数
    dV2 = 0; %dが2次の項を入れる変数

    for i = iLQR.N:-1:1
        x = X(:,i);
        u = U(:,i);
        Ak = CalcA(x,u,iLQR.dt,car);
        Bk = CalcB(x,iLQR.dt,car);

        dEdx = evasion_dx1(x,iLQR.xe,iLQR) + evasion_dx1(x,iLQR.xe2,iLQR) + evasion_dx1(x,iLQR.xe3,iLQR);
        dEdxx = evasion_hes_x1(x,iLQR.xe,iLQR) + evasion_hes_x1(x,iLQR.xe2,iLQR) + evasion_hes_x1(x,iLQR.xe3,iLQR);
        Qx = ((x-iLQR.x_ob)'*Qk)*iLQR.dt + sk*Ak + iLQR.e*dEdx*iLQR.dt;
        Qxx = Qk*iLQR.dt + Ak'*Sk*Ak + iLQR.e*dEdxx*iLQR.dt;

        dBdu = barrier_du(u,iLQR);
        dBduu = barrier_hes_u(u,iLQR);
        Qu = ((u-iLQR.u_ob)'*Rk)*iLQR.dt + sk*Bk + iLQR.b*dBdu*iLQR.dt;
        Quu = (Rk)*iLQR.dt + Bk'*Sk*Bk + iLQR.b*dBduu*iLQR.dt;

        Qux = Bk'*Sk*Ak;
        
        K_= -inv(Quu)*Qux; %閉ループゲインの計算
        K(:,:,i) = K_;
        d_ = -inv(Quu)*Qu'; %開ループフィードバックの計算
        d(:,i) = d_;
        sk = Qx + d_'*Quu*K_ + Qu*K_ + d_'*Qux; %Vxの更新
        Sk = Qxx + K_'*Quu*K_ + K_'*Qux + Qux'*K_; %Vxxの更新
        dV1 = dV1 + Qu*d_;
        dV2 = dV2 + (1/2)*d_'*Quu*d_;
    end
end


%Forward
function [X,U,J] = Forward(X,U,K,d,dV1,dV2,J,iLQR,car)
    alpha = 1; %直線探索の係数を初期化

    X_ = zeros(iLQR.len_x,iLQR.N+1); %新しいxの値を入れていく変数
    X_(:,1) = X(:,1); %xの初期値は変化しない

    U_ = zeros(iLQR.len_u,iLQR.N); %新しいuの値を入れていく変数

    %直線探索を終わらせるためのカウント
    count = 0;
    %評価関数の最少記録とその周辺
    J_min = J;
    U_min = U;
    X_min = X;
    
    while true
        for i = 1:1:iLQR.N
            U_(:,i) = U(:,i) + K(:,:,i)*(X_(:,i)-X(:,i)) + alpha*d(:,i);
            X_(:,i+1) = func_risan(X_(:,i),U_(:,i),iLQR.dt,car);
        end

        J_new = CalcJ(X_,U_,iLQR);
        dV1_ = alpha*dV1;
        dV2_ = (alpha^2)*dV2;
        z = (J_new-J)/(dV1_+dV2_);

        if 1e-4 <= z && z <= 10 %直線探索が条件を満たしていれば
            J = J_new;
            U = U_;
            X = X_;
            break
        end

        if J_min > J_new %評価関数の最少記録を更新したら
            J_min = J_new;
            U_min = U_;
            X_min = X_;
        end

        if count == 20 %20回やっても直線探索が上手く行かなければ
            J = J_min;
            U = U_min;
            X = X_min;
            break
        end
       
        alpha = (1/2)*alpha;

        count = count + 1;

    end
end


%Akの計算
function Ak = CalcA(x,u,dt,car)
    Ak = eye(3) + ...
        [0 0 -car.r*sin(x(3))*(u(1)+u(2))
        0 0 car.r*cos(x(3))*(u(1)+u(2))
        0 0 0]*dt;
end


%Bkの計算
function Bk = CalcB(x,dt,car)
    cos_ = cos(x(3));
    sin_ = sin(x(3));
    Bk = [car.r*cos_ car.r*cos_
        car.r*sin_ car.r*sin_
        car.rT -car.rT]*dt;
end


%差分駆動型二輪車のモデル(ode用)
function dxi = two_wheel_car(t,xi,u,car)
    dxi = zeros(3,1); %dxiの型を定義
    r = car.r;
    rT = car.rT;
    cos_ = cos(xi(3));
    sin_ = sin(xi(3));
    dxi(1) = r*cos_*u(1) + r*cos_*u(2);
    dxi(2) = r*sin_*u(1) + r*sin_*u(2);
    dxi(3) = rT*u(1) - rT*u(2);
end


%二輪車の離散状態方程式
function xk1 = func_risan(xk,u,dt,car)
    xk1 = zeros(3,1); %xk1の型を定義
    r = car.r;
    rT = car.rT;
    cos_ = cos(xk(3));
    sin_ = sin(xk(3));
    xk1(1) = xk(1) + (r * cos_ * (u(1) + u(2)))*dt;
    xk1(2) = xk(2) + (r * sin_ * (u(1) + u(2)))*dt;
    xk1(3) = xk(3) + (rT * (u(1) - u(2)))*dt;
end


%バリア関数全体
function Bar = barrier(u,control)
    zs = control.bar_d - control.bar_C*u;
    
    Bar = 0;

    %拘束条件の数だけ
    for i = 1:control.con_dim
        Bar = Bar + barrier_z(zs(i),control.del_bar); %値を足していく
    end
end

%バリア関数のu微分（答えは行ベクトルになる）
function dBdu = barrier_du(u,control)
    zs = control.bar_d - control.bar_C*u;

    dBdu = zeros(1,control.len_u);

    for i = 1:control.con_dim
        dBdz = barrier_dz(zs(i),control.del_bar);
        dBdu = dBdu + dBdz*(-control.bar_C(i,:));
    end
end

%バリア関数のuヘッシアン（答えは行列になる）
function dBduu = barrier_hes_u(u,control)
    zs = control.bar_d - control.bar_C*u;

    dBduu = zeros(control.len_u,control.len_u);

    for i = 1:control.con_dim
        dBdzz = barrier_hes_z(zs(i),control.del_bar);
        dBduu = dBduu + control.bar_C(i,:)'*dBdzz*control.bar_C(i,:);
    end
end

%バリア関数（-logか二次関数かを勝手に切り替えてくれる）
function value = barrier_z(z,delta)
    if z > delta
        value = -log(z);
    else
        value = (1/2)*( ((z-2*delta)/delta)^2 -1) - log(delta);
    end
end

%バリア関数の微分値（B(z)のz微分）
function value = barrier_dz(z,delta)
    if z > delta
        value = -(1/z);
    else
        value = (z-2*delta)/delta;
    end
end

%バリア関数のz二階微分（B(z)のz二階微分）
function value = barrier_hes_z(z,delta)
    if z > delta
        value = 1/(z^2);
    else
        value = 1/delta;
    end
end


%
%logバリア関数型回避関数
function eva = evasion(x1,x2,control)
    z = (x1-x2)'*control.E*(x1-x2) - control.radius^2;
    eva = barrier_z(z,control.del_bar);
end

%logバリア関数型回避関数のx1微分
function dBdx1 = evasion_dx1(x1,x2,control)
    z = (x1-x2)'*control.E*(x1-x2) - control.radius^2;
    dzdx1 = 2*(x1-x2)'*control.E;
    dBdz = barrier_dz(z,control.del_bar);
    dBdx1 = dzdx1*dBdz;
end

%logバリア関数型回避関数のx全体微分
function dBdx = evasion_dx(x1,x2,control)
    z = (x1-x2)'*control.E*(x1-x2) - control.radius^2;
    dzdx1 = 2*(x1-x2)'*control.E;
    dzdx2 = 2*(x2-x1)'*control.E;
    dBdz = barrier_dz(z,control.del_bar);
    dBdx = [dzdx1*dBdz dzdx2*dBdz];
end

%logバリア関数型回避関数のヘッシアン（x1のみ）
function dBdx1x1 = evasion_hes_x1(x1,x2,control)
    z = (x1-x2)'*control.E*(x1-x2) - control.radius^2;
    dzdx1 = 2*(x1-x2)'*control.E;
    dzdx1x1 = 2*control.E;
    dBdz = barrier_dz(z,control.del_bar);
    dBdzz = barrier_hes_z(z,control.del_bar);
    dBdx1x1 = dBdz*dzdx1x1 + dzdx1'*dBdzz*dzdx1;
end

%logバリア関数型回避関数のヘッシアン（x全体）
function dBdxx = evasion_hes_x(x1,x2,control)
    z = (x1-x2)'*control.E*(x1-x2) - control.radius^2;
    dzdx1 = 2*(x1-x2)'*control.E;
    dzdx2 = 2*(x2-x1)'*control.E;
    dzdx1x1 = 2*control.E;
    dzdx1x2 = -2*control.E;
    dzdx2x2 = 2*control.E;
    dzdx = [dzdx1 dzdx2];
    dzdxx = [dzdx1x1 dzdx1x2; dzdx1x2 dzdx2x2];
    dBdz = barrier_dz(z,control.del_bar);
    dBdzz = barrier_hes_z(z,control.del_bar);
    dBdxx = dBdz*dzdxx + dzdx'*dBdzz*dzdx;
end
%}


%{
%指数関数型回避関数
function value = evasion(x1,x2,control)
    value = 1/(1-exp(norm(x1-x2)^2-control.radius^2));
end

%回避関数の微分
function value = evasion_dx1(x1,x2,control)
    value = 2*(x1-x2)' * exp(norm(x1-x2)^2-control.radius^2)/... 
        ((1-exp(norm(x1-x2)^2-control.radius^2))^2);
end

%回避関数の微分だが、どっちで微分するかが違う
function value = evasion_dx2(x1,x2,control)
    value = 2*(x2-x1)' * exp(norm(x1-x2)^2-control.radius^2)/... 
        ((1-exp(norm(x1-x2)^2-control.radius^2))^2);
end

%回避関数のx全体での微分
function value = evasion_dx(x1,x2,control)
    val = 2*exp(norm(x1-x2)^2-control.radius^2)/... 
        ((1-exp(norm(x1-x2)^2-control.radius^2))^2);
    value = [(x1-x2)'*val 0 (x2-x1)'*val 0];
end

%回避関数のヘッシアン（近似）
function value = evasion_hes_x1(x1,x2,control)
    val = 2*exp(norm(x1-x2)^2-control.radius^2)/... 
        ((1-exp(norm(x1-x2)^2-control.radius^2))^2);
    value = val*diag([1 1 0]);
end

%回避関数のヘッシアン（近似）
function value = evasion_hes_x(x1,x2,control)
    val = 2*exp(norm(x1-x2)^2-control.radius^2)/... 
        ((1-exp(norm(x1-x2)^2-control.radius^2))^2);
    value = [val*diag([1 1 0]) -val*diag([1 1 0])
        -val*diag([1 1 0]) val*diag([1 1 0])];
end
%}

%{
%二次関数型回避関数
function value = evasion(x1,x2,control)
    value = sigmoid(-norm(x1(1:2)-x2(1:2))^2);
end

function value = evasion_dx1(x1,x2,control)
    value = -2*norm(x1(1:2)-x2(1:2))*sigmoid(-norm(x1(1:2)-x2(1:2))^2)*(1-sigmoid(-norm(x1(1:2)-x2(1:2))^2));
end

function value = evasion_hes_x1(x1,x2,control)
    value = 
end
%}