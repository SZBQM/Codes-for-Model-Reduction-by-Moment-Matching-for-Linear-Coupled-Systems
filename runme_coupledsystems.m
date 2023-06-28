clear all;
clc;


A10=xlsread('matrix_data.xlsx','Sheet1');
A20=xlsread('matrix_data.xlsx','Sheet2');
A30=xlsread('matrix_data.xlsx','Sheet3');
B1=xlsread('matrix_data.xlsx','Sheet4');
B2=xlsread('matrix_data.xlsx','Sheet5');
B3=xlsread('matrix_data.xlsx','Sheet6');
C1=xlsread('matrix_data.xlsx','Sheet7');
C2=xlsread('matrix_data.xlsx','Sheet8');
C3=xlsread('matrix_data.xlsx','Sheet9');

%  -------------------------------------------------------x1-6 order
D10 = 0;

K11 = 2;
A1 = A10+B1*K11*C1;
[num10,den10] = ss2tf(A1,B1,C1,D10);
M10 = tf(num10,den10);%----------------------------------tf of x1

eig_A1 = eig(A1);
% for i=1:length(eig_A1)
%     if eig_A1(i)>=0
%         disp('x1不稳定');
%     end
% end

% CONT1 = ctrb(A1,B1);
% if rank(CONT1) ~= size(A1)
%     disp('x1不可控');
% else
%     ;
% end
% OBSER1 = obsv(A1,C1);
% if rank(OBSER1) ~= size(A1)
%     disp('x1不可观');
% else
%     ;
% end
 
%  -------------------------------------------------------x2-8 order 
D20 = 0;

K22 = 5;
A2 = A20+B2*K22*C2;
[num20,den20] = ss2tf(A2,B2,C2,D20);
% [z20,p20,k20] = ss2zp(A2,B2,C2,D20);
M20 = tf(num20,den20);%----------------------------------tf of x2

eig_A2 = eig(A2);
% for i=1:length(eig_A2)
%     if eig_A2(i)>=0
%         disp('x2不稳定');
%     end
% end

% CONT2 = ctrb(A2,B2);
% if rank(CONT2) ~= size(A2)
%     disp('x2不可控');
% else
%     ;
% end
% OBSER2 = obsv(A2,C2);
% if rank(OBSER2) ~= size(A2)
%     disp('x2不可观');
% else
%     ;
% end

%  -------------------------------------------------------x3-30 order
D30 = 0;

K33 = 1;
A3 = A30+B3*K33*C3;
[num30,den30] = ss2tf(A3,B3,C3,D30);
M30 = tf(num30,den30);%---------------------------------- tf of x3

eig_A3 = eig(A3);
% for i=1:length(eig_A3)
%     if eig_A3(i)>=0
%         disp('x3不稳定');
%     end
% end
% 
% CONT3 = ctrb(A3,B3);
% if rank(CONT3) ~= size(A3)
%     disp('x3不可控');
% else
%     ;
% end
% OBSER3 = obsv(A3,C3);
% if rank(OBSER3) ~= size(A3)
%     disp('x3不可观');
% else
%     ;
% end

% -------------------------------------------------------S and L 
S = [3 1;0 4];
eig_S = eig(S);
L = [1,0];
OBSER1 = obsv(S,L);
if rank(OBSER1) ~= size(S)
    disp('L1 and S1 error');
else
    ;
end
% -------------------------------------------------------input
w0 = [1;1];
w = [];
for t = 0:0.01:1;
    w = [w expm(t*S)*w0];
end
figure;
subplot(2,2,1);
plot(w(1,:),w(2,:),'-');
title('Input');
grid on
hold on
ceshi_1 = expm(t*S)*w0;
u = L*expm(t*S)*w0;

K12 = 2;
K13 = 1;
G1 = 1;
G2 = 3;
G3 = 1;
R1 = 1;
R2 = 1;
R3 = 2;
% -------------------------------------------------------coupled relation
K = [2,2,1;0,5,0;0,0,1];
G = [1;3;1];
R = [1,1,2];
Ac = blkdiag(A1,A2,A3);
Bc = blkdiag(B1,B2,B3);
Cc = blkdiag(C1,C2,C3);
A = Ac+Bc*K*Cc;
B = Bc*G;
C = R*Cc;
D = zeros(1,1);
eig_A = eig(A);
[z,p,k] = ss2zp(A,B,C,D);
[num,den] = ss2tf(A,B,C,D);
M = tf(num,den);%----------------------------------tf of OS
% assignin('base','M.Variable',3)

% for i = 1:length(p)
%     if p(i) >=0
%         disp('系统不稳定');
%     else
%         ;
%     end
% end
% -----------------------------------------------------------------------OS
% figure;
[n m] = size(A);
x0 = zeros(n,1);
tic
[t,x] = ode15s(@(t,x)A*x+B*u,0:0.01:1,x0);
toc
y = C*x';
% plot(t,y,'-k');
% hold on;

% ---------------------------------------------------solution of the sylvester equation

L2 = G2*L;
% u2 = L2*expm(t*S1)*w0;
X2 = sylvester(A2,-S,-B2*L2);

L3 = G3*L;
X3 = sylvester(A3,-S,-B3*L3);

L1 = K12*C2*X2+K13*C3*X3+G1*L;
X1 = sylvester(A1,-S,-B1*L1);


% % ---------------------------------------------------------------------先对子系统降阶，再耦合
% 
% Y1 = sylvester(A10,-S,-B1*L);
% [Q1 T1] = qr(Y1,0);
% Ar1 = Q1'*A10*Q1;
% Br1 = Q1'*B1;
% Cr1 = C1*Q1;
% Dr1 = D10;
% 
% Y2 = sylvester(A20,-S,-B2*L);
% [Q2 T2] = qr(Y2,0);
% Ar2 = Q2'*A20*Q2;
% Br2 = Q2'*B2;
% Cr2 = C2*Q2;
% Dr2 = D20;
% 
% Y3 = sylvester(A30,-S,-B3*L);
% [Q3 T3] = qr(Y3,0);
% Ar3 = Q3'*A30*Q3;
% Br3 = Q3'*B3;
% Cr3 = C3*Q3;
% Dr3 = D30;
% 
% Ac_r = blkdiag(Ar1,Ar2,Ar3);
% Bc_r = blkdiag(Br1,Br2,Br3);
% Cc_r = blkdiag(Cr1,Cr2,Cr3);
% Ar_r = Ac_r+Bc_r*K*Cc_r;
% Br_r = Bc_r*G;
% Cr_r = R*Cc_r;
% Dr_r = zeros(1,1);
% 
% [nr_2 mr_2] = size(Ac_r);
% xcr0 = zeros(nr_2,1);
% tic
% [tr_2,xr_2] = ode15s(@(tr_2,xr_2)Ar_r*xr_2+Br_r*u,0:0.01:1,xcr0);
% toc
% yr_r = Cr_r*xr_2';
% 
% % ----------------------------------------------------------------------先对子系统降阶，再耦合



% -----------------------------------------------------------------------RS
X = blkdiag(X1,X2,X3);
[Q W] = qr(X,0);
Ar = Q'*A*Q;
Br = Q'*B;
Cr = C*Q;
Dr = D;
eig_Ar = eig(Ar);
[num_r,den_r] = ss2tf(Ar,Br,Cr,Dr);
Mr_1 = tf(num_r,den_r);
[zr,pr,kr]=ss2zp(Ar,Br,Cr,Dr);
Mr_2=zpk(zr,pr,kr);%----------------------------------tf of RS
% for i=1:length(eig_Ar)
%     if eig_Ar(i)>=0
%         disp('RS不稳定');
%     end
% end
% 
% CONT_r = ctrb(Ar,Br);
% if rank(CONT_r) ~= size(Ar)
%     disp('RS不可控');
% else
%     ;
% end
% OBSER_r = obsv(Ar,Cr);
% if rank(OBSER_r) ~= size(Ar)
%     disp('RS不可观');
% else
%     ;
% end


[nr mr] = size(Ar);
xr0 = zeros(nr,1);
tic
[tr,xr] = ode15s(@(tr,xr)Ar*xr+Br*u,0:0.01:1,xr0);
toc
yr = Cr*xr';
% figure;
subplot(2,2,2);
plot(t,y,'-k',tr,yr,'-+r');
legend('y','yr');
title('Output of OS and RS');
hold on;
%------------------------------------------------------------Relative error
% figure;
subplot(2,2,3);
for j = 1:101;
    Relative_error(j) = abs(y(j)-yr(j))/abs(y(j));
    
end
semilogy(t,Relative_error,'--k');
hold on;
% -----------------------------------------------------------Bode plots of OS and RS
[num,den]=ss2tf(A,B,C,D);
TF=tf(num,den);
[num_r,den_r]=ss2tf(Ar,Br,Cr,Dr);
TFr=tf(num_r,den_r);
% [num_r2,den_r2]=ss2tf(Ar_r,Br_r,Cr_r,Dr_r);
% TFr2=tf(num_r2,den_r2);

[mag,phase,w]=bode(TF);
[mag_r,phase_r,w_r]=bode(TFr);
% [mag_r2,phase_r2,w_r2]=bode(TFr2);


mag_data=squeeze(mag(1,1,:));
magr_data=squeeze(mag_r(1,1,:));
phase_data=squeeze(phase(1,1,:));
phaser_data=squeeze(phase_r(1,1,:));

xlswrite('D:\Program Files\Polyspace\R2021a\bin\ChugaoData\data_bode.xlsx',w,1);
xlswrite('D:\Program Files\Polyspace\R2021a\bin\ChugaoData\data_bode.xlsx',mag_data,1,'B');
xlswrite('D:\Program Files\Polyspace\R2021a\bin\ChugaoData\data_bode.xlsx',phase_data,1,'C');

xlswrite('D:\Program Files\Polyspace\R2021a\bin\ChugaoData\data_bode.xlsx',w_r,1,'D');
xlswrite('D:\Program Files\Polyspace\R2021a\bin\ChugaoData\data_bode.xlsx',magr_data,1,'E');
xlswrite('D:\Program Files\Polyspace\R2021a\bin\ChugaoData\data_bode.xlsx',phaser_data,1,'F');



% figure;
subplot(2,2,4);
margin(TF);
grid on;
hold on;
margin(TFr,'r');
grid on;
hold on;
% margin(TFr2,'b');
% grid on;
% hold on;
legend('OS','RS');

%------------------------------------------------------------提取bode图数据
% figure;
% margin(TFr,'r');
% grid on;
% hold on;
% lh=findall(gca,'type','line'); % 从当前图(gca)中取出曲线的handle,
% xc=get(lh,'xdata'); % 取出x轴数据，注意，这个x和y是以cell的数据结构保存的
% yc=get(lh,'ydata'); % 取出y轴数据x=xc{1};从cell转换成矩阵，可以这样写y=yc{1};
% 
% frequency_r = xc{3};
% mag_r = yc{3};

% ---------------------------------------------------------the memoents of OS and RS
eig_S = [eig_S,eig_S,eig_S];
for i = 1:2
    for l = 1:nr
        moments_os(i,l) = [C*((inv(eig_S(l)*eye(n)-A))^(i))*B];
        moments_rs(i,l) = [Cr*((inv(eig_S(l)*eye(nr)-Ar))^(i))*Br];
%         moments_rs2(i,l) = [Cr_r*((inv(eig_S(l)*eye(nr)-Ar_r))^(i))*Br_r];
    end
    
end
 disp('The First and Second Moments of OS and RS')
 moments_os
 moments_rs
