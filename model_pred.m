clear variables


t=1;

x=0;
y=2;
z=2;

X=[x;y;z];

v_x=5;
v_y=-0.1;
v_z=2;
a_x=0;
a_y=0;
a_z=0;
wx=1;
wy=5;
wz=2;
g=9.81;
sampling_time=0.001;
%sample_time=0;

k_d=0.1;
k_m=0.1;

V=[v_x;v_y;v_z];

p=0;

while(t<2)

for i=1:2000
    mod_v=sqrt(v_x^2+v_y^2+v_z^2);
    mod_a=sqrt(a_x^2+a_y^2+a_z^2);
    mod_w=sqrt(wx^2+wy^2+wz^2);

    H=[-k_d*mod_v,-k_m*wz, k_m*wy; k_m*wz, -k_d*mod_v, -k_m*wx; -k_m*wy, k_m*wx, -k_d*mod_v];
    G=[0;0;-g];
    A(:,p+i+1)=H*V(:,p+i)+G;
    A=horzcat(A,[a_x;a_y;a_z]);
    a_x=A(1,p+i);
    a_y=A(2,p+i);
    a_z=A(3,p+i);
    V(:,p+i+1)= V(:,p+i)+A(:,p+i+1)*sampling_time;
    v_x=V(1,p+i);
    v_y=V(2,p+i);
    v_z=V(3,p+i);
    V=horzcat(V,[v_x;v_y;v_z]);
    X=horzcat(X,[x;y;z]);
    X(:,p+i+1)=X(:,p+i)+V(:,p+i+1)*sampling_time;
    if(X(3,p+i)<0)
        disp(p+i);
        break;
    end
end

i=p+i;

V=V(:,1:i+1);
A=A(:,1:i+1);
X=X(:,1:i+1);

ep=0.72;
mu=0.4;
r=0.002;

mod_V_T=sqrt((V(1,i+1)-r*wy)^2+(V(2,i+1)+r*wx)^2);

alpha=mu*(1+ep)*V(3,i+1)/mod_V_T;

B_t=[0 alpha*r 0;-alpha*r 0 0;0 0 0];
A_t=[1-alpha 0 0;0 1-alpha 0;0 0 -ep];

V(:,i+2)=A_t*V(:,i+1)+B_t*[wx;wy;wz];
X(:,i+2)=X(:,i+1);

disp(X(:,i+1));

v_x=V(1,i+2);
v_y=V(2,i+2);
v_z=V(3,i+2);
a_x=0;
a_y=0;
a_z=0;
j=i+2;

for i=0:2000
    mod_v=sqrt(v_x^2+v_y^2+v_z^2);
    mod_a=sqrt(a_x^2+a_y^2+a_z^2);
    mod_w=sqrt(wx^2+wy^2+wz^2);

    H=[-k_d*mod_v,-k_m*wz, k_m*wy; k_m*wz, -k_d*mod_v, -k_m*wx; -k_m*wy, k_m*wx, -k_d*mod_v];
    G=[0;0;-g];
    A(:,j+i+1)=H*V(:,j+i)+G;
    A=horzcat(A,[a_x;a_y;a_z]);
    a_x=A(1,j+i);
    a_y=A(2,j+i);
    a_z=A(3,j+i);
    V(:,j+i+1)= V(:,j+i)+A(:,j+i+1)*sampling_time;
    v_x=V(1,j+i);
    v_y=V(2,j+i);
    v_z=V(3,j+i);
    V=horzcat(V,[v_x;v_y;v_z]);
    X=horzcat(X,[x;y;z]);
    X(:,j+i+1)=X(:,j+i)+V(:,j+i+1)*sampling_time;

        if(X(1,j+i)>6)
            disp(j+i);
            break;
        end
end

ep=0.52;
mu=0.6;

a=1;
b=1;
c=2;

%orientation of racket
n=(1/sqrt(a^2+b^2+c^2))*[a;b;c];


R=n/[V(1,j+i+1);V(2,j+i+1);V(3,j+i+1)];

v_r=V(:,j+i+1)'*n;
mod_v_t=sqrt((V(1,j+i+1)-r*wy)^2+(V(2,j+i+1)+r*wx)^2);

k=mu*(1+ep)*v_r/mod_v_t;

A_r=[1-k 0 0;0 1-k 0;0 0 -ep];
B_r=[0 k*r 0;-k*r 0 0;0 0 0];

A_r_1= R*A_r*R';

V_r=[-5;-0.5;4];
mod_V_r=sqrt(V_r(1)^2+V_r(2)^2+V_r(3)^2);

V(:,j+i+2)=(eye(3)-A_r_1)*V_r + A_r_1*V(:,j+i+1)+R*B_r*[wx;wy;wz];

X(:,j+i+2)=X(:,j+i+1);

disp(X(:,j+i+2));


v_x=V(1,j+i+2);
v_y=V(2,j+i+2);
v_z=V(3,j+i+2);

a_x=0;
a_y=0;
a_z=0;
k1=j+i+2;


for i=0:2000
    mod_v=sqrt(v_x^2+v_y^2+v_z^2);
    mod_a=sqrt(a_x^2+a_y^2+a_z^2);
    mod_w=sqrt(wx^2+wy^2+wz^2);

    H=[-k_d*mod_v,-k_m*wz, k_m*wy; k_m*wz, -k_d*mod_v, -k_m*wx; -k_m*wy, k_m*wx, -k_d*mod_v];
    G=[0;0;-g];
    A(:,k1+i+1)=H*V(:,k1+i)+G;
    A=horzcat(A,[a_x;a_y;a_z]);
    a_x=A(1,k1+i);
    a_y=A(2,k1+i);
    a_z=A(3,k1+i);
    V(:,k1+i+1)= V(:,k1+i)+A(:,k1+i+1)*sampling_time;
    v_x=V(1,k1+i);
    v_y=V(2,k1+i);
    v_z=V(3,k1+i);
    V=horzcat(V,[v_x;v_y;v_z]);
    X=horzcat(X,[x;y;z]);
    X(:,k1+i+1)=X(:,k1+i)+V(:,k1+i+1)*sampling_time;
        
        if(X(1,k1+i)<0 || X(3,k1+i)<0)
            disp(k1+i);
            break;
        end
end


i=k1+i;

V=V(:,1:i+1);
A=A(:,1:i+1);
X=X(:,1:i+1);

ep=0.72;
mu=0.4;
r=0.002;

mod_V_T=sqrt((V(1,i+1)-r*wy)^2+(V(2,i+1)+r*wx)^2);

alpha=mu*(1+ep)*V(3,i+1)/mod_V_T;

B_t=[0 alpha*r 0;-alpha*r 0 0;0 0 0];
A_t=[1-alpha 0 0;0 1-alpha 0;0 0 -ep];

V(:,i+2)=A_t*V(:,i+1)+B_t*[wx;wy;wz];
X(:,i+2)=X(:,i+1);

disp(X(:,i+1));

v_x=V(1,i+2);
v_y=V(2,i+2);
v_z=V(3,i+2);
a_x=0;
a_y=0;
a_z=0;
j=i+2;

for i=0:2000
    mod_v=sqrt(v_x^2+v_y^2+v_z^2);
    mod_a=sqrt(a_x^2+a_y^2+a_z^2);
    mod_w=sqrt(wx^2+wy^2+wz^2);

    H=[-k_d*mod_v,-k_m*wz, k_m*wy; k_m*wz, -k_d*mod_v, -k_m*wx; -k_m*wy, k_m*wx, -k_d*mod_v];
    G=[0;0;-g];
    A(:,j+i+1)=H*V(:,j+i)+G;
    A=horzcat(A,[a_x;a_y;a_z]);
    a_x=A(1,j+i);
    a_y=A(2,j+i);
    a_z=A(3,j+i);
    V(:,j+i+1)= V(:,j+i)+A(:,j+i+1)*sampling_time;
    v_x=V(1,j+i);
    v_y=V(2,j+i);
    v_z=V(3,j+i);
    V=horzcat(V,[v_x;v_y;v_z]);
    X=horzcat(X,[x;y;z]);
    X(:,j+i+1)=X(:,j+i)+V(:,j+i+1)*sampling_time;

        if(X(1,j+i)>6 || X(1,j+i)<0)
            disp(j+i);
            break;
        end
end


ep=0.52;
mu=0.6;

a=1;
b=1;
c=2;

%orientation of racket
n=(1/sqrt(a^2+b^2+c^2))*[a;b;c];


R=n/[V(1,j+i+1);V(2,j+i+1);V(3,j+i+1)];

v_r=V(:,j+i+1)'*n;
mod_v_t=sqrt((V(1,j+i+1)-r*wy)^2+(V(2,j+i+1)+r*wx)^2);

k=mu*(1+ep)*v_r/mod_v_t;

A_r=[1-k 0 0;0 1-k 0;0 0 -ep];
B_r=[0 k*r 0;-k*r 0 0;0 0 0];

A_r_1= R*A_r*R';

V_r=[5+t;-0.5;4+t];
mod_V_r=sqrt(V_r(1)^2+V_r(2)^2+V_r(3)^2);

V(:,j+i+2)=(eye(3)-A_r_1)*V_r + A_r_1*V(:,j+i+1)+R*B_r*[wx;wy;wz];

X(:,j+i+2)=X(:,j+i+1);

disp(X(:,j+i+2));

v_x=V(1,j+i+2);
v_y=V(2,j+i+2);
v_z=V(3,j+i+2);

a_x=0;
a_y=0;
a_z=0;


p=j+i+1;

t=t+1;
end

figure;
plot3(X(1,:),X(2,:),X(3,:));