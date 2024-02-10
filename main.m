%% Read the image
clear
close all 
clc
img = imread('PalazzoTe.jpg');
%Rotating the image to make it vertical
img = imrotate(img, -90);
FNT_SZ = 28;
plotw = length(img)*10;
%% Preprocessing
% Canny Edge Detection
%converting the original image to a gray scale image
gray_im = rgb2gray(img);
resized_gray_im = imresize(gray_im, 1);
curves = uint8([repelem(255,10) linspace(255,0,96) zeros(1,150)]);
preprocessed_im = intlut(resized_gray_im, curves);
%applying the Canny filter to the image
img_canny = edge(preprocessed_im, 'Canny');
subplot(1, 2, 1);  % First subplot
imshow(img);
title('Original Image');
subplot(1, 2, 2);  % Second subplot
imshow(img_canny);
title("Image with Canny edge detection");
%% 1-From 𝐶1, 𝐶2 find the horizon (vanishing) line ℎ of the plane 
% orthogonal to the cone axis 
figure();
imshow(img_canny);
title('Pick points for C1');
hold on;
%selecting 6 points in order to calculate parameters
[x, y]=getpts;
scatter(x,y,100,'filled');
%% Estimation of conic parameters
A1=[x.^2 x.*y y.^2 x y ones(size(x))];
[~,~,V1] = svd(A1);
N1 = V1(:,end);
cc1 = N1(:, 1);
% changing the name of variables
[a1 b1 c1 d1 e1 f1] = deal(cc1(1),cc1(2),cc1(3),cc1(4),cc1(5),cc1(6));
% matrix of the conic: off-diagonal elements must be divided by two
cc1=[a1 b1/2 d1/2; b1/2 c1 e1/2; d1/2 e1/2 f1];
%% selecting 6 points in order to calculate parameters
title('Pick points for C2')
hold on;
[c, z]=getpts;
scatter(c,z,100,'filled');
%% Estimation of conic parameters

A2=[c.^2 c.*z z.^2 c z ones(size(c))];
[~,~,V2] = svd(A2);
N2 = V2(:,end);
cc2 = N2(:, 1);
% changing the name of variables
[a2 b2 c2 d2 e2 f2] = deal(cc2(1),cc2(2),cc2(3),cc2(4),cc2(5),cc2(6));
% matrix of the conic: off-diagonal elements must be divided by two
cc2=[a2 b2/2 d2/2; b2/2 c2 e2/2; d2/2 e2/2 f2];
%% lines
A=[x(1), y(1), 1];
B=[x(6), y(6), 1];

C=[c(1), z(1), 1];
D=[c(6), z(6), 1];
%generatrix lines 
l1 = cross(A,C);
l2 = cross(B,D);

%print all
plot([A(1), C(1)], [A(2), C(2)], 'g');
plot([B(1), D(1)], [B(2), D(2)], 'g');

%% Problem 1 - horizone line h of the plane orthogonal to the cylinder axis
%Computation of image of circular points by computing the intersection
%between the two conics
syms 'g';
syms 's';
eq1 = a1*g^2 + b1*g*s + c1*s^2 + d1*g + e1*s + f1;
eq2 = a2*g^2 + b2*g*s + c2*s^2 + d2*g + e2*s + f2;
eqns= [eq1==0, eq2==0];
fimplicit(eqns);
S = solve(eqns, [g,s]);
s1 = [double(S.g(1));double(S.s(1));1];
s2 = [double(S.g(2));double(S.s(2));1];
s3 = [double(S.g(3));double(S.s(3));1];
s4 = [double(S.g(4));double(S.s(4));1];
% images of circular points
II = s1;
JJ = s2;

%Computation and plot of h
h=cross(II, JJ);
h=h/h(3);

x = linspace(1,100000,1000000);
y=((JJ(2)-II(2))/(JJ(1)-II(1)))*(x-II(1))+II(2);
plot(x,y,'linewidth',2,'color','b')

%% Problem 2 - image projection 𝑎 of the cylinder axis, and its 
% vanishing point V
%centers of cc1 and cc2 by computing their diameters
d1=cc1*II;
d2=cc1*JJ;
d3=cc2*II;
d4=cc2*JJ;
d1= d1/d1(3);
d2= d2/d2(3);
d3= d3/d3(3);
d4= d4/d4(3);
%Computation and plot of centres of C1 and C2
O1=cross(d1,d2);
O2=cross(d3,d4);

O1=O1./O1(3);
O2=O2./O2(3);


scatter(O1(1), O1(2), 100,'filled');
scatter(O2(1), O2(2), 100,'filled');

%Computation of the image of the axis
a=cross(O1,O2);
x = linspace(1,100000,1000000);
y=((O2(2)-O1(2))/(O2(1)-O1(1)))*(x-O1(1))+O1(2);
plot(x,y,'linewidth',2,'color','r')

V = cross(l1,l2);
%% Problem 3 - find the calibration matrix 𝐾
syms aa U0 V0 fy;
omega = [aa^2 0 -U0*aa^2; 0 1 -V0; 0 0 (fy^2+aa^2*U0^2+V0^2)];
%computing vanishing points of diameters which are parallel to x and y
v1 = cross(d1,d3);
v2 = cross(d2,d4);
% 4 equations for 4 unknowns
eq1 = v1.' * omega * v2 == 0;
eq2 = V * omega *  v1 == 0;
eq3 = II.' * omega * II == 0;
eq4 = JJ.' * omega * JJ == 0;
S = solve([eq1, eq2, eq3, eq4],[aa U0 V0 fy]);
U0 = double(S.U0);
V0 = double(S.V0);
aa = double(S.aa);
fy = double(S.fy);
U0 = abs(U0(1));
V0 = abs(V0(1));
aa = abs(aa(1));
fy = abs(fy(1));
fx = aa * fy;
K = [fx 0 U0; 0 fy V0; 0 0 1];

%% Problem 4 - determine the orientation of the cylinder axis wrt the 
% camera reference
%I am computing the angle between the horizon and the axis of the cylinder
angle = acos((h(1) * a(1) + h(2) * a(2)) / ( sqrt((h(1)^2 + h(2)^2) * (a(1)^2 + a(2)^2))));
angle = 180*angle/pi;

%% Problem 5 - Compute the ratio between the radius of the circular cross
% sections and their distance -- NOT FINISHED
% I am computing the image of the dual absolute conic
DIAC = II*JJ' + JJ*II';
[U, S, Vt] = svd(DIAC);
% I am computing two pairs of diametrically opposed points on the conics 
syms 'g';
syms 's';
eq1 = a1*g^2 + b1*g*s + c1*s^2 + d1*g + e1*s + f1;
eq2 = d1(1)*g + d1(2)*s == d1(3);
S1 = solve(eqns, [g,s]);
p1_1 = [abs(double(S1.g(1)));abs(double(S1.s(1)));1];
p2_1 = [abs(double(S1.g(2)));abs(double(S1.s(2)));1];

syms 'g';
syms 's';
eq3 = a2*g^2 + b2*g*s + c2*s^2 + d2*g + e2*s + f2;
eq4 = d3(1)*g + d3(2)*s == d3(3);
S2 = solve(eqns, [g,s]);
p1_2 = [abs(double(S2.g(1)));abs(double(S2.s(1)));1];
p2_2 = [abs(double(S2.g(2)));abs(double(S2.s(2)));1];
Hr = inv(U);
hold on;
scatter(p1_1(1), p1_1(2), 300,'filled');
scatter(p2_1(1), p2_1(2), 300,'filled');
scatter(p1_2(1), p1_2(2), 300,'filled');
scatter(p2_2(1), p2_2(2), 300,'filled'); 

%then I should rectify the points and calculate the ratio with them
