%% Read the image
clear
close all 
clc
img = imread('PalazzoTe.jpg');
img = imrotate(img, -90);
img = imresize(img, 0.2);
%% corners
AA = [250*0.2 1873.8*0.2 1];
BB = [3205.7*0.2 2533.4*0.2 1];
CC = [267.3*0.2 4332.3*0.2 1];
DD = [1850.4*0.2 4338.3*0.2 1];

% showing the vertical plane to be rectified in the image scene
imshow(img), hold on;
% lines
xy = [AA(1:2); BB(1:2)];
plot(xy(:,1),xy(:,2),'LineWidth',4,'Color', 'g');
xy = [BB(1:2); DD(1:2)];
plot(xy(:,1),xy(:,2),'LineWidth',4,'Color', 'g');
xy = [CC(1:2); DD(1:2)];
plot(xy(:,1),xy(:,2),'LineWidth',4,'Color', 'g');
xy = [CC(1:2); AA(1:2)];
plot(xy(:,1),xy(:,2),'LineWidth',4,'Color', 'g');
%%
% points
plot(AA(1), AA(2),'.r','MarkerSize',12);
text(AA(1), AA(2), 'A', 'FontSize', 24, 'Color', 'r');
plot(BB(1), BB(2),'.r','MarkerSize',12);
text(BB(1), BB(2), 'B', 'FontSize', 24, 'Color', 'r');
plot(CC(1), CC(2),'.r','MarkerSize',12);
text(CC(1), CC(2), 'C', 'FontSize', 24, 'Color', 'r');
plot(DD(1), DD(2),'.r','MarkerSize',12);
text(DD(1), DD(2), 'D', 'FontSize', 24, 'Color', 'r');
hold off;
%%
% first vanishing point
line_ab = cross(AA, BB);
line_ab = line_ab / line_ab(3);
line_cd = cross(CC, DD);
line_cd = line_cd / line_cd(3);
VV1 = cross(line_ab, line_cd);
VV1 = VV1 / VV1(3);

% second vanishing point
line_bd = cross(DD, BB);
line_bd = line_bd / line_bd(3);
line_ac = cross(AA, CC);
line_ac = line_ac / line_ac(3);
VV2 = cross(line_bd, line_ac);
VV2 = VV2 / VV2(3);

% image of the line at the infinity
inf_line = cross(VV1, VV2);
inf_line = inf_line./inf_line(3);
hold on;
plot(VV1(1),VV1(2),'r.','MarkerSize',100);
plot(VV2(1),VV2(2),'b.','MarkerSize',100);
xy = [VV1(1:2); VV2(1:2)];
plot(xy(:,1),xy(:,2),'LineWidth',4,'Color', 'g');
%%
H = [eye(2),zeros(2,1); inf_line(:)'];

% applying the homography to the image
tform = projective2d(H');
J = imwarp(img,tform);

% plotting the final outcome
figure();
subplot(1, 2, 1);  % First subplot
imshow(img);
title('Original Image');
subplot(1, 2, 2);  % Second subplot
imshow(J);
title("Image rectified");
imwrite(J,'images/affRect.JPG');