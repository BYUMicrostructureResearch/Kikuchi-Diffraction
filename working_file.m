clear;
clc;

A = Kikuchi_Simulation();
A.Readfile('test_Ni.txt',[],4000);
A.Kikuchi_Calculation([1 0 0], 256, [], 3);
set(gcf(), 'name', 'Ni Sphere')
I = getimage(gcf());
G = fspecial('gaussian',[5 5],5);
Ig = imfilter(I,G,'same');
imshow(Ig);
% 
% B = Kikuchi_Simulation();
% B.Readfile('test_Al.txt',[],4000);
% B.Kikuchi_Calculation([0 1 1], 512, [], 2, 0);
% set(gcf(), 'name', 'Al_Fundamental Zone')
% 
% C = Kikuchi_Simulation();
% C.Readfile('test_Zr.txt',[],4000);
% C.Kikuchi_Calculation([1 0 10], 512, [], 2, 0);
% set(gcf(), 'name', 'Zr_Fundamental Zone')
% 
% D = Kikuchi_Simulation();
% D.Readfile('test_Fe.txt',[],4000);
% D.Kikuchi_Calculation([0 1 1], 512, [], 2, 0);
% set(gcf(), 'name', 'Fe_Fundamental Zone')
% 
% imsave();