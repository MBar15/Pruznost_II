clc;
clear all;
close all;


syms

b = 80;
h = 160;
t = 3;

xt = (t*b*b/2 + (h-2*t)*t*t/2 + b*t*b/2) / (b*t + (h-2*t)*t + b*t)
yt = (t*b*t/2 + (h-2*t)*t*h/2 + b*t*(h-t/2)) / (b*t + (h-2*t)*t + b*t)


U = b*t

