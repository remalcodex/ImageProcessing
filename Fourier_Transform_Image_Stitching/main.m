clc;
close all;
imtool close all;
addpath(genpath('Functions/'));

%phaseCorrelation('Inputs/Test/');
phaseCorrelation('Inputs/Cells1/');
%phaseCorrelation('Inputs/Cells2/');

%Remove close all if you wish to see intermediate outputs.
close all;