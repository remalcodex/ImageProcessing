clc;
close all;
imtool close all;
addpath(genpath('Functions/'));
%addpath(genpath('Inputs/'));
ceilfix = @(x)ceil(abs(x)).*sign(x);

%Step 1: Morphing the image.
%morph('Inputs/morph_harvey.txt'); %Morphing of joker.
%morph('Inputs/morph_harvey_wrong.txt'); %Morphing of joker.
morph('Inputs/morph_joker.txt'); %Morphing of harvey.
%morph('Inputs/params_morph.txt'); %Morphing of the shape.

%Step 2: Atlas creation.
%atlas('Inputs/params_atlas.txt'); %Atlas of shapes.
%atlas('Inputs/params_atlas_wrong.txt'); %Atlas of shapes.
%atlas('Inputs/brain/params_atlas.txt'); %Atlas of brain.
%atlas('Inputs/pumpkin/params_atlas.txt'); %Atlas of pumpkin.

%Remove close all if you wish to see intermediate outputs.
close all;