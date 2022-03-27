global arcs = [1 8; 1 25; 1 27; 2 27; 2 34; 3 2; 3 7; 3 21; 4 13; 4 23; 4 28; 4 31; 5 7; 5 20; 5 30; 5 39; 6 4; 6 17; 6 18; 6 19; 6 36; 7 3; 7 35; 8 5; 8 6; 8 11; 8 18; 8 27; 8 28; 8 35; 9 26; 10 12; 10 16; 10 18; 10 26; 11 8; 11 9; 11 13; 11 19; 11 26; 11 28; 11 38; 12 2; 12 7; 12 9; 12 29; 12 30; 12 31; 12 33; 12 36; 13 14; 13 25; 13 28; 13 29; 13 40; 14 5; 14 9; 14 27; 14 30; 15 2; 15 6; 15 7; 15 16; 15 34; 16 14; 16 29; 16 32; 16 40; 17 11; 17 12; 17 27; 18 2; 18 5; 18 31; 18 35; 19 2; 19 24; 19 27; 20 34; 21 11; 21 32; 22 31; 23 8; 23 16; 23 18; 23 26; 23 27; 23 28; 23 32; 23 37; 24 4; 24 17; 25 4; 25 5; 25 32; 25 36; 26 29; 27 8; 27 14; 27 19; 27 37; 27 38; 28 8; 28 27; 29 3; 29 13; 29 25; 29 33; 29 34; 30 13; 31 6; 31 13; 31 15; 31 18; 31 26; 31 38; 32 21; 33 6; 33 7; 33 8; 33 10; 33 11; 33 26; 33 29; 33 35; 34 2; 34 5; 34 8; 34 16; 34 18; 34 27; 34 30; 35 14; 35 16; 35 22; 35 24; 35 26; 35 32; 35 37; 36 3; 36 14; 36 16; 36 22; 36 25; 37 2; 37 15; 37 22; 37 27; 37 31; 38 2; 38 7; 38 12; 38 19; 38 35; 39 20; 39 31; 39 32]
global d_x = [10.0, 9.0, 8.0, 8.0, 8.0, 4.0, 3.0, 2.0, 2.0, 3.0, 8.0, 4.0, 3.0, 5.0, 6.0, 5.0, 8.0, 9.0, 9.0, 2.0, 2.0, 7.0, 3.0, 1.0, 1.0, 8.0, 7.0, 3.0, 6.0, 9.0, 8.0, 8.0, 6.0, 5.0, 5.0, 3.0, 10.0, 1.0, 3.0, 4.0, 2.0, 1.0, 3.0, 5.0, 8.0, 2.0, 2.0, 3.0, 1.0, 3.0, 2.0, 7.0, 3.0, 4.0, 1.0, 2.0, 1.0, 5.0, 5.0, 7.0, 10.0, 8.0, 6.0, 2.0, 3.0, 4.0, 9.0, 9.0, 7.0, 6.0, 5.0, 10.0, 5.0, 9.0, 8.0, 5.0, 10.0, 7.0, 1.0, 8.0, 9.0, 4.0, 6.0, 7.0, 4.0, 5.0, 3.0, 9.0, 3.0, 2.0, 9.0, 6.0, 1.0, 1.0, 8.0, 10.0, 1.0, 8.0, 10.0, 4.0, 10.0, 9.0, 1.0, 5.0, 3.0, 6.0, 9.0, 2.0, 2.0, 1.0, 9.0, 6.0, 9.0, 7.0, 4.0, 7.0, 9.0, 4.0, 9.0, 6.0, 10.0, 1.0, 6.0, 5.0, 6.0, 6.0, 6.0, 1.0, 6.0, 5.0, 10.0, 2.0, 6.0, 8.0, 1.0, 7.0, 3.0, 2.0, 8.0, 2.0, 4.0, 10.0, 7.0, 7.0, 5.0, 2.0, 8.0, 7.0, 5.0, 8.0, 3.0, 2.0, 3.0, 1.0, 10.0, 10.0, 6.0]
global b_x = 5
global d_y = [5.0, 8.0, 8.0, 4.0, 6.0, 3.0, 1.0, 3.0, 3.0, 6.0, 2.0, 9.0, 1.0, 6.0, 3.0, 1.0, 5.0, 8.0, 2.0, 3.0, 4.0, 7.0, 1.0, 4.0, 5.0, 8.0, 4.0, 2.0, 5.0, 5.0, 3.0, 1.0, 3.0, 8.0, 5.0, 4.0, 8.0, 3.0, 4.0, 3.0, 2.0, 7.0, 8.0, 5.0, 6.0, 3.0, 5.0, 4.0, 10.0, 4.0, 9.0, 5.0, 10.0, 8.0, 3.0, 10.0, 2.0, 8.0, 1.0, 7.0, 2.0, 6.0, 3.0, 1.0, 3.0, 9.0, 9.0, 1.0, 1.0, 9.0, 3.0, 3.0, 5.0, 4.0, 9.0, 4.0, 3.0, 2.0, 10.0, 6.0, 1.0, 4.0, 7.0, 9.0, 10.0, 6.0, 8.0, 6.0, 1.0, 7.0, 10.0, 1.0, 8.0, 1.0, 1.0, 10.0, 3.0, 2.0, 4.0, 4.0, 9.0, 7.0, 10.0, 2.0, 2.0, 9.0, 7.0, 5.0, 3.0, 10.0, 5.0, 2.0, 4.0, 1.0, 3.0, 4.0, 6.0, 8.0, 5.0, 1.0, 7.0, 3.0, 10.0, 1.0, 5.0, 10.0, 10.0, 10.0, 10.0, 7.0, 8.0, 3.0, 10.0, 2.0, 7.0, 4.0, 10.0, 4.0, 10.0, 10.0, 6.0, 4.0, 5.0, 3.0, 4.0, 8.0, 2.0, 4.0, 10.0, 7.0, 7.0, 9.0, 1.0, 5.0, 4.0, 7.0, 5.0]
global b_y = 10
global p = [0.599, 0.117, 0.313, 0.812, 0.604, 0.532, 0.109, 0.897, 0.659, 0.687, 0.973, 0.61, 0.89, 0.506, 0.212, 0.93, 0.764, 0.473, 0.088, 0.056, 0.776, 0.523, 0.696, 0.532, 0.106, 0.483, 0.917, 0.89, 0.349, 0.302, 0.366, 0.172, 0.786, 0.715, 0.792, 0.481, 0.273, 0.628, 0.007, 0.933, 0.922, 0.888, 0.087, 0.978, 0.576, 0.667, 0.446, 0.312, 0.981, 0.698, 0.529, 0.223, 0.692, 0.476, 0.945, 0.15, 0.711, 0.572, 0.675, 0.486, 0.384, 0.641, 0.535, 0.391, 0.067, 0.432, 0.621, 0.2, 0.87, 0.241, 0.54, 0.061, 0.942, 0.161, 0.656, 0.299, 0.504, 0.235, 0.686, 0.826, 0.677, 0.4, 0.193, 0.085, 0.985, 0.649, 0.326, 0.132, 0.569, 0.72, 0.623, 0.747, 0.217, 0.466, 0.844, 0.237, 0.96, 0.131, 0.941, 0.83, 0.447, 0.541, 0.457, 0.078, 0.512, 0.349, 0.29, 0.246, 0.683, 0.784, 0.018, 0.247, 0.259, 0.029, 0.464, 0.932, 0.576, 0.614, 0.63, 0.098, 0.055, 0.397, 0.718, 0.882, 0.932, 0.954, 0.808, 0.867, 0.587, 0.415, 0.021, 0.257, 0.768, 0.159, 0.423, 0.767, 0.596, 0.448, 0.857, 0.66, 0.303, 0.46, 0.59, 0.878, 0.706, 0.139, 0.43, 0.087, 0.186, 0.235, 0.128, 0.041, 0.005, 0.567, 0.807, 0.78, 0.522]
global q = [0.637, 0.574, 0.7, 0.935, 0.829, 0.672, 0.793, 0.955, 0.767, 0.828, 0.976, 0.921, 0.913, 0.792, 0.313, 0.997, 0.775, 0.694, 0.359, 0.367, 0.897, 0.644, 0.742, 0.713, 0.474, 0.65, 0.962, 0.924, 0.489, 0.375, 0.972, 0.393, 0.964, 0.78, 0.854, 0.523, 0.366, 0.706, 0.864, 0.984, 0.942, 0.914, 0.728, 0.998, 0.812, 0.96, 0.976, 0.458, 0.982, 0.837, 0.633, 0.78, 0.725, 0.641, 0.974, 0.361, 0.918, 0.643, 0.858, 0.689, 0.515, 0.831, 0.842, 0.403, 0.96, 0.965, 0.961, 0.783, 0.981, 0.392, 0.618, 0.648, 0.974, 0.245, 0.67, 0.318, 0.68, 0.565, 0.768, 0.995, 0.802, 0.659, 0.87, 0.143, 0.99, 0.969, 0.7, 0.669, 0.613, 0.808, 0.987, 0.944, 0.233, 0.969, 0.901, 0.855, 0.982, 0.409, 0.994, 0.869, 0.693, 0.951, 0.934, 0.41, 0.974, 0.398, 0.639, 0.771, 0.978, 0.801, 0.064, 0.857, 0.653, 0.261, 0.984, 0.94, 0.677, 0.917, 0.992, 0.946, 0.798, 0.799, 0.817, 0.978, 0.974, 0.958, 0.894, 0.917, 0.946, 0.518, 0.172, 0.961, 0.861, 0.724, 0.478, 0.884, 0.598, 0.814, 0.969, 0.707, 0.676, 0.877, 0.657, 0.94, 0.828, 0.756, 0.446, 0.695, 0.734, 0.484, 0.366, 0.178, 0.478, 0.826, 0.92, 0.859, 0.921]
global origin = 1
global destination = 40