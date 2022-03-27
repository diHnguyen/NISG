global arcs = [1 6; 1 9; 1 10; 1 21; 1 24; 1 28; 1 29; 2 6; 2 8; 2 12; 2 16; 2 21; 2 26; 2 35; 3 6; 3 16; 3 32; 4 10; 4 28; 4 36; 5 7; 5 22; 5 33; 6 5; 6 19; 6 21; 6 29; 6 32; 6 40; 7 4; 7 36; 8 30; 9 3; 9 13; 9 19; 9 21; 10 7; 10 29; 10 39; 10 40; 11 28; 11 36; 12 3; 12 13; 12 17; 13 5; 13 7; 13 14; 13 22; 13 29; 13 34; 14 6; 14 7; 14 15; 14 28; 14 30; 14 32; 14 40; 15 11; 15 31; 15 35; 15 38; 16 5; 16 8; 17 7; 17 14; 17 22; 17 24; 17 39; 18 15; 18 22; 18 23; 19 9; 19 13; 19 33; 20 15; 20 19; 20 23; 21 2; 21 25; 22 7; 22 9; 23 7; 23 27; 23 34; 23 39; 24 12; 24 19; 25 24; 25 28; 26 6; 26 22; 26 30; 26 34; 26 39; 27 5; 27 13; 27 17; 27 28; 27 29; 27 32; 27 35; 27 39; 27 40; 28 8; 28 19; 28 22; 28 37; 28 40; 29 2; 29 13; 29 15; 29 18; 30 11; 30 12; 30 26; 31 9; 31 23; 31 34; 32 13; 32 15; 33 19; 33 40; 34 21; 34 23; 34 25; 35 8; 35 11; 35 39; 36 8; 36 19; 36 24; 36 38; 36 40; 37 18; 37 20; 37 24; 37 40; 38 12; 38 15; 38 17; 38 27; 38 36; 38 40; 39 35]
global d_x = [1.0, 4.0, 4.0, 6.0, 4.0, 9.0, 8.0, 5.0, 3.0, 8.0, 5.0, 8.0, 7.0, 8.0, 10.0, 3.0, 6.0, 7.0, 9.0, 6.0, 9.0, 7.0, 5.0, 4.0, 6.0, 2.0, 5.0, 3.0, 2.0, 4.0, 3.0, 7.0, 6.0, 3.0, 9.0, 5.0, 2.0, 7.0, 9.0, 10.0, 1.0, 8.0, 5.0, 5.0, 1.0, 7.0, 9.0, 4.0, 3.0, 3.0, 8.0, 6.0, 9.0, 9.0, 4.0, 3.0, 3.0, 5.0, 2.0, 1.0, 3.0, 4.0, 4.0, 2.0, 5.0, 1.0, 2.0, 10.0, 5.0, 1.0, 3.0, 9.0, 8.0, 6.0, 7.0, 4.0, 4.0, 5.0, 8.0, 6.0, 1.0, 5.0, 5.0, 9.0, 2.0, 2.0, 4.0, 7.0, 2.0, 5.0, 2.0, 1.0, 4.0, 10.0, 1.0, 2.0, 9.0, 3.0, 6.0, 2.0, 5.0, 9.0, 9.0, 5.0, 9.0, 6.0, 8.0, 1.0, 2.0, 3.0, 5.0, 4.0, 4.0, 9.0, 8.0, 6.0, 7.0, 10.0, 4.0, 10.0, 5.0, 2.0, 6.0, 1.0, 10.0, 5.0, 8.0, 10.0, 6.0, 4.0, 1.0, 6.0, 6.0, 1.0, 6.0, 7.0, 5.0, 6.0, 2.0, 6.0, 10.0, 2.0, 6.0, 2.0, 5.0]
global b_x = 5
global d_y = [9.0, 2.0, 4.0, 7.0, 6.0, 6.0, 9.0, 8.0, 5.0, 5.0, 7.0, 5.0, 7.0, 8.0, 3.0, 6.0, 8.0, 1.0, 10.0, 9.0, 1.0, 9.0, 6.0, 3.0, 1.0, 8.0, 5.0, 3.0, 7.0, 2.0, 8.0, 9.0, 4.0, 4.0, 3.0, 6.0, 10.0, 9.0, 3.0, 6.0, 9.0, 6.0, 6.0, 10.0, 2.0, 1.0, 2.0, 5.0, 6.0, 4.0, 8.0, 3.0, 10.0, 2.0, 7.0, 3.0, 10.0, 5.0, 7.0, 8.0, 8.0, 1.0, 1.0, 3.0, 1.0, 7.0, 2.0, 9.0, 3.0, 10.0, 3.0, 1.0, 8.0, 1.0, 9.0, 5.0, 3.0, 8.0, 4.0, 6.0, 1.0, 10.0, 4.0, 3.0, 10.0, 10.0, 9.0, 1.0, 6.0, 6.0, 3.0, 10.0, 7.0, 1.0, 4.0, 4.0, 8.0, 7.0, 1.0, 1.0, 1.0, 10.0, 5.0, 5.0, 7.0, 10.0, 6.0, 8.0, 1.0, 2.0, 7.0, 9.0, 5.0, 5.0, 1.0, 7.0, 10.0, 7.0, 9.0, 4.0, 4.0, 1.0, 3.0, 2.0, 5.0, 2.0, 1.0, 5.0, 10.0, 10.0, 6.0, 1.0, 9.0, 1.0, 10.0, 6.0, 4.0, 1.0, 6.0, 10.0, 2.0, 6.0, 5.0, 6.0, 1.0]
global b_y = 10
global p = [0.175, 0.308, 0.247, 0.042, 0.97, 0.684, 0.947, 0.411, 0.006, 0.326, 0.701, 0.531, 0.884, 0.599, 0.339, 0.969, 0.088, 0.357, 0.767, 0.171, 0.798, 0.838, 0.07, 0.741, 0.969, 0.081, 0.315, 0.072, 0.733, 0.744, 0.17, 0.675, 0.394, 0.074, 0.085, 0.919, 0.165, 0.42, 0.407, 0.371, 0.239, 0.878, 0.695, 0.7, 0.375, 0.1, 0.663, 0.431, 0.415, 0.129, 0.073, 0.555, 0.979, 0.784, 0.369, 0.174, 0.622, 0.873, 0.777, 0.254, 0.756, 0.685, 0.057, 0.979, 0.029, 0.001, 0.093, 0.345, 0.596, 0.769, 0.376, 0.485, 0.024, 0.807, 0.406, 0.248, 0.832, 0.012, 0.544, 0.146, 0.354, 0.072, 0.525, 0.681, 0.275, 0.682, 0.207, 0.798, 0.901, 0.312, 0.19, 0.088, 0.809, 0.384, 0.495, 0.744, 0.795, 0.894, 0.247, 0.902, 0.816, 0.457, 0.334, 0.445, 0.516, 0.353, 0.014, 0.023, 0.798, 0.07, 0.722, 0.538, 0.35, 0.278, 0.675, 0.653, 0.796, 0.238, 0.52, 0.002, 0.358, 0.483, 0.536, 0.047, 0.524, 0.422, 0.474, 0.558, 0.062, 0.95, 0.425, 0.405, 0.699, 0.029, 0.168, 0.2, 0.18, 0.422, 0.923, 0.638, 0.968, 0.272, 0.999, 0.912, 0.112]
global q = [0.892, 0.885, 0.671, 0.448, 0.999, 0.828, 0.999, 0.844, 0.595, 0.955, 0.999, 0.762, 0.964, 0.924, 0.467, 0.984, 0.58, 0.693, 0.808, 0.228, 0.88, 0.936, 0.817, 0.876, 0.998, 0.671, 0.681, 0.08, 0.864, 0.95, 0.447, 0.704, 0.435, 0.974, 0.23, 0.979, 0.868, 0.623, 0.536, 0.739, 0.53, 0.933, 0.934, 0.799, 0.696, 0.119, 0.957, 0.54, 0.574, 0.43, 0.248, 0.908, 0.982, 0.889, 0.981, 0.388, 0.82, 0.884, 0.893, 0.532, 0.868, 0.871, 0.69, 0.998, 0.37, 0.833, 0.793, 0.897, 0.796, 0.868, 0.467, 0.975, 0.188, 0.863, 0.628, 0.621, 0.845, 0.977, 0.798, 0.513, 0.585, 0.217, 0.789, 0.999, 0.925, 0.763, 0.871, 0.906, 0.974, 0.633, 0.401, 0.951, 0.812, 0.389, 0.587, 0.838, 0.994, 0.917, 0.813, 0.978, 0.867, 0.777, 0.485, 0.941, 0.822, 0.929, 0.675, 0.98, 0.994, 0.353, 0.957, 0.982, 0.885, 0.843, 0.686, 0.683, 0.859, 0.337, 0.833, 0.008, 0.747, 0.637, 0.908, 0.882, 0.95, 0.525, 0.825, 0.795, 0.178, 0.963, 0.939, 0.854, 0.742, 0.712, 0.95, 0.923, 0.55, 0.72, 0.967, 0.96, 0.981, 0.887, 0.999, 0.963, 0.921]
global origin = 1
global destination = 40