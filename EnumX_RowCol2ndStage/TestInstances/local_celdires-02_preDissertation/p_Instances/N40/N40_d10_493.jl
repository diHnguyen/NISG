global arcs = [1 2; 1 6; 1 10; 1 23; 1 30; 2 3; 2 16; 2 27; 2 29; 2 38; 3 7; 3 10; 3 36; 3 37; 4 12; 4 29; 4 33; 5 14; 5 40; 6 4; 6 20; 6 22; 7 10; 7 17; 7 28; 7 33; 8 6; 9 15; 9 16; 9 18; 9 28; 9 31; 9 39; 10 2; 10 18; 11 3; 11 13; 11 23; 12 14; 12 15; 12 28; 12 30; 12 40; 13 2; 13 3; 13 31; 14 11; 14 19; 14 23; 14 29; 14 32; 15 36; 16 3; 16 17; 16 23; 16 26; 16 31; 17 9; 17 11; 17 14; 18 15; 18 25; 18 32; 18 34; 18 39; 19 15; 19 28; 19 35; 19 38; 20 13; 20 22; 20 28; 20 35; 21 3; 21 5; 21 35; 21 37; 21 38; 22 11; 22 12; 22 30; 22 34; 23 24; 23 29; 24 6; 24 17; 24 22; 24 33; 25 3; 25 8; 25 12; 26 2; 26 13; 26 27; 26 32; 26 39; 27 4; 27 19; 27 21; 27 22; 28 13; 28 29; 28 36; 29 4; 29 18; 29 25; 29 31; 29 32; 29 37; 30 2; 30 7; 30 8; 30 24; 30 27; 31 7; 31 16; 31 17; 31 24; 32 18; 32 27; 32 35; 32 36; 33 9; 33 11; 33 15; 33 22; 33 26; 33 28; 33 30; 33 36; 34 2; 34 10; 34 15; 34 25; 34 33; 35 19; 35 27; 35 29; 36 5; 37 10; 37 11; 38 9; 38 12; 39 32; 39 35]
global d_x = [1.0, 4.0, 1.0, 2.0, 7.0, 2.0, 10.0, 5.0, 4.0, 9.0, 8.0, 5.0, 8.0, 3.0, 3.0, 7.0, 8.0, 2.0, 2.0, 6.0, 6.0, 8.0, 9.0, 5.0, 5.0, 6.0, 5.0, 6.0, 4.0, 7.0, 3.0, 7.0, 4.0, 3.0, 3.0, 4.0, 6.0, 6.0, 10.0, 3.0, 6.0, 4.0, 4.0, 7.0, 10.0, 1.0, 1.0, 6.0, 6.0, 5.0, 8.0, 3.0, 7.0, 1.0, 9.0, 2.0, 4.0, 8.0, 5.0, 1.0, 1.0, 9.0, 7.0, 5.0, 10.0, 5.0, 7.0, 1.0, 1.0, 8.0, 7.0, 7.0, 6.0, 10.0, 4.0, 9.0, 10.0, 2.0, 7.0, 3.0, 8.0, 1.0, 10.0, 4.0, 10.0, 2.0, 4.0, 1.0, 2.0, 1.0, 4.0, 1.0, 7.0, 9.0, 3.0, 7.0, 5.0, 9.0, 7.0, 8.0, 3.0, 8.0, 6.0, 10.0, 10.0, 9.0, 7.0, 10.0, 7.0, 5.0, 8.0, 5.0, 10.0, 10.0, 4.0, 10.0, 4.0, 5.0, 5.0, 4.0, 3.0, 9.0, 10.0, 8.0, 4.0, 5.0, 2.0, 9.0, 8.0, 3.0, 4.0, 8.0, 6.0, 6.0, 8.0, 9.0, 5.0, 1.0, 3.0, 5.0, 10.0, 1.0, 5.0, 7.0, 1.0]
global b_x = 5
global d_y = [5.0, 4.0, 9.0, 6.0, 10.0, 2.0, 9.0, 10.0, 8.0, 4.0, 10.0, 7.0, 5.0, 7.0, 8.0, 1.0, 6.0, 6.0, 7.0, 2.0, 9.0, 6.0, 2.0, 5.0, 1.0, 3.0, 8.0, 4.0, 6.0, 10.0, 2.0, 10.0, 8.0, 1.0, 8.0, 2.0, 1.0, 4.0, 5.0, 7.0, 3.0, 2.0, 5.0, 1.0, 2.0, 3.0, 3.0, 3.0, 9.0, 6.0, 10.0, 3.0, 4.0, 5.0, 10.0, 8.0, 2.0, 5.0, 9.0, 5.0, 5.0, 1.0, 9.0, 10.0, 6.0, 10.0, 10.0, 3.0, 1.0, 10.0, 2.0, 10.0, 8.0, 3.0, 6.0, 9.0, 9.0, 4.0, 5.0, 3.0, 7.0, 1.0, 2.0, 6.0, 7.0, 5.0, 1.0, 3.0, 10.0, 9.0, 2.0, 2.0, 5.0, 7.0, 7.0, 8.0, 3.0, 4.0, 9.0, 4.0, 7.0, 2.0, 10.0, 6.0, 2.0, 10.0, 7.0, 10.0, 10.0, 3.0, 3.0, 5.0, 9.0, 6.0, 2.0, 1.0, 2.0, 8.0, 3.0, 8.0, 2.0, 5.0, 3.0, 10.0, 1.0, 4.0, 9.0, 8.0, 3.0, 6.0, 2.0, 8.0, 10.0, 7.0, 5.0, 8.0, 2.0, 5.0, 6.0, 5.0, 9.0, 5.0, 9.0, 6.0, 9.0]
global b_y = 10
global p = [0.533, 0.653, 0.208, 0.333, 0.221, 0.769, 0.941, 0.258, 0.7, 0.174, 0.091, 0.308, 0.867, 0.334, 0.435, 0.685, 0.823, 0.306, 0.335, 0.28, 0.865, 0.529, 0.962, 0.979, 0.626, 0.923, 0.613, 0.962, 0.886, 0.896, 0.834, 0.759, 0.609, 0.28, 0.833, 0.331, 0.484, 0.466, 0.274, 0.706, 0.162, 0.111, 0.935, 0.161, 0.152, 0.488, 0.638, 0.921, 0.12, 0.871, 0.869, 0.08, 0.32, 0.875, 0.721, 0.327, 0.695, 0.642, 0.078, 0.059, 0.546, 0.773, 0.509, 0.695, 0.432, 0.767, 0.584, 0.109, 0.253, 0.707, 0.505, 0.443, 0.215, 0.838, 0.116, 0.85, 0.414, 0.802, 0.693, 0.821, 0.556, 0.816, 0.72, 0.235, 0.686, 0.394, 0.595, 0.748, 0.881, 0.782, 0.58, 0.81, 0.212, 0.483, 0.548, 0.891, 0.326, 0.679, 0.643, 0.503, 0.597, 0.186, 0.723, 0.147, 0.17, 0.356, 0.163, 0.3, 0.712, 0.499, 0.933, 0.125, 0.338, 0.016, 0.88, 0.885, 0.96, 0.114, 0.234, 0.603, 0.082, 0.864, 0.294, 0.02, 0.534, 0.096, 0.125, 0.096, 0.776, 0.093, 0.568, 0.232, 0.121, 0.861, 0.538, 0.088, 0.812, 0.177, 0.194, 0.131, 0.437, 0.931, 0.208, 0.731, 0.658]
global q = [0.797, 0.89, 0.398, 0.827, 0.335, 0.981, 0.999, 0.28, 0.829, 0.701, 0.204, 0.937, 0.891, 0.665, 0.639, 0.967, 0.844, 0.439, 0.424, 0.992, 0.976, 0.529, 0.965, 0.985, 0.731, 0.971, 0.749, 0.988, 0.907, 0.933, 0.934, 0.974, 0.874, 0.887, 0.895, 0.658, 0.889, 0.716, 0.618, 0.792, 0.922, 0.507, 0.961, 0.243, 0.596, 0.832, 0.654, 0.956, 0.593, 0.931, 0.91, 0.85, 0.797, 0.935, 0.746, 0.339, 0.962, 0.726, 0.174, 0.472, 0.747, 0.899, 0.802, 0.737, 0.437, 0.99, 0.693, 0.577, 0.279, 0.862, 0.551, 0.813, 0.337, 0.933, 0.135, 0.997, 0.815, 0.866, 0.955, 0.916, 0.958, 0.913, 0.838, 0.502, 0.705, 0.983, 0.728, 0.946, 0.988, 0.934, 0.67, 0.903, 0.59, 0.921, 0.788, 0.971, 0.674, 0.983, 0.927, 0.917, 0.861, 0.613, 0.899, 0.592, 0.795, 0.788, 0.328, 0.556, 0.777, 0.946, 0.971, 0.826, 0.442, 0.972, 0.969, 0.98, 0.967, 0.527, 0.978, 0.718, 0.987, 0.91, 0.564, 0.305, 0.927, 0.47, 0.194, 0.53, 0.974, 0.565, 0.605, 0.924, 0.667, 0.901, 0.731, 0.131, 0.911, 0.896, 0.374, 0.976, 0.463, 0.965, 0.623, 0.838, 0.695]
global origin = 1
global destination = 40