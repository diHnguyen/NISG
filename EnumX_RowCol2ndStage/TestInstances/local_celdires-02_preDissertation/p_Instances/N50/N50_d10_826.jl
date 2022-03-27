global arcs = [1 4; 1 10; 1 21; 1 32; 1 48; 2 13; 2 25; 2 26; 2 37; 3 21; 3 23; 3 35; 4 6; 4 40; 4 43; 4 49; 4 50; 5 3; 5 12; 5 13; 5 27; 5 39; 5 48; 6 11; 6 47; 6 48; 7 5; 7 35; 8 3; 8 4; 8 6; 8 7; 8 33; 9 13; 9 16; 9 19; 9 39; 9 41; 9 46; 10 12; 10 23; 10 40; 10 50; 11 48; 12 2; 12 46; 13 8; 13 24; 13 33; 13 35; 13 38; 13 40; 13 50; 14 5; 14 10; 14 19; 14 38; 14 39; 14 43; 15 4; 15 28; 15 39; 15 42; 16 3; 16 4; 16 15; 16 20; 16 39; 17 19; 17 39; 18 9; 18 20; 18 26; 18 36; 18 39; 18 46; 19 5; 19 30; 19 31; 19 32; 19 49; 20 7; 20 17; 20 19; 20 21; 20 22; 20 27; 20 32; 20 36; 20 46; 21 17; 21 47; 22 10; 22 12; 22 49; 22 50; 23 10; 23 21; 23 34; 24 2; 24 26; 24 40; 25 10; 25 12; 25 22; 25 26; 25 28; 25 33; 25 36; 25 37; 26 13; 27 34; 28 5; 28 12; 28 14; 28 17; 28 19; 28 20; 28 25; 28 50; 29 12; 29 18; 29 21; 29 22; 29 31; 29 38; 30 5; 30 14; 30 15; 31 2; 31 29; 31 45; 32 5; 32 15; 32 22; 32 24; 32 28; 32 29; 32 35; 32 43; 33 41; 33 44; 34 2; 34 9; 34 24; 34 29; 34 35; 34 37; 34 39; 34 48; 35 22; 35 29; 35 30; 35 45; 36 13; 36 27; 36 42; 37 3; 37 7; 37 8; 37 24; 37 28; 37 45; 38 22; 38 33; 38 37; 39 15; 39 18; 39 19; 39 29; 39 32; 39 50; 40 8; 40 32; 40 34; 40 35; 40 41; 41 10; 41 13; 41 19; 41 22; 42 13; 42 20; 42 21; 42 24; 42 27; 42 29; 42 39; 43 11; 43 14; 43 18; 43 23; 43 31; 43 42; 43 46; 44 8; 44 23; 44 31; 44 34; 44 37; 44 50; 45 23; 45 25; 45 34; 45 35; 45 38; 46 5; 46 14; 46 16; 46 17; 46 28; 46 50; 47 16; 47 22; 47 31; 47 38; 48 4; 48 21; 48 35; 48 38; 48 42; 49 10; 49 17; 49 18; 49 31; 49 39; 49 40; 49 44]
global d_x = [9.0, 10.0, 5.0, 9.0, 6.0, 3.0, 8.0, 7.0, 3.0, 7.0, 6.0, 1.0, 10.0, 6.0, 2.0, 10.0, 4.0, 6.0, 10.0, 2.0, 9.0, 1.0, 2.0, 5.0, 9.0, 9.0, 5.0, 8.0, 1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 9.0, 2.0, 7.0, 3.0, 9.0, 4.0, 4.0, 10.0, 2.0, 9.0, 9.0, 3.0, 6.0, 1.0, 1.0, 7.0, 4.0, 5.0, 9.0, 4.0, 1.0, 8.0, 10.0, 6.0, 7.0, 1.0, 7.0, 10.0, 3.0, 1.0, 6.0, 7.0, 9.0, 8.0, 8.0, 6.0, 8.0, 7.0, 8.0, 6.0, 2.0, 9.0, 5.0, 4.0, 8.0, 2.0, 4.0, 10.0, 1.0, 4.0, 4.0, 2.0, 6.0, 6.0, 10.0, 1.0, 2.0, 6.0, 2.0, 4.0, 4.0, 5.0, 8.0, 2.0, 1.0, 10.0, 7.0, 8.0, 5.0, 5.0, 1.0, 5.0, 6.0, 4.0, 6.0, 3.0, 3.0, 4.0, 4.0, 5.0, 3.0, 3.0, 10.0, 4.0, 6.0, 6.0, 10.0, 7.0, 8.0, 4.0, 7.0, 9.0, 4.0, 4.0, 3.0, 3.0, 1.0, 3.0, 4.0, 6.0, 10.0, 5.0, 6.0, 8.0, 6.0, 2.0, 8.0, 8.0, 9.0, 10.0, 3.0, 6.0, 4.0, 5.0, 8.0, 2.0, 10.0, 6.0, 9.0, 3.0, 4.0, 5.0, 3.0, 5.0, 1.0, 7.0, 4.0, 8.0, 8.0, 3.0, 1.0, 9.0, 6.0, 5.0, 9.0, 2.0, 6.0, 8.0, 9.0, 5.0, 6.0, 4.0, 10.0, 5.0, 9.0, 4.0, 1.0, 5.0, 3.0, 5.0, 2.0, 2.0, 8.0, 9.0, 3.0, 6.0, 10.0, 3.0, 7.0, 9.0, 4.0, 7.0, 1.0, 6.0, 8.0, 3.0, 1.0, 10.0, 10.0, 2.0, 6.0, 2.0, 4.0, 9.0, 3.0, 5.0, 2.0, 8.0, 10.0, 9.0, 8.0, 8.0, 6.0, 8.0, 6.0, 7.0, 8.0, 9.0, 5.0, 3.0, 5.0, 5.0, 8.0, 5.0]
global b_x = 5
global d_y = [7.0, 8.0, 5.0, 9.0, 5.0, 6.0, 3.0, 5.0, 9.0, 2.0, 9.0, 3.0, 6.0, 1.0, 10.0, 4.0, 2.0, 5.0, 9.0, 4.0, 9.0, 6.0, 7.0, 7.0, 9.0, 3.0, 6.0, 2.0, 10.0, 1.0, 5.0, 3.0, 9.0, 9.0, 5.0, 8.0, 5.0, 2.0, 3.0, 2.0, 6.0, 5.0, 9.0, 3.0, 3.0, 7.0, 3.0, 6.0, 9.0, 7.0, 8.0, 7.0, 6.0, 4.0, 6.0, 2.0, 8.0, 10.0, 10.0, 5.0, 5.0, 6.0, 2.0, 10.0, 2.0, 5.0, 7.0, 10.0, 2.0, 5.0, 10.0, 2.0, 2.0, 4.0, 6.0, 9.0, 1.0, 3.0, 8.0, 4.0, 9.0, 4.0, 10.0, 8.0, 1.0, 8.0, 5.0, 7.0, 10.0, 9.0, 8.0, 5.0, 4.0, 1.0, 5.0, 1.0, 3.0, 3.0, 8.0, 5.0, 1.0, 2.0, 5.0, 2.0, 8.0, 9.0, 9.0, 7.0, 4.0, 5.0, 5.0, 7.0, 4.0, 7.0, 8.0, 4.0, 2.0, 5.0, 1.0, 8.0, 7.0, 6.0, 1.0, 1.0, 1.0, 4.0, 9.0, 4.0, 9.0, 6.0, 9.0, 2.0, 10.0, 7.0, 6.0, 3.0, 1.0, 9.0, 8.0, 4.0, 10.0, 8.0, 6.0, 2.0, 10.0, 6.0, 2.0, 8.0, 9.0, 8.0, 7.0, 3.0, 8.0, 4.0, 5.0, 6.0, 9.0, 10.0, 2.0, 4.0, 3.0, 6.0, 5.0, 7.0, 9.0, 1.0, 3.0, 3.0, 8.0, 8.0, 1.0, 6.0, 6.0, 7.0, 10.0, 3.0, 4.0, 5.0, 10.0, 9.0, 3.0, 7.0, 5.0, 1.0, 8.0, 7.0, 2.0, 10.0, 9.0, 9.0, 9.0, 1.0, 3.0, 8.0, 7.0, 7.0, 10.0, 1.0, 9.0, 4.0, 9.0, 2.0, 2.0, 1.0, 1.0, 1.0, 6.0, 9.0, 9.0, 6.0, 6.0, 9.0, 10.0, 5.0, 10.0, 5.0, 1.0, 4.0, 8.0, 6.0, 9.0, 3.0, 3.0, 8.0, 6.0, 3.0, 4.0, 6.0]
global b_y = 10
global p = [0.431, 0.079, 0.868, 0.335, 0.615, 0.877, 0.677, 0.481, 0.958, 0.905, 0.22, 0.534, 0.903, 0.926, 0.667, 0.16, 0.563, 0.434, 0.827, 0.116, 0.904, 0.89, 0.617, 0.558, 0.12, 0.36, 0.362, 0.539, 0.082, 0.704, 0.891, 0.981, 0.832, 0.353, 0.032, 0.725, 0.924, 0.985, 0.351, 0.088, 0.484, 0.739, 0.993, 0.46, 0.143, 0.898, 0.783, 0.324, 0.139, 0.521, 0.276, 0.728, 0.829, 0.953, 0.593, 0.935, 0.425, 0.679, 0.083, 0.125, 0.603, 0.031, 0.085, 0.658, 0.811, 0.113, 0.761, 0.747, 0.158, 0.216, 0.277, 0.205, 0.6, 0.111, 0.288, 0.03, 0.727, 0.995, 0.525, 0.524, 0.45, 0.603, 0.554, 0.797, 0.569, 0.425, 0.944, 0.064, 0.328, 0.351, 0.278, 0.099, 0.359, 0.424, 0.079, 0.492, 0.882, 0.182, 0.557, 0.42, 0.61, 0.189, 0.541, 0.76, 0.871, 0.405, 0.864, 0.428, 0.193, 0.449, 0.722, 0.403, 0.679, 0.477, 0.722, 0.757, 0.803, 0.962, 0.326, 0.656, 0.76, 0.26, 0.356, 0.463, 0.611, 0.303, 0.02, 0.099, 0.327, 0.958, 0.608, 0.69, 0.261, 0.201, 0.779, 0.334, 0.535, 0.264, 0.253, 0.745, 0.771, 0.894, 0.786, 0.819, 0.088, 0.321, 0.561, 0.928, 0.553, 0.859, 0.78, 0.615, 0.065, 0.118, 0.154, 0.802, 0.289, 0.145, 0.488, 0.873, 0.606, 0.104, 0.029, 0.075, 0.197, 0.294, 0.677, 0.476, 0.644, 0.666, 0.379, 0.914, 0.258, 0.344, 0.594, 0.47, 0.638, 0.581, 0.846, 0.054, 0.048, 0.358, 0.411, 0.443, 0.99, 0.156, 0.926, 0.514, 0.644, 0.992, 0.608, 0.52, 0.802, 0.399, 0.206, 0.761, 0.931, 0.032, 0.693, 0.532, 0.96, 0.569, 0.362, 0.858, 0.427, 0.119, 0.168, 0.448, 0.287, 0.007, 0.334, 0.7, 0.247, 0.196, 0.014, 0.811, 0.703, 0.597, 0.856, 0.161, 0.502, 0.486, 0.744, 0.637, 0.353, 0.428, 0.549, 0.339]
global q = [0.649, 0.348, 0.874, 0.903, 0.744, 0.889, 0.943, 0.754, 0.981, 0.992, 0.458, 0.686, 0.928, 0.989, 0.694, 0.176, 0.746, 0.722, 0.943, 0.996, 0.922, 0.983, 0.635, 0.738, 0.24, 0.375, 0.697, 0.685, 0.8, 0.985, 0.959, 0.985, 0.86, 0.692, 0.773, 0.926, 0.965, 0.987, 0.973, 0.11, 0.719, 0.771, 0.998, 0.96, 0.737, 0.953, 0.901, 0.363, 0.89, 0.96, 0.678, 0.964, 0.955, 0.953, 0.601, 0.974, 0.881, 0.856, 0.323, 0.135, 0.71, 0.071, 0.695, 0.973, 0.911, 0.395, 0.913, 0.894, 0.36, 0.583, 0.47, 0.426, 0.788, 0.685, 0.571, 0.686, 0.825, 0.997, 0.546, 0.872, 0.456, 0.635, 0.893, 0.939, 0.581, 0.484, 0.961, 0.749, 0.847, 0.469, 0.902, 0.875, 0.822, 0.536, 0.418, 0.581, 0.922, 0.511, 0.622, 0.554, 0.629, 0.603, 0.758, 0.832, 0.875, 0.485, 0.989, 0.773, 0.328, 0.773, 0.874, 0.594, 0.718, 0.751, 0.778, 0.84, 0.934, 0.962, 0.411, 0.739, 0.782, 0.723, 0.506, 0.756, 0.681, 0.667, 0.768, 0.55, 0.897, 0.989, 0.963, 0.917, 0.741, 0.339, 0.872, 0.815, 0.857, 0.958, 0.427, 0.993, 0.785, 0.993, 0.788, 0.933, 0.384, 0.419, 0.677, 0.978, 0.979, 0.945, 0.995, 0.824, 0.153, 0.122, 0.39, 0.923, 0.418, 0.535, 0.799, 0.892, 0.906, 0.136, 0.808, 0.657, 0.37, 0.715, 0.91, 0.65, 0.801, 0.851, 0.534, 0.977, 0.895, 0.46, 0.625, 0.591, 0.991, 0.885, 0.852, 0.227, 0.536, 0.937, 0.942, 0.851, 0.998, 0.831, 0.944, 0.68, 0.77, 0.994, 0.76, 0.874, 0.844, 0.634, 0.422, 0.968, 0.951, 0.063, 0.765, 0.994, 0.969, 0.624, 0.944, 0.985, 0.96, 0.231, 0.337, 0.607, 0.402, 0.25, 0.426, 0.973, 0.606, 0.496, 0.551, 0.964, 0.829, 0.677, 0.875, 0.325, 0.83, 0.662, 0.89, 0.839, 0.538, 0.441, 0.809, 0.613]
global origin = 1
global destination = 50