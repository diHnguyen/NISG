global arcs = [1 6; 1 16; 1 43; 1 58; 2 16; 2 36; 2 38; 2 48; 2 59; 3 6; 3 8; 3 16; 3 19; 3 38; 3 39; 3 41; 3 44; 3 57; 4 13; 4 26; 4 35; 5 11; 5 13; 5 34; 5 35; 5 59; 6 9; 6 11; 6 12; 6 13; 6 23; 6 42; 6 47; 6 50; 6 52; 6 54; 6 55; 6 60; 7 28; 7 31; 8 5; 8 11; 8 16; 8 25; 8 43; 8 54; 9 11; 9 14; 9 23; 9 25; 9 40; 9 45; 9 58; 10 23; 10 34; 10 35; 10 43; 10 51; 11 4; 11 7; 11 8; 11 26; 11 49; 11 50; 11 53; 11 57; 12 11; 12 18; 12 49; 12 52; 13 11; 13 14; 13 21; 13 24; 13 27; 13 33; 13 54; 13 57; 14 5; 14 19; 14 20; 14 40; 15 3; 15 13; 15 19; 15 48; 15 51; 15 54; 16 20; 16 21; 16 32; 17 2; 17 34; 17 37; 17 47; 18 19; 18 24; 18 31; 18 32; 18 42; 19 22; 19 26; 19 33; 19 50; 19 52; 20 2; 20 8; 20 12; 21 27; 21 33; 21 34; 21 60; 22 11; 22 12; 22 47; 22 51; 23 2; 23 9; 23 15; 23 24; 23 41; 24 7; 24 11; 24 13; 24 16; 24 41; 24 43; 24 48; 24 49; 24 51; 25 2; 25 4; 25 12; 25 20; 25 21; 25 43; 26 9; 26 13; 26 56; 26 60; 27 10; 27 13; 27 28; 27 34; 27 37; 27 38; 27 49; 27 51; 27 60; 28 3; 28 37; 28 43; 28 49; 28 53; 28 54; 28 55; 28 56; 29 13; 29 19; 29 30; 29 36; 29 38; 29 42; 29 43; 29 51; 30 5; 30 9; 30 13; 30 46; 30 50; 30 56; 31 6; 32 7; 32 9; 32 25; 32 55; 33 15; 33 16; 33 19; 33 20; 33 27; 33 36; 33 39; 34 3; 34 5; 34 18; 34 40; 34 53; 34 54; 34 59; 34 60; 35 18; 35 21; 35 49; 36 10; 36 17; 36 25; 36 44; 37 2; 37 6; 37 8; 37 13; 37 18; 37 55; 38 8; 38 15; 38 22; 38 37; 38 55; 38 56; 39 8; 39 17; 39 45; 39 51; 39 53; 39 56; 40 8; 40 11; 40 19; 40 24; 41 5; 41 26; 41 32; 41 39; 41 43; 41 47; 42 24; 42 33; 42 37; 42 39; 42 45; 42 60; 43 17; 43 21; 43 23; 43 34; 44 5; 44 28; 44 32; 45 9; 45 11; 45 20; 45 22; 45 25; 45 31; 45 37; 45 54; 46 8; 46 28; 46 58; 47 9; 47 32; 47 33; 47 44; 47 45; 47 46; 47 58; 48 4; 48 20; 48 34; 48 35; 48 36; 48 44; 49 6; 49 8; 49 14; 49 34; 49 39; 49 42; 49 60; 50 2; 50 19; 50 27; 50 37; 50 38; 51 3; 51 4; 51 5; 51 14; 51 42; 51 52; 52 5; 52 13; 52 17; 52 18; 52 27; 52 36; 52 38; 52 41; 53 3; 53 5; 53 23; 53 31; 53 47; 53 48; 54 7; 54 51; 54 60; 55 4; 55 9; 55 12; 55 16; 55 22; 55 27; 55 29; 55 53; 56 15; 56 25; 56 60; 57 6; 57 8; 57 26; 57 37; 57 52; 58 7; 58 8; 58 9; 58 12; 58 16; 58 27; 58 30; 58 42; 58 46; 58 49; 58 50; 59 6; 59 16; 59 25; 59 37; 59 51; 59 54; 59 57]
global d_x = [6.0, 3.0, 8.0, 4.0, 1.0, 10.0, 2.0, 2.0, 7.0, 10.0, 1.0, 2.0, 3.0, 4.0, 8.0, 10.0, 1.0, 10.0, 6.0, 7.0, 6.0, 10.0, 3.0, 7.0, 3.0, 7.0, 7.0, 8.0, 4.0, 2.0, 9.0, 3.0, 7.0, 3.0, 2.0, 5.0, 9.0, 3.0, 9.0, 1.0, 8.0, 5.0, 5.0, 9.0, 2.0, 10.0, 4.0, 1.0, 9.0, 8.0, 6.0, 2.0, 2.0, 4.0, 5.0, 4.0, 8.0, 9.0, 8.0, 8.0, 2.0, 5.0, 5.0, 2.0, 1.0, 3.0, 7.0, 2.0, 4.0, 5.0, 2.0, 10.0, 7.0, 4.0, 1.0, 5.0, 1.0, 8.0, 7.0, 2.0, 2.0, 4.0, 2.0, 6.0, 10.0, 4.0, 5.0, 4.0, 3.0, 10.0, 10.0, 8.0, 6.0, 9.0, 3.0, 3.0, 4.0, 5.0, 6.0, 2.0, 5.0, 9.0, 7.0, 4.0, 4.0, 8.0, 4.0, 5.0, 10.0, 10.0, 8.0, 4.0, 3.0, 1.0, 10.0, 9.0, 9.0, 2.0, 4.0, 8.0, 2.0, 1.0, 5.0, 5.0, 3.0, 9.0, 8.0, 10.0, 7.0, 5.0, 9.0, 4.0, 8.0, 9.0, 8.0, 8.0, 8.0, 2.0, 1.0, 8.0, 3.0, 4.0, 10.0, 10.0, 8.0, 5.0, 1.0, 10.0, 1.0, 5.0, 8.0, 6.0, 5.0, 4.0, 3.0, 1.0, 3.0, 10.0, 6.0, 6.0, 1.0, 3.0, 8.0, 5.0, 3.0, 5.0, 7.0, 6.0, 7.0, 10.0, 1.0, 8.0, 5.0, 4.0, 9.0, 5.0, 6.0, 6.0, 1.0, 5.0, 3.0, 3.0, 6.0, 7.0, 5.0, 10.0, 7.0, 6.0, 10.0, 4.0, 1.0, 4.0, 3.0, 8.0, 5.0, 8.0, 1.0, 3.0, 8.0, 3.0, 5.0, 9.0, 3.0, 5.0, 1.0, 7.0, 9.0, 7.0, 8.0, 3.0, 10.0, 2.0, 8.0, 1.0, 6.0, 6.0, 2.0, 1.0, 7.0, 5.0, 3.0, 3.0, 5.0, 1.0, 7.0, 8.0, 6.0, 4.0, 4.0, 5.0, 2.0, 8.0, 1.0, 2.0, 5.0, 2.0, 10.0, 10.0, 7.0, 5.0, 3.0, 8.0, 10.0, 4.0, 5.0, 3.0, 6.0, 7.0, 9.0, 3.0, 4.0, 3.0, 2.0, 6.0, 6.0, 1.0, 7.0, 5.0, 7.0, 2.0, 3.0, 5.0, 7.0, 9.0, 3.0, 8.0, 2.0, 6.0, 8.0, 2.0, 10.0, 1.0, 8.0, 5.0, 7.0, 4.0, 7.0, 6.0, 4.0, 8.0, 6.0, 7.0, 2.0, 8.0, 3.0, 9.0, 2.0, 2.0, 9.0, 7.0, 1.0, 1.0, 8.0, 9.0, 4.0, 10.0, 4.0, 7.0, 3.0, 4.0, 6.0, 9.0, 1.0, 4.0, 9.0, 4.0, 9.0, 6.0, 6.0, 5.0, 3.0, 10.0, 5.0, 2.0, 7.0, 1.0, 10.0, 9.0, 10.0, 8.0, 9.0, 1.0, 10.0, 8.0, 3.0, 2.0, 2.0, 10.0, 4.0, 1.0, 10.0, 8.0]
global b_x = 5
global d_y = [1.0, 3.0, 3.0, 7.0, 4.0, 1.0, 4.0, 2.0, 6.0, 6.0, 7.0, 8.0, 2.0, 1.0, 2.0, 1.0, 8.0, 2.0, 2.0, 7.0, 2.0, 3.0, 4.0, 3.0, 4.0, 4.0, 7.0, 5.0, 4.0, 6.0, 4.0, 7.0, 5.0, 8.0, 6.0, 3.0, 3.0, 3.0, 5.0, 9.0, 8.0, 1.0, 2.0, 5.0, 1.0, 2.0, 3.0, 8.0, 5.0, 2.0, 7.0, 10.0, 3.0, 7.0, 7.0, 3.0, 10.0, 5.0, 2.0, 8.0, 6.0, 3.0, 10.0, 3.0, 8.0, 9.0, 8.0, 4.0, 10.0, 7.0, 6.0, 1.0, 3.0, 10.0, 10.0, 8.0, 1.0, 3.0, 3.0, 5.0, 5.0, 9.0, 2.0, 9.0, 4.0, 7.0, 9.0, 6.0, 3.0, 7.0, 1.0, 8.0, 6.0, 6.0, 5.0, 9.0, 8.0, 8.0, 2.0, 10.0, 10.0, 8.0, 4.0, 10.0, 2.0, 5.0, 9.0, 2.0, 9.0, 2.0, 2.0, 6.0, 10.0, 8.0, 10.0, 7.0, 8.0, 6.0, 6.0, 5.0, 8.0, 2.0, 8.0, 7.0, 2.0, 7.0, 3.0, 9.0, 4.0, 4.0, 9.0, 5.0, 2.0, 10.0, 3.0, 6.0, 7.0, 2.0, 10.0, 10.0, 1.0, 9.0, 6.0, 4.0, 10.0, 5.0, 1.0, 9.0, 4.0, 8.0, 1.0, 7.0, 8.0, 7.0, 9.0, 5.0, 6.0, 8.0, 3.0, 6.0, 2.0, 6.0, 10.0, 10.0, 8.0, 3.0, 8.0, 1.0, 9.0, 5.0, 5.0, 2.0, 5.0, 4.0, 10.0, 8.0, 2.0, 5.0, 1.0, 7.0, 6.0, 6.0, 2.0, 5.0, 9.0, 10.0, 9.0, 6.0, 5.0, 5.0, 10.0, 7.0, 10.0, 7.0, 8.0, 10.0, 6.0, 9.0, 2.0, 7.0, 6.0, 1.0, 10.0, 6.0, 6.0, 6.0, 6.0, 3.0, 6.0, 2.0, 4.0, 9.0, 5.0, 10.0, 2.0, 10.0, 10.0, 3.0, 3.0, 4.0, 4.0, 1.0, 2.0, 3.0, 10.0, 2.0, 7.0, 9.0, 8.0, 10.0, 4.0, 6.0, 4.0, 2.0, 9.0, 4.0, 7.0, 4.0, 8.0, 3.0, 2.0, 5.0, 7.0, 10.0, 1.0, 7.0, 3.0, 1.0, 4.0, 4.0, 7.0, 2.0, 5.0, 5.0, 2.0, 5.0, 8.0, 4.0, 8.0, 9.0, 10.0, 1.0, 2.0, 9.0, 2.0, 9.0, 4.0, 2.0, 6.0, 6.0, 7.0, 6.0, 3.0, 5.0, 6.0, 6.0, 3.0, 4.0, 6.0, 7.0, 5.0, 5.0, 10.0, 9.0, 9.0, 2.0, 7.0, 1.0, 5.0, 10.0, 3.0, 3.0, 1.0, 6.0, 1.0, 4.0, 9.0, 8.0, 8.0, 2.0, 1.0, 7.0, 2.0, 1.0, 8.0, 6.0, 2.0, 6.0, 9.0, 1.0, 7.0, 6.0, 2.0, 10.0, 4.0, 5.0, 3.0, 4.0, 10.0, 10.0, 8.0, 8.0, 1.0, 9.0, 1.0, 2.0, 1.0, 6.0, 6.0, 7.0, 8.0, 9.0]
global b_y = 10
global p = [0.055, 0.3, 0.076, 0.958, 0.672, 0.456, 0.957, 0.426, 0.197, 0.495, 0.003, 0.248, 0.938, 0.976, 0.273, 0.451, 0.583, 0.244, 0.034, 0.423, 0.323, 0.411, 0.362, 0.216, 0.828, 0.881, 0.609, 0.318, 0.995, 0.028, 0.376, 0.448, 0.619, 0.475, 0.28, 0.692, 0.461, 0.928, 0.535, 0.517, 0.089, 0.205, 0.215, 0.943, 0.607, 0.114, 0.426, 0.847, 0.416, 0.874, 0.533, 0.078, 0.142, 0.517, 0.343, 0.245, 0.104, 0.986, 0.258, 0.634, 0.88, 0.741, 0.038, 0.24, 0.24, 0.41, 0.403, 0.176, 0.598, 0.351, 0.145, 0.242, 0.832, 0.951, 0.203, 0.696, 0.048, 0.971, 0.854, 0.38, 0.391, 0.048, 0.613, 0.511, 0.445, 0.109, 0.734, 0.521, 0.211, 0.215, 0.092, 0.958, 0.709, 0.233, 0.216, 0.667, 0.53, 0.05, 0.413, 0.683, 0.184, 0.476, 0.827, 0.754, 0.881, 0.626, 0.143, 0.138, 0.802, 0.061, 0.24, 0.837, 0.086, 0.914, 0.259, 0.946, 0.644, 0.714, 0.96, 0.843, 0.812, 0.335, 0.092, 0.46, 0.426, 0.813, 0.605, 0.758, 0.463, 0.748, 0.682, 0.907, 0.205, 0.284, 0.104, 0.641, 0.144, 0.73, 0.204, 0.623, 0.329, 0.151, 0.537, 0.689, 0.272, 0.33, 0.615, 0.561, 0.544, 0.772, 0.249, 0.379, 0.076, 0.75, 0.866, 0.899, 0.297, 0.758, 0.033, 0.216, 0.923, 0.996, 0.235, 0.44, 0.587, 0.027, 0.512, 0.978, 0.71, 0.349, 0.81, 0.249, 0.25, 0.78, 0.103, 0.854, 0.043, 0.743, 0.843, 0.287, 0.373, 0.875, 0.869, 0.474, 0.942, 0.88, 0.329, 0.882, 0.537, 0.085, 0.544, 0.244, 0.219, 0.335, 0.411, 0.573, 0.372, 0.78, 0.746, 0.175, 0.522, 0.897, 0.4, 0.658, 0.343, 0.489, 0.433, 0.963, 0.272, 0.608, 0.538, 0.61, 0.215, 0.157, 0.095, 0.853, 0.393, 0.402, 0.391, 0.601, 0.027, 0.599, 0.175, 0.212, 0.9, 0.478, 0.188, 0.485, 0.622, 0.001, 0.883, 0.012, 0.582, 0.831, 0.295, 0.917, 0.151, 0.444, 0.013, 0.606, 0.426, 0.14, 0.115, 0.301, 0.191, 0.088, 0.715, 0.477, 0.865, 0.092, 0.537, 0.68, 0.29, 0.929, 0.324, 0.468, 0.946, 0.71, 0.811, 0.253, 0.966, 0.431, 0.874, 0.693, 0.544, 0.297, 0.794, 0.142, 0.529, 0.191, 0.242, 0.279, 0.214, 0.1, 0.712, 0.398, 0.21, 0.079, 0.796, 0.74, 0.685, 0.13, 0.912, 0.112, 0.266, 0.636, 0.78, 0.642, 0.664, 0.897, 0.375, 0.701, 0.283, 0.977, 0.636, 0.214, 0.109, 0.274, 0.124, 0.88, 0.897, 0.572, 0.866, 0.861, 0.118, 0.216, 0.66, 0.501, 0.58, 0.962, 0.948, 0.492, 0.751, 0.489, 0.143, 0.937, 0.815, 0.395, 0.163, 0.928, 0.419, 0.571, 0.387, 0.845, 0.583, 0.384, 0.66, 0.909, 0.115, 0.422, 0.009, 0.442]
global q = [0.729, 0.747, 0.992, 0.995, 0.96, 0.859, 0.963, 0.92, 0.611, 0.54, 0.961, 0.769, 0.97, 0.984, 0.593, 0.794, 0.696, 0.883, 0.3, 0.59, 0.92, 0.817, 0.72, 0.43, 0.862, 0.965, 0.711, 0.922, 0.995, 0.453, 0.94, 0.625, 0.776, 0.78, 0.807, 0.942, 0.486, 0.989, 0.752, 0.568, 0.333, 0.38, 0.913, 0.986, 0.658, 0.263, 0.505, 0.869, 0.854, 0.897, 0.762, 0.107, 0.988, 0.696, 0.495, 0.303, 0.87, 0.995, 0.284, 0.819, 0.969, 0.956, 0.097, 0.703, 0.76, 0.555, 0.822, 0.654, 0.721, 0.529, 0.578, 0.831, 0.875, 0.981, 0.355, 0.924, 0.594, 0.977, 0.95, 0.698, 0.393, 0.061, 0.974, 0.79, 0.715, 0.696, 0.814, 0.972, 0.266, 0.982, 0.573, 0.974, 0.823, 0.638, 0.788, 0.687, 0.674, 0.662, 0.796, 0.804, 0.247, 0.504, 0.988, 0.894, 0.896, 0.725, 0.411, 0.534, 0.854, 0.248, 0.527, 0.902, 0.695, 0.925, 0.663, 0.979, 0.989, 0.905, 0.988, 0.937, 0.815, 0.771, 0.585, 0.658, 0.993, 0.827, 0.739, 0.845, 0.5, 0.818, 0.985, 0.947, 0.547, 0.928, 0.526, 0.655, 0.589, 0.775, 0.773, 0.737, 0.742, 0.687, 0.723, 0.77, 0.407, 0.77, 0.905, 0.764, 0.98, 0.985, 0.77, 0.761, 0.824, 0.935, 0.883, 0.997, 0.801, 0.762, 0.035, 0.356, 0.958, 0.998, 0.947, 0.98, 0.884, 0.064, 0.697, 0.978, 0.937, 0.715, 0.91, 0.504, 0.949, 0.838, 0.447, 0.985, 0.84, 0.955, 0.859, 0.514, 0.682, 0.913, 0.973, 0.653, 0.965, 0.969, 0.727, 0.984, 0.955, 0.674, 0.635, 0.992, 0.952, 0.796, 0.647, 0.749, 0.59, 0.976, 0.76, 0.618, 0.988, 0.999, 0.462, 0.871, 0.587, 0.806, 0.507, 0.983, 0.757, 0.889, 0.565, 0.792, 0.617, 0.654, 0.823, 0.945, 0.994, 0.807, 0.911, 0.709, 0.771, 0.671, 0.589, 0.785, 0.989, 0.75, 0.563, 0.625, 0.866, 0.881, 0.919, 0.86, 0.858, 0.869, 0.382, 0.963, 0.765, 0.736, 0.853, 0.901, 0.447, 0.33, 0.573, 0.578, 0.523, 0.837, 0.941, 0.797, 0.941, 0.567, 0.886, 0.695, 0.306, 0.983, 0.588, 0.965, 0.992, 0.745, 0.94, 0.861, 0.981, 0.633, 0.964, 0.949, 0.544, 0.939, 0.842, 0.687, 0.795, 0.406, 0.845, 0.81, 0.528, 0.157, 0.927, 0.543, 0.407, 0.772, 0.975, 0.762, 0.867, 0.624, 0.944, 0.139, 0.451, 0.805, 0.903, 0.642, 0.963, 0.976, 0.803, 0.779, 0.879, 0.988, 0.667, 0.389, 0.994, 0.84, 0.65, 0.898, 0.907, 0.824, 0.876, 0.877, 0.629, 0.986, 0.988, 0.761, 0.902, 0.985, 0.953, 0.777, 0.838, 0.52, 0.451, 0.995, 0.826, 0.693, 0.643, 0.968, 0.486, 0.836, 0.987, 0.901, 0.965, 0.932, 0.696, 0.984, 0.618, 0.462, 0.921, 0.445]
global origin = 1
global destination = 60