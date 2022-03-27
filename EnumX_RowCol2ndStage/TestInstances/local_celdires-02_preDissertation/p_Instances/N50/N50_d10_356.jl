global arcs = [1 10; 1 19; 1 25; 1 26; 1 27; 1 28; 1 37; 1 46; 1 49; 2 5; 2 10; 2 35; 2 39; 2 45; 3 6; 3 27; 3 29; 3 47; 4 2; 4 27; 5 12; 5 34; 5 47; 6 11; 6 13; 6 19; 6 36; 7 14; 7 15; 7 27; 7 29; 7 38; 7 43; 7 46; 7 49; 8 3; 8 15; 8 16; 8 38; 8 44; 9 8; 9 10; 9 16; 9 21; 9 30; 9 34; 10 17; 10 28; 10 30; 10 43; 10 46; 11 13; 11 34; 11 43; 12 5; 12 8; 12 9; 12 10; 12 17; 12 22; 12 28; 12 35; 12 36; 13 4; 13 6; 13 11; 13 15; 13 17; 13 23; 13 39; 14 26; 14 30; 14 31; 14 38; 14 39; 15 4; 15 40; 15 41; 15 47; 16 8; 16 25; 16 30; 16 34; 16 37; 16 45; 16 49; 17 3; 17 9; 17 14; 17 34; 17 40; 18 8; 18 28; 18 33; 19 6; 19 11; 19 15; 19 21; 19 27; 19 49; 20 4; 20 6; 20 8; 20 15; 20 21; 20 28; 20 32; 20 42; 20 44; 20 46; 21 6; 21 10; 21 22; 21 38; 22 16; 22 36; 22 46; 22 50; 23 10; 23 24; 23 32; 23 43; 23 45; 24 8; 24 15; 24 17; 24 21; 24 25; 24 32; 25 3; 25 35; 25 40; 25 49; 26 7; 26 8; 26 9; 27 6; 27 8; 27 14; 27 23; 27 29; 27 39; 27 49; 28 9; 28 16; 28 17; 28 25; 28 37; 28 41; 28 44; 28 48; 29 5; 29 6; 29 9; 29 15; 29 22; 29 36; 29 44; 30 5; 30 10; 30 19; 30 21; 30 43; 30 48; 31 2; 31 18; 31 24; 31 50; 32 9; 32 12; 32 22; 32 27; 32 38; 32 45; 33 2; 33 6; 33 11; 33 14; 33 34; 34 3; 34 10; 34 20; 34 28; 34 47; 35 24; 35 32; 35 34; 35 46; 36 16; 36 19; 36 41; 37 3; 37 14; 37 16; 37 39; 37 47; 37 48; 37 49; 38 5; 38 7; 38 13; 38 15; 38 18; 38 39; 38 44; 39 49; 40 10; 40 29; 40 35; 40 37; 40 45; 40 47; 41 8; 41 17; 41 32; 41 43; 42 4; 42 5; 42 29; 42 32; 42 40; 42 45; 42 46; 43 15; 43 24; 43 31; 43 38; 43 41; 43 45; 43 46; 44 11; 44 24; 44 35; 44 39; 45 2; 45 4; 45 9; 45 11; 45 21; 45 24; 46 17; 46 21; 46 27; 46 29; 47 3; 47 11; 47 13; 47 36; 47 41; 47 43; 48 28; 48 37; 48 41; 48 49; 49 7; 49 28; 49 37; 49 39; 49 46; 49 47]
global d_x = [4.0, 8.0, 5.0, 2.0, 6.0, 7.0, 7.0, 8.0, 5.0, 7.0, 7.0, 3.0, 7.0, 3.0, 1.0, 1.0, 10.0, 10.0, 6.0, 5.0, 4.0, 8.0, 4.0, 6.0, 6.0, 6.0, 6.0, 6.0, 9.0, 5.0, 2.0, 1.0, 4.0, 8.0, 9.0, 4.0, 8.0, 8.0, 9.0, 10.0, 10.0, 2.0, 9.0, 3.0, 2.0, 1.0, 10.0, 4.0, 4.0, 6.0, 5.0, 8.0, 1.0, 10.0, 2.0, 10.0, 3.0, 6.0, 6.0, 10.0, 7.0, 7.0, 5.0, 8.0, 9.0, 1.0, 2.0, 8.0, 6.0, 10.0, 1.0, 7.0, 4.0, 4.0, 7.0, 5.0, 10.0, 7.0, 2.0, 9.0, 4.0, 5.0, 4.0, 3.0, 10.0, 2.0, 1.0, 6.0, 2.0, 10.0, 6.0, 5.0, 10.0, 1.0, 9.0, 9.0, 1.0, 2.0, 8.0, 1.0, 8.0, 2.0, 6.0, 3.0, 9.0, 5.0, 3.0, 8.0, 7.0, 10.0, 1.0, 2.0, 2.0, 7.0, 5.0, 3.0, 9.0, 1.0, 8.0, 5.0, 1.0, 1.0, 3.0, 1.0, 3.0, 8.0, 3.0, 1.0, 2.0, 1.0, 5.0, 2.0, 1.0, 9.0, 2.0, 7.0, 7.0, 4.0, 5.0, 1.0, 5.0, 10.0, 1.0, 6.0, 2.0, 9.0, 4.0, 4.0, 9.0, 10.0, 3.0, 8.0, 5.0, 3.0, 5.0, 10.0, 4.0, 7.0, 8.0, 8.0, 4.0, 4.0, 9.0, 3.0, 5.0, 6.0, 2.0, 6.0, 10.0, 10.0, 1.0, 7.0, 5.0, 5.0, 4.0, 8.0, 1.0, 10.0, 4.0, 6.0, 5.0, 4.0, 10.0, 4.0, 2.0, 4.0, 7.0, 10.0, 10.0, 9.0, 6.0, 1.0, 4.0, 3.0, 6.0, 4.0, 2.0, 4.0, 6.0, 3.0, 3.0, 8.0, 7.0, 2.0, 9.0, 8.0, 10.0, 8.0, 10.0, 3.0, 1.0, 1.0, 10.0, 1.0, 3.0, 7.0, 9.0, 7.0, 1.0, 7.0, 9.0, 4.0, 10.0, 1.0, 4.0, 8.0, 6.0, 2.0, 9.0, 9.0, 7.0, 6.0, 6.0, 1.0, 7.0, 7.0, 10.0, 4.0, 10.0, 6.0, 5.0, 6.0, 3.0, 9.0, 5.0, 7.0, 2.0, 10.0, 10.0, 6.0, 1.0, 6.0, 5.0, 3.0, 3.0, 6.0, 1.0, 7.0, 3.0, 3.0]
global b_x = 5
global d_y = [10.0, 2.0, 7.0, 2.0, 2.0, 2.0, 10.0, 4.0, 6.0, 8.0, 10.0, 3.0, 2.0, 2.0, 6.0, 6.0, 1.0, 7.0, 7.0, 6.0, 8.0, 4.0, 9.0, 7.0, 10.0, 3.0, 5.0, 3.0, 8.0, 7.0, 2.0, 5.0, 9.0, 9.0, 7.0, 5.0, 7.0, 8.0, 2.0, 5.0, 5.0, 10.0, 4.0, 1.0, 7.0, 6.0, 1.0, 4.0, 6.0, 6.0, 10.0, 9.0, 7.0, 2.0, 9.0, 10.0, 6.0, 2.0, 5.0, 1.0, 1.0, 9.0, 1.0, 2.0, 10.0, 4.0, 4.0, 10.0, 5.0, 6.0, 10.0, 6.0, 1.0, 6.0, 4.0, 4.0, 3.0, 10.0, 1.0, 7.0, 5.0, 1.0, 2.0, 8.0, 7.0, 6.0, 5.0, 9.0, 2.0, 7.0, 6.0, 7.0, 4.0, 3.0, 10.0, 2.0, 7.0, 8.0, 6.0, 9.0, 10.0, 8.0, 7.0, 5.0, 5.0, 6.0, 8.0, 3.0, 6.0, 10.0, 2.0, 9.0, 9.0, 3.0, 9.0, 2.0, 5.0, 3.0, 4.0, 7.0, 6.0, 8.0, 2.0, 9.0, 5.0, 6.0, 7.0, 7.0, 8.0, 9.0, 2.0, 9.0, 9.0, 10.0, 2.0, 2.0, 9.0, 4.0, 5.0, 4.0, 4.0, 9.0, 2.0, 6.0, 7.0, 8.0, 10.0, 9.0, 1.0, 9.0, 5.0, 2.0, 1.0, 2.0, 4.0, 3.0, 3.0, 3.0, 5.0, 10.0, 5.0, 2.0, 5.0, 8.0, 3.0, 1.0, 9.0, 4.0, 4.0, 5.0, 5.0, 6.0, 6.0, 6.0, 3.0, 1.0, 5.0, 2.0, 8.0, 8.0, 7.0, 5.0, 4.0, 7.0, 6.0, 2.0, 2.0, 2.0, 2.0, 7.0, 10.0, 10.0, 5.0, 6.0, 7.0, 6.0, 6.0, 1.0, 3.0, 6.0, 7.0, 1.0, 3.0, 10.0, 1.0, 5.0, 6.0, 2.0, 6.0, 8.0, 8.0, 5.0, 8.0, 2.0, 1.0, 4.0, 4.0, 4.0, 8.0, 10.0, 5.0, 3.0, 4.0, 7.0, 8.0, 4.0, 7.0, 1.0, 10.0, 9.0, 9.0, 2.0, 10.0, 5.0, 10.0, 4.0, 8.0, 10.0, 9.0, 8.0, 10.0, 1.0, 5.0, 4.0, 2.0, 9.0, 10.0, 9.0, 3.0, 5.0, 6.0, 3.0, 3.0, 2.0, 10.0, 3.0, 5.0, 6.0, 2.0, 5.0]
global b_y = 10
global p = [0.558, 0.121, 0.721, 0.81, 0.013, 0.149, 0.926, 0.174, 0.871, 0.52, 0.236, 0.523, 0.345, 0.448, 0.531, 0.617, 0.512, 0.949, 0.51, 0.984, 0.139, 0.393, 0.562, 0.738, 0.498, 0.817, 0.43, 0.859, 0.038, 0.089, 0.716, 0.011, 0.661, 0.069, 0.65, 0.037, 0.133, 0.066, 0.617, 0.313, 0.765, 0.839, 0.166, 0.365, 0.968, 0.474, 0.557, 0.296, 0.461, 0.295, 0.386, 0.814, 0.24, 0.814, 0.575, 0.585, 0.413, 0.513, 0.107, 0.78, 0.559, 0.199, 0.36, 0.379, 0.983, 0.502, 0.046, 0.466, 0.719, 0.278, 0.422, 0.413, 0.127, 0.711, 0.902, 0.559, 0.629, 0.961, 0.706, 0.539, 0.527, 0.673, 0.567, 0.902, 0.177, 0.228, 0.758, 0.767, 0.941, 0.378, 0.918, 0.45, 0.479, 0.04, 0.175, 0.096, 0.198, 0.126, 0.014, 0.618, 0.178, 0.067, 0.642, 0.828, 0.027, 0.215, 0.228, 0.786, 0.427, 0.983, 0.469, 0.937, 0.83, 0.762, 0.552, 0.701, 0.403, 0.793, 0.956, 0.988, 0.633, 0.329, 0.544, 0.4, 0.945, 0.503, 0.6, 0.024, 0.144, 0.124, 0.327, 0.786, 0.896, 0.533, 0.49, 0.356, 0.443, 0.242, 0.308, 0.86, 0.846, 0.01, 0.011, 0.938, 0.935, 0.01, 0.532, 0.683, 0.3, 0.642, 0.219, 0.78, 0.419, 0.186, 0.961, 0.546, 0.368, 0.717, 0.09, 0.529, 0.389, 0.67, 0.635, 0.585, 0.67, 0.277, 0.272, 0.99, 0.139, 0.382, 0.633, 0.363, 0.623, 0.435, 0.429, 0.346, 0.038, 0.224, 0.583, 0.318, 0.561, 0.869, 0.083, 0.168, 0.95, 0.102, 0.11, 0.694, 0.226, 0.221, 0.173, 0.72, 0.37, 0.488, 0.601, 0.976, 0.232, 0.232, 0.886, 0.989, 0.036, 0.948, 0.861, 0.698, 0.614, 0.389, 0.057, 0.921, 0.263, 0.606, 0.742, 0.809, 0.726, 0.705, 0.975, 0.605, 0.93, 0.058, 0.519, 0.012, 0.369, 0.462, 0.078, 0.151, 0.092, 0.934, 0.173, 0.127, 0.196, 0.804, 0.305, 0.178, 0.759, 0.74, 0.859, 0.028, 0.221, 0.911, 0.342, 0.376, 0.32, 0.963, 0.479, 0.397, 0.905, 0.461, 0.468, 0.442, 0.87, 0.774, 0.036, 0.475, 0.361, 0.21, 0.041, 0.875, 0.704, 0.583, 0.615, 0.21]
global q = [0.694, 0.539, 0.891, 0.842, 0.973, 0.614, 0.963, 0.432, 0.872, 0.926, 0.924, 0.798, 0.764, 0.459, 0.532, 0.993, 0.664, 0.967, 0.961, 0.995, 0.181, 0.756, 0.716, 0.82, 0.502, 0.95, 0.726, 0.998, 0.197, 0.2, 0.976, 0.602, 0.835, 0.649, 0.683, 0.652, 0.695, 0.736, 0.914, 0.722, 0.995, 0.903, 0.776, 0.556, 0.988, 0.681, 0.619, 0.312, 0.556, 0.604, 0.912, 0.859, 0.482, 0.853, 0.936, 0.933, 0.749, 0.838, 0.212, 0.961, 0.744, 0.997, 0.501, 0.398, 0.991, 0.683, 0.16, 0.573, 0.842, 0.555, 0.632, 0.461, 0.553, 0.748, 0.946, 0.936, 0.776, 0.992, 0.779, 0.643, 0.82, 0.91, 0.591, 0.925, 0.619, 0.805, 0.893, 0.788, 0.971, 0.617, 0.933, 0.825, 0.648, 0.447, 0.942, 0.366, 0.804, 0.572, 0.244, 0.788, 0.679, 0.845, 0.897, 0.94, 0.739, 0.565, 0.596, 0.853, 0.441, 0.988, 0.974, 0.947, 0.985, 0.925, 0.868, 0.95, 0.859, 0.922, 0.967, 0.998, 0.909, 0.885, 0.594, 0.518, 0.968, 0.991, 0.659, 0.203, 0.724, 0.503, 0.966, 0.804, 0.939, 0.533, 0.871, 0.8, 0.92, 0.703, 0.756, 0.899, 0.942, 0.067, 0.294, 0.951, 0.935, 0.846, 0.612, 0.963, 0.698, 0.904, 0.763, 0.891, 0.558, 0.503, 0.994, 0.623, 0.746, 0.981, 0.41, 0.595, 0.533, 0.909, 0.915, 0.847, 0.99, 0.64, 0.586, 0.999, 0.786, 0.593, 0.883, 0.696, 0.743, 0.898, 0.596, 0.902, 0.575, 0.996, 0.673, 0.699, 0.565, 0.912, 0.089, 0.45, 0.953, 0.468, 0.605, 0.941, 0.516, 0.498, 0.748, 0.832, 0.39, 0.778, 0.983, 0.998, 0.62, 0.427, 0.886, 0.995, 0.403, 0.979, 0.903, 0.843, 0.817, 0.509, 0.147, 0.922, 0.911, 0.783, 0.768, 0.949, 0.792, 0.873, 0.983, 0.971, 0.981, 0.275, 0.709, 0.867, 0.933, 0.826, 0.09, 0.209, 0.438, 0.955, 0.412, 0.316, 0.543, 0.973, 0.948, 0.316, 0.783, 0.947, 0.923, 0.563, 0.958, 0.988, 0.809, 0.833, 0.669, 0.97, 0.664, 0.467, 0.921, 0.875, 0.787, 0.98, 0.875, 0.86, 0.836, 0.622, 0.803, 0.404, 0.419, 0.9, 0.915, 0.596, 0.795, 0.26]
global origin = 1
global destination = 50