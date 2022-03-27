global arcs = [1 9; 1 19; 1 29; 1 37; 1 50; 2 3; 2 10; 2 14; 2 18; 2 34; 2 42; 2 50; 3 6; 3 16; 3 19; 3 30; 3 34; 3 38; 3 46; 4 40; 5 6; 5 13; 5 19; 5 29; 5 31; 5 44; 5 45; 6 13; 6 34; 6 40; 7 8; 7 15; 7 35; 7 41; 8 5; 8 27; 8 44; 8 47; 8 49; 9 18; 9 32; 9 37; 10 8; 10 20; 10 22; 10 28; 10 49; 11 21; 11 26; 11 28; 11 34; 11 35; 11 41; 11 42; 12 5; 12 15; 12 34; 13 8; 13 14; 13 28; 13 35; 13 44; 13 49; 13 50; 14 13; 14 22; 14 33; 14 41; 15 2; 15 28; 15 35; 15 40; 16 21; 16 24; 16 40; 16 45; 17 5; 17 12; 17 37; 17 38; 17 49; 18 4; 18 5; 18 10; 18 15; 18 20; 18 29; 18 34; 18 38; 19 21; 19 22; 19 26; 19 31; 19 33; 19 37; 20 10; 20 14; 20 29; 20 38; 20 45; 20 49; 21 6; 21 20; 21 25; 21 29; 21 30; 21 43; 21 46; 21 48; 21 50; 22 7; 22 23; 22 30; 22 31; 22 40; 22 43; 22 49; 22 50; 23 4; 23 14; 23 18; 23 19; 23 29; 24 5; 24 9; 25 2; 25 30; 25 34; 25 36; 25 40; 26 10; 26 22; 26 25; 26 38; 27 11; 27 19; 27 33; 27 38; 27 40; 27 42; 28 6; 28 11; 28 16; 28 22; 28 24; 28 43; 28 47; 28 49; 29 12; 29 17; 29 26; 29 38; 29 43; 29 46; 30 4; 30 16; 30 24; 31 3; 31 12; 31 18; 31 32; 31 49; 32 15; 32 23; 32 36; 32 49; 33 7; 33 13; 34 15; 34 17; 34 36; 34 37; 34 47; 35 9; 35 11; 35 25; 35 34; 36 17; 36 26; 36 39; 37 13; 37 14; 37 19; 38 2; 38 24; 38 32; 38 33; 39 20; 39 25; 39 33; 39 38; 39 46; 39 49; 40 5; 40 6; 40 7; 40 26; 40 39; 41 22; 41 24; 41 31; 41 37; 42 17; 42 23; 42 25; 42 27; 42 31; 42 35; 43 2; 43 15; 43 35; 43 48; 44 17; 44 28; 44 48; 45 2; 45 7; 45 10; 45 24; 45 34; 45 37; 45 40; 45 42; 46 2; 46 19; 46 22; 46 45; 47 31; 47 33; 47 41; 48 14; 48 42; 48 43; 48 44; 48 49; 49 21; 49 24; 49 28]
global d_x = [10.0, 8.0, 7.0, 1.0, 4.0, 4.0, 1.0, 10.0, 4.0, 3.0, 1.0, 3.0, 6.0, 8.0, 1.0, 8.0, 4.0, 1.0, 8.0, 5.0, 1.0, 9.0, 5.0, 7.0, 5.0, 3.0, 7.0, 7.0, 10.0, 1.0, 8.0, 8.0, 7.0, 6.0, 1.0, 2.0, 4.0, 3.0, 10.0, 9.0, 10.0, 8.0, 9.0, 2.0, 6.0, 9.0, 9.0, 9.0, 6.0, 7.0, 1.0, 9.0, 6.0, 4.0, 4.0, 3.0, 6.0, 3.0, 6.0, 4.0, 6.0, 5.0, 10.0, 10.0, 8.0, 10.0, 1.0, 6.0, 6.0, 2.0, 9.0, 5.0, 10.0, 1.0, 8.0, 4.0, 2.0, 7.0, 9.0, 7.0, 5.0, 2.0, 3.0, 4.0, 1.0, 1.0, 7.0, 10.0, 4.0, 2.0, 3.0, 7.0, 1.0, 10.0, 4.0, 2.0, 6.0, 7.0, 5.0, 6.0, 4.0, 1.0, 8.0, 9.0, 7.0, 8.0, 7.0, 1.0, 7.0, 10.0, 6.0, 4.0, 3.0, 6.0, 4.0, 6.0, 2.0, 1.0, 9.0, 2.0, 1.0, 10.0, 7.0, 6.0, 2.0, 3.0, 9.0, 3.0, 4.0, 8.0, 10.0, 3.0, 5.0, 5.0, 6.0, 7.0, 7.0, 9.0, 1.0, 6.0, 4.0, 4.0, 5.0, 1.0, 3.0, 8.0, 9.0, 8.0, 1.0, 10.0, 9.0, 2.0, 2.0, 1.0, 9.0, 2.0, 5.0, 3.0, 5.0, 1.0, 6.0, 5.0, 3.0, 1.0, 10.0, 10.0, 10.0, 3.0, 3.0, 3.0, 1.0, 4.0, 2.0, 8.0, 1.0, 6.0, 7.0, 8.0, 6.0, 7.0, 7.0, 3.0, 3.0, 3.0, 3.0, 7.0, 7.0, 8.0, 5.0, 5.0, 9.0, 8.0, 7.0, 7.0, 1.0, 4.0, 7.0, 4.0, 10.0, 4.0, 7.0, 9.0, 4.0, 8.0, 4.0, 6.0, 6.0, 3.0, 6.0, 5.0, 9.0, 1.0, 9.0, 4.0, 7.0, 2.0, 5.0, 2.0, 6.0, 7.0, 2.0, 5.0, 4.0, 8.0, 9.0, 6.0, 9.0, 1.0, 7.0, 2.0, 7.0, 6.0, 3.0, 9.0, 6.0, 10.0, 3.0, 5.0]
global b_x = 5
global d_y = [3.0, 9.0, 4.0, 4.0, 4.0, 10.0, 5.0, 1.0, 8.0, 4.0, 10.0, 9.0, 7.0, 10.0, 4.0, 7.0, 10.0, 2.0, 2.0, 10.0, 10.0, 7.0, 8.0, 8.0, 6.0, 3.0, 10.0, 4.0, 1.0, 8.0, 10.0, 7.0, 9.0, 4.0, 2.0, 1.0, 2.0, 3.0, 9.0, 9.0, 4.0, 4.0, 1.0, 2.0, 6.0, 1.0, 4.0, 2.0, 5.0, 8.0, 1.0, 2.0, 9.0, 8.0, 9.0, 3.0, 2.0, 2.0, 3.0, 9.0, 5.0, 4.0, 2.0, 1.0, 10.0, 2.0, 1.0, 4.0, 5.0, 8.0, 1.0, 3.0, 1.0, 10.0, 1.0, 2.0, 10.0, 8.0, 7.0, 3.0, 1.0, 9.0, 7.0, 10.0, 2.0, 5.0, 10.0, 4.0, 8.0, 2.0, 3.0, 2.0, 7.0, 5.0, 6.0, 7.0, 7.0, 10.0, 4.0, 1.0, 2.0, 1.0, 5.0, 5.0, 6.0, 6.0, 8.0, 8.0, 6.0, 5.0, 3.0, 10.0, 10.0, 2.0, 6.0, 8.0, 6.0, 9.0, 8.0, 5.0, 2.0, 6.0, 5.0, 2.0, 1.0, 9.0, 5.0, 1.0, 1.0, 1.0, 9.0, 1.0, 3.0, 6.0, 1.0, 2.0, 5.0, 1.0, 3.0, 9.0, 10.0, 9.0, 9.0, 4.0, 3.0, 2.0, 8.0, 5.0, 7.0, 2.0, 10.0, 6.0, 6.0, 3.0, 5.0, 9.0, 8.0, 6.0, 1.0, 10.0, 4.0, 9.0, 5.0, 7.0, 1.0, 8.0, 2.0, 8.0, 8.0, 10.0, 10.0, 2.0, 10.0, 7.0, 5.0, 1.0, 6.0, 4.0, 5.0, 2.0, 4.0, 1.0, 6.0, 9.0, 3.0, 6.0, 2.0, 6.0, 5.0, 8.0, 5.0, 9.0, 9.0, 1.0, 10.0, 2.0, 4.0, 10.0, 3.0, 4.0, 10.0, 5.0, 3.0, 6.0, 1.0, 6.0, 2.0, 3.0, 5.0, 8.0, 7.0, 8.0, 10.0, 4.0, 2.0, 10.0, 6.0, 4.0, 10.0, 10.0, 6.0, 5.0, 7.0, 2.0, 6.0, 9.0, 10.0, 3.0, 9.0, 10.0, 5.0, 1.0, 6.0, 10.0, 8.0, 6.0, 3.0, 5.0]
global b_y = 10
global p = [0.091, 0.139, 0.28, 0.434, 0.893, 0.896, 0.329, 0.629, 0.025, 0.859, 0.176, 0.067, 0.23, 0.002, 0.86, 0.787, 0.364, 0.78, 0.559, 0.147, 0.375, 0.313, 0.239, 0.896, 0.003, 0.67, 0.142, 0.98, 0.144, 0.98, 0.343, 0.799, 0.289, 0.907, 0.523, 0.555, 0.625, 0.77, 0.998, 0.748, 0.696, 0.063, 0.589, 0.368, 0.398, 0.292, 0.642, 0.108, 0.236, 0.545, 0.325, 0.564, 0.681, 0.029, 0.538, 0.601, 0.576, 0.203, 0.556, 0.838, 0.509, 0.928, 0.247, 0.146, 0.208, 0.782, 0.727, 0.203, 0.904, 0.972, 0.698, 0.317, 0.921, 0.655, 0.675, 0.321, 0.092, 0.336, 0.506, 0.405, 0.66, 0.157, 0.204, 0.023, 0.063, 0.404, 0.156, 0.301, 0.373, 0.76, 0.57, 0.198, 0.258, 0.64, 0.635, 0.953, 0.764, 0.053, 0.586, 0.145, 0.601, 0.587, 0.745, 0.196, 0.952, 0.1, 0.633, 0.395, 0.975, 0.452, 0.813, 0.738, 0.419, 0.938, 0.695, 0.692, 0.656, 0.786, 0.848, 0.242, 0.93, 0.103, 0.76, 0.948, 0.096, 0.467, 0.361, 0.488, 0.166, 0.899, 0.48, 0.974, 0.498, 0.968, 0.108, 0.761, 0.968, 0.095, 0.906, 0.563, 0.188, 0.256, 0.263, 0.356, 0.755, 0.706, 0.964, 0.301, 0.654, 0.435, 0.248, 0.091, 0.01, 0.189, 0.614, 0.903, 0.033, 0.084, 0.41, 0.33, 0.716, 0.906, 0.938, 0.153, 0.257, 0.397, 0.839, 0.838, 0.308, 0.293, 0.879, 0.107, 0.744, 0.321, 0.463, 0.19, 0.474, 0.234, 0.403, 0.76, 0.532, 0.899, 0.168, 0.83, 0.243, 0.055, 0.005, 0.163, 0.426, 0.672, 0.126, 0.662, 0.55, 0.991, 0.315, 0.005, 0.17, 0.607, 0.013, 0.948, 0.881, 0.107, 0.445, 0.93, 0.249, 0.019, 0.31, 0.156, 0.872, 0.593, 0.315, 0.75, 0.547, 0.117, 0.653, 0.786, 0.332, 0.533, 0.803, 0.687, 0.325, 0.501, 0.174, 0.425, 0.347, 0.243, 0.849, 0.835, 0.286, 0.742, 0.115, 0.607, 0.532, 0.614, 0.187, 0.049, 0.43, 0.924]
global q = [0.127, 0.732, 0.629, 0.953, 0.996, 0.914, 0.643, 0.652, 0.592, 0.972, 0.463, 0.249, 0.868, 0.68, 0.997, 0.802, 0.562, 0.904, 0.629, 0.703, 0.642, 0.651, 0.656, 0.94, 0.287, 0.946, 0.459, 0.995, 0.797, 0.982, 0.953, 0.858, 0.867, 0.971, 0.62, 0.582, 0.742, 0.874, 0.998, 0.88, 0.762, 0.285, 0.71, 0.405, 0.435, 0.741, 0.869, 0.725, 0.281, 0.678, 0.729, 0.767, 0.906, 0.091, 0.824, 0.765, 0.793, 0.656, 0.883, 0.853, 0.777, 0.93, 0.557, 0.905, 0.646, 0.887, 0.89, 0.38, 0.923, 0.975, 0.941, 0.488, 0.987, 0.751, 0.863, 0.918, 0.94, 0.865, 0.699, 0.811, 0.939, 0.52, 0.672, 0.973, 0.762, 0.537, 0.506, 0.671, 0.723, 0.942, 0.719, 0.561, 0.691, 0.901, 0.895, 0.956, 0.98, 0.158, 0.827, 0.722, 0.614, 0.733, 0.768, 0.204, 0.974, 0.588, 0.954, 0.984, 0.99, 0.661, 0.816, 0.822, 0.89, 0.941, 0.706, 0.891, 0.786, 0.822, 0.998, 0.713, 0.979, 0.563, 0.901, 0.992, 0.815, 0.935, 0.991, 0.854, 0.32, 0.969, 0.578, 0.975, 0.701, 0.981, 0.794, 0.923, 0.982, 0.297, 0.938, 0.563, 0.391, 0.576, 0.461, 0.861, 0.967, 0.996, 0.964, 0.807, 0.939, 0.753, 0.577, 0.566, 0.169, 0.59, 0.694, 0.993, 0.096, 0.21, 0.811, 0.969, 0.833, 0.939, 0.961, 0.86, 0.899, 0.488, 0.851, 0.957, 0.411, 0.558, 0.959, 0.78, 0.863, 0.969, 0.484, 0.925, 0.822, 0.843, 0.646, 0.979, 0.638, 0.976, 0.284, 0.91, 0.618, 0.068, 0.038, 0.697, 0.721, 0.706, 0.233, 0.863, 0.722, 0.998, 0.575, 0.872, 0.46, 0.833, 0.798, 0.996, 0.901, 0.775, 0.564, 0.993, 0.577, 0.847, 0.578, 0.446, 0.987, 0.794, 0.457, 0.878, 0.987, 0.661, 0.656, 0.828, 0.463, 0.793, 0.903, 0.827, 0.953, 0.727, 0.86, 0.982, 0.809, 0.749, 0.992, 0.876, 0.917, 0.96, 0.257, 0.916, 0.731, 0.654, 0.621, 0.76, 0.971, 0.958]
global origin = 1
global destination = 50