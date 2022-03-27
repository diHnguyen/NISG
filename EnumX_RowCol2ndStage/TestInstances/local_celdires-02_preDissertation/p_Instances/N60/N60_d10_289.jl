global arcs = [1 35; 1 38; 1 40; 1 42; 1 44; 2 3; 2 6; 2 9; 2 10; 2 25; 2 46; 2 47; 2 55; 3 18; 3 19; 3 22; 3 23; 3 26; 3 31; 3 44; 3 51; 4 17; 4 21; 4 25; 4 28; 4 33; 4 42; 4 55; 4 60; 5 13; 5 21; 5 23; 5 29; 5 32; 5 45; 6 5; 6 21; 6 37; 6 60; 7 21; 7 27; 7 31; 7 44; 7 49; 7 55; 8 30; 8 35; 8 37; 8 41; 9 22; 9 57; 10 44; 10 52; 11 16; 11 36; 11 39; 11 44; 11 53; 12 9; 12 14; 12 18; 12 22; 12 29; 12 38; 12 47; 13 8; 13 17; 13 22; 13 35; 13 37; 13 40; 13 45; 13 56; 13 57; 14 20; 14 23; 14 28; 14 46; 14 49; 14 50; 15 12; 15 14; 15 24; 15 37; 15 43; 15 48; 15 55; 16 28; 16 34; 16 47; 17 16; 17 31; 17 32; 17 41; 17 60; 18 54; 18 59; 19 29; 19 50; 19 53; 20 19; 20 24; 20 38; 21 2; 21 19; 21 26; 21 35; 21 44; 21 54; 22 2; 22 11; 22 26; 22 30; 22 40; 22 49; 23 15; 23 16; 23 29; 23 42; 23 58; 24 6; 24 36; 25 2; 25 16; 25 30; 25 32; 25 33; 26 7; 26 21; 26 25; 26 38; 26 45; 26 46; 27 6; 27 22; 27 33; 27 43; 27 55; 28 10; 28 22; 28 29; 28 38; 28 39; 28 47; 28 54; 29 9; 29 22; 29 33; 29 35; 29 47; 29 54; 29 56; 29 57; 30 5; 30 7; 30 10; 30 22; 30 28; 30 32; 30 45; 30 48; 31 11; 31 29; 31 33; 31 40; 31 49; 32 15; 32 16; 32 19; 32 38; 32 43; 32 57; 32 59; 33 7; 33 20; 33 24; 33 36; 33 57; 34 8; 34 22; 34 32; 34 53; 34 57; 35 21; 35 54; 36 20; 36 34; 36 45; 36 49; 37 3; 37 25; 37 33; 37 42; 37 58; 38 24; 38 34; 38 35; 38 40; 38 47; 38 52; 39 7; 39 35; 39 46; 39 56; 40 28; 40 49; 41 16; 41 28; 41 48; 41 59; 42 7; 42 11; 42 24; 42 27; 42 40; 42 48; 42 53; 42 56; 42 60; 43 7; 43 8; 43 18; 43 38; 43 40; 43 47; 43 49; 43 51; 43 56; 44 2; 44 4; 44 18; 44 47; 45 13; 45 19; 45 26; 45 33; 45 35; 45 38; 45 55; 45 56; 46 3; 46 5; 46 8; 46 34; 46 37; 46 44; 46 45; 47 6; 47 10; 47 14; 47 15; 47 16; 47 19; 47 28; 47 35; 47 49; 48 13; 48 18; 48 26; 48 47; 49 6; 49 35; 49 40; 49 46; 50 21; 50 32; 50 34; 50 53; 50 54; 50 58; 51 8; 51 10; 51 20; 51 21; 51 28; 51 36; 51 43; 51 45; 51 49; 51 53; 52 8; 52 9; 52 12; 52 18; 52 24; 52 39; 53 3; 53 16; 53 22; 53 34; 53 36; 53 46; 54 13; 54 18; 54 23; 55 13; 55 21; 55 29; 55 30; 55 36; 56 13; 56 19; 56 28; 56 43; 56 46; 56 50; 56 52; 57 9; 57 20; 57 26; 57 30; 57 39; 57 47; 58 11; 58 28; 58 29; 58 33; 58 38; 58 42; 58 45; 58 59; 59 17; 59 21; 59 24; 59 32; 59 33; 59 36]
global d_x = [1.0, 7.0, 6.0, 6.0, 5.0, 4.0, 9.0, 10.0, 7.0, 1.0, 4.0, 5.0, 10.0, 2.0, 5.0, 9.0, 8.0, 5.0, 5.0, 1.0, 8.0, 4.0, 9.0, 2.0, 7.0, 7.0, 8.0, 4.0, 2.0, 1.0, 8.0, 10.0, 10.0, 6.0, 2.0, 1.0, 4.0, 9.0, 10.0, 9.0, 10.0, 9.0, 3.0, 9.0, 4.0, 3.0, 10.0, 5.0, 9.0, 2.0, 9.0, 7.0, 10.0, 2.0, 10.0, 2.0, 9.0, 4.0, 10.0, 6.0, 3.0, 2.0, 7.0, 2.0, 1.0, 3.0, 3.0, 10.0, 10.0, 6.0, 5.0, 1.0, 7.0, 3.0, 7.0, 7.0, 10.0, 5.0, 2.0, 9.0, 1.0, 6.0, 9.0, 10.0, 6.0, 2.0, 4.0, 5.0, 4.0, 1.0, 8.0, 8.0, 4.0, 10.0, 10.0, 8.0, 3.0, 2.0, 1.0, 1.0, 7.0, 7.0, 8.0, 2.0, 4.0, 8.0, 4.0, 10.0, 1.0, 1.0, 3.0, 10.0, 9.0, 3.0, 5.0, 3.0, 6.0, 7.0, 5.0, 4.0, 10.0, 8.0, 10.0, 2.0, 10.0, 9.0, 2.0, 1.0, 5.0, 7.0, 2.0, 10.0, 6.0, 3.0, 3.0, 3.0, 1.0, 1.0, 5.0, 5.0, 6.0, 9.0, 7.0, 5.0, 1.0, 9.0, 9.0, 2.0, 8.0, 6.0, 3.0, 6.0, 2.0, 3.0, 1.0, 3.0, 2.0, 9.0, 7.0, 1.0, 4.0, 9.0, 9.0, 9.0, 2.0, 1.0, 2.0, 8.0, 10.0, 1.0, 5.0, 7.0, 3.0, 10.0, 3.0, 2.0, 10.0, 8.0, 4.0, 6.0, 8.0, 1.0, 3.0, 6.0, 7.0, 9.0, 9.0, 7.0, 10.0, 5.0, 10.0, 6.0, 4.0, 2.0, 5.0, 1.0, 3.0, 7.0, 3.0, 4.0, 9.0, 4.0, 9.0, 1.0, 2.0, 5.0, 2.0, 8.0, 2.0, 2.0, 9.0, 6.0, 10.0, 4.0, 8.0, 7.0, 3.0, 10.0, 3.0, 2.0, 10.0, 3.0, 10.0, 10.0, 9.0, 10.0, 8.0, 9.0, 4.0, 6.0, 8.0, 8.0, 4.0, 2.0, 6.0, 6.0, 8.0, 7.0, 8.0, 10.0, 3.0, 9.0, 8.0, 6.0, 7.0, 1.0, 1.0, 1.0, 2.0, 5.0, 8.0, 10.0, 4.0, 4.0, 6.0, 6.0, 3.0, 5.0, 9.0, 8.0, 7.0, 3.0, 8.0, 10.0, 3.0, 3.0, 9.0, 9.0, 10.0, 5.0, 10.0, 1.0, 6.0, 5.0, 3.0, 8.0, 1.0, 6.0, 7.0, 1.0, 7.0, 2.0, 2.0, 2.0, 6.0, 1.0, 10.0, 6.0, 3.0, 6.0, 2.0, 1.0, 7.0, 3.0, 1.0, 8.0, 5.0, 1.0, 7.0, 2.0, 2.0, 1.0, 6.0, 4.0, 4.0, 9.0, 3.0, 7.0, 4.0, 10.0, 1.0, 2.0, 4.0, 3.0, 9.0, 6.0, 5.0, 2.0, 10.0, 6.0, 2.0, 4.0, 2.0, 9.0, 10.0, 6.0, 4.0]
global b_x = 5
global d_y = [2.0, 4.0, 8.0, 1.0, 5.0, 9.0, 10.0, 5.0, 9.0, 2.0, 6.0, 1.0, 4.0, 4.0, 2.0, 4.0, 4.0, 10.0, 3.0, 1.0, 10.0, 6.0, 9.0, 4.0, 3.0, 5.0, 7.0, 10.0, 10.0, 3.0, 6.0, 4.0, 8.0, 6.0, 10.0, 6.0, 2.0, 5.0, 6.0, 7.0, 4.0, 4.0, 6.0, 1.0, 7.0, 3.0, 7.0, 5.0, 6.0, 5.0, 5.0, 8.0, 10.0, 9.0, 2.0, 3.0, 8.0, 10.0, 4.0, 6.0, 2.0, 2.0, 4.0, 10.0, 6.0, 1.0, 10.0, 9.0, 5.0, 7.0, 1.0, 2.0, 2.0, 8.0, 7.0, 5.0, 6.0, 6.0, 8.0, 5.0, 2.0, 1.0, 4.0, 8.0, 7.0, 5.0, 9.0, 4.0, 2.0, 7.0, 7.0, 3.0, 2.0, 2.0, 5.0, 7.0, 6.0, 6.0, 2.0, 3.0, 9.0, 5.0, 6.0, 3.0, 10.0, 7.0, 2.0, 9.0, 1.0, 8.0, 2.0, 7.0, 10.0, 6.0, 10.0, 8.0, 2.0, 3.0, 1.0, 4.0, 6.0, 6.0, 4.0, 8.0, 5.0, 7.0, 5.0, 9.0, 9.0, 10.0, 5.0, 9.0, 8.0, 7.0, 6.0, 3.0, 4.0, 7.0, 1.0, 3.0, 2.0, 4.0, 8.0, 7.0, 10.0, 9.0, 1.0, 7.0, 3.0, 10.0, 9.0, 1.0, 4.0, 1.0, 1.0, 5.0, 4.0, 9.0, 5.0, 9.0, 6.0, 3.0, 2.0, 9.0, 3.0, 5.0, 3.0, 10.0, 7.0, 6.0, 10.0, 6.0, 2.0, 2.0, 7.0, 6.0, 9.0, 5.0, 6.0, 3.0, 6.0, 4.0, 8.0, 8.0, 4.0, 2.0, 9.0, 2.0, 10.0, 4.0, 9.0, 2.0, 5.0, 6.0, 10.0, 1.0, 3.0, 6.0, 3.0, 7.0, 7.0, 5.0, 5.0, 10.0, 3.0, 10.0, 2.0, 10.0, 5.0, 4.0, 1.0, 5.0, 7.0, 1.0, 1.0, 7.0, 7.0, 3.0, 1.0, 10.0, 7.0, 6.0, 2.0, 9.0, 9.0, 5.0, 3.0, 7.0, 6.0, 10.0, 3.0, 6.0, 5.0, 6.0, 7.0, 2.0, 5.0, 7.0, 10.0, 2.0, 4.0, 3.0, 4.0, 9.0, 3.0, 8.0, 10.0, 7.0, 3.0, 8.0, 4.0, 6.0, 8.0, 6.0, 4.0, 1.0, 6.0, 7.0, 3.0, 8.0, 5.0, 3.0, 7.0, 4.0, 3.0, 7.0, 2.0, 5.0, 10.0, 6.0, 6.0, 2.0, 2.0, 4.0, 6.0, 6.0, 2.0, 7.0, 9.0, 3.0, 1.0, 1.0, 6.0, 5.0, 2.0, 8.0, 3.0, 6.0, 6.0, 4.0, 1.0, 5.0, 8.0, 2.0, 8.0, 6.0, 2.0, 6.0, 7.0, 9.0, 4.0, 6.0, 8.0, 8.0, 2.0, 7.0, 4.0, 2.0, 2.0, 9.0, 10.0, 4.0, 6.0, 1.0, 3.0, 10.0, 3.0, 1.0, 8.0, 2.0, 3.0, 4.0, 3.0, 1.0, 6.0, 8.0, 6.0]
global b_y = 10
global p = [0.713, 0.731, 0.744, 0.725, 0.801, 0.753, 0.83, 0.452, 0.162, 0.131, 0.296, 0.553, 0.66, 0.644, 0.485, 0.932, 0.176, 0.426, 0.583, 0.746, 0.887, 0.537, 0.195, 0.59, 0.317, 0.842, 0.056, 0.069, 0.823, 0.209, 0.935, 0.379, 0.57, 0.825, 0.099, 0.636, 0.228, 0.285, 0.634, 0.517, 0.434, 0.417, 0.929, 0.983, 0.254, 0.428, 0.832, 0.839, 0.427, 0.814, 0.722, 0.503, 0.341, 0.622, 0.785, 0.033, 0.072, 0.365, 0.097, 0.815, 0.161, 0.102, 0.558, 0.126, 0.667, 0.208, 0.514, 0.718, 0.901, 0.052, 0.314, 0.806, 0.779, 0.122, 0.91, 0.216, 0.699, 0.535, 0.536, 0.613, 0.511, 0.725, 0.258, 0.449, 0.398, 0.554, 0.806, 0.162, 0.331, 0.172, 0.185, 0.11, 0.397, 0.532, 0.021, 0.805, 0.831, 0.237, 0.103, 0.98, 0.738, 0.681, 0.358, 0.216, 0.902, 0.56, 0.516, 0.542, 0.273, 0.729, 0.319, 0.425, 0.31, 0.276, 0.529, 0.144, 0.29, 0.27, 0.7, 0.012, 0.793, 0.918, 0.756, 0.542, 0.03, 0.994, 0.4, 0.856, 0.093, 0.99, 0.683, 0.008, 0.566, 0.86, 0.436, 0.552, 0.183, 0.066, 0.203, 0.881, 0.247, 0.505, 0.45, 0.887, 0.984, 0.068, 0.162, 0.705, 0.041, 0.814, 0.927, 0.304, 0.428, 0.416, 0.086, 0.736, 0.685, 0.586, 0.853, 0.223, 0.547, 0.997, 0.103, 0.894, 0.4, 0.859, 0.028, 0.544, 0.862, 0.382, 0.49, 0.161, 0.155, 0.22, 0.895, 0.197, 0.549, 0.017, 0.544, 0.918, 0.305, 0.887, 0.961, 0.086, 0.055, 0.728, 0.629, 0.003, 0.438, 0.847, 0.319, 0.877, 0.632, 0.148, 0.648, 0.169, 0.443, 0.74, 0.021, 0.832, 0.565, 0.914, 0.282, 0.259, 0.943, 0.172, 0.996, 0.575, 0.697, 0.973, 0.558, 0.196, 0.444, 0.016, 0.267, 0.621, 0.777, 0.524, 0.344, 0.378, 0.061, 0.974, 0.007, 0.591, 0.363, 0.602, 0.699, 0.179, 0.366, 0.363, 0.503, 0.478, 0.692, 0.308, 0.838, 0.341, 0.104, 0.812, 0.015, 0.606, 0.805, 0.102, 0.675, 0.325, 0.797, 0.244, 0.614, 0.57, 0.372, 0.473, 0.953, 0.018, 0.614, 0.511, 0.09, 0.943, 0.67, 0.524, 0.739, 0.671, 0.983, 0.362, 0.985, 0.758, 0.681, 0.016, 0.3, 0.427, 0.047, 0.962, 0.704, 0.778, 0.647, 0.189, 0.724, 0.443, 0.606, 0.991, 0.171, 0.391, 0.98, 0.879, 0.043, 0.374, 0.766, 0.618, 0.988, 0.458, 0.28, 0.274, 0.115, 0.247, 0.231, 0.121, 0.556, 0.78, 0.397, 0.416, 0.786, 0.012, 0.963, 0.669, 0.44, 0.965, 0.58, 0.965, 0.995, 0.149, 0.777, 0.194, 0.136, 0.537, 0.015, 0.432, 0.939, 0.849, 0.318, 0.619, 0.099, 0.958, 0.814, 0.756, 0.511, 0.915, 0.393, 0.646, 0.51]
global q = [0.911, 0.732, 0.92, 0.987, 0.831, 0.881, 0.933, 0.803, 0.272, 0.976, 0.384, 0.92, 0.808, 0.692, 0.55, 0.996, 0.357, 0.906, 0.749, 0.845, 0.994, 0.795, 0.539, 0.961, 0.699, 0.901, 0.602, 0.625, 0.971, 0.583, 0.97, 0.748, 0.694, 0.947, 0.899, 0.724, 0.482, 0.933, 0.875, 0.713, 0.687, 0.438, 0.936, 0.984, 0.545, 0.825, 0.865, 0.956, 0.562, 0.908, 0.972, 0.717, 0.963, 0.794, 0.96, 0.151, 0.593, 0.843, 0.203, 0.987, 0.866, 0.186, 0.997, 0.147, 0.945, 0.504, 0.689, 0.79, 0.984, 0.478, 0.857, 0.912, 0.997, 0.284, 0.998, 0.401, 0.947, 0.568, 0.927, 0.813, 0.915, 0.957, 0.968, 0.8, 0.632, 0.702, 0.861, 0.212, 0.779, 0.737, 0.41, 0.248, 0.847, 0.771, 0.66, 0.963, 0.873, 0.78, 0.152, 0.993, 0.766, 0.685, 0.811, 0.779, 0.902, 0.87, 0.94, 0.895, 0.273, 0.773, 0.425, 0.968, 0.522, 0.935, 0.562, 0.551, 0.754, 0.48, 0.74, 0.837, 0.885, 0.935, 0.868, 0.644, 0.36, 0.996, 0.97, 0.948, 0.951, 0.995, 0.694, 0.357, 0.871, 0.906, 0.964, 0.736, 0.95, 0.097, 0.642, 0.977, 0.55, 0.564, 0.588, 0.959, 0.994, 0.087, 0.372, 0.755, 0.241, 0.87, 0.979, 0.617, 0.478, 0.894, 0.787, 0.83, 0.878, 0.884, 0.996, 0.496, 0.936, 0.997, 0.365, 0.912, 0.873, 0.896, 0.457, 0.616, 0.944, 0.482, 0.969, 0.379, 0.551, 0.449, 0.954, 0.451, 0.777, 0.361, 0.954, 0.955, 0.93, 0.971, 0.968, 0.726, 0.48, 0.793, 0.881, 0.169, 0.765, 0.914, 0.737, 0.979, 0.853, 0.185, 0.744, 0.202, 0.748, 0.76, 0.156, 0.923, 0.721, 0.914, 0.667, 0.387, 0.986, 0.777, 0.996, 0.678, 0.792, 0.98, 0.949, 0.688, 0.871, 0.716, 0.694, 0.986, 0.91, 0.807, 0.966, 0.907, 0.783, 0.978, 0.256, 0.817, 0.824, 0.666, 0.802, 0.642, 0.639, 0.446, 0.67, 0.708, 0.711, 0.817, 0.978, 0.583, 0.803, 0.82, 0.2, 0.681, 0.822, 0.261, 0.973, 0.583, 0.924, 0.785, 0.967, 0.982, 0.716, 0.604, 0.989, 0.703, 0.909, 0.617, 0.346, 0.965, 0.671, 0.595, 0.857, 0.717, 0.989, 0.575, 0.995, 0.95, 0.827, 0.816, 0.659, 0.68, 0.555, 0.989, 0.816, 0.805, 0.702, 0.959, 0.792, 0.934, 0.9, 0.995, 0.853, 0.513, 0.999, 0.895, 0.36, 0.696, 0.918, 0.734, 0.991, 0.577, 0.632, 0.805, 0.572, 0.346, 0.313, 0.841, 0.576, 0.786, 0.936, 0.925, 0.903, 0.164, 0.983, 0.738, 0.523, 0.968, 0.761, 0.999, 0.998, 0.867, 0.915, 0.738, 0.549, 0.542, 0.793, 0.765, 0.951, 0.856, 0.328, 0.837, 0.857, 0.961, 0.899, 0.866, 0.904, 0.99, 0.874, 0.934, 0.947]
global origin = 1
global destination = 60