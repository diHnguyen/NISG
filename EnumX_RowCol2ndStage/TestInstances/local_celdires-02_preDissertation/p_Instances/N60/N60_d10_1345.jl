global arcs = [1 12; 1 18; 1 27; 1 28; 1 31; 1 34; 1 37; 1 49; 1 51; 1 53; 1 57; 2 6; 2 41; 3 19; 3 23; 3 30; 3 32; 4 31; 4 34; 4 43; 4 55; 5 24; 5 60; 6 27; 6 28; 6 56; 7 3; 7 10; 7 30; 7 35; 7 46; 7 54; 7 57; 7 59; 8 4; 8 14; 8 16; 8 27; 8 31; 8 34; 8 45; 8 54; 8 55; 9 11; 9 33; 9 43; 9 55; 9 58; 10 6; 10 26; 10 27; 10 29; 10 44; 10 47; 11 8; 11 13; 11 17; 11 26; 11 30; 11 41; 11 48; 12 13; 12 20; 12 49; 12 52; 12 54; 12 60; 13 8; 13 25; 13 34; 13 36; 13 50; 13 55; 13 56; 14 6; 14 10; 14 23; 14 31; 14 35; 14 40; 14 48; 15 6; 15 19; 15 43; 15 45; 16 8; 16 14; 16 24; 16 39; 16 42; 16 52; 16 58; 17 9; 17 15; 17 20; 18 2; 18 14; 19 5; 19 6; 19 22; 19 24; 19 30; 19 48; 20 7; 20 10; 20 29; 20 33; 20 51; 20 55; 20 59; 21 5; 21 24; 22 13; 22 47; 23 27; 23 28; 23 40; 23 45; 23 51; 23 57; 24 8; 24 40; 24 48; 24 52; 24 53; 24 58; 25 2; 25 3; 25 5; 25 19; 25 29; 25 34; 25 35; 26 14; 26 15; 26 16; 26 23; 26 24; 26 27; 26 38; 26 44; 26 47; 26 49; 26 55; 27 2; 27 7; 27 15; 27 21; 27 24; 27 42; 27 53; 27 60; 28 11; 28 20; 28 26; 28 29; 28 39; 28 51; 28 57; 29 8; 29 17; 29 18; 29 22; 29 40; 29 58; 29 59; 30 32; 30 60; 31 6; 31 8; 31 33; 31 38; 31 56; 32 35; 33 2; 33 13; 33 24; 33 32; 33 44; 33 47; 33 51; 33 55; 34 7; 34 20; 34 28; 34 36; 34 44; 34 50; 35 26; 35 36; 36 5; 36 16; 36 26; 36 40; 36 41; 36 45; 36 49; 36 53; 36 56; 36 58; 37 2; 37 22; 37 59; 38 14; 38 18; 38 31; 38 41; 38 51; 38 60; 39 14; 39 21; 39 25; 39 33; 39 36; 39 45; 39 55; 39 57; 40 9; 40 11; 40 23; 40 33; 40 39; 40 50; 40 55; 41 11; 41 13; 41 21; 41 22; 41 30; 41 33; 41 37; 41 43; 41 52; 41 56; 41 58; 42 8; 42 11; 42 18; 42 26; 42 34; 42 41; 42 46; 42 51; 42 56; 42 57; 42 59; 43 32; 43 40; 43 42; 43 47; 43 50; 44 4; 44 14; 44 21; 44 37; 44 52; 44 53; 45 7; 45 22; 45 23; 45 31; 45 37; 45 39; 46 7; 46 11; 46 27; 46 38; 46 45; 46 47; 46 56; 47 3; 47 15; 47 27; 47 29; 47 44; 47 53; 48 2; 48 24; 48 35; 48 42; 49 8; 49 13; 49 14; 49 19; 49 34; 49 40; 49 41; 49 59; 50 9; 50 22; 50 25; 50 46; 51 4; 51 17; 51 36; 51 37; 51 43; 51 53; 52 2; 52 25; 52 30; 53 8; 53 21; 53 24; 53 52; 54 11; 54 35; 54 37; 54 59; 55 2; 55 3; 55 5; 55 30; 55 44; 55 49; 56 26; 57 3; 57 4; 57 38; 57 39; 58 4; 58 12; 58 16; 58 21; 58 26; 58 39; 58 50; 59 27]
global d_x = [3.0, 8.0, 6.0, 5.0, 7.0, 5.0, 8.0, 10.0, 2.0, 5.0, 9.0, 1.0, 3.0, 6.0, 2.0, 9.0, 10.0, 6.0, 2.0, 2.0, 2.0, 7.0, 8.0, 3.0, 10.0, 6.0, 2.0, 4.0, 1.0, 6.0, 10.0, 1.0, 1.0, 1.0, 1.0, 10.0, 8.0, 7.0, 3.0, 5.0, 3.0, 4.0, 7.0, 2.0, 2.0, 5.0, 10.0, 7.0, 5.0, 5.0, 4.0, 7.0, 7.0, 4.0, 5.0, 4.0, 9.0, 9.0, 5.0, 9.0, 6.0, 4.0, 6.0, 3.0, 4.0, 2.0, 4.0, 2.0, 8.0, 4.0, 7.0, 5.0, 7.0, 9.0, 4.0, 5.0, 2.0, 10.0, 5.0, 6.0, 2.0, 9.0, 9.0, 2.0, 8.0, 6.0, 6.0, 4.0, 5.0, 1.0, 7.0, 5.0, 7.0, 9.0, 8.0, 5.0, 5.0, 1.0, 6.0, 5.0, 9.0, 5.0, 8.0, 4.0, 2.0, 7.0, 8.0, 2.0, 7.0, 2.0, 5.0, 7.0, 6.0, 8.0, 3.0, 3.0, 7.0, 1.0, 1.0, 4.0, 8.0, 4.0, 8.0, 5.0, 9.0, 5.0, 4.0, 8.0, 10.0, 3.0, 10.0, 6.0, 4.0, 3.0, 5.0, 8.0, 4.0, 9.0, 9.0, 6.0, 10.0, 7.0, 9.0, 3.0, 10.0, 9.0, 5.0, 10.0, 4.0, 10.0, 9.0, 6.0, 1.0, 2.0, 2.0, 3.0, 9.0, 5.0, 4.0, 3.0, 7.0, 4.0, 7.0, 10.0, 2.0, 2.0, 3.0, 10.0, 6.0, 2.0, 2.0, 4.0, 2.0, 8.0, 2.0, 9.0, 5.0, 10.0, 7.0, 2.0, 5.0, 10.0, 10.0, 5.0, 1.0, 1.0, 9.0, 8.0, 3.0, 9.0, 4.0, 10.0, 2.0, 1.0, 6.0, 7.0, 2.0, 10.0, 8.0, 7.0, 5.0, 2.0, 9.0, 6.0, 6.0, 1.0, 2.0, 10.0, 6.0, 9.0, 2.0, 7.0, 7.0, 10.0, 9.0, 3.0, 7.0, 8.0, 9.0, 9.0, 3.0, 9.0, 10.0, 9.0, 8.0, 3.0, 7.0, 4.0, 6.0, 5.0, 9.0, 4.0, 8.0, 4.0, 10.0, 4.0, 9.0, 6.0, 2.0, 3.0, 9.0, 7.0, 8.0, 5.0, 2.0, 10.0, 5.0, 2.0, 2.0, 8.0, 5.0, 8.0, 9.0, 1.0, 8.0, 6.0, 3.0, 10.0, 8.0, 3.0, 6.0, 8.0, 4.0, 5.0, 6.0, 1.0, 7.0, 6.0, 3.0, 5.0, 10.0, 2.0, 4.0, 3.0, 6.0, 9.0, 6.0, 7.0, 3.0, 9.0, 4.0, 1.0, 7.0, 7.0, 6.0, 5.0, 8.0, 4.0, 2.0, 1.0, 4.0, 10.0, 10.0, 8.0, 9.0, 2.0, 7.0, 1.0, 2.0, 6.0, 1.0, 5.0, 9.0, 10.0, 9.0, 1.0, 5.0, 4.0, 7.0, 6.0, 4.0, 3.0, 10.0, 1.0, 3.0, 2.0, 10.0, 1.0, 1.0, 1.0, 6.0, 6.0, 4.0, 1.0, 9.0, 9.0, 3.0, 6.0]
global b_x = 5
global d_y = [9.0, 2.0, 2.0, 1.0, 5.0, 10.0, 9.0, 8.0, 2.0, 9.0, 6.0, 5.0, 4.0, 9.0, 3.0, 8.0, 9.0, 2.0, 2.0, 10.0, 1.0, 1.0, 7.0, 9.0, 9.0, 2.0, 4.0, 7.0, 7.0, 9.0, 7.0, 7.0, 3.0, 8.0, 7.0, 6.0, 6.0, 9.0, 5.0, 2.0, 7.0, 7.0, 4.0, 5.0, 1.0, 2.0, 8.0, 5.0, 2.0, 5.0, 4.0, 2.0, 3.0, 10.0, 4.0, 6.0, 2.0, 2.0, 1.0, 8.0, 5.0, 7.0, 6.0, 4.0, 5.0, 8.0, 7.0, 1.0, 10.0, 4.0, 6.0, 8.0, 8.0, 9.0, 9.0, 4.0, 8.0, 10.0, 7.0, 6.0, 8.0, 1.0, 6.0, 7.0, 6.0, 3.0, 2.0, 10.0, 5.0, 1.0, 10.0, 7.0, 8.0, 5.0, 10.0, 2.0, 1.0, 1.0, 10.0, 1.0, 3.0, 4.0, 6.0, 9.0, 1.0, 8.0, 8.0, 9.0, 7.0, 3.0, 9.0, 5.0, 1.0, 10.0, 1.0, 9.0, 10.0, 3.0, 2.0, 4.0, 9.0, 1.0, 4.0, 6.0, 1.0, 4.0, 7.0, 3.0, 10.0, 9.0, 4.0, 10.0, 7.0, 10.0, 4.0, 2.0, 1.0, 6.0, 3.0, 4.0, 8.0, 5.0, 1.0, 6.0, 7.0, 9.0, 3.0, 6.0, 4.0, 6.0, 6.0, 9.0, 8.0, 1.0, 1.0, 5.0, 5.0, 1.0, 6.0, 4.0, 5.0, 5.0, 8.0, 2.0, 8.0, 6.0, 7.0, 8.0, 9.0, 6.0, 5.0, 3.0, 10.0, 1.0, 7.0, 1.0, 3.0, 6.0, 3.0, 2.0, 3.0, 4.0, 7.0, 4.0, 4.0, 9.0, 5.0, 8.0, 2.0, 5.0, 4.0, 10.0, 5.0, 2.0, 10.0, 5.0, 7.0, 5.0, 1.0, 3.0, 7.0, 6.0, 5.0, 6.0, 7.0, 6.0, 10.0, 9.0, 1.0, 1.0, 9.0, 5.0, 1.0, 5.0, 3.0, 4.0, 4.0, 5.0, 2.0, 3.0, 7.0, 1.0, 9.0, 9.0, 5.0, 6.0, 6.0, 5.0, 2.0, 9.0, 7.0, 7.0, 2.0, 9.0, 8.0, 4.0, 9.0, 6.0, 7.0, 10.0, 3.0, 4.0, 7.0, 3.0, 1.0, 6.0, 8.0, 2.0, 6.0, 7.0, 4.0, 5.0, 7.0, 6.0, 6.0, 3.0, 9.0, 2.0, 1.0, 9.0, 9.0, 3.0, 6.0, 7.0, 5.0, 2.0, 10.0, 5.0, 9.0, 6.0, 3.0, 10.0, 4.0, 1.0, 4.0, 8.0, 9.0, 7.0, 1.0, 6.0, 7.0, 3.0, 1.0, 10.0, 9.0, 5.0, 5.0, 4.0, 2.0, 2.0, 3.0, 7.0, 7.0, 8.0, 6.0, 2.0, 6.0, 6.0, 3.0, 6.0, 8.0, 4.0, 1.0, 7.0, 9.0, 5.0, 5.0, 5.0, 4.0, 7.0, 6.0, 10.0, 7.0, 3.0, 5.0, 7.0, 6.0, 3.0, 7.0, 7.0, 1.0, 5.0, 9.0, 7.0, 5.0, 3.0, 6.0, 6.0]
global b_y = 10
global p = [0.229, 0.309, 0.377, 0.468, 0.129, 0.71, 0.92, 0.968, 0.834, 0.518, 0.173, 0.363, 0.726, 0.472, 0.442, 0.348, 0.721, 0.513, 0.888, 0.457, 0.804, 0.29, 0.463, 0.907, 0.972, 0.707, 0.597, 0.603, 0.511, 0.412, 0.521, 0.302, 0.345, 0.163, 0.983, 0.86, 0.514, 0.899, 0.92, 0.624, 0.605, 0.004, 0.958, 0.111, 0.989, 0.025, 0.93, 0.355, 0.727, 0.69, 0.576, 0.084, 0.508, 0.254, 0.292, 0.644, 0.307, 0.358, 0.665, 0.05, 0.693, 0.128, 0.344, 0.506, 0.493, 0.381, 0.678, 0.433, 0.908, 0.877, 0.643, 0.592, 0.707, 0.143, 0.355, 0.424, 0.559, 0.357, 0.079, 0.212, 0.368, 0.474, 0.239, 0.552, 0.149, 0.909, 0.873, 0.914, 0.451, 0.259, 0.834, 0.507, 0.659, 0.066, 0.555, 0.758, 0.481, 0.241, 0.593, 0.592, 0.67, 0.055, 0.53, 0.213, 0.713, 0.668, 0.104, 0.484, 0.824, 0.793, 0.787, 0.473, 0.057, 0.151, 0.634, 0.571, 0.601, 0.171, 0.764, 0.014, 0.312, 0.657, 0.519, 0.977, 0.034, 0.603, 0.105, 0.963, 0.585, 0.794, 0.622, 0.565, 0.903, 0.522, 0.312, 0.552, 0.341, 0.452, 0.03, 0.975, 0.942, 0.741, 0.292, 0.325, 0.857, 0.222, 0.82, 0.526, 0.95, 0.908, 0.285, 0.158, 0.93, 0.15, 0.643, 0.88, 0.894, 0.129, 0.537, 0.49, 0.316, 0.039, 0.848, 0.606, 0.782, 0.957, 0.462, 0.788, 0.225, 0.558, 0.384, 0.702, 0.102, 0.818, 0.11, 0.887, 0.101, 0.95, 0.555, 0.286, 0.422, 0.504, 0.589, 0.106, 0.29, 0.689, 0.619, 0.575, 0.272, 0.14, 0.382, 0.105, 0.462, 0.81, 0.713, 0.547, 0.154, 0.733, 0.136, 0.106, 0.117, 0.338, 0.999, 0.426, 0.778, 0.282, 0.198, 0.989, 0.482, 0.296, 0.089, 0.378, 0.477, 0.263, 0.44, 0.051, 0.472, 0.057, 0.012, 0.345, 0.027, 0.491, 0.737, 0.89, 0.738, 0.365, 0.602, 0.187, 0.205, 0.8, 0.538, 0.843, 0.156, 0.322, 0.473, 0.618, 0.672, 0.543, 0.907, 0.313, 0.393, 0.821, 0.404, 0.027, 0.548, 0.292, 0.12, 0.791, 0.081, 0.251, 0.008, 0.969, 0.205, 0.162, 0.415, 0.164, 0.596, 0.818, 0.342, 0.269, 0.244, 0.784, 0.633, 0.302, 0.296, 0.052, 0.466, 0.426, 0.319, 0.006, 0.488, 0.938, 0.854, 0.688, 0.931, 0.913, 0.997, 0.676, 0.855, 0.149, 0.709, 0.803, 0.548, 0.364, 0.648, 0.764, 0.154, 0.346, 0.646, 0.522, 0.138, 0.157, 0.044, 0.374, 0.446, 0.413, 0.176, 0.421, 0.848, 0.47, 0.127, 0.173, 0.048, 0.065, 0.712, 0.2, 0.768, 0.727, 0.713, 0.228, 0.675, 0.061, 0.412, 0.328, 0.557, 0.467, 0.778, 0.481, 0.622, 0.517, 0.265, 0.126, 0.816, 0.327, 0.821, 0.175, 0.888, 0.583]
global q = [0.746, 0.847, 0.869, 0.85, 0.255, 0.755, 0.985, 0.985, 0.907, 0.602, 0.881, 0.592, 0.863, 0.626, 0.947, 0.539, 0.76, 0.697, 0.957, 0.806, 0.907, 0.89, 0.53, 0.993, 0.978, 0.771, 0.986, 0.802, 0.984, 0.859, 0.522, 0.407, 0.496, 0.532, 0.998, 0.902, 0.606, 0.994, 0.932, 0.767, 0.783, 0.303, 0.992, 0.904, 0.995, 0.928, 0.939, 0.759, 0.842, 0.887, 0.776, 0.911, 0.6, 0.451, 0.328, 0.743, 0.571, 0.838, 0.947, 0.355, 0.838, 0.732, 0.535, 0.591, 0.645, 0.693, 0.974, 0.7, 0.924, 0.94, 0.766, 0.929, 0.962, 0.377, 0.86, 0.436, 0.818, 0.944, 0.95, 0.214, 0.561, 0.8, 0.426, 0.556, 0.64, 0.925, 0.948, 0.96, 0.576, 0.342, 0.954, 0.805, 0.7, 0.73, 0.863, 0.907, 0.704, 0.772, 0.799, 0.836, 0.735, 0.126, 0.633, 0.43, 0.836, 0.735, 0.38, 0.629, 0.961, 0.944, 0.836, 0.822, 0.797, 0.842, 0.903, 0.889, 0.924, 0.798, 0.937, 0.656, 0.859, 0.828, 0.875, 0.988, 0.734, 0.974, 0.931, 0.972, 0.758, 0.918, 0.736, 0.63, 0.951, 0.68, 0.85, 0.987, 0.983, 0.688, 0.798, 0.985, 0.964, 0.77, 0.369, 0.952, 0.894, 0.385, 0.907, 0.73, 0.97, 0.912, 0.926, 0.311, 0.993, 0.513, 0.79, 0.985, 0.96, 0.622, 0.753, 0.576, 0.87, 0.787, 0.858, 0.644, 0.942, 0.963, 0.466, 0.869, 0.637, 0.962, 0.634, 0.789, 0.59, 0.9, 0.81, 0.926, 0.237, 0.988, 0.771, 0.718, 0.477, 0.588, 0.812, 0.416, 0.979, 0.757, 0.825, 0.982, 0.712, 0.689, 0.691, 0.11, 0.623, 0.959, 0.749, 0.701, 0.886, 0.753, 0.317, 0.287, 0.463, 0.62, 0.999, 0.995, 0.88, 0.889, 0.256, 0.994, 0.955, 0.667, 0.602, 0.576, 0.738, 0.505, 0.795, 0.223, 0.53, 0.124, 0.475, 0.402, 0.665, 0.55, 0.769, 0.927, 0.856, 0.681, 0.648, 0.232, 0.806, 0.923, 0.565, 0.987, 0.303, 0.893, 0.967, 0.896, 0.747, 0.677, 0.965, 0.491, 0.76, 0.981, 0.769, 0.814, 0.594, 0.976, 0.669, 0.847, 0.801, 0.277, 0.203, 0.973, 0.572, 0.679, 0.8, 0.487, 0.888, 0.833, 0.42, 0.79, 0.823, 0.934, 0.938, 0.843, 0.709, 0.761, 0.901, 0.955, 0.754, 0.422, 0.637, 0.973, 0.974, 0.776, 0.942, 0.979, 0.999, 0.713, 0.918, 0.29, 0.9, 0.999, 0.642, 0.958, 0.86, 0.973, 0.954, 0.424, 0.986, 0.975, 0.443, 0.328, 0.548, 0.916, 0.673, 0.683, 0.605, 0.535, 0.932, 0.494, 0.682, 0.649, 0.493, 0.217, 0.99, 0.619, 0.852, 0.881, 0.956, 0.376, 0.77, 0.074, 0.415, 0.776, 0.834, 0.814, 0.829, 0.644, 0.759, 0.783, 0.872, 0.782, 0.951, 0.593, 0.995, 0.889, 0.939, 0.704]
global origin = 1
global destination = 60