global arcs = [1 42; 1 55; 2 15; 2 21; 2 24; 2 36; 2 38; 2 51; 3 4; 3 17; 3 21; 3 23; 3 27; 3 28; 3 29; 3 50; 4 13; 4 20; 4 21; 4 27; 4 43; 4 49; 5 36; 5 42; 6 8; 6 11; 6 15; 6 22; 6 24; 6 25; 6 29; 6 33; 6 44; 6 46; 6 55; 7 2; 7 10; 7 12; 7 21; 7 34; 7 37; 7 38; 7 42; 8 2; 8 3; 8 4; 8 15; 8 25; 8 30; 8 36; 8 44; 8 46; 8 49; 8 51; 8 54; 9 19; 9 28; 9 33; 9 34; 9 35; 9 37; 9 47; 9 48; 9 51; 9 58; 10 4; 10 8; 10 13; 10 28; 10 30; 10 34; 10 41; 10 53; 10 55; 11 3; 11 5; 11 10; 11 14; 11 21; 11 25; 11 27; 11 50; 11 54; 11 56; 11 59; 12 8; 12 23; 12 29; 12 51; 12 57; 13 2; 13 10; 13 17; 13 25; 13 51; 13 58; 13 59; 14 21; 14 49; 14 58; 15 10; 15 47; 15 60; 16 8; 16 13; 16 14; 16 30; 16 32; 16 35; 16 41; 16 45; 16 48; 17 14; 17 20; 17 23; 17 30; 17 31; 17 46; 17 47; 17 54; 18 12; 18 27; 18 28; 18 44; 18 59; 19 16; 20 10; 20 19; 20 26; 20 48; 20 53; 20 60; 21 23; 21 34; 21 37; 22 13; 22 39; 23 17; 23 33; 24 3; 24 7; 24 11; 24 13; 24 14; 24 17; 24 28; 24 30; 25 37; 25 40; 25 46; 25 60; 26 12; 26 14; 26 20; 26 21; 26 22; 26 36; 26 60; 27 11; 27 32; 28 4; 28 23; 28 27; 28 37; 28 39; 28 51; 28 55; 28 59; 29 9; 29 13; 29 17; 29 27; 29 35; 29 39; 29 40; 29 42; 29 43; 29 53; 29 56; 30 15; 30 16; 30 20; 30 33; 30 35; 30 52; 30 54; 30 55; 31 43; 31 58; 32 7; 32 56; 32 58; 33 8; 33 13; 33 18; 33 24; 33 36; 33 42; 33 51; 34 30; 34 31; 34 38; 34 55; 35 2; 35 17; 35 28; 35 32; 35 33; 35 36; 35 56; 35 58; 36 3; 36 14; 36 19; 36 32; 36 33; 36 45; 36 51; 36 58; 37 17; 38 8; 38 10; 38 24; 38 30; 38 41; 38 47; 38 54; 38 55; 39 4; 39 12; 39 31; 39 33; 39 34; 39 37; 39 44; 40 16; 40 27; 40 41; 40 60; 41 6; 41 10; 41 14; 41 37; 41 43; 42 10; 42 32; 42 45; 42 54; 43 18; 43 27; 43 33; 43 37; 43 46; 43 49; 43 55; 44 9; 44 28; 44 39; 44 40; 44 47; 44 50; 45 2; 45 7; 45 30; 45 34; 45 49; 45 54; 45 59; 46 5; 46 60; 47 2; 47 11; 47 18; 47 26; 47 37; 47 48; 47 51; 47 59; 48 10; 48 27; 48 29; 48 31; 48 35; 48 42; 48 57; 49 9; 49 16; 49 33; 50 5; 50 10; 50 13; 50 19; 50 30; 50 32; 50 54; 51 3; 51 16; 51 28; 51 50; 52 4; 52 25; 52 35; 52 57; 53 6; 53 30; 53 31; 53 42; 54 21; 54 28; 54 35; 55 12; 55 22; 56 47; 57 19; 57 20; 57 45; 57 50; 57 56; 58 3; 58 47; 59 11; 59 15; 59 33]
global d_x = [8.0, 9.0, 4.0, 9.0, 3.0, 4.0, 4.0, 2.0, 10.0, 10.0, 1.0, 4.0, 3.0, 3.0, 4.0, 3.0, 9.0, 7.0, 10.0, 4.0, 8.0, 5.0, 6.0, 8.0, 2.0, 4.0, 1.0, 9.0, 8.0, 1.0, 10.0, 6.0, 1.0, 4.0, 4.0, 3.0, 1.0, 9.0, 2.0, 9.0, 3.0, 6.0, 8.0, 1.0, 6.0, 1.0, 8.0, 3.0, 3.0, 9.0, 8.0, 1.0, 1.0, 9.0, 3.0, 2.0, 5.0, 5.0, 6.0, 9.0, 7.0, 2.0, 1.0, 4.0, 8.0, 5.0, 5.0, 10.0, 10.0, 6.0, 4.0, 10.0, 4.0, 4.0, 10.0, 7.0, 7.0, 7.0, 1.0, 2.0, 6.0, 5.0, 5.0, 7.0, 8.0, 9.0, 5.0, 4.0, 3.0, 1.0, 7.0, 7.0, 5.0, 3.0, 1.0, 5.0, 3.0, 3.0, 2.0, 2.0, 10.0, 10.0, 2.0, 5.0, 6.0, 6.0, 5.0, 9.0, 7.0, 10.0, 9.0, 5.0, 8.0, 4.0, 7.0, 3.0, 2.0, 10.0, 10.0, 9.0, 5.0, 1.0, 8.0, 8.0, 3.0, 5.0, 8.0, 8.0, 8.0, 2.0, 5.0, 9.0, 2.0, 9.0, 8.0, 8.0, 6.0, 10.0, 3.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 1.0, 7.0, 5.0, 9.0, 8.0, 2.0, 8.0, 8.0, 6.0, 8.0, 4.0, 8.0, 3.0, 7.0, 3.0, 1.0, 9.0, 10.0, 10.0, 5.0, 1.0, 8.0, 8.0, 3.0, 8.0, 5.0, 4.0, 4.0, 10.0, 6.0, 1.0, 2.0, 8.0, 1.0, 2.0, 10.0, 8.0, 10.0, 9.0, 10.0, 4.0, 10.0, 8.0, 6.0, 4.0, 2.0, 7.0, 5.0, 8.0, 7.0, 10.0, 2.0, 1.0, 8.0, 8.0, 10.0, 2.0, 1.0, 10.0, 7.0, 10.0, 10.0, 8.0, 4.0, 4.0, 5.0, 3.0, 5.0, 10.0, 6.0, 10.0, 7.0, 7.0, 3.0, 10.0, 4.0, 7.0, 3.0, 5.0, 4.0, 9.0, 5.0, 1.0, 1.0, 4.0, 10.0, 7.0, 1.0, 7.0, 10.0, 8.0, 7.0, 6.0, 10.0, 3.0, 7.0, 8.0, 10.0, 10.0, 10.0, 1.0, 6.0, 3.0, 5.0, 8.0, 9.0, 2.0, 10.0, 2.0, 9.0, 4.0, 5.0, 4.0, 7.0, 1.0, 4.0, 10.0, 4.0, 5.0, 4.0, 8.0, 8.0, 1.0, 5.0, 3.0, 10.0, 1.0, 5.0, 9.0, 8.0, 8.0, 10.0, 4.0, 3.0, 3.0, 8.0, 2.0, 10.0, 7.0, 3.0, 2.0, 1.0, 2.0, 8.0, 2.0, 3.0, 8.0, 10.0, 4.0, 10.0, 3.0, 5.0, 9.0, 1.0, 8.0, 2.0, 1.0, 1.0, 6.0, 9.0, 4.0, 8.0, 1.0, 9.0, 7.0, 1.0, 1.0, 5.0, 6.0, 10.0, 8.0, 10.0, 10.0, 7.0, 1.0, 3.0, 5.0, 10.0]
global b_x = 5
global d_y = [2.0, 5.0, 8.0, 6.0, 7.0, 2.0, 7.0, 9.0, 9.0, 10.0, 9.0, 7.0, 7.0, 3.0, 2.0, 7.0, 9.0, 6.0, 9.0, 5.0, 9.0, 9.0, 9.0, 2.0, 10.0, 8.0, 5.0, 8.0, 9.0, 2.0, 10.0, 1.0, 10.0, 6.0, 6.0, 10.0, 10.0, 5.0, 6.0, 1.0, 2.0, 9.0, 2.0, 5.0, 8.0, 3.0, 9.0, 1.0, 9.0, 8.0, 4.0, 4.0, 3.0, 9.0, 6.0, 7.0, 7.0, 7.0, 6.0, 2.0, 4.0, 2.0, 3.0, 6.0, 9.0, 9.0, 3.0, 4.0, 7.0, 7.0, 3.0, 4.0, 2.0, 4.0, 5.0, 10.0, 7.0, 1.0, 8.0, 3.0, 9.0, 1.0, 10.0, 2.0, 7.0, 2.0, 1.0, 3.0, 3.0, 3.0, 9.0, 8.0, 9.0, 5.0, 3.0, 6.0, 8.0, 2.0, 4.0, 5.0, 8.0, 1.0, 7.0, 3.0, 9.0, 3.0, 8.0, 4.0, 1.0, 1.0, 4.0, 6.0, 3.0, 4.0, 2.0, 7.0, 6.0, 10.0, 4.0, 4.0, 9.0, 1.0, 6.0, 10.0, 6.0, 8.0, 8.0, 8.0, 3.0, 1.0, 2.0, 4.0, 2.0, 1.0, 7.0, 7.0, 6.0, 9.0, 4.0, 1.0, 5.0, 8.0, 9.0, 9.0, 4.0, 2.0, 2.0, 7.0, 6.0, 7.0, 1.0, 3.0, 7.0, 8.0, 8.0, 1.0, 3.0, 6.0, 6.0, 5.0, 5.0, 3.0, 8.0, 9.0, 3.0, 5.0, 9.0, 6.0, 4.0, 5.0, 2.0, 7.0, 9.0, 5.0, 5.0, 8.0, 3.0, 10.0, 2.0, 5.0, 10.0, 8.0, 7.0, 2.0, 2.0, 9.0, 4.0, 9.0, 6.0, 9.0, 5.0, 3.0, 5.0, 7.0, 7.0, 10.0, 4.0, 7.0, 3.0, 1.0, 2.0, 1.0, 2.0, 6.0, 8.0, 3.0, 8.0, 9.0, 2.0, 10.0, 8.0, 10.0, 6.0, 6.0, 2.0, 2.0, 3.0, 3.0, 7.0, 3.0, 7.0, 2.0, 4.0, 1.0, 6.0, 4.0, 1.0, 8.0, 6.0, 10.0, 6.0, 7.0, 7.0, 6.0, 3.0, 5.0, 2.0, 8.0, 3.0, 7.0, 5.0, 4.0, 9.0, 2.0, 3.0, 8.0, 5.0, 8.0, 3.0, 2.0, 4.0, 5.0, 9.0, 8.0, 10.0, 7.0, 10.0, 3.0, 1.0, 3.0, 3.0, 10.0, 3.0, 2.0, 5.0, 10.0, 5.0, 9.0, 5.0, 6.0, 9.0, 9.0, 5.0, 3.0, 7.0, 8.0, 6.0, 6.0, 3.0, 7.0, 3.0, 3.0, 3.0, 10.0, 6.0, 3.0, 1.0, 9.0, 5.0, 7.0, 7.0, 5.0, 3.0, 7.0, 2.0, 2.0, 2.0, 7.0, 4.0, 6.0, 1.0, 9.0, 3.0, 3.0, 3.0, 6.0, 7.0, 9.0, 1.0, 1.0, 1.0, 10.0, 4.0, 1.0, 9.0, 8.0, 10.0, 6.0, 2.0, 8.0, 3.0, 3.0, 2.0]
global b_y = 10
global p = [0.295, 0.892, 0.57, 0.775, 0.486, 0.206, 0.843, 0.077, 0.733, 0.22, 0.952, 0.193, 0.295, 0.293, 0.995, 0.12, 0.607, 0.63, 0.911, 0.665, 0.602, 0.758, 0.214, 0.013, 0.539, 0.014, 0.736, 0.557, 0.179, 0.766, 0.129, 0.137, 0.977, 0.569, 0.046, 0.809, 0.426, 0.313, 0.981, 0.745, 0.865, 0.562, 0.777, 0.114, 0.649, 0.438, 0.898, 0.873, 0.165, 0.087, 0.597, 0.264, 0.335, 0.171, 0.534, 0.358, 0.997, 0.403, 0.054, 0.252, 0.737, 0.351, 0.487, 0.082, 0.387, 0.824, 0.529, 0.73, 0.361, 0.849, 0.552, 0.877, 0.278, 0.294, 0.219, 0.72, 0.654, 0.065, 0.239, 0.131, 0.53, 0.413, 0.608, 0.671, 0.873, 0.301, 0.051, 0.982, 0.608, 0.214, 0.09, 0.523, 0.836, 0.889, 0.644, 0.344, 0.09, 0.561, 0.375, 0.456, 0.688, 0.624, 0.775, 0.033, 0.209, 0.721, 0.879, 0.329, 0.328, 0.034, 0.556, 0.022, 0.261, 0.306, 0.985, 0.498, 0.189, 0.532, 0.49, 0.647, 0.137, 0.384, 0.908, 0.134, 0.815, 0.465, 0.714, 0.352, 0.461, 0.84, 0.589, 0.347, 0.749, 0.688, 0.284, 0.21, 0.001, 0.061, 0.375, 0.315, 0.473, 0.655, 0.058, 0.761, 0.779, 0.601, 0.789, 0.37, 0.941, 0.224, 0.312, 0.724, 0.545, 0.694, 0.746, 0.108, 0.52, 0.61, 0.658, 0.504, 0.443, 0.434, 0.516, 0.679, 0.912, 0.441, 0.049, 0.711, 0.836, 0.771, 0.529, 0.914, 0.2, 0.642, 0.982, 0.159, 0.487, 0.935, 0.186, 0.026, 0.976, 0.295, 0.995, 0.15, 0.377, 0.018, 0.072, 0.167, 0.885, 0.531, 0.169, 0.338, 0.931, 0.016, 0.757, 0.67, 0.814, 0.673, 0.391, 0.117, 0.76, 0.771, 0.971, 0.031, 0.499, 0.359, 0.693, 0.839, 0.982, 0.561, 0.183, 0.179, 0.296, 0.802, 0.927, 0.688, 0.672, 0.736, 0.725, 0.284, 0.969, 0.948, 0.792, 0.9, 0.563, 0.619, 0.368, 0.258, 0.846, 0.713, 0.656, 0.877, 0.37, 0.037, 0.825, 0.29, 0.754, 0.594, 0.852, 0.108, 0.406, 0.521, 0.208, 0.065, 0.952, 0.563, 0.046, 0.363, 0.009, 0.636, 0.843, 0.739, 0.581, 0.02, 0.72, 0.769, 0.65, 0.602, 0.46, 0.6, 0.286, 0.787, 0.961, 0.077, 0.944, 0.818, 0.539, 0.004, 0.512, 0.167, 0.558, 0.313, 0.964, 0.868, 0.522, 0.598, 0.797, 0.913, 0.109, 0.263, 0.909, 0.206, 0.911, 0.404, 0.451, 0.219, 0.795, 0.867, 0.917, 0.677, 0.404, 0.808, 0.235, 0.197, 0.025, 0.842, 0.286, 0.336, 0.473, 0.814, 0.16, 0.951, 0.716, 0.857, 0.21, 0.522, 0.711, 0.597, 0.617, 0.604, 0.615, 0.676, 0.417, 0.74, 0.611, 0.531, 0.855, 0.758, 0.94, 0.019, 0.622, 0.149, 0.523]
global q = [0.374, 0.972, 0.767, 0.859, 0.596, 0.545, 0.953, 0.711, 0.771, 0.502, 0.957, 0.725, 0.518, 0.632, 0.998, 0.4, 0.645, 0.702, 0.922, 0.694, 0.848, 0.938, 0.893, 0.283, 0.697, 0.857, 0.781, 0.674, 0.84, 0.805, 0.695, 0.677, 0.995, 0.996, 0.765, 0.878, 0.854, 0.418, 0.995, 0.881, 0.981, 0.586, 0.78, 0.725, 0.763, 0.552, 0.957, 0.922, 0.876, 0.744, 0.81, 0.993, 0.781, 0.792, 0.546, 0.364, 0.997, 0.426, 0.205, 0.757, 0.97, 0.687, 0.53, 0.243, 0.821, 0.837, 0.652, 0.908, 0.984, 0.97, 0.693, 0.971, 0.702, 0.873, 0.693, 0.72, 0.694, 0.779, 0.516, 0.9, 0.746, 0.611, 0.688, 0.986, 0.998, 0.88, 0.489, 0.984, 0.806, 0.675, 0.618, 0.841, 0.999, 0.952, 0.872, 0.857, 0.144, 0.977, 0.397, 0.518, 0.883, 0.972, 0.923, 0.449, 0.856, 0.939, 0.965, 0.873, 0.826, 0.761, 0.804, 0.491, 0.351, 0.861, 0.991, 0.751, 0.977, 0.687, 0.502, 0.809, 0.97, 0.673, 0.969, 0.285, 0.901, 0.846, 0.857, 0.783, 0.789, 0.888, 0.955, 0.867, 0.839, 0.89, 0.446, 0.916, 0.115, 0.997, 0.829, 0.348, 0.573, 0.964, 0.925, 0.853, 0.984, 0.76, 0.839, 0.842, 0.97, 0.931, 0.882, 0.792, 0.737, 0.805, 0.854, 0.677, 0.941, 0.901, 0.776, 0.996, 0.715, 0.755, 0.827, 0.725, 0.971, 0.982, 0.796, 0.893, 0.918, 0.962, 0.895, 0.962, 0.705, 0.797, 0.994, 0.836, 0.554, 0.988, 0.511, 0.453, 0.99, 0.911, 0.999, 0.37, 0.606, 0.954, 0.113, 0.995, 0.913, 0.577, 0.376, 0.498, 0.971, 0.786, 0.909, 0.974, 0.928, 0.851, 0.53, 0.487, 0.937, 0.774, 0.989, 0.366, 0.814, 0.788, 0.893, 0.913, 0.998, 0.935, 0.349, 0.363, 0.333, 0.926, 0.952, 0.958, 0.99, 0.888, 0.933, 0.656, 0.998, 0.974, 0.945, 0.927, 0.968, 0.767, 0.727, 0.639, 0.99, 0.806, 0.751, 0.974, 0.811, 0.433, 0.927, 0.808, 0.986, 0.862, 0.867, 0.92, 0.833, 0.994, 0.893, 0.955, 0.969, 0.752, 0.161, 0.891, 0.49, 0.649, 0.859, 0.897, 0.592, 0.701, 0.797, 0.875, 0.925, 0.89, 0.47, 0.797, 0.574, 0.859, 0.994, 0.513, 0.992, 0.875, 0.559, 0.801, 0.853, 0.256, 0.829, 0.396, 0.972, 0.882, 0.728, 0.797, 0.857, 0.99, 0.413, 0.486, 0.919, 0.311, 0.995, 0.7, 0.979, 0.962, 0.943, 0.998, 0.947, 0.94, 0.672, 0.937, 0.668, 0.709, 0.662, 0.909, 0.984, 0.585, 0.772, 0.843, 0.613, 0.969, 0.864, 0.857, 0.835, 0.974, 0.869, 0.641, 0.675, 0.802, 0.658, 0.744, 0.595, 0.989, 0.653, 0.712, 0.885, 0.983, 0.978, 0.998, 0.75, 0.621, 0.951]
global origin = 1
global destination = 60