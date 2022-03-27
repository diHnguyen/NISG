global arcs = [1 12; 1 13; 1 36; 1 47; 2 11; 2 13; 2 26; 2 46; 3 13; 3 15; 3 28; 3 41; 3 49; 4 18; 4 43; 4 46; 5 2; 5 10; 5 13; 5 15; 6 3; 6 8; 6 9; 6 20; 6 22; 6 27; 6 35; 6 37; 6 42; 6 50; 7 10; 7 21; 7 30; 7 34; 7 46; 7 50; 8 14; 8 21; 8 22; 8 23; 8 31; 8 40; 8 45; 9 22; 9 26; 9 37; 10 11; 10 44; 11 5; 11 17; 11 18; 11 25; 11 36; 11 41; 11 44; 11 46; 11 48; 12 8; 12 36; 12 38; 12 48; 13 3; 13 9; 13 14; 13 22; 13 25; 13 31; 13 34; 13 42; 14 15; 14 20; 14 35; 14 36; 14 41; 14 42; 14 43; 15 3; 15 4; 15 6; 15 44; 15 46; 15 48; 15 49; 15 50; 16 5; 16 15; 16 17; 16 30; 17 32; 17 36; 17 44; 17 46; 18 6; 18 8; 18 13; 18 16; 19 8; 19 11; 19 14; 19 22; 19 30; 19 36; 19 45; 20 6; 20 10; 20 29; 20 34; 20 37; 20 39; 20 41; 21 3; 21 34; 21 41; 21 48; 22 7; 22 20; 22 30; 22 37; 23 40; 24 6; 25 3; 25 12; 25 17; 25 32; 25 33; 25 34; 25 45; 26 10; 26 21; 26 44; 27 8; 27 15; 27 25; 27 37; 27 47; 28 9; 28 29; 28 31; 28 42; 28 49; 29 2; 29 14; 29 28; 29 31; 29 36; 29 42; 29 47; 30 8; 30 23; 30 38; 30 49; 31 9; 31 18; 31 19; 31 28; 31 34; 31 39; 31 40; 31 41; 32 6; 32 7; 32 8; 32 14; 32 17; 32 22; 32 34; 32 36; 32 37; 33 22; 34 7; 34 26; 34 33; 34 37; 34 50; 35 19; 35 50; 36 4; 36 6; 36 8; 36 26; 36 27; 36 29; 36 46; 36 49; 37 11; 37 17; 37 22; 38 13; 38 24; 38 35; 38 42; 39 42; 39 50; 40 4; 40 6; 40 17; 40 19; 40 24; 40 26; 40 48; 41 13; 41 20; 41 23; 41 30; 42 6; 42 25; 42 29; 42 46; 43 4; 44 5; 44 22; 44 40; 45 14; 45 17; 45 28; 45 31; 45 46; 46 2; 46 8; 46 23; 46 31; 47 19; 47 33; 47 40; 47 50; 48 2; 48 30; 48 36; 49 3; 49 10; 49 21; 49 24; 49 35; 49 48; 49 50]
global d_x = [1.0, 6.0, 5.0, 10.0, 6.0, 2.0, 4.0, 7.0, 7.0, 2.0, 10.0, 5.0, 2.0, 4.0, 3.0, 5.0, 6.0, 3.0, 3.0, 2.0, 9.0, 10.0, 10.0, 5.0, 3.0, 1.0, 3.0, 9.0, 9.0, 3.0, 1.0, 9.0, 3.0, 6.0, 9.0, 5.0, 5.0, 2.0, 3.0, 1.0, 5.0, 9.0, 10.0, 5.0, 3.0, 10.0, 4.0, 7.0, 1.0, 6.0, 7.0, 1.0, 7.0, 8.0, 3.0, 5.0, 2.0, 7.0, 3.0, 1.0, 9.0, 9.0, 8.0, 4.0, 2.0, 5.0, 6.0, 10.0, 9.0, 2.0, 9.0, 8.0, 8.0, 8.0, 2.0, 6.0, 1.0, 3.0, 1.0, 6.0, 8.0, 10.0, 7.0, 7.0, 7.0, 6.0, 8.0, 10.0, 8.0, 2.0, 3.0, 10.0, 8.0, 3.0, 5.0, 7.0, 7.0, 1.0, 9.0, 4.0, 7.0, 1.0, 10.0, 10.0, 2.0, 1.0, 5.0, 10.0, 5.0, 5.0, 5.0, 2.0, 1.0, 10.0, 7.0, 3.0, 8.0, 10.0, 4.0, 7.0, 8.0, 1.0, 1.0, 7.0, 9.0, 6.0, 7.0, 9.0, 2.0, 7.0, 2.0, 7.0, 1.0, 3.0, 5.0, 8.0, 6.0, 8.0, 1.0, 6.0, 4.0, 1.0, 2.0, 2.0, 10.0, 2.0, 4.0, 4.0, 4.0, 1.0, 7.0, 10.0, 8.0, 7.0, 7.0, 10.0, 3.0, 5.0, 6.0, 1.0, 5.0, 4.0, 10.0, 1.0, 4.0, 3.0, 5.0, 1.0, 3.0, 8.0, 6.0, 1.0, 1.0, 1.0, 9.0, 9.0, 7.0, 3.0, 7.0, 7.0, 6.0, 1.0, 5.0, 2.0, 8.0, 3.0, 5.0, 9.0, 10.0, 3.0, 2.0, 4.0, 7.0, 7.0, 6.0, 2.0, 8.0, 10.0, 7.0, 1.0, 5.0, 8.0, 4.0, 2.0, 6.0, 3.0, 4.0, 5.0, 10.0, 1.0, 1.0, 10.0, 3.0, 2.0, 5.0, 6.0, 8.0, 2.0, 3.0, 4.0, 10.0, 1.0, 4.0, 8.0, 9.0, 9.0, 2.0, 1.0, 4.0, 1.0, 7.0, 10.0, 6.0, 9.0, 2.0]
global b_x = 5
global d_y = [6.0, 6.0, 1.0, 6.0, 2.0, 6.0, 7.0, 1.0, 7.0, 3.0, 3.0, 4.0, 8.0, 4.0, 10.0, 7.0, 8.0, 1.0, 6.0, 5.0, 2.0, 8.0, 10.0, 1.0, 3.0, 4.0, 9.0, 7.0, 7.0, 2.0, 3.0, 10.0, 10.0, 7.0, 2.0, 5.0, 3.0, 8.0, 2.0, 3.0, 9.0, 6.0, 7.0, 5.0, 7.0, 5.0, 4.0, 3.0, 2.0, 6.0, 7.0, 7.0, 3.0, 9.0, 6.0, 4.0, 5.0, 3.0, 4.0, 1.0, 9.0, 5.0, 7.0, 2.0, 1.0, 3.0, 2.0, 7.0, 3.0, 4.0, 8.0, 5.0, 10.0, 6.0, 1.0, 6.0, 1.0, 9.0, 9.0, 4.0, 3.0, 3.0, 3.0, 10.0, 3.0, 3.0, 6.0, 3.0, 5.0, 5.0, 3.0, 5.0, 4.0, 4.0, 1.0, 6.0, 1.0, 3.0, 10.0, 3.0, 10.0, 1.0, 8.0, 6.0, 10.0, 10.0, 8.0, 4.0, 6.0, 7.0, 7.0, 2.0, 1.0, 7.0, 8.0, 10.0, 2.0, 2.0, 5.0, 9.0, 1.0, 10.0, 6.0, 4.0, 5.0, 7.0, 6.0, 6.0, 10.0, 5.0, 3.0, 2.0, 9.0, 1.0, 1.0, 6.0, 10.0, 8.0, 8.0, 3.0, 1.0, 1.0, 2.0, 10.0, 10.0, 10.0, 1.0, 3.0, 1.0, 1.0, 2.0, 1.0, 4.0, 7.0, 4.0, 5.0, 5.0, 2.0, 8.0, 3.0, 1.0, 1.0, 7.0, 8.0, 3.0, 4.0, 9.0, 7.0, 5.0, 10.0, 10.0, 10.0, 9.0, 6.0, 8.0, 8.0, 9.0, 2.0, 10.0, 8.0, 6.0, 4.0, 1.0, 7.0, 3.0, 10.0, 2.0, 9.0, 1.0, 6.0, 5.0, 7.0, 4.0, 3.0, 5.0, 1.0, 10.0, 7.0, 9.0, 4.0, 9.0, 5.0, 1.0, 3.0, 7.0, 4.0, 1.0, 8.0, 9.0, 5.0, 3.0, 6.0, 6.0, 9.0, 8.0, 7.0, 8.0, 6.0, 4.0, 9.0, 9.0, 7.0, 6.0, 1.0, 10.0, 9.0, 2.0, 2.0, 7.0, 5.0, 6.0, 9.0, 1.0, 9.0, 4.0]
global b_y = 10
global p = [0.376, 0.336, 0.787, 0.038, 0.234, 0.254, 0.722, 0.144, 0.549, 0.545, 0.244, 0.193, 0.227, 0.463, 0.928, 0.628, 0.933, 0.303, 0.556, 0.313, 0.554, 0.577, 0.975, 0.629, 0.638, 0.731, 0.34, 0.097, 0.26, 0.195, 0.92, 0.407, 0.24, 0.638, 0.508, 0.986, 0.57, 0.237, 0.309, 0.295, 0.435, 0.596, 0.788, 0.052, 0.47, 0.242, 0.379, 0.643, 0.636, 0.839, 0.595, 0.087, 0.79, 0.143, 0.291, 0.815, 0.64, 0.18, 0.918, 0.916, 0.838, 0.166, 0.065, 0.697, 0.32, 0.632, 0.185, 0.59, 0.667, 0.584, 0.46, 0.569, 0.47, 0.694, 0.466, 0.196, 0.892, 0.912, 0.327, 0.496, 0.781, 0.443, 0.878, 0.499, 0.872, 0.673, 0.291, 0.462, 0.002, 0.972, 0.245, 0.444, 0.13, 0.656, 0.052, 0.924, 0.148, 0.818, 0.382, 0.666, 0.933, 0.087, 0.53, 0.724, 0.532, 0.945, 0.361, 0.887, 0.711, 0.479, 0.969, 0.065, 0.448, 0.836, 0.4, 0.948, 0.995, 0.101, 0.496, 0.316, 0.987, 0.213, 0.574, 0.195, 0.541, 0.809, 0.231, 0.038, 0.961, 0.198, 0.243, 0.432, 0.268, 0.192, 0.74, 0.821, 0.577, 0.343, 0.849, 0.766, 0.337, 0.243, 0.001, 0.177, 0.574, 0.623, 0.519, 0.298, 0.325, 0.427, 0.401, 0.617, 0.186, 0.657, 0.011, 0.285, 0.035, 0.838, 0.734, 0.064, 0.457, 0.315, 0.897, 0.904, 0.597, 0.796, 0.353, 0.025, 0.195, 0.902, 0.92, 0.548, 0.537, 0.888, 0.166, 0.397, 0.612, 0.538, 0.489, 0.575, 0.655, 0.419, 0.679, 0.265, 0.638, 0.647, 0.011, 0.962, 0.271, 0.963, 0.865, 0.577, 0.238, 0.113, 0.239, 0.202, 0.592, 0.869, 0.147, 0.239, 0.466, 0.775, 0.887, 0.654, 0.678, 0.565, 0.852, 0.671, 0.933, 0.737, 0.034, 0.975, 0.425, 0.7, 0.438, 0.999, 0.959, 0.31, 0.612, 0.407, 0.795, 0.662, 0.234, 0.195, 0.707, 0.175, 0.813, 0.638, 0.556, 0.453, 0.765, 0.541, 0.246, 0.167, 0.959]
global q = [0.654, 0.562, 0.883, 0.679, 0.345, 0.959, 0.784, 0.235, 0.757, 0.914, 0.953, 0.972, 0.41, 0.604, 0.987, 0.809, 0.98, 0.802, 0.81, 0.676, 0.59, 0.594, 0.991, 0.8, 0.892, 0.911, 0.957, 0.929, 0.302, 0.438, 0.997, 0.466, 0.728, 0.656, 0.622, 0.987, 0.951, 0.629, 0.844, 0.967, 0.545, 0.985, 0.889, 0.109, 0.647, 0.438, 0.647, 0.891, 0.935, 0.908, 0.816, 0.848, 0.963, 0.734, 0.836, 0.821, 0.706, 0.819, 0.945, 0.992, 0.892, 0.421, 0.249, 0.701, 0.569, 0.854, 0.791, 0.964, 0.972, 0.601, 0.723, 0.936, 0.483, 0.812, 0.978, 0.425, 0.929, 0.98, 0.6, 0.872, 0.924, 0.837, 0.958, 0.904, 0.98, 0.971, 0.41, 0.891, 0.818, 0.98, 0.367, 0.813, 0.843, 0.703, 0.47, 0.962, 0.738, 0.899, 0.704, 0.691, 0.959, 0.586, 0.658, 0.912, 0.682, 0.95, 0.433, 0.961, 0.833, 0.818, 0.976, 0.571, 0.944, 0.847, 0.867, 0.981, 0.995, 0.89, 0.647, 0.708, 0.99, 0.874, 0.731, 0.473, 0.72, 0.985, 0.938, 0.937, 0.964, 0.983, 0.773, 0.524, 0.279, 0.894, 0.82, 0.83, 0.822, 0.495, 0.94, 0.864, 0.479, 0.471, 0.318, 0.47, 0.619, 0.995, 0.713, 0.374, 0.499, 0.706, 0.822, 0.957, 0.409, 0.678, 0.203, 0.969, 0.603, 0.881, 0.885, 0.455, 0.82, 0.64, 0.951, 0.917, 0.699, 0.973, 0.369, 0.058, 0.421, 0.982, 0.943, 0.924, 0.996, 0.944, 0.289, 0.924, 0.693, 0.922, 0.643, 0.915, 0.751, 0.95, 0.776, 0.508, 0.964, 0.848, 0.017, 0.968, 0.42, 0.965, 0.942, 0.788, 0.894, 0.648, 0.738, 0.869, 0.705, 0.91, 0.152, 0.478, 0.491, 0.817, 0.928, 0.708, 0.916, 0.671, 0.946, 0.882, 0.974, 0.812, 0.784, 0.989, 0.46, 0.8, 0.578, 0.999, 0.96, 0.963, 0.619, 0.613, 0.836, 0.715, 0.982, 0.98, 0.833, 0.79, 0.842, 0.813, 0.911, 0.551, 0.937, 0.899, 0.774, 0.687, 0.984]
global origin = 1
global destination = 50