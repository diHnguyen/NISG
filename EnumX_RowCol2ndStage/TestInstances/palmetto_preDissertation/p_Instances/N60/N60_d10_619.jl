global arcs = [1 7; 1 13; 1 14; 1 19; 1 29; 1 32; 1 48; 1 51; 2 4; 2 6; 2 15; 2 17; 2 34; 2 44; 2 56; 3 2; 3 6; 3 14; 3 15; 3 16; 3 26; 3 57; 4 5; 4 9; 4 31; 4 34; 4 37; 4 39; 4 50; 4 55; 4 57; 5 12; 5 14; 5 20; 5 23; 5 24; 5 35; 5 43; 6 12; 6 18; 6 21; 6 32; 7 29; 7 34; 7 56; 8 6; 8 20; 8 24; 8 26; 8 49; 8 59; 9 5; 9 10; 9 14; 9 19; 9 32; 9 36; 9 56; 10 5; 10 12; 10 45; 11 13; 11 21; 11 40; 11 47; 12 13; 12 36; 12 40; 12 46; 12 55; 12 57; 13 6; 13 17; 13 21; 13 22; 13 39; 13 41; 13 46; 14 13; 14 23; 15 17; 15 31; 15 46; 16 5; 16 6; 17 10; 17 19; 17 21; 17 24; 17 25; 17 26; 17 28; 17 43; 17 46; 17 53; 18 11; 18 16; 18 25; 18 33; 18 43; 18 53; 19 15; 19 26; 19 35; 19 58; 20 9; 21 13; 21 17; 21 22; 21 23; 21 27; 21 34; 21 36; 22 14; 22 15; 22 33; 22 35; 22 50; 22 51; 23 18; 23 22; 23 27; 23 43; 23 53; 23 55; 24 11; 24 51; 25 2; 25 14; 25 22; 25 26; 25 35; 25 36; 25 53; 26 7; 26 17; 26 22; 26 30; 27 45; 27 46; 27 48; 27 53; 28 3; 28 7; 28 9; 28 13; 28 22; 28 29; 28 32; 28 40; 28 53; 28 55; 29 22; 29 34; 29 39; 29 58; 30 2; 30 4; 30 17; 30 31; 30 40; 30 41; 30 60; 31 22; 31 32; 31 47; 31 51; 31 54; 31 56; 31 60; 32 23; 32 26; 32 35; 32 37; 32 45; 32 56; 33 25; 33 28; 33 46; 34 16; 34 31; 34 42; 34 48; 35 6; 35 7; 35 10; 35 38; 35 50; 35 60; 36 6; 36 16; 36 24; 36 30; 36 31; 36 55; 37 13; 37 14; 37 28; 37 35; 37 46; 37 60; 38 8; 38 17; 38 26; 38 33; 38 44; 38 56; 39 8; 39 13; 39 20; 39 32; 39 34; 40 9; 40 11; 40 17; 40 20; 40 33; 40 47; 40 48; 40 59; 41 4; 41 29; 41 44; 41 51; 42 26; 42 47; 42 51; 43 18; 43 33; 43 34; 43 36; 43 42; 43 55; 44 5; 44 10; 44 12; 44 18; 44 21; 44 22; 44 33; 44 42; 44 52; 44 54; 45 17; 45 27; 45 30; 45 32; 45 34; 45 43; 45 50; 45 54; 46 18; 46 20; 46 23; 46 41; 46 50; 46 51; 47 2; 47 5; 47 9; 47 13; 47 22; 47 25; 47 27; 47 51; 47 57; 48 13; 48 19; 48 21; 48 47; 49 3; 49 5; 49 9; 49 10; 49 16; 49 18; 49 28; 49 47; 49 55; 50 16; 50 20; 50 34; 50 37; 50 48; 50 58; 50 60; 51 14; 51 25; 51 32; 51 48; 51 55; 51 57; 52 31; 52 43; 52 51; 53 2; 53 11; 53 16; 53 19; 53 42; 53 44; 54 10; 54 20; 54 21; 54 22; 54 42; 54 49; 54 51; 54 59; 55 15; 55 36; 55 41; 56 8; 56 22; 56 27; 56 38; 56 46; 56 57; 57 2; 57 7; 57 11; 57 41; 58 5; 58 11; 58 15; 58 22; 58 38; 59 17; 59 20; 59 23; 59 27; 59 41]
global d_x = [8.0, 5.0, 2.0, 6.0, 6.0, 3.0, 8.0, 4.0, 8.0, 3.0, 3.0, 4.0, 7.0, 5.0, 7.0, 3.0, 7.0, 7.0, 1.0, 10.0, 4.0, 9.0, 10.0, 6.0, 5.0, 7.0, 7.0, 1.0, 7.0, 8.0, 10.0, 5.0, 9.0, 8.0, 3.0, 1.0, 6.0, 8.0, 6.0, 4.0, 5.0, 9.0, 5.0, 7.0, 6.0, 10.0, 8.0, 5.0, 2.0, 3.0, 4.0, 1.0, 9.0, 8.0, 1.0, 10.0, 7.0, 5.0, 3.0, 3.0, 9.0, 8.0, 7.0, 6.0, 9.0, 2.0, 4.0, 7.0, 2.0, 10.0, 6.0, 5.0, 10.0, 8.0, 1.0, 8.0, 10.0, 3.0, 9.0, 4.0, 5.0, 7.0, 2.0, 7.0, 5.0, 6.0, 5.0, 8.0, 8.0, 5.0, 6.0, 5.0, 6.0, 2.0, 4.0, 3.0, 1.0, 3.0, 1.0, 10.0, 7.0, 9.0, 7.0, 4.0, 7.0, 4.0, 8.0, 7.0, 6.0, 10.0, 2.0, 9.0, 5.0, 4.0, 3.0, 4.0, 8.0, 4.0, 4.0, 9.0, 10.0, 6.0, 5.0, 9.0, 1.0, 4.0, 4.0, 5.0, 1.0, 1.0, 8.0, 1.0, 3.0, 2.0, 9.0, 4.0, 7.0, 9.0, 5.0, 5.0, 3.0, 7.0, 2.0, 8.0, 10.0, 6.0, 2.0, 7.0, 2.0, 9.0, 6.0, 7.0, 2.0, 10.0, 8.0, 9.0, 6.0, 2.0, 6.0, 7.0, 1.0, 4.0, 1.0, 10.0, 3.0, 8.0, 10.0, 7.0, 7.0, 6.0, 3.0, 8.0, 1.0, 2.0, 1.0, 4.0, 1.0, 4.0, 10.0, 2.0, 5.0, 8.0, 5.0, 3.0, 7.0, 3.0, 8.0, 8.0, 5.0, 10.0, 9.0, 7.0, 8.0, 9.0, 5.0, 6.0, 2.0, 3.0, 7.0, 8.0, 4.0, 10.0, 9.0, 9.0, 9.0, 8.0, 7.0, 10.0, 10.0, 3.0, 7.0, 5.0, 8.0, 6.0, 1.0, 4.0, 7.0, 7.0, 10.0, 7.0, 5.0, 1.0, 3.0, 7.0, 7.0, 5.0, 1.0, 8.0, 1.0, 8.0, 9.0, 1.0, 7.0, 5.0, 10.0, 6.0, 1.0, 7.0, 3.0, 10.0, 8.0, 4.0, 5.0, 9.0, 2.0, 6.0, 5.0, 9.0, 8.0, 4.0, 8.0, 3.0, 8.0, 1.0, 1.0, 4.0, 9.0, 4.0, 3.0, 9.0, 10.0, 9.0, 8.0, 8.0, 9.0, 9.0, 6.0, 7.0, 7.0, 4.0, 7.0, 9.0, 2.0, 2.0, 3.0, 1.0, 9.0, 5.0, 7.0, 5.0, 7.0, 4.0, 5.0, 1.0, 9.0, 5.0, 7.0, 2.0, 7.0, 10.0, 5.0, 10.0, 2.0, 6.0, 7.0, 4.0, 8.0, 3.0, 3.0, 2.0, 2.0, 5.0, 10.0, 4.0, 5.0, 9.0, 10.0, 10.0, 3.0, 8.0, 10.0, 8.0, 5.0, 9.0, 7.0, 7.0, 6.0, 8.0, 1.0, 9.0, 9.0, 5.0, 3.0, 4.0, 8.0, 8.0, 5.0, 4.0, 5.0, 1.0, 7.0, 5.0]
global b_x = 5
global d_y = [9.0, 6.0, 5.0, 2.0, 8.0, 10.0, 8.0, 10.0, 7.0, 10.0, 8.0, 10.0, 2.0, 10.0, 6.0, 4.0, 7.0, 9.0, 2.0, 9.0, 2.0, 4.0, 2.0, 4.0, 8.0, 6.0, 6.0, 7.0, 7.0, 6.0, 10.0, 5.0, 5.0, 6.0, 5.0, 7.0, 2.0, 1.0, 3.0, 4.0, 7.0, 5.0, 3.0, 10.0, 10.0, 2.0, 8.0, 8.0, 5.0, 7.0, 6.0, 1.0, 2.0, 4.0, 9.0, 5.0, 8.0, 1.0, 3.0, 4.0, 2.0, 3.0, 3.0, 4.0, 4.0, 3.0, 9.0, 10.0, 7.0, 2.0, 8.0, 3.0, 2.0, 10.0, 5.0, 5.0, 9.0, 10.0, 9.0, 6.0, 7.0, 5.0, 10.0, 6.0, 3.0, 8.0, 2.0, 3.0, 2.0, 10.0, 10.0, 1.0, 3.0, 4.0, 9.0, 7.0, 8.0, 2.0, 9.0, 4.0, 5.0, 9.0, 6.0, 1.0, 7.0, 2.0, 5.0, 6.0, 9.0, 8.0, 9.0, 9.0, 9.0, 2.0, 2.0, 4.0, 10.0, 2.0, 4.0, 9.0, 1.0, 5.0, 8.0, 1.0, 1.0, 5.0, 8.0, 5.0, 4.0, 5.0, 2.0, 4.0, 4.0, 2.0, 7.0, 1.0, 4.0, 9.0, 3.0, 2.0, 5.0, 5.0, 4.0, 7.0, 5.0, 10.0, 7.0, 6.0, 3.0, 4.0, 5.0, 2.0, 10.0, 1.0, 8.0, 2.0, 8.0, 10.0, 10.0, 6.0, 6.0, 5.0, 3.0, 8.0, 8.0, 2.0, 6.0, 1.0, 4.0, 6.0, 1.0, 8.0, 2.0, 1.0, 2.0, 3.0, 3.0, 9.0, 2.0, 9.0, 10.0, 8.0, 5.0, 1.0, 8.0, 4.0, 7.0, 5.0, 10.0, 8.0, 3.0, 7.0, 7.0, 4.0, 10.0, 10.0, 4.0, 7.0, 6.0, 2.0, 4.0, 4.0, 5.0, 9.0, 9.0, 6.0, 4.0, 5.0, 9.0, 3.0, 3.0, 6.0, 1.0, 8.0, 8.0, 7.0, 8.0, 10.0, 2.0, 7.0, 8.0, 10.0, 2.0, 8.0, 2.0, 2.0, 9.0, 3.0, 3.0, 6.0, 4.0, 10.0, 5.0, 2.0, 3.0, 9.0, 1.0, 4.0, 1.0, 2.0, 5.0, 5.0, 5.0, 8.0, 1.0, 10.0, 9.0, 2.0, 8.0, 7.0, 1.0, 3.0, 3.0, 5.0, 10.0, 4.0, 3.0, 9.0, 2.0, 6.0, 8.0, 3.0, 10.0, 8.0, 9.0, 8.0, 10.0, 6.0, 1.0, 10.0, 3.0, 5.0, 4.0, 1.0, 5.0, 7.0, 2.0, 3.0, 9.0, 2.0, 3.0, 5.0, 2.0, 5.0, 1.0, 1.0, 3.0, 5.0, 5.0, 6.0, 5.0, 9.0, 6.0, 9.0, 9.0, 8.0, 2.0, 2.0, 7.0, 5.0, 3.0, 3.0, 6.0, 4.0, 3.0, 2.0, 1.0, 10.0, 6.0, 5.0, 1.0, 8.0, 4.0, 9.0, 2.0, 6.0, 1.0, 6.0, 4.0, 5.0, 6.0, 8.0, 2.0, 7.0, 1.0, 9.0, 9.0, 6.0, 8.0, 4.0, 8.0, 5.0]
global b_y = 10
global p = [0.329, 0.362, 0.542, 0.672, 0.283, 0.584, 0.129, 0.801, 0.688, 0.501, 0.76, 0.293, 0.09, 0.192, 0.471, 0.909, 0.196, 0.168, 0.572, 0.325, 0.1, 0.671, 0.838, 0.048, 0.662, 0.322, 0.104, 0.221, 0.374, 0.722, 0.402, 0.882, 0.923, 0.468, 0.764, 0.602, 0.22, 0.23, 0.922, 0.893, 0.718, 0.506, 0.693, 0.313, 0.04, 0.576, 0.579, 0.192, 0.63, 0.877, 0.193, 0.962, 0.78, 0.345, 0.592, 0.989, 0.46, 0.995, 0.982, 0.962, 0.273, 0.058, 0.247, 0.517, 0.708, 0.695, 0.275, 0.795, 0.922, 0.365, 0.394, 0.367, 0.96, 0.971, 0.189, 0.156, 0.781, 0.926, 0.73, 0.656, 0.214, 0.631, 0.617, 0.655, 0.179, 0.992, 0.067, 0.044, 0.767, 0.932, 0.615, 0.869, 0.321, 0.934, 0.744, 0.84, 0.692, 0.68, 0.091, 0.868, 0.933, 0.137, 0.832, 0.507, 0.046, 0.838, 0.596, 0.205, 0.698, 0.056, 0.507, 0.41, 0.146, 0.568, 0.614, 0.309, 0.57, 0.365, 0.103, 0.629, 0.248, 0.924, 0.676, 0.922, 0.464, 0.723, 0.804, 0.71, 0.613, 0.707, 0.414, 0.619, 0.449, 0.852, 0.496, 0.403, 0.275, 0.254, 0.123, 0.307, 0.008, 0.282, 0.906, 0.221, 0.34, 0.631, 0.83, 0.879, 0.809, 0.826, 0.988, 0.258, 0.86, 0.734, 0.672, 0.341, 0.795, 0.344, 0.376, 0.022, 0.132, 0.103, 0.794, 0.196, 0.28, 0.97, 0.587, 0.808, 0.099, 0.236, 0.215, 0.642, 0.745, 0.209, 0.028, 0.905, 0.152, 0.087, 0.173, 0.002, 0.724, 0.077, 0.874, 0.748, 0.613, 0.425, 0.432, 0.848, 0.877, 0.877, 0.341, 0.339, 0.03, 0.028, 0.935, 0.15, 0.639, 0.539, 0.267, 0.142, 0.895, 0.56, 0.96, 0.979, 0.267, 0.527, 0.057, 0.11, 0.561, 0.347, 0.468, 0.473, 0.613, 0.946, 0.127, 0.834, 0.906, 0.177, 0.807, 0.889, 0.609, 0.556, 0.127, 0.82, 0.991, 0.857, 0.264, 0.208, 0.025, 0.436, 0.846, 0.889, 0.881, 0.146, 0.75, 0.046, 0.169, 0.54, 0.604, 0.605, 0.016, 0.664, 0.863, 0.965, 0.391, 0.968, 0.084, 0.346, 0.599, 0.309, 0.004, 0.051, 0.268, 0.659, 0.594, 0.822, 0.095, 0.932, 0.585, 0.378, 0.521, 0.204, 0.982, 0.898, 0.741, 0.377, 0.696, 0.59, 0.9, 0.08, 0.014, 0.493, 0.065, 0.598, 0.581, 0.932, 0.617, 0.878, 0.189, 0.511, 0.65, 0.246, 0.335, 0.682, 0.38, 0.294, 0.444, 0.283, 0.874, 0.578, 0.27, 0.305, 0.573, 0.789, 0.701, 0.864, 0.167, 0.124, 0.291, 0.617, 0.416, 0.479, 0.543, 0.499, 0.701, 0.771, 0.257, 0.159, 0.459, 0.994, 0.718, 0.525, 0.379, 0.952, 0.757, 0.118, 0.451, 0.372, 0.333, 0.319, 0.855, 0.155, 0.928, 0.864, 0.371, 0.905, 0.978, 0.268, 0.815, 0.41, 0.801, 0.031]
global q = [0.928, 0.59, 0.984, 0.997, 0.284, 0.609, 0.878, 0.811, 0.733, 0.522, 0.938, 0.617, 0.848, 0.329, 0.554, 0.932, 0.956, 0.914, 0.836, 0.69, 0.876, 0.812, 0.928, 0.105, 0.729, 0.914, 0.88, 0.804, 0.388, 0.913, 0.615, 0.908, 0.959, 0.807, 0.925, 0.67, 0.634, 0.448, 0.943, 0.94, 0.867, 0.556, 0.792, 0.803, 0.986, 0.991, 0.709, 0.393, 0.887, 0.878, 0.334, 0.963, 0.916, 0.463, 0.753, 0.996, 0.723, 0.997, 0.985, 0.988, 0.777, 0.736, 0.613, 0.76, 0.712, 0.988, 0.97, 0.821, 0.956, 0.999, 0.913, 0.77, 0.993, 0.999, 0.994, 0.874, 0.798, 0.95, 0.879, 0.996, 0.874, 0.658, 0.989, 0.775, 0.83, 0.995, 0.111, 0.859, 0.972, 0.962, 0.811, 0.96, 0.476, 0.976, 0.754, 0.89, 0.834, 0.859, 0.685, 0.892, 0.969, 0.195, 0.913, 0.746, 0.231, 0.891, 0.625, 0.334, 0.817, 0.366, 0.836, 0.95, 0.946, 0.985, 0.824, 0.673, 0.66, 0.783, 0.314, 0.839, 0.761, 0.978, 0.711, 0.939, 0.901, 0.831, 0.974, 0.816, 0.891, 0.837, 0.565, 0.769, 0.964, 0.91, 0.674, 0.752, 0.296, 0.295, 0.235, 0.602, 0.808, 0.778, 0.96, 0.827, 0.37, 0.872, 0.983, 0.879, 0.843, 0.845, 0.992, 0.696, 0.927, 0.777, 0.847, 0.654, 0.922, 0.349, 0.487, 0.181, 0.961, 0.27, 0.946, 0.603, 0.851, 0.999, 0.957, 0.959, 0.739, 0.913, 0.334, 0.928, 0.957, 0.719, 0.149, 0.925, 0.908, 0.538, 0.782, 0.4, 0.755, 0.723, 0.952, 0.983, 0.885, 0.862, 0.614, 0.946, 0.944, 0.927, 0.903, 0.574, 0.521, 0.464, 0.94, 0.927, 0.675, 0.629, 0.999, 0.4, 0.942, 0.89, 0.965, 0.986, 0.89, 0.885, 0.453, 0.294, 0.833, 0.517, 0.876, 0.56, 0.666, 0.985, 0.736, 0.998, 0.994, 0.515, 0.814, 0.958, 0.866, 0.596, 0.48, 0.995, 0.996, 0.954, 0.784, 0.552, 0.397, 0.987, 0.963, 0.904, 0.91, 0.779, 0.909, 0.874, 0.189, 0.751, 0.623, 0.816, 0.194, 0.734, 0.997, 0.974, 0.413, 0.979, 0.655, 0.791, 0.626, 0.464, 0.01, 0.606, 0.66, 0.846, 0.843, 0.873, 0.182, 0.983, 0.605, 0.632, 0.927, 0.632, 0.995, 0.959, 0.926, 0.891, 0.986, 0.785, 0.919, 0.52, 0.305, 0.992, 0.435, 0.896, 0.767, 0.934, 0.836, 0.943, 0.573, 0.953, 0.902, 0.325, 0.581, 0.714, 0.513, 0.554, 0.778, 0.358, 0.943, 0.994, 0.479, 0.727, 0.7, 0.834, 0.889, 0.942, 0.694, 0.691, 0.495, 0.845, 0.668, 0.786, 0.981, 0.77, 0.979, 0.93, 0.478, 0.81, 0.959, 0.998, 0.98, 0.649, 0.666, 0.966, 0.87, 0.686, 0.515, 0.944, 0.809, 0.912, 0.922, 0.845, 0.952, 0.954, 0.669, 0.973, 0.986, 0.884, 0.978, 0.681, 0.974, 0.491]
global origin = 1
global destination = 60