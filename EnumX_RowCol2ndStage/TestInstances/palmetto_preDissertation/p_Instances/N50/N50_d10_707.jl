global arcs = [1 7; 1 12; 1 25; 1 37; 1 49; 2 10; 2 15; 2 25; 3 12; 3 34; 3 36; 3 41; 4 20; 4 22; 4 24; 4 30; 4 32; 4 33; 4 37; 5 4; 5 6; 5 13; 5 15; 5 21; 5 35; 6 50; 7 9; 7 22; 7 30; 7 35; 7 50; 8 6; 8 18; 8 24; 8 39; 8 48; 9 10; 9 12; 9 37; 9 47; 10 42; 11 3; 11 6; 11 14; 11 24; 11 26; 11 36; 11 39; 11 42; 12 15; 12 16; 12 27; 12 41; 13 3; 13 4; 13 28; 13 33; 13 36; 13 38; 13 39; 13 43; 14 9; 14 27; 14 34; 14 37; 14 49; 14 50; 15 2; 16 18; 16 20; 16 41; 16 44; 17 30; 17 40; 17 43; 17 44; 17 50; 18 2; 18 6; 18 11; 18 22; 18 31; 18 38; 19 36; 20 12; 20 17; 20 38; 20 49; 21 17; 21 20; 21 22; 21 30; 21 38; 21 39; 21 40; 21 41; 22 2; 22 18; 22 30; 22 40; 22 48; 23 4; 23 33; 23 37; 23 42; 24 26; 24 35; 24 39; 24 43; 24 47; 25 19; 25 22; 25 23; 25 30; 25 32; 25 37; 25 42; 25 50; 26 5; 26 8; 26 10; 26 32; 26 42; 26 47; 27 5; 27 13; 27 21; 27 29; 27 43; 28 7; 28 45; 28 47; 29 2; 29 7; 29 26; 29 28; 29 33; 29 38; 29 39; 29 47; 30 32; 31 2; 31 5; 31 25; 32 9; 32 10; 32 12; 32 33; 32 36; 33 4; 33 17; 33 47; 34 6; 34 17; 34 19; 34 25; 34 28; 35 14; 35 16; 35 17; 35 21; 35 38; 36 16; 36 19; 36 20; 36 30; 36 34; 36 47; 37 7; 37 8; 37 17; 37 18; 37 25; 37 27; 37 47; 38 7; 38 11; 38 14; 38 26; 38 36; 38 43; 39 24; 39 35; 39 47; 40 2; 40 7; 40 8; 40 11; 40 20; 41 8; 41 33; 41 39; 42 28; 42 31; 42 39; 43 7; 43 16; 43 19; 43 25; 43 32; 43 34; 44 9; 45 15; 45 20; 45 35; 45 38; 45 48; 45 49; 45 50; 46 3; 46 7; 46 9; 47 3; 47 4; 47 32; 47 35; 47 36; 48 2; 48 8; 48 9; 48 13; 48 30; 48 44; 48 45; 49 21; 49 22; 49 46; 49 50]
global d_x = [3.0, 9.0, 7.0, 8.0, 5.0, 2.0, 7.0, 3.0, 3.0, 9.0, 10.0, 4.0, 5.0, 3.0, 5.0, 1.0, 3.0, 6.0, 7.0, 5.0, 6.0, 10.0, 3.0, 10.0, 7.0, 5.0, 8.0, 10.0, 8.0, 10.0, 9.0, 7.0, 7.0, 3.0, 10.0, 4.0, 1.0, 2.0, 6.0, 5.0, 5.0, 9.0, 5.0, 4.0, 6.0, 4.0, 4.0, 1.0, 4.0, 10.0, 10.0, 8.0, 4.0, 6.0, 1.0, 3.0, 10.0, 4.0, 10.0, 8.0, 5.0, 2.0, 9.0, 3.0, 8.0, 10.0, 9.0, 7.0, 2.0, 4.0, 6.0, 3.0, 4.0, 5.0, 5.0, 2.0, 10.0, 10.0, 4.0, 1.0, 8.0, 4.0, 5.0, 6.0, 5.0, 3.0, 3.0, 6.0, 9.0, 7.0, 2.0, 2.0, 1.0, 9.0, 10.0, 6.0, 5.0, 4.0, 6.0, 3.0, 2.0, 2.0, 10.0, 8.0, 10.0, 10.0, 2.0, 10.0, 3.0, 6.0, 1.0, 3.0, 5.0, 1.0, 4.0, 7.0, 2.0, 1.0, 9.0, 1.0, 10.0, 4.0, 4.0, 5.0, 5.0, 1.0, 7.0, 8.0, 8.0, 6.0, 7.0, 7.0, 3.0, 5.0, 9.0, 5.0, 2.0, 2.0, 1.0, 7.0, 8.0, 5.0, 10.0, 9.0, 4.0, 2.0, 7.0, 10.0, 1.0, 5.0, 1.0, 3.0, 10.0, 7.0, 8.0, 1.0, 6.0, 4.0, 4.0, 10.0, 1.0, 4.0, 2.0, 8.0, 8.0, 2.0, 1.0, 5.0, 8.0, 2.0, 7.0, 5.0, 5.0, 5.0, 1.0, 3.0, 7.0, 3.0, 4.0, 6.0, 2.0, 6.0, 8.0, 3.0, 6.0, 2.0, 4.0, 8.0, 9.0, 1.0, 3.0, 2.0, 9.0, 8.0, 9.0, 6.0, 3.0, 10.0, 4.0, 7.0, 7.0, 1.0, 10.0, 10.0, 3.0, 2.0, 7.0, 1.0, 4.0, 5.0, 4.0, 9.0, 4.0, 6.0, 9.0, 2.0, 7.0, 2.0, 6.0, 3.0, 9.0, 7.0, 6.0, 2.0, 7.0, 8.0, 8.0, 8.0]
global b_x = 5
global d_y = [6.0, 2.0, 3.0, 8.0, 3.0, 8.0, 9.0, 3.0, 3.0, 2.0, 8.0, 1.0, 8.0, 3.0, 8.0, 7.0, 9.0, 1.0, 8.0, 7.0, 5.0, 5.0, 6.0, 9.0, 4.0, 1.0, 8.0, 7.0, 4.0, 2.0, 6.0, 7.0, 1.0, 10.0, 1.0, 10.0, 5.0, 1.0, 6.0, 8.0, 6.0, 8.0, 8.0, 8.0, 4.0, 6.0, 6.0, 6.0, 7.0, 4.0, 4.0, 10.0, 4.0, 5.0, 2.0, 1.0, 5.0, 4.0, 5.0, 1.0, 3.0, 7.0, 1.0, 4.0, 5.0, 2.0, 9.0, 4.0, 10.0, 2.0, 4.0, 6.0, 7.0, 10.0, 4.0, 2.0, 6.0, 7.0, 7.0, 7.0, 9.0, 7.0, 4.0, 10.0, 1.0, 9.0, 1.0, 5.0, 7.0, 6.0, 10.0, 9.0, 3.0, 6.0, 9.0, 4.0, 4.0, 6.0, 3.0, 7.0, 5.0, 10.0, 9.0, 5.0, 5.0, 7.0, 9.0, 6.0, 8.0, 4.0, 5.0, 5.0, 7.0, 1.0, 1.0, 8.0, 10.0, 1.0, 4.0, 9.0, 8.0, 4.0, 6.0, 6.0, 4.0, 2.0, 9.0, 9.0, 7.0, 1.0, 10.0, 6.0, 3.0, 1.0, 3.0, 7.0, 6.0, 10.0, 8.0, 9.0, 10.0, 3.0, 3.0, 9.0, 6.0, 7.0, 9.0, 3.0, 7.0, 2.0, 4.0, 8.0, 2.0, 4.0, 4.0, 7.0, 4.0, 7.0, 7.0, 9.0, 9.0, 8.0, 3.0, 8.0, 6.0, 3.0, 6.0, 3.0, 7.0, 7.0, 10.0, 2.0, 10.0, 7.0, 10.0, 2.0, 6.0, 6.0, 10.0, 2.0, 5.0, 5.0, 6.0, 9.0, 2.0, 2.0, 3.0, 7.0, 4.0, 1.0, 2.0, 1.0, 7.0, 6.0, 8.0, 10.0, 1.0, 8.0, 1.0, 5.0, 7.0, 8.0, 7.0, 4.0, 3.0, 5.0, 5.0, 9.0, 5.0, 10.0, 1.0, 1.0, 1.0, 5.0, 3.0, 10.0, 8.0, 1.0, 4.0, 7.0, 5.0, 6.0, 8.0, 5.0, 3.0, 7.0, 4.0, 1.0]
global b_y = 10
global p = [0.961, 0.282, 0.575, 0.155, 0.278, 0.603, 0.434, 0.229, 0.259, 0.864, 0.585, 0.366, 0.141, 0.309, 0.801, 0.759, 0.481, 0.715, 0.042, 0.137, 0.414, 0.379, 0.697, 0.621, 0.545, 0.492, 0.016, 0.617, 0.673, 0.54, 0.632, 0.256, 0.48, 0.27, 0.009, 0.858, 0.749, 0.302, 0.31, 0.757, 0.432, 0.789, 0.835, 0.827, 0.755, 0.936, 0.859, 0.156, 0.846, 0.024, 0.798, 0.261, 0.981, 0.653, 0.644, 0.577, 0.898, 0.965, 0.425, 0.025, 0.147, 0.01, 0.537, 0.684, 0.314, 0.343, 0.889, 0.496, 0.268, 0.713, 0.28, 0.037, 0.48, 0.887, 0.415, 0.043, 0.094, 0.671, 0.129, 0.516, 0.923, 0.487, 0.26, 0.378, 0.67, 0.213, 0.367, 0.089, 0.261, 0.929, 0.639, 0.837, 0.34, 0.495, 0.631, 0.802, 0.335, 0.925, 0.356, 0.361, 0.379, 0.99, 0.708, 0.015, 0.311, 0.482, 0.992, 0.876, 0.649, 0.691, 0.752, 0.862, 0.025, 0.952, 0.952, 0.311, 0.968, 0.785, 0.242, 0.359, 0.443, 0.193, 0.239, 0.214, 0.949, 0.421, 0.476, 0.8, 0.737, 0.581, 0.582, 0.283, 0.309, 0.398, 0.741, 0.21, 0.019, 0.658, 0.391, 0.271, 0.181, 0.317, 0.291, 0.751, 0.72, 0.526, 0.398, 0.775, 0.775, 0.57, 0.545, 0.245, 0.901, 0.157, 0.898, 0.584, 0.625, 0.447, 0.344, 0.197, 0.695, 0.31, 0.59, 0.148, 0.655, 0.657, 0.907, 0.767, 0.588, 0.031, 0.536, 0.481, 0.384, 0.025, 0.425, 0.448, 0.945, 0.13, 0.989, 0.382, 0.642, 0.818, 0.992, 0.913, 0.275, 0.82, 0.316, 0.68, 0.178, 0.496, 0.241, 0.987, 0.643, 0.181, 0.165, 0.752, 0.166, 0.252, 0.757, 0.991, 0.056, 0.164, 0.152, 0.5, 0.668, 0.73, 0.264, 0.077, 0.609, 0.508, 0.049, 0.885, 0.13, 0.431, 0.708, 0.238, 0.444, 0.595, 0.386, 0.683, 0.714, 0.791, 0.932, 0.26, 0.196, 0.704, 0.662, 0.244]
global q = [0.993, 0.588, 0.959, 0.465, 0.768, 0.761, 0.915, 0.332, 0.326, 0.979, 0.615, 0.504, 0.576, 0.345, 0.88, 0.873, 0.822, 0.949, 0.311, 0.81, 0.558, 0.414, 0.866, 0.667, 0.604, 0.606, 0.749, 0.683, 0.879, 0.611, 0.82, 0.836, 0.7, 0.371, 0.962, 0.932, 0.938, 0.404, 0.941, 0.995, 0.452, 0.99, 0.987, 0.98, 0.971, 0.993, 0.958, 0.562, 0.968, 0.922, 0.967, 0.263, 0.998, 0.962, 0.93, 0.712, 0.934, 0.997, 0.474, 0.496, 0.347, 0.916, 0.577, 0.784, 0.626, 0.843, 0.989, 0.691, 0.39, 0.907, 0.487, 0.318, 0.986, 0.946, 0.594, 0.464, 0.299, 0.737, 0.718, 0.556, 0.936, 0.756, 0.669, 0.485, 0.987, 0.517, 0.408, 0.369, 0.558, 0.982, 0.839, 0.875, 0.849, 0.648, 0.876, 0.817, 0.954, 0.933, 0.381, 0.73, 0.986, 0.999, 0.738, 0.934, 0.925, 0.806, 0.999, 0.958, 0.724, 0.942, 0.792, 0.896, 0.573, 0.997, 0.966, 0.514, 0.978, 0.942, 0.449, 0.787, 0.967, 0.746, 0.918, 0.424, 0.964, 0.71, 0.707, 0.95, 0.84, 0.749, 0.782, 0.819, 0.879, 0.732, 0.792, 0.784, 0.485, 0.842, 0.409, 0.929, 0.974, 0.946, 0.711, 0.865, 0.898, 0.794, 0.727, 0.912, 0.859, 0.976, 0.994, 0.603, 0.968, 0.564, 0.96, 0.742, 0.824, 0.701, 0.58, 0.361, 0.793, 0.985, 0.667, 0.594, 0.916, 0.922, 0.909, 0.972, 0.658, 0.665, 0.886, 0.481, 0.516, 0.444, 0.92, 0.619, 0.994, 0.155, 0.991, 0.945, 0.893, 0.9, 0.999, 0.931, 0.787, 0.824, 0.762, 0.764, 0.694, 0.5, 0.347, 0.995, 0.979, 0.427, 0.319, 0.994, 0.802, 0.411, 0.881, 0.996, 0.528, 0.333, 0.997, 0.993, 0.706, 0.897, 0.927, 0.446, 0.923, 0.679, 0.228, 0.981, 0.329, 0.486, 0.716, 0.275, 0.701, 0.704, 0.726, 0.939, 0.785, 0.984, 0.941, 0.561, 0.915, 0.728, 0.746, 0.872]
global origin = 1
global destination = 50