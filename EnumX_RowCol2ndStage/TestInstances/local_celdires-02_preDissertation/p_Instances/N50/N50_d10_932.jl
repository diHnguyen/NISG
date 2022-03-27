global arcs = [1 3; 1 7; 1 17; 1 33; 1 38; 1 45; 1 48; 2 7; 2 28; 2 31; 2 37; 2 41; 3 7; 3 9; 3 15; 3 28; 3 35; 3 37; 4 2; 4 5; 4 9; 4 14; 4 46; 4 49; 5 4; 5 7; 5 10; 5 18; 5 27; 5 41; 5 47; 6 32; 7 8; 7 9; 7 32; 7 40; 7 49; 8 6; 8 13; 8 23; 8 29; 8 33; 8 40; 8 43; 8 46; 9 14; 9 16; 9 28; 9 37; 9 43; 9 46; 10 21; 10 26; 10 28; 10 30; 10 45; 10 48; 11 2; 11 15; 11 17; 11 18; 11 25; 11 28; 11 29; 11 31; 11 32; 11 48; 11 49; 12 10; 12 18; 12 19; 12 25; 12 26; 12 30; 12 33; 13 12; 13 18; 13 28; 13 29; 13 37; 14 9; 14 30; 14 34; 14 38; 14 43; 15 7; 15 14; 15 17; 16 8; 16 45; 16 48; 17 3; 17 11; 17 27; 17 31; 17 40; 17 49; 18 16; 18 23; 18 36; 18 44; 18 50; 19 11; 19 13; 20 6; 20 13; 20 25; 20 48; 21 9; 21 14; 21 17; 21 18; 21 26; 22 15; 23 9; 23 39; 24 13; 24 14; 24 21; 24 31; 24 33; 24 38; 24 42; 24 46; 25 6; 25 27; 25 32; 25 36; 25 49; 26 19; 26 32; 26 36; 27 9; 27 24; 27 29; 27 30; 28 8; 28 12; 28 14; 28 19; 28 24; 28 48; 29 43; 30 4; 30 6; 30 11; 30 25; 30 45; 30 49; 31 9; 31 12; 31 37; 31 40; 31 42; 31 45; 32 12; 32 21; 32 26; 32 34; 32 35; 33 22; 33 26; 33 31; 33 34; 34 20; 34 36; 34 40; 35 7; 35 21; 35 26; 35 27; 35 30; 35 31; 35 36; 35 49; 36 8; 36 18; 36 29; 36 31; 36 35; 36 39; 36 41; 36 47; 37 13; 37 17; 37 29; 37 30; 37 35; 37 38; 37 39; 37 49; 38 8; 38 29; 38 34; 38 40; 38 42; 38 48; 39 17; 39 18; 39 33; 39 36; 39 50; 40 9; 40 25; 40 30; 40 39; 40 46; 41 4; 41 8; 41 22; 41 24; 41 34; 41 40; 41 48; 42 2; 42 13; 42 16; 42 24; 42 36; 42 46; 42 47; 43 8; 43 10; 43 22; 43 28; 43 40; 43 45; 44 13; 44 15; 44 48; 44 49; 45 8; 45 14; 45 50; 46 10; 46 17; 46 20; 46 25; 46 34; 46 41; 46 48; 47 13; 47 19; 47 28; 47 41; 47 43; 48 3; 48 18; 48 37; 48 39; 48 44; 48 49; 49 11; 49 25; 49 32; 49 40]
global d_x = [5.0, 2.0, 10.0, 5.0, 1.0, 5.0, 9.0, 9.0, 6.0, 10.0, 10.0, 6.0, 9.0, 3.0, 7.0, 6.0, 2.0, 7.0, 9.0, 6.0, 1.0, 10.0, 9.0, 4.0, 6.0, 8.0, 9.0, 6.0, 9.0, 4.0, 3.0, 4.0, 3.0, 4.0, 7.0, 7.0, 1.0, 1.0, 10.0, 6.0, 1.0, 9.0, 8.0, 3.0, 3.0, 10.0, 2.0, 3.0, 7.0, 5.0, 6.0, 5.0, 5.0, 8.0, 1.0, 6.0, 4.0, 4.0, 10.0, 9.0, 5.0, 1.0, 4.0, 10.0, 1.0, 4.0, 6.0, 1.0, 4.0, 10.0, 6.0, 9.0, 7.0, 10.0, 4.0, 8.0, 2.0, 4.0, 1.0, 6.0, 5.0, 1.0, 1.0, 4.0, 6.0, 1.0, 10.0, 6.0, 6.0, 2.0, 1.0, 3.0, 6.0, 10.0, 3.0, 9.0, 2.0, 8.0, 3.0, 8.0, 8.0, 2.0, 9.0, 9.0, 8.0, 3.0, 10.0, 4.0, 7.0, 8.0, 10.0, 2.0, 2.0, 8.0, 2.0, 2.0, 10.0, 3.0, 8.0, 4.0, 10.0, 7.0, 6.0, 1.0, 4.0, 7.0, 10.0, 2.0, 10.0, 4.0, 10.0, 5.0, 5.0, 1.0, 8.0, 1.0, 9.0, 5.0, 8.0, 8.0, 2.0, 3.0, 6.0, 4.0, 3.0, 5.0, 7.0, 8.0, 2.0, 4.0, 3.0, 6.0, 3.0, 5.0, 9.0, 1.0, 3.0, 7.0, 6.0, 7.0, 2.0, 3.0, 3.0, 9.0, 7.0, 1.0, 7.0, 2.0, 9.0, 9.0, 9.0, 3.0, 5.0, 10.0, 1.0, 3.0, 9.0, 3.0, 9.0, 2.0, 10.0, 6.0, 5.0, 1.0, 1.0, 1.0, 4.0, 5.0, 3.0, 6.0, 1.0, 5.0, 3.0, 4.0, 6.0, 9.0, 9.0, 5.0, 6.0, 1.0, 1.0, 9.0, 1.0, 7.0, 2.0, 10.0, 7.0, 1.0, 1.0, 2.0, 5.0, 6.0, 2.0, 7.0, 4.0, 7.0, 8.0, 5.0, 1.0, 6.0, 4.0, 1.0, 2.0, 5.0, 5.0, 2.0, 1.0, 2.0, 10.0, 10.0, 10.0, 3.0, 3.0, 3.0, 4.0, 9.0, 2.0, 8.0, 4.0, 8.0, 6.0, 4.0, 1.0, 2.0, 6.0, 5.0, 1.0, 9.0, 9.0, 2.0, 5.0, 5.0, 2.0, 4.0, 9.0, 5.0]
global b_x = 5
global d_y = [3.0, 5.0, 4.0, 5.0, 10.0, 5.0, 3.0, 8.0, 6.0, 5.0, 3.0, 7.0, 6.0, 6.0, 6.0, 5.0, 4.0, 9.0, 10.0, 4.0, 4.0, 5.0, 6.0, 7.0, 2.0, 3.0, 10.0, 7.0, 8.0, 8.0, 9.0, 9.0, 10.0, 6.0, 5.0, 7.0, 7.0, 5.0, 6.0, 6.0, 7.0, 3.0, 8.0, 4.0, 7.0, 3.0, 9.0, 9.0, 4.0, 10.0, 1.0, 2.0, 8.0, 2.0, 3.0, 3.0, 4.0, 7.0, 5.0, 9.0, 8.0, 7.0, 3.0, 6.0, 2.0, 8.0, 1.0, 1.0, 1.0, 4.0, 6.0, 3.0, 10.0, 5.0, 9.0, 1.0, 1.0, 2.0, 3.0, 5.0, 3.0, 1.0, 10.0, 1.0, 6.0, 9.0, 4.0, 4.0, 10.0, 1.0, 9.0, 9.0, 9.0, 6.0, 1.0, 4.0, 1.0, 5.0, 7.0, 9.0, 1.0, 8.0, 9.0, 2.0, 5.0, 9.0, 2.0, 7.0, 4.0, 1.0, 2.0, 8.0, 9.0, 7.0, 2.0, 9.0, 7.0, 2.0, 7.0, 8.0, 8.0, 2.0, 6.0, 7.0, 6.0, 4.0, 7.0, 5.0, 3.0, 9.0, 9.0, 5.0, 8.0, 7.0, 3.0, 1.0, 8.0, 1.0, 6.0, 1.0, 8.0, 6.0, 6.0, 10.0, 2.0, 6.0, 3.0, 1.0, 4.0, 3.0, 8.0, 6.0, 10.0, 2.0, 8.0, 5.0, 6.0, 3.0, 8.0, 2.0, 3.0, 8.0, 9.0, 2.0, 7.0, 8.0, 5.0, 5.0, 3.0, 10.0, 5.0, 4.0, 2.0, 4.0, 9.0, 5.0, 9.0, 9.0, 5.0, 10.0, 5.0, 8.0, 1.0, 3.0, 4.0, 9.0, 9.0, 7.0, 6.0, 6.0, 7.0, 2.0, 6.0, 10.0, 7.0, 2.0, 7.0, 2.0, 3.0, 9.0, 7.0, 5.0, 7.0, 1.0, 6.0, 10.0, 8.0, 10.0, 3.0, 7.0, 4.0, 1.0, 6.0, 4.0, 5.0, 5.0, 10.0, 1.0, 4.0, 2.0, 4.0, 3.0, 5.0, 4.0, 7.0, 9.0, 4.0, 2.0, 4.0, 10.0, 1.0, 5.0, 9.0, 1.0, 6.0, 6.0, 5.0, 1.0, 3.0, 9.0, 2.0, 8.0, 2.0, 1.0, 5.0, 1.0, 2.0, 9.0, 9.0, 6.0, 4.0, 10.0, 2.0, 7.0, 7.0, 1.0]
global b_y = 10
global p = [0.165, 0.704, 0.67, 0.717, 0.515, 0.386, 0.237, 0.381, 0.602, 0.107, 0.392, 0.54, 0.462, 0.213, 0.945, 0.63, 0.888, 0.587, 0.754, 0.787, 0.049, 0.924, 0.824, 0.788, 0.802, 0.733, 0.924, 0.543, 0.343, 0.103, 0.692, 0.903, 0.273, 0.824, 0.977, 0.266, 0.091, 0.563, 0.544, 0.381, 0.953, 0.41, 0.286, 0.19, 0.864, 0.344, 0.329, 0.246, 0.821, 0.619, 0.672, 0.268, 0.238, 0.085, 0.994, 0.339, 0.954, 0.903, 0.111, 0.832, 0.558, 0.935, 0.593, 0.153, 0.55, 0.143, 0.306, 0.11, 0.66, 0.499, 0.772, 0.578, 0.974, 0.153, 0.059, 0.258, 0.24, 0.944, 0.801, 0.124, 0.247, 0.933, 0.274, 0.323, 0.51, 0.371, 0.183, 0.322, 0.283, 0.695, 0.875, 0.189, 0.744, 0.392, 0.022, 0.296, 0.648, 0.499, 0.845, 0.217, 0.931, 0.899, 0.695, 0.016, 0.367, 0.316, 0.69, 0.938, 0.193, 0.33, 0.706, 0.079, 0.998, 0.095, 0.693, 0.433, 0.789, 0.977, 0.284, 0.167, 0.208, 0.234, 0.261, 0.134, 0.127, 0.005, 0.972, 0.709, 0.495, 0.967, 0.228, 0.715, 0.446, 0.098, 0.172, 0.443, 0.872, 0.39, 0.28, 0.965, 0.274, 0.327, 0.352, 0.848, 0.412, 0.274, 0.018, 0.356, 0.719, 0.938, 0.898, 0.991, 0.624, 0.789, 0.189, 0.656, 0.771, 0.092, 0.119, 0.757, 0.264, 0.294, 0.842, 0.712, 0.836, 0.599, 0.256, 0.748, 0.492, 0.949, 0.844, 0.562, 0.623, 0.476, 0.41, 0.95, 0.386, 0.063, 0.433, 0.299, 0.312, 0.625, 0.816, 0.498, 0.5, 0.192, 0.398, 0.262, 0.041, 0.101, 0.729, 0.618, 0.47, 0.127, 0.734, 0.063, 0.387, 0.935, 0.425, 0.485, 0.071, 0.173, 0.621, 0.82, 0.668, 0.358, 0.361, 0.04, 0.49, 0.95, 0.778, 0.677, 0.323, 0.367, 0.724, 0.186, 0.859, 0.405, 0.766, 0.152, 0.937, 0.595, 0.968, 0.941, 0.834, 0.375, 0.269, 0.244, 0.636, 0.666, 0.654, 0.285, 0.261, 0.215, 0.613, 0.932, 0.239, 0.622, 0.049, 0.674, 0.455, 0.022, 0.42, 0.222, 0.782, 0.043, 0.891, 0.992, 0.198, 0.759, 0.645, 0.957, 0.599, 0.51, 0.548, 0.332]
global q = [0.795, 0.74, 0.701, 0.872, 0.733, 0.776, 0.321, 0.517, 0.651, 0.734, 0.873, 0.946, 0.743, 0.701, 0.964, 0.67, 0.92, 0.965, 0.812, 0.825, 0.819, 0.967, 0.999, 0.914, 0.924, 0.752, 0.977, 0.655, 0.709, 0.652, 0.911, 0.977, 0.284, 0.87, 0.986, 0.484, 0.817, 0.855, 0.912, 0.599, 0.988, 0.435, 0.46, 0.847, 0.951, 0.391, 0.669, 0.692, 0.973, 0.878, 0.865, 0.272, 0.688, 0.449, 0.997, 0.796, 0.981, 0.905, 0.952, 0.926, 0.619, 0.94, 0.876, 0.614, 0.917, 0.333, 0.549, 0.387, 0.73, 0.728, 0.88, 0.746, 0.987, 0.212, 0.27, 0.434, 0.287, 0.972, 0.859, 0.343, 0.25, 0.954, 0.949, 0.596, 0.58, 0.604, 0.625, 0.64, 0.804, 0.958, 0.938, 0.42, 0.843, 0.747, 0.696, 0.707, 0.921, 0.8, 0.969, 0.865, 0.967, 0.945, 0.714, 0.469, 0.761, 0.449, 0.782, 0.941, 0.301, 0.738, 0.992, 0.762, 0.999, 0.692, 0.698, 0.817, 0.855, 0.998, 0.842, 0.608, 0.721, 0.658, 0.601, 0.245, 0.687, 0.856, 0.996, 0.977, 0.597, 0.995, 0.567, 0.729, 0.859, 0.912, 0.496, 0.85, 0.91, 0.923, 0.76, 0.992, 0.824, 0.633, 0.403, 0.9, 0.877, 0.679, 0.475, 0.563, 0.742, 0.954, 0.955, 0.999, 0.731, 0.872, 0.634, 0.885, 0.829, 0.79, 0.304, 0.997, 0.509, 0.543, 0.952, 0.883, 0.971, 0.665, 0.324, 0.758, 0.769, 0.967, 0.902, 0.861, 0.918, 0.594, 0.579, 0.993, 0.686, 0.205, 0.659, 0.683, 0.714, 0.945, 0.928, 0.865, 0.798, 0.551, 0.955, 0.411, 0.461, 0.989, 0.897, 0.884, 0.966, 0.717, 0.836, 0.714, 0.582, 0.97, 0.802, 0.938, 0.46, 0.232, 0.76, 0.867, 0.955, 0.68, 0.759, 0.792, 0.639, 0.983, 0.957, 0.889, 0.85, 0.865, 0.986, 0.28, 0.956, 0.777, 0.912, 0.203, 0.942, 0.941, 0.968, 0.987, 0.88, 0.414, 0.9, 0.697, 0.849, 0.68, 0.972, 0.603, 0.345, 0.87, 0.881, 0.968, 0.735, 0.991, 0.331, 0.791, 0.687, 0.518, 0.996, 0.261, 0.806, 0.71, 0.927, 0.996, 0.992, 0.855, 0.964, 0.984, 0.728, 0.817, 0.927, 0.969]
global origin = 1
global destination = 50