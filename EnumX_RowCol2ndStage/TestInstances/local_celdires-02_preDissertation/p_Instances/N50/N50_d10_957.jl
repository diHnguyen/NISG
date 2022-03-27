global arcs = [1 6; 1 7; 1 15; 1 25; 1 28; 1 29; 1 30; 1 36; 1 45; 2 9; 2 18; 2 19; 2 37; 3 7; 3 16; 3 23; 3 28; 4 2; 4 10; 4 16; 4 18; 4 21; 4 30; 4 32; 4 39; 5 4; 5 15; 5 16; 5 25; 5 34; 5 46; 6 4; 6 7; 6 10; 6 12; 6 25; 6 41; 7 3; 7 22; 7 23; 7 42; 7 49; 8 30; 9 15; 9 22; 9 31; 9 47; 10 7; 11 3; 11 9; 11 10; 11 16; 11 23; 11 29; 11 32; 11 45; 12 18; 12 20; 13 5; 13 10; 13 14; 13 17; 13 36; 14 4; 14 17; 15 8; 15 11; 15 17; 15 35; 15 36; 15 37; 15 42; 16 4; 16 7; 16 11; 16 42; 17 5; 17 8; 17 16; 17 21; 17 24; 17 33; 18 4; 18 8; 18 11; 18 33; 18 37; 19 8; 20 3; 20 5; 20 15; 20 35; 20 45; 21 9; 21 16; 21 30; 21 32; 21 37; 21 38; 21 49; 22 13; 22 36; 22 37; 22 40; 23 21; 23 27; 23 47; 24 3; 24 8; 24 23; 24 32; 24 50; 25 5; 25 10; 25 14; 25 26; 25 36; 26 2; 26 3; 26 45; 27 16; 27 25; 27 26; 27 40; 27 41; 28 26; 28 36; 28 43; 28 48; 29 17; 29 22; 29 24; 29 42; 30 14; 30 28; 30 32; 30 47; 31 2; 31 11; 31 19; 31 34; 31 39; 32 4; 32 5; 32 20; 32 28; 32 37; 33 4; 33 11; 33 19; 33 21; 33 26; 34 2; 34 4; 34 11; 34 14; 34 29; 34 32; 34 50; 35 14; 35 23; 35 36; 35 46; 36 9; 36 26; 36 46; 37 9; 37 21; 37 23; 37 29; 37 48; 38 2; 38 9; 38 14; 38 29; 38 46; 38 48; 39 7; 39 44; 39 48; 40 2; 40 29; 40 38; 40 41; 40 50; 41 17; 41 20; 41 34; 41 43; 42 13; 42 39; 43 34; 43 40; 43 47; 44 3; 44 8; 44 9; 44 10; 44 18; 44 45; 45 2; 45 6; 45 8; 45 10; 45 12; 45 15; 45 28; 45 30; 45 36; 46 2; 46 11; 46 47; 47 5; 47 7; 47 8; 47 9; 47 31; 48 14; 48 16; 48 27; 48 29; 48 37; 48 38; 48 40; 48 42; 48 44; 49 17; 49 18; 49 24; 49 43]
global d_x = [2.0, 5.0, 2.0, 2.0, 9.0, 4.0, 10.0, 7.0, 3.0, 4.0, 7.0, 7.0, 5.0, 6.0, 8.0, 3.0, 2.0, 6.0, 10.0, 5.0, 3.0, 6.0, 1.0, 8.0, 10.0, 6.0, 7.0, 2.0, 8.0, 5.0, 2.0, 5.0, 10.0, 1.0, 4.0, 4.0, 10.0, 3.0, 1.0, 1.0, 7.0, 10.0, 8.0, 10.0, 7.0, 2.0, 6.0, 10.0, 6.0, 3.0, 4.0, 8.0, 7.0, 7.0, 2.0, 1.0, 4.0, 7.0, 2.0, 3.0, 4.0, 5.0, 3.0, 7.0, 4.0, 8.0, 3.0, 5.0, 7.0, 9.0, 6.0, 9.0, 6.0, 5.0, 5.0, 9.0, 3.0, 5.0, 5.0, 3.0, 9.0, 7.0, 8.0, 5.0, 10.0, 5.0, 3.0, 10.0, 2.0, 5.0, 1.0, 5.0, 5.0, 6.0, 2.0, 3.0, 3.0, 4.0, 1.0, 6.0, 2.0, 6.0, 2.0, 9.0, 5.0, 1.0, 9.0, 10.0, 1.0, 2.0, 7.0, 5.0, 1.0, 6.0, 9.0, 3.0, 9.0, 2.0, 3.0, 6.0, 9.0, 3.0, 3.0, 8.0, 4.0, 8.0, 9.0, 1.0, 2.0, 2.0, 5.0, 10.0, 8.0, 1.0, 9.0, 7.0, 6.0, 7.0, 5.0, 7.0, 10.0, 3.0, 8.0, 8.0, 10.0, 2.0, 3.0, 1.0, 3.0, 5.0, 3.0, 5.0, 8.0, 2.0, 6.0, 2.0, 2.0, 7.0, 9.0, 2.0, 3.0, 4.0, 1.0, 1.0, 3.0, 2.0, 5.0, 1.0, 5.0, 6.0, 10.0, 4.0, 5.0, 1.0, 3.0, 3.0, 1.0, 2.0, 1.0, 8.0, 6.0, 8.0, 10.0, 5.0, 5.0, 9.0, 4.0, 4.0, 5.0, 1.0, 10.0, 5.0, 8.0, 6.0, 9.0, 3.0, 4.0, 8.0, 4.0, 10.0, 6.0, 3.0, 3.0, 4.0, 7.0, 8.0, 8.0, 3.0, 10.0, 1.0, 6.0, 4.0, 1.0, 3.0, 1.0, 5.0, 6.0, 9.0, 1.0, 6.0, 1.0, 10.0, 6.0, 2.0, 2.0, 10.0, 10.0, 10.0, 1.0, 1.0]
global b_x = 5
global d_y = [3.0, 6.0, 6.0, 8.0, 1.0, 4.0, 4.0, 5.0, 9.0, 10.0, 8.0, 2.0, 7.0, 2.0, 5.0, 6.0, 6.0, 7.0, 8.0, 7.0, 3.0, 7.0, 5.0, 4.0, 4.0, 8.0, 7.0, 1.0, 9.0, 8.0, 4.0, 2.0, 8.0, 1.0, 10.0, 4.0, 10.0, 2.0, 8.0, 9.0, 7.0, 9.0, 1.0, 3.0, 7.0, 9.0, 8.0, 6.0, 1.0, 6.0, 7.0, 8.0, 1.0, 3.0, 8.0, 1.0, 4.0, 7.0, 3.0, 8.0, 6.0, 4.0, 2.0, 4.0, 4.0, 9.0, 2.0, 5.0, 5.0, 5.0, 9.0, 3.0, 2.0, 1.0, 4.0, 9.0, 4.0, 6.0, 10.0, 4.0, 6.0, 9.0, 7.0, 4.0, 5.0, 3.0, 7.0, 1.0, 7.0, 10.0, 5.0, 4.0, 10.0, 7.0, 4.0, 8.0, 9.0, 9.0, 9.0, 6.0, 7.0, 4.0, 8.0, 6.0, 6.0, 6.0, 6.0, 5.0, 6.0, 9.0, 4.0, 7.0, 7.0, 6.0, 2.0, 8.0, 3.0, 9.0, 8.0, 9.0, 3.0, 7.0, 5.0, 3.0, 2.0, 5.0, 7.0, 7.0, 7.0, 4.0, 1.0, 1.0, 4.0, 3.0, 5.0, 8.0, 3.0, 9.0, 10.0, 5.0, 4.0, 4.0, 5.0, 10.0, 3.0, 3.0, 4.0, 6.0, 7.0, 2.0, 2.0, 3.0, 4.0, 7.0, 10.0, 4.0, 6.0, 6.0, 2.0, 1.0, 4.0, 2.0, 6.0, 7.0, 6.0, 9.0, 10.0, 8.0, 8.0, 3.0, 8.0, 3.0, 5.0, 10.0, 2.0, 2.0, 10.0, 9.0, 1.0, 9.0, 8.0, 4.0, 8.0, 3.0, 5.0, 5.0, 5.0, 6.0, 4.0, 4.0, 9.0, 4.0, 2.0, 4.0, 7.0, 10.0, 8.0, 10.0, 2.0, 9.0, 2.0, 9.0, 9.0, 9.0, 9.0, 8.0, 5.0, 6.0, 6.0, 7.0, 7.0, 5.0, 5.0, 2.0, 3.0, 6.0, 7.0, 5.0, 6.0, 7.0, 8.0, 6.0, 5.0, 6.0, 9.0, 2.0, 10.0, 5.0, 5.0, 8.0]
global b_y = 10
global p = [0.376, 0.311, 0.494, 0.055, 0.837, 0.217, 0.714, 0.923, 0.417, 0.875, 0.167, 0.749, 0.003, 0.135, 0.172, 0.791, 0.274, 0.437, 0.099, 0.374, 0.271, 0.904, 0.043, 0.57, 0.323, 0.429, 0.233, 0.273, 0.3, 0.739, 0.835, 0.712, 0.278, 0.452, 0.815, 0.867, 0.25, 0.268, 0.47, 0.536, 0.259, 0.865, 0.349, 0.893, 0.05, 0.415, 0.392, 0.349, 0.241, 0.765, 0.752, 0.421, 0.885, 0.997, 0.044, 0.543, 0.11, 0.593, 0.078, 0.08, 0.376, 0.686, 0.673, 0.747, 0.669, 0.804, 0.086, 0.652, 0.818, 0.546, 0.223, 0.213, 0.894, 0.461, 0.751, 0.223, 0.672, 0.041, 0.827, 0.607, 0.859, 0.284, 0.178, 0.711, 0.432, 0.641, 0.238, 0.656, 0.057, 0.356, 0.598, 0.971, 0.192, 0.869, 0.41, 0.494, 0.004, 0.957, 0.144, 0.571, 0.477, 0.543, 0.974, 0.4, 0.336, 0.588, 0.994, 0.992, 0.73, 0.949, 0.18, 0.989, 0.091, 0.176, 0.984, 0.519, 0.333, 0.11, 0.086, 0.147, 0.472, 0.816, 0.29, 0.967, 0.285, 0.05, 0.862, 0.843, 0.064, 0.159, 0.246, 0.671, 0.429, 0.254, 0.445, 0.632, 0.925, 0.602, 0.071, 0.267, 0.319, 0.58, 0.556, 0.057, 0.466, 0.039, 0.976, 0.542, 0.875, 0.446, 0.839, 0.716, 0.624, 0.557, 0.991, 0.654, 0.583, 0.513, 0.542, 0.725, 0.632, 0.922, 0.122, 0.049, 0.248, 0.163, 0.214, 0.864, 0.066, 0.043, 0.193, 0.357, 0.848, 0.409, 0.925, 0.854, 0.158, 0.742, 0.679, 0.579, 0.306, 0.67, 0.734, 0.404, 0.291, 0.023, 0.205, 0.666, 0.802, 0.38, 0.105, 0.633, 0.43, 0.223, 0.979, 0.534, 0.264, 0.791, 0.669, 0.808, 0.242, 0.104, 0.14, 0.991, 0.067, 0.562, 0.675, 0.927, 0.677, 0.281, 0.372, 0.094, 0.05, 0.956, 0.044, 0.957, 0.423, 0.641, 0.088, 0.829, 0.188, 0.938, 0.898, 0.453, 0.547, 0.179, 0.003, 0.181, 0.577, 0.166]
global q = [0.65, 0.678, 0.898, 0.756, 0.946, 0.332, 0.762, 0.931, 0.599, 0.991, 0.697, 0.89, 0.375, 0.655, 0.928, 0.837, 0.603, 0.76, 0.276, 0.426, 0.882, 0.97, 0.471, 0.735, 0.469, 0.555, 0.926, 0.684, 0.899, 0.993, 0.935, 0.805, 0.817, 0.64, 0.927, 0.982, 0.706, 0.708, 0.595, 0.974, 0.429, 0.889, 0.916, 0.909, 0.687, 0.787, 0.74, 0.811, 0.535, 0.961, 0.972, 0.885, 0.999, 0.997, 0.775, 0.8, 0.272, 0.708, 0.673, 0.957, 0.74, 0.913, 0.687, 0.793, 0.957, 0.983, 0.907, 0.717, 0.863, 0.967, 0.791, 0.748, 0.919, 0.872, 0.845, 0.225, 0.912, 0.486, 0.978, 0.718, 0.864, 0.639, 0.706, 0.846, 0.826, 0.749, 0.884, 0.663, 0.464, 0.451, 0.852, 0.972, 0.607, 0.919, 0.912, 0.719, 0.359, 0.971, 0.511, 0.625, 0.49, 0.929, 0.975, 0.423, 0.397, 0.689, 0.995, 0.995, 0.987, 0.967, 0.842, 0.995, 0.688, 0.338, 0.997, 0.753, 0.963, 0.901, 0.442, 0.659, 0.959, 0.964, 0.433, 0.987, 0.302, 0.723, 0.888, 0.897, 0.781, 0.845, 0.701, 0.842, 0.817, 0.681, 0.986, 0.691, 0.994, 0.778, 0.257, 0.513, 0.508, 0.837, 0.896, 0.344, 0.634, 0.536, 0.983, 0.83, 0.896, 0.971, 0.858, 0.908, 0.92, 0.827, 0.993, 0.798, 0.642, 0.659, 0.858, 0.969, 0.95, 0.958, 0.773, 0.101, 0.462, 0.809, 0.408, 0.943, 0.385, 0.656, 0.196, 0.465, 0.951, 0.855, 0.988, 0.996, 0.39, 0.798, 0.766, 0.661, 0.409, 0.978, 0.832, 0.796, 0.687, 0.594, 0.376, 0.876, 0.942, 0.728, 0.356, 0.649, 0.699, 0.645, 0.989, 0.833, 0.74, 0.811, 0.769, 0.943, 0.878, 0.6, 0.362, 0.996, 0.788, 0.797, 0.879, 0.956, 0.711, 0.603, 0.815, 0.701, 0.118, 0.966, 0.947, 0.974, 0.796, 0.903, 0.547, 0.952, 0.89, 0.945, 0.929, 0.496, 0.958, 0.74, 0.383, 0.205, 0.772, 0.34]
global origin = 1
global destination = 50