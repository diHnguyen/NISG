global arcs = [1 20; 1 28; 1 31; 1 38; 1 43; 1 47; 2 4; 2 7; 2 27; 2 40; 3 23; 4 6; 4 12; 4 13; 4 19; 5 9; 5 18; 6 2; 6 11; 6 15; 6 21; 6 37; 6 44; 7 33; 7 39; 7 42; 7 46; 8 17; 8 21; 9 3; 9 8; 9 24; 9 35; 9 46; 9 50; 10 18; 10 45; 10 48; 11 18; 11 24; 11 41; 12 14; 12 18; 12 30; 12 36; 12 38; 13 2; 13 4; 13 19; 13 25; 13 30; 13 31; 13 36; 13 38; 13 47; 14 16; 14 20; 14 27; 14 29; 14 45; 15 4; 15 7; 15 19; 15 26; 15 34; 15 35; 15 37; 15 45; 15 47; 16 11; 16 40; 16 43; 16 47; 16 50; 17 21; 17 28; 17 44; 17 45; 17 46; 18 4; 18 5; 18 7; 18 11; 18 43; 18 47; 19 16; 19 18; 19 21; 19 29; 19 31; 19 32; 19 37; 20 3; 20 13; 20 17; 20 27; 20 31; 21 5; 21 25; 21 34; 21 46; 22 10; 22 24; 22 28; 22 37; 22 42; 23 7; 23 16; 23 22; 23 30; 23 31; 23 48; 23 49; 24 16; 24 19; 24 26; 24 31; 24 49; 25 10; 25 17; 25 38; 25 45; 26 2; 26 8; 26 13; 26 15; 26 21; 26 27; 26 37; 26 46; 27 22; 27 46; 27 50; 28 18; 28 37; 28 49; 29 9; 29 20; 29 23; 29 37; 29 45; 30 17; 30 21; 30 22; 30 27; 30 33; 30 42; 30 48; 31 12; 31 14; 31 23; 31 29; 31 36; 31 39; 32 6; 32 31; 32 38; 32 42; 33 43; 34 2; 34 3; 34 11; 34 14; 34 30; 34 31; 34 44; 34 48; 35 19; 35 34; 35 42; 35 46; 36 24; 36 34; 36 35; 36 46; 36 47; 37 7; 37 8; 37 16; 37 19; 37 28; 37 40; 37 45; 37 47; 38 6; 38 17; 38 29; 38 39; 38 44; 38 47; 39 6; 39 11; 39 12; 39 43; 39 49; 40 18; 40 25; 40 28; 40 33; 40 35; 40 36; 41 8; 41 23; 42 2; 42 10; 42 11; 42 18; 42 28; 42 48; 42 49; 43 8; 43 18; 43 20; 43 35; 43 45; 44 26; 44 35; 44 38; 44 43; 44 49; 45 8; 45 9; 45 15; 45 17; 45 23; 45 24; 45 27; 45 39; 45 41; 45 48; 46 6; 46 7; 46 15; 46 31; 46 34; 46 37; 47 10; 47 14; 47 15; 47 18; 47 25; 47 39; 47 45; 47 46; 48 15; 48 32; 48 35; 48 42; 48 50; 49 6; 49 18; 49 24; 49 35; 49 42; 49 45]
global d_x = [7.0, 7.0, 4.0, 5.0, 5.0, 1.0, 4.0, 9.0, 2.0, 5.0, 7.0, 3.0, 5.0, 8.0, 4.0, 10.0, 3.0, 3.0, 7.0, 8.0, 7.0, 9.0, 4.0, 1.0, 1.0, 9.0, 1.0, 9.0, 6.0, 9.0, 6.0, 1.0, 3.0, 4.0, 7.0, 2.0, 9.0, 7.0, 9.0, 3.0, 10.0, 3.0, 9.0, 4.0, 6.0, 3.0, 1.0, 2.0, 1.0, 7.0, 4.0, 4.0, 1.0, 5.0, 9.0, 4.0, 1.0, 9.0, 9.0, 2.0, 1.0, 4.0, 10.0, 8.0, 9.0, 5.0, 6.0, 6.0, 10.0, 6.0, 4.0, 4.0, 1.0, 7.0, 5.0, 5.0, 9.0, 8.0, 9.0, 8.0, 8.0, 2.0, 4.0, 7.0, 1.0, 3.0, 2.0, 4.0, 2.0, 3.0, 4.0, 9.0, 4.0, 7.0, 8.0, 1.0, 3.0, 2.0, 4.0, 2.0, 9.0, 6.0, 4.0, 6.0, 10.0, 2.0, 5.0, 2.0, 2.0, 6.0, 7.0, 8.0, 6.0, 8.0, 2.0, 5.0, 6.0, 5.0, 3.0, 6.0, 3.0, 4.0, 3.0, 4.0, 5.0, 5.0, 3.0, 7.0, 3.0, 3.0, 3.0, 9.0, 1.0, 10.0, 9.0, 9.0, 8.0, 7.0, 9.0, 5.0, 7.0, 10.0, 4.0, 2.0, 2.0, 3.0, 8.0, 3.0, 8.0, 3.0, 8.0, 8.0, 1.0, 8.0, 1.0, 10.0, 7.0, 9.0, 1.0, 3.0, 7.0, 7.0, 8.0, 1.0, 9.0, 10.0, 3.0, 1.0, 2.0, 3.0, 10.0, 6.0, 5.0, 8.0, 5.0, 1.0, 1.0, 10.0, 7.0, 5.0, 9.0, 1.0, 8.0, 1.0, 6.0, 9.0, 8.0, 5.0, 2.0, 9.0, 1.0, 8.0, 6.0, 9.0, 9.0, 1.0, 7.0, 7.0, 2.0, 4.0, 10.0, 6.0, 2.0, 3.0, 10.0, 1.0, 3.0, 7.0, 1.0, 1.0, 3.0, 8.0, 4.0, 3.0, 3.0, 2.0, 7.0, 3.0, 7.0, 9.0, 8.0, 2.0, 4.0, 10.0, 5.0, 7.0, 4.0, 4.0, 10.0, 8.0, 1.0, 5.0, 4.0, 8.0, 5.0, 10.0, 7.0, 5.0, 6.0, 1.0, 1.0, 7.0, 1.0, 7.0, 10.0, 4.0, 6.0, 8.0, 3.0, 4.0, 3.0, 1.0, 6.0, 10.0, 1.0]
global b_x = 5
global d_y = [2.0, 9.0, 9.0, 1.0, 9.0, 2.0, 8.0, 7.0, 8.0, 7.0, 4.0, 4.0, 9.0, 2.0, 6.0, 1.0, 6.0, 9.0, 1.0, 2.0, 8.0, 5.0, 4.0, 7.0, 6.0, 10.0, 7.0, 4.0, 9.0, 3.0, 6.0, 9.0, 3.0, 4.0, 6.0, 1.0, 3.0, 10.0, 4.0, 3.0, 5.0, 9.0, 3.0, 7.0, 1.0, 7.0, 6.0, 2.0, 8.0, 1.0, 3.0, 7.0, 9.0, 7.0, 3.0, 4.0, 7.0, 9.0, 7.0, 3.0, 3.0, 2.0, 6.0, 9.0, 6.0, 5.0, 10.0, 4.0, 8.0, 4.0, 8.0, 5.0, 9.0, 5.0, 9.0, 6.0, 7.0, 5.0, 9.0, 7.0, 6.0, 3.0, 4.0, 8.0, 7.0, 9.0, 2.0, 1.0, 10.0, 10.0, 10.0, 8.0, 10.0, 6.0, 9.0, 9.0, 5.0, 6.0, 1.0, 3.0, 2.0, 8.0, 6.0, 9.0, 8.0, 1.0, 7.0, 9.0, 8.0, 6.0, 4.0, 1.0, 7.0, 8.0, 3.0, 9.0, 4.0, 7.0, 10.0, 6.0, 9.0, 5.0, 1.0, 9.0, 7.0, 3.0, 3.0, 2.0, 10.0, 2.0, 2.0, 3.0, 8.0, 4.0, 1.0, 9.0, 6.0, 6.0, 3.0, 5.0, 10.0, 2.0, 6.0, 6.0, 10.0, 1.0, 9.0, 5.0, 5.0, 5.0, 1.0, 10.0, 3.0, 7.0, 10.0, 1.0, 2.0, 9.0, 10.0, 2.0, 8.0, 4.0, 6.0, 10.0, 5.0, 7.0, 6.0, 10.0, 8.0, 5.0, 10.0, 9.0, 5.0, 7.0, 4.0, 7.0, 9.0, 10.0, 1.0, 10.0, 6.0, 2.0, 5.0, 7.0, 1.0, 1.0, 5.0, 9.0, 6.0, 10.0, 7.0, 6.0, 6.0, 3.0, 10.0, 8.0, 7.0, 8.0, 5.0, 6.0, 5.0, 9.0, 6.0, 3.0, 4.0, 7.0, 9.0, 1.0, 10.0, 1.0, 4.0, 3.0, 6.0, 2.0, 1.0, 3.0, 10.0, 3.0, 7.0, 1.0, 7.0, 7.0, 2.0, 8.0, 4.0, 2.0, 2.0, 5.0, 3.0, 5.0, 1.0, 9.0, 4.0, 3.0, 8.0, 1.0, 2.0, 9.0, 9.0, 3.0, 10.0, 3.0, 9.0, 10.0, 7.0, 5.0, 10.0, 8.0, 7.0, 5.0, 6.0, 2.0, 3.0, 7.0, 10.0]
global b_y = 10
global p = [0.127, 0.455, 0.828, 0.265, 0.19, 0.408, 0.408, 0.931, 0.067, 0.941, 0.99, 0.077, 0.845, 0.803, 0.201, 0.483, 0.076, 0.037, 0.55, 0.574, 0.226, 0.9, 0.34, 0.844, 0.954, 0.199, 0.929, 0.186, 0.363, 0.752, 0.756, 0.718, 0.29, 0.978, 0.839, 0.597, 0.53, 0.207, 0.87, 0.629, 0.048, 0.512, 0.316, 0.944, 0.324, 0.295, 0.47, 0.151, 0.195, 0.85, 0.222, 0.741, 0.518, 0.814, 0.364, 0.69, 0.018, 0.069, 0.509, 0.634, 0.667, 0.734, 0.594, 0.185, 0.56, 0.759, 0.12, 0.728, 0.218, 0.933, 0.97, 0.785, 0.235, 0.925, 0.706, 0.41, 0.443, 0.5, 0.064, 0.897, 0.976, 0.516, 0.193, 0.484, 0.317, 0.692, 0.224, 0.895, 0.426, 0.942, 0.938, 0.77, 0.22, 0.351, 0.969, 0.606, 0.339, 0.298, 0.875, 0.537, 0.837, 0.029, 0.942, 0.294, 0.135, 0.71, 0.65, 0.471, 0.622, 0.743, 0.198, 0.72, 0.473, 0.791, 0.211, 0.893, 0.642, 0.753, 0.217, 0.522, 0.235, 0.145, 0.985, 0.096, 0.843, 0.036, 0.002, 0.88, 0.91, 0.052, 0.633, 0.439, 0.272, 0.207, 0.508, 0.743, 0.225, 0.71, 0.406, 0.493, 0.005, 0.476, 0.647, 0.566, 0.374, 0.085, 0.874, 0.724, 0.696, 0.875, 0.188, 0.535, 0.156, 0.966, 0.423, 0.472, 0.656, 0.054, 0.662, 0.259, 0.194, 0.473, 0.674, 0.444, 0.163, 0.291, 0.731, 0.934, 0.489, 0.734, 0.414, 0.311, 0.961, 0.354, 0.92, 0.499, 0.47, 0.004, 0.39, 0.006, 0.563, 0.6, 0.766, 0.344, 0.647, 0.089, 0.826, 0.213, 0.102, 0.82, 0.325, 0.573, 0.529, 0.405, 0.656, 0.76, 0.659, 0.709, 0.833, 0.612, 0.673, 0.114, 0.719, 0.072, 0.894, 0.281, 0.272, 0.145, 0.204, 0.675, 0.59, 0.898, 0.89, 0.209, 0.791, 0.46, 0.393, 0.751, 0.134, 0.379, 0.362, 0.447, 0.044, 0.015, 0.329, 0.851, 0.373, 0.344, 0.313, 0.303, 0.861, 0.814, 0.574, 0.375, 0.165, 0.678, 0.727, 0.667, 0.369, 0.562, 0.928, 0.767, 0.435, 0.383, 0.857, 0.508, 0.044, 0.59, 0.343, 0.433, 0.934, 0.45, 0.02, 0.271, 0.733]
global q = [0.683, 0.668, 0.877, 0.362, 0.512, 0.566, 0.866, 0.976, 0.661, 0.97, 0.993, 0.207, 0.875, 0.842, 0.42, 0.581, 0.894, 0.505, 0.927, 0.972, 0.653, 0.977, 0.567, 0.913, 0.991, 0.859, 0.964, 0.326, 0.392, 0.758, 0.996, 0.801, 0.625, 0.999, 0.889, 0.664, 0.57, 0.427, 0.944, 0.726, 0.21, 0.683, 0.513, 0.95, 0.693, 0.736, 0.488, 0.414, 0.729, 0.942, 0.788, 0.78, 0.742, 0.924, 0.727, 0.752, 0.639, 0.916, 0.802, 0.919, 0.881, 0.87, 0.828, 0.23, 0.679, 0.822, 0.849, 0.734, 0.755, 0.97, 0.974, 0.831, 0.808, 0.973, 0.758, 0.583, 0.57, 0.834, 0.964, 0.901, 0.986, 0.829, 0.616, 0.867, 0.972, 0.809, 0.755, 0.941, 0.562, 0.944, 0.944, 0.894, 0.566, 0.685, 0.998, 0.831, 0.372, 0.63, 0.913, 0.706, 0.915, 0.078, 0.976, 0.41, 0.163, 0.968, 0.963, 0.677, 0.921, 0.82, 0.47, 0.834, 0.57, 0.816, 0.845, 0.952, 0.657, 0.88, 0.747, 0.759, 0.657, 0.488, 0.992, 0.576, 0.849, 0.065, 0.881, 0.993, 0.912, 0.608, 0.808, 0.702, 0.894, 0.979, 0.787, 0.767, 0.384, 0.79, 0.851, 0.592, 0.52, 0.642, 0.793, 0.885, 0.677, 0.429, 0.951, 0.835, 0.754, 0.938, 0.242, 0.909, 0.185, 0.97, 0.529, 0.911, 0.689, 0.661, 0.746, 0.405, 0.716, 0.577, 0.942, 0.785, 0.458, 0.933, 0.773, 0.979, 0.715, 0.745, 0.607, 0.874, 0.973, 0.546, 0.951, 0.653, 0.601, 0.192, 0.57, 0.31, 0.618, 0.614, 0.932, 0.612, 0.781, 0.952, 0.93, 0.565, 0.356, 0.918, 0.62, 0.952, 0.631, 0.647, 0.835, 0.957, 0.82, 0.95, 0.932, 0.748, 0.933, 0.517, 0.907, 0.562, 0.903, 0.795, 0.807, 0.67, 0.973, 0.874, 0.818, 0.989, 0.908, 0.65, 0.903, 0.722, 0.467, 0.948, 0.73, 0.857, 0.41, 0.502, 0.202, 0.689, 0.883, 0.951, 0.463, 0.969, 0.955, 0.982, 0.942, 0.897, 0.934, 0.591, 0.432, 0.907, 0.927, 0.941, 0.703, 0.932, 0.982, 0.956, 0.787, 0.807, 0.925, 0.869, 0.58, 0.775, 0.494, 0.813, 0.965, 0.722, 0.319, 0.393, 0.931]
global origin = 1
global destination = 50