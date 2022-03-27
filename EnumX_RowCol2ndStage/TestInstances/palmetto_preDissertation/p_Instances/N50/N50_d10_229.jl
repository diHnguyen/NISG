global arcs = [1 5; 1 8; 1 10; 1 16; 1 36; 1 39; 2 10; 2 13; 2 20; 2 34; 2 39; 2 43; 3 2; 3 6; 3 21; 3 22; 3 27; 3 30; 3 34; 3 49; 4 9; 4 15; 4 27; 4 36; 4 38; 5 8; 5 9; 5 29; 6 11; 6 15; 6 23; 6 44; 7 8; 7 25; 7 31; 7 32; 8 15; 8 34; 8 37; 8 49; 9 7; 9 12; 9 15; 9 26; 10 8; 10 42; 11 5; 11 10; 11 18; 11 28; 11 32; 11 41; 12 28; 12 29; 12 32; 13 7; 13 9; 13 22; 13 26; 13 27; 13 37; 13 45; 14 7; 14 26; 14 33; 14 43; 15 2; 15 11; 15 21; 15 22; 15 49; 16 38; 16 48; 17 18; 17 19; 17 31; 18 7; 18 16; 18 20; 18 21; 18 38; 19 4; 19 16; 19 27; 19 39; 19 41; 19 43; 19 45; 20 7; 20 29; 20 49; 20 50; 21 45; 22 10; 22 15; 22 25; 22 40; 22 41; 22 47; 23 11; 23 25; 23 29; 23 32; 24 20; 25 5; 25 10; 25 16; 26 8; 26 9; 26 21; 26 28; 26 30; 26 47; 26 49; 27 8; 27 32; 27 42; 28 4; 28 9; 28 10; 28 14; 28 16; 28 39; 28 44; 29 6; 29 25; 29 31; 29 35; 29 37; 29 43; 29 46; 30 9; 30 14; 30 19; 30 26; 30 43; 30 47; 30 50; 31 2; 31 10; 31 16; 31 37; 31 45; 32 7; 32 11; 32 24; 32 36; 32 37; 33 4; 33 11; 33 24; 33 48; 34 7; 34 16; 34 19; 34 39; 35 14; 35 17; 35 21; 35 24; 35 26; 35 30; 36 17; 36 18; 36 21; 36 32; 36 39; 37 8; 37 9; 37 13; 37 16; 37 40; 37 43; 37 47; 38 18; 38 23; 38 30; 39 18; 39 21; 39 25; 39 40; 39 44; 39 46; 40 2; 40 11; 40 20; 40 22; 40 41; 40 47; 40 49; 41 23; 41 44; 41 46; 42 14; 42 18; 42 28; 42 46; 43 5; 43 29; 44 25; 44 40; 44 48; 45 10; 45 34; 45 36; 45 37; 46 3; 46 6; 46 9; 46 23; 46 38; 46 50; 47 28; 47 34; 48 14; 48 17; 48 20; 48 23; 48 33; 48 36; 48 38; 48 41; 49 9; 49 27; 49 34; 49 35; 49 38; 49 46]
global d_x = [8.0, 10.0, 10.0, 3.0, 8.0, 1.0, 3.0, 7.0, 7.0, 5.0, 1.0, 2.0, 3.0, 6.0, 10.0, 3.0, 1.0, 1.0, 1.0, 5.0, 3.0, 1.0, 2.0, 4.0, 7.0, 5.0, 4.0, 4.0, 6.0, 8.0, 8.0, 9.0, 8.0, 2.0, 8.0, 1.0, 3.0, 8.0, 1.0, 6.0, 8.0, 4.0, 5.0, 9.0, 8.0, 2.0, 7.0, 10.0, 8.0, 1.0, 4.0, 6.0, 9.0, 5.0, 4.0, 2.0, 6.0, 2.0, 5.0, 4.0, 5.0, 3.0, 5.0, 10.0, 7.0, 2.0, 2.0, 9.0, 4.0, 2.0, 3.0, 10.0, 3.0, 8.0, 3.0, 3.0, 10.0, 6.0, 6.0, 3.0, 9.0, 9.0, 6.0, 4.0, 2.0, 8.0, 3.0, 4.0, 6.0, 6.0, 3.0, 5.0, 8.0, 2.0, 5.0, 5.0, 5.0, 7.0, 5.0, 8.0, 7.0, 5.0, 3.0, 6.0, 4.0, 3.0, 8.0, 1.0, 10.0, 8.0, 7.0, 10.0, 8.0, 10.0, 5.0, 4.0, 9.0, 4.0, 2.0, 5.0, 7.0, 10.0, 5.0, 10.0, 10.0, 9.0, 8.0, 3.0, 8.0, 6.0, 5.0, 6.0, 2.0, 9.0, 5.0, 8.0, 2.0, 9.0, 9.0, 7.0, 3.0, 5.0, 3.0, 9.0, 7.0, 8.0, 4.0, 7.0, 9.0, 1.0, 7.0, 4.0, 6.0, 10.0, 3.0, 4.0, 6.0, 8.0, 1.0, 3.0, 2.0, 2.0, 10.0, 8.0, 3.0, 4.0, 10.0, 2.0, 5.0, 1.0, 1.0, 10.0, 8.0, 6.0, 10.0, 2.0, 8.0, 10.0, 8.0, 8.0, 2.0, 1.0, 4.0, 3.0, 6.0, 6.0, 2.0, 7.0, 6.0, 1.0, 9.0, 2.0, 7.0, 7.0, 9.0, 1.0, 1.0, 6.0, 7.0, 2.0, 9.0, 1.0, 6.0, 3.0, 1.0, 4.0, 3.0, 9.0, 1.0, 8.0, 8.0, 10.0, 5.0, 6.0, 1.0, 3.0, 7.0, 5.0, 4.0, 8.0, 4.0, 3.0, 1.0, 10.0, 6.0, 3.0, 1.0, 9.0]
global b_x = 5
global d_y = [8.0, 9.0, 1.0, 2.0, 10.0, 10.0, 8.0, 4.0, 6.0, 7.0, 1.0, 2.0, 9.0, 5.0, 1.0, 8.0, 8.0, 10.0, 7.0, 10.0, 6.0, 7.0, 10.0, 4.0, 6.0, 6.0, 7.0, 5.0, 6.0, 9.0, 6.0, 2.0, 5.0, 4.0, 9.0, 6.0, 4.0, 10.0, 3.0, 10.0, 3.0, 7.0, 7.0, 9.0, 3.0, 10.0, 3.0, 1.0, 10.0, 5.0, 6.0, 3.0, 6.0, 1.0, 5.0, 4.0, 8.0, 10.0, 2.0, 10.0, 8.0, 2.0, 6.0, 9.0, 8.0, 2.0, 8.0, 6.0, 2.0, 8.0, 7.0, 1.0, 10.0, 7.0, 4.0, 9.0, 10.0, 7.0, 10.0, 7.0, 6.0, 10.0, 2.0, 3.0, 10.0, 10.0, 10.0, 3.0, 7.0, 3.0, 9.0, 10.0, 10.0, 7.0, 4.0, 3.0, 8.0, 2.0, 6.0, 5.0, 1.0, 7.0, 4.0, 3.0, 3.0, 9.0, 6.0, 2.0, 8.0, 1.0, 5.0, 8.0, 10.0, 10.0, 6.0, 1.0, 5.0, 2.0, 7.0, 4.0, 7.0, 10.0, 2.0, 3.0, 10.0, 6.0, 4.0, 10.0, 9.0, 8.0, 2.0, 9.0, 6.0, 4.0, 3.0, 5.0, 7.0, 2.0, 10.0, 9.0, 3.0, 6.0, 4.0, 1.0, 2.0, 7.0, 1.0, 5.0, 3.0, 10.0, 5.0, 9.0, 7.0, 10.0, 7.0, 5.0, 6.0, 1.0, 4.0, 5.0, 10.0, 6.0, 3.0, 2.0, 3.0, 5.0, 4.0, 5.0, 8.0, 6.0, 4.0, 6.0, 8.0, 1.0, 5.0, 4.0, 7.0, 5.0, 7.0, 8.0, 4.0, 9.0, 6.0, 6.0, 4.0, 7.0, 1.0, 7.0, 5.0, 4.0, 9.0, 8.0, 1.0, 8.0, 8.0, 6.0, 2.0, 8.0, 2.0, 6.0, 10.0, 1.0, 9.0, 2.0, 8.0, 6.0, 5.0, 1.0, 1.0, 2.0, 10.0, 3.0, 5.0, 9.0, 2.0, 3.0, 1.0, 5.0, 8.0, 3.0, 6.0, 10.0, 1.0, 5.0, 6.0, 6.0, 2.0, 1.0]
global b_y = 10
global p = [0.153, 0.531, 0.176, 0.377, 0.75, 0.099, 0.685, 0.26, 0.808, 0.926, 0.324, 0.318, 0.328, 0.341, 0.634, 0.364, 0.416, 0.363, 0.146, 0.449, 0.146, 0.368, 0.054, 0.237, 0.119, 0.134, 0.091, 0.854, 0.717, 0.243, 0.65, 0.476, 0.317, 0.073, 0.971, 0.258, 0.971, 0.799, 0.951, 0.081, 0.443, 0.485, 0.36, 0.238, 0.343, 0.071, 0.485, 0.156, 0.554, 0.547, 0.126, 0.785, 0.847, 0.129, 0.44, 0.003, 0.115, 0.366, 0.428, 0.954, 0.155, 0.935, 0.122, 0.21, 0.967, 0.812, 0.538, 0.845, 0.965, 0.77, 0.21, 0.446, 0.449, 0.542, 0.452, 0.488, 0.8, 0.449, 0.973, 0.727, 0.854, 0.537, 0.794, 0.807, 0.494, 0.535, 0.147, 0.439, 0.271, 0.689, 0.227, 0.336, 0.983, 0.764, 0.905, 0.907, 0.067, 0.083, 0.374, 0.977, 0.662, 0.734, 0.803, 0.828, 0.246, 0.232, 0.155, 0.277, 0.16, 0.194, 0.043, 0.83, 0.25, 0.875, 0.356, 0.657, 0.857, 0.431, 0.4, 0.245, 0.171, 0.382, 0.93, 0.837, 0.051, 0.528, 0.883, 0.444, 0.179, 0.516, 0.36, 0.978, 0.787, 0.099, 0.499, 0.545, 0.644, 0.119, 0.736, 0.547, 0.823, 0.281, 0.907, 0.263, 0.363, 0.324, 0.169, 0.238, 0.191, 0.301, 0.893, 0.24, 0.422, 0.839, 0.358, 0.647, 0.202, 0.747, 0.851, 0.759, 0.467, 0.671, 0.904, 0.353, 0.452, 0.807, 0.183, 0.495, 0.251, 0.191, 0.559, 0.657, 0.7, 0.189, 0.84, 0.06, 0.823, 0.256, 0.747, 0.75, 0.876, 0.407, 0.197, 0.787, 0.714, 0.313, 0.105, 0.827, 0.901, 0.09, 0.813, 0.548, 0.228, 0.075, 0.89, 0.96, 0.872, 0.395, 0.982, 0.393, 0.861, 0.723, 0.33, 0.388, 0.59, 0.532, 0.366, 0.144, 0.144, 0.075, 0.814, 0.88, 0.474, 0.53, 0.728, 0.038, 0.634, 0.41, 0.415, 0.755, 0.467, 0.844, 0.614, 0.322, 0.538, 0.658, 0.815, 0.386]
global q = [0.453, 0.951, 0.803, 0.703, 0.824, 0.967, 0.806, 0.514, 0.907, 0.977, 0.765, 0.571, 0.572, 0.651, 0.669, 0.745, 0.611, 0.383, 0.639, 0.578, 0.568, 0.808, 0.145, 0.24, 0.414, 0.562, 0.283, 0.959, 0.854, 0.543, 0.707, 0.608, 0.62, 0.982, 0.994, 0.554, 0.974, 0.822, 0.968, 0.578, 0.55, 0.601, 0.99, 0.936, 0.523, 0.608, 0.489, 0.693, 0.871, 0.692, 0.395, 0.789, 0.855, 0.199, 0.642, 0.316, 0.821, 0.511, 0.432, 0.983, 0.266, 0.99, 0.541, 0.718, 0.994, 0.944, 0.696, 0.857, 0.971, 0.988, 0.385, 0.653, 0.565, 0.67, 0.671, 0.916, 0.879, 0.539, 0.997, 0.894, 0.957, 0.705, 0.855, 0.968, 0.877, 0.644, 0.84, 0.947, 0.844, 0.942, 0.469, 0.821, 0.992, 0.879, 0.928, 0.949, 0.882, 0.926, 0.441, 0.986, 0.948, 0.958, 0.925, 0.867, 0.64, 0.434, 0.679, 0.508, 0.917, 0.757, 0.733, 0.937, 0.737, 0.916, 0.359, 0.74, 0.979, 0.604, 0.892, 0.585, 0.634, 0.865, 0.993, 0.898, 0.092, 0.834, 0.969, 0.9, 0.996, 0.583, 0.795, 0.997, 0.961, 0.954, 0.657, 0.687, 0.859, 0.388, 0.744, 0.625, 0.94, 0.805, 0.999, 0.618, 0.478, 0.888, 0.865, 0.326, 0.857, 0.85, 0.919, 0.79, 0.596, 0.893, 0.466, 0.781, 0.406, 0.927, 0.925, 0.915, 0.83, 0.758, 0.979, 0.494, 0.579, 0.96, 0.33, 0.757, 0.875, 0.863, 0.97, 0.719, 0.78, 0.466, 0.875, 0.441, 0.827, 0.414, 0.972, 0.944, 0.962, 0.937, 0.702, 0.864, 0.847, 0.442, 0.555, 0.947, 0.905, 0.23, 0.915, 0.99, 0.508, 0.198, 0.997, 0.971, 0.99, 0.692, 0.997, 0.905, 0.896, 0.922, 0.824, 0.813, 0.858, 0.671, 0.743, 0.825, 0.319, 0.302, 0.821, 0.922, 0.909, 0.762, 0.824, 0.792, 0.918, 0.879, 0.918, 0.984, 0.96, 0.918, 0.85, 0.657, 0.734, 0.754, 0.919, 0.72]
global origin = 1
global destination = 50