global arcs = [1 10; 1 14; 1 15; 1 32; 1 34; 1 40; 1 44; 1 50; 2 4; 2 23; 2 26; 2 43; 2 49; 3 2; 3 24; 3 28; 3 31; 3 35; 4 18; 4 30; 4 50; 5 11; 5 15; 5 29; 6 5; 6 24; 6 29; 6 42; 6 43; 7 2; 7 13; 7 16; 7 23; 7 34; 7 44; 7 46; 8 15; 8 31; 8 33; 9 14; 9 33; 9 46; 10 27; 10 34; 11 4; 11 9; 11 20; 11 21; 11 41; 11 46; 12 22; 12 24; 12 39; 13 7; 13 18; 13 24; 14 8; 14 12; 14 17; 14 24; 14 45; 14 46; 14 50; 15 3; 15 4; 15 5; 15 21; 15 23; 15 47; 15 48; 16 5; 16 20; 16 26; 16 37; 17 2; 17 11; 17 12; 17 36; 17 38; 18 20; 18 23; 18 41; 18 46; 19 5; 19 29; 19 31; 19 35; 19 37; 20 8; 20 10; 20 36; 21 7; 21 13; 21 15; 21 23; 21 25; 21 45; 21 47; 22 10; 22 21; 22 45; 23 26; 23 27; 23 42; 23 44; 23 47; 24 23; 24 37; 24 47; 25 2; 25 15; 25 28; 25 33; 25 34; 25 40; 25 49; 26 8; 26 9; 26 46; 27 4; 27 17; 27 24; 27 28; 27 40; 28 2; 28 12; 28 21; 28 23; 28 44; 29 3; 29 4; 29 8; 29 11; 30 2; 30 3; 30 24; 30 28; 30 37; 30 43; 30 50; 31 3; 31 6; 31 16; 31 27; 31 28; 31 33; 31 43; 32 3; 32 5; 32 17; 32 25; 32 38; 33 8; 33 46; 34 7; 34 14; 34 26; 35 9; 35 40; 36 2; 36 8; 36 37; 36 50; 37 4; 37 14; 37 17; 37 19; 37 20; 37 31; 37 41; 38 7; 38 26; 38 42; 38 48; 39 17; 39 31; 40 12; 40 20; 40 29; 40 30; 40 31; 40 50; 41 21; 41 30; 41 33; 41 44; 41 46; 42 7; 42 8; 42 30; 43 26; 43 28; 43 29; 44 7; 44 13; 44 18; 44 30; 45 21; 45 27; 45 35; 45 43; 45 50; 46 2; 46 11; 46 32; 46 37; 46 41; 47 4; 47 10; 47 11; 47 14; 47 39; 47 45; 48 19; 48 27; 48 37; 48 40; 49 4; 49 13; 49 14; 49 17; 49 20; 49 25; 49 31; 49 38; 49 42]
global d_x = [9.0, 3.0, 4.0, 5.0, 4.0, 3.0, 4.0, 7.0, 1.0, 5.0, 6.0, 7.0, 3.0, 4.0, 1.0, 9.0, 5.0, 7.0, 5.0, 10.0, 10.0, 5.0, 4.0, 4.0, 1.0, 6.0, 6.0, 6.0, 2.0, 5.0, 10.0, 4.0, 2.0, 5.0, 8.0, 4.0, 6.0, 8.0, 6.0, 1.0, 2.0, 2.0, 3.0, 4.0, 8.0, 1.0, 6.0, 2.0, 5.0, 8.0, 1.0, 8.0, 4.0, 5.0, 8.0, 4.0, 5.0, 7.0, 3.0, 4.0, 8.0, 7.0, 4.0, 7.0, 1.0, 4.0, 6.0, 9.0, 8.0, 9.0, 10.0, 1.0, 9.0, 9.0, 9.0, 1.0, 8.0, 7.0, 6.0, 6.0, 7.0, 2.0, 5.0, 9.0, 9.0, 8.0, 7.0, 3.0, 5.0, 8.0, 8.0, 6.0, 10.0, 2.0, 6.0, 5.0, 5.0, 1.0, 9.0, 7.0, 4.0, 2.0, 6.0, 9.0, 7.0, 3.0, 5.0, 9.0, 5.0, 1.0, 10.0, 8.0, 8.0, 6.0, 9.0, 1.0, 5.0, 4.0, 4.0, 7.0, 6.0, 5.0, 6.0, 6.0, 9.0, 5.0, 6.0, 6.0, 4.0, 8.0, 7.0, 5.0, 6.0, 1.0, 6.0, 9.0, 8.0, 7.0, 3.0, 8.0, 4.0, 8.0, 3.0, 3.0, 7.0, 2.0, 10.0, 5.0, 1.0, 7.0, 9.0, 8.0, 1.0, 5.0, 9.0, 7.0, 7.0, 4.0, 7.0, 2.0, 6.0, 2.0, 4.0, 5.0, 2.0, 9.0, 8.0, 10.0, 9.0, 7.0, 3.0, 6.0, 10.0, 2.0, 9.0, 9.0, 10.0, 9.0, 10.0, 9.0, 7.0, 7.0, 4.0, 6.0, 9.0, 6.0, 2.0, 7.0, 7.0, 9.0, 7.0, 10.0, 6.0, 3.0, 8.0, 10.0, 4.0, 3.0, 2.0, 6.0, 4.0, 5.0, 1.0, 6.0, 4.0, 8.0, 8.0, 1.0, 3.0, 4.0, 1.0, 4.0, 3.0, 4.0, 6.0, 8.0, 8.0, 7.0, 4.0, 8.0, 10.0, 4.0, 8.0, 10.0, 2.0, 7.0]
global b_x = 5
global d_y = [5.0, 10.0, 4.0, 1.0, 7.0, 8.0, 8.0, 1.0, 3.0, 5.0, 10.0, 3.0, 2.0, 8.0, 1.0, 8.0, 1.0, 9.0, 7.0, 10.0, 5.0, 6.0, 3.0, 1.0, 10.0, 4.0, 5.0, 3.0, 2.0, 1.0, 8.0, 2.0, 6.0, 5.0, 1.0, 7.0, 6.0, 6.0, 1.0, 10.0, 2.0, 8.0, 3.0, 2.0, 7.0, 10.0, 6.0, 6.0, 10.0, 5.0, 8.0, 3.0, 5.0, 6.0, 10.0, 8.0, 4.0, 9.0, 9.0, 9.0, 5.0, 9.0, 2.0, 2.0, 7.0, 9.0, 10.0, 10.0, 4.0, 9.0, 6.0, 10.0, 2.0, 1.0, 1.0, 8.0, 6.0, 3.0, 8.0, 6.0, 7.0, 6.0, 5.0, 2.0, 6.0, 8.0, 8.0, 6.0, 7.0, 2.0, 4.0, 8.0, 9.0, 10.0, 5.0, 6.0, 1.0, 1.0, 1.0, 2.0, 6.0, 8.0, 5.0, 9.0, 8.0, 1.0, 9.0, 10.0, 10.0, 10.0, 2.0, 3.0, 9.0, 3.0, 7.0, 7.0, 4.0, 10.0, 7.0, 9.0, 3.0, 1.0, 10.0, 10.0, 5.0, 5.0, 7.0, 9.0, 2.0, 5.0, 6.0, 1.0, 3.0, 4.0, 8.0, 2.0, 3.0, 5.0, 7.0, 10.0, 3.0, 6.0, 2.0, 5.0, 1.0, 9.0, 1.0, 2.0, 7.0, 9.0, 6.0, 8.0, 6.0, 10.0, 4.0, 2.0, 10.0, 4.0, 4.0, 6.0, 10.0, 5.0, 7.0, 9.0, 1.0, 7.0, 8.0, 1.0, 8.0, 10.0, 2.0, 7.0, 1.0, 3.0, 9.0, 10.0, 2.0, 7.0, 2.0, 3.0, 10.0, 3.0, 8.0, 4.0, 8.0, 3.0, 9.0, 6.0, 1.0, 7.0, 4.0, 7.0, 5.0, 7.0, 5.0, 9.0, 8.0, 8.0, 10.0, 2.0, 4.0, 5.0, 7.0, 3.0, 5.0, 8.0, 1.0, 6.0, 7.0, 3.0, 10.0, 2.0, 2.0, 2.0, 8.0, 2.0, 7.0, 5.0, 5.0, 3.0, 8.0, 3.0, 2.0, 6.0, 10.0, 2.0]
global b_y = 10
global p = [0.391, 0.405, 0.345, 0.182, 0.789, 0.458, 0.165, 0.029, 0.539, 0.54, 0.082, 0.976, 0.929, 0.39, 0.448, 0.3, 0.251, 0.328, 0.202, 0.495, 0.231, 0.144, 0.829, 0.432, 0.489, 0.422, 0.755, 0.092, 0.737, 0.095, 0.862, 0.595, 0.749, 0.327, 0.581, 0.952, 0.625, 0.112, 0.952, 0.133, 0.014, 0.298, 0.804, 0.035, 0.162, 0.689, 0.738, 0.78, 0.002, 0.149, 0.654, 0.501, 0.936, 0.75, 0.295, 0.779, 0.059, 0.604, 0.847, 0.025, 0.952, 0.794, 0.079, 0.15, 0.917, 0.357, 0.619, 0.556, 0.676, 0.43, 0.59, 0.804, 0.502, 0.378, 0.801, 0.325, 0.018, 0.314, 0.455, 0.086, 0.817, 0.435, 0.587, 0.952, 0.891, 0.475, 0.608, 0.464, 0.668, 0.062, 0.52, 0.906, 0.216, 0.403, 0.359, 0.462, 0.12, 0.872, 0.714, 0.318, 0.846, 0.668, 0.042, 0.554, 0.557, 0.135, 0.703, 0.237, 0.991, 0.911, 0.492, 0.298, 0.88, 0.691, 0.432, 0.203, 0.685, 0.441, 0.572, 0.764, 0.13, 0.424, 0.343, 0.698, 0.028, 0.34, 0.704, 0.968, 0.204, 0.398, 0.615, 0.807, 0.085, 0.586, 0.464, 0.931, 0.886, 0.408, 0.004, 0.499, 0.932, 0.904, 0.514, 0.988, 0.734, 0.121, 0.899, 0.533, 0.082, 0.811, 0.342, 0.676, 0.362, 0.564, 0.922, 0.245, 0.675, 0.027, 0.663, 0.805, 0.119, 0.18, 0.187, 0.854, 0.157, 0.515, 0.26, 0.008, 0.818, 0.504, 0.182, 0.35, 0.071, 0.847, 0.059, 0.114, 0.884, 0.452, 0.976, 0.785, 0.728, 0.107, 0.24, 0.973, 0.63, 0.741, 0.994, 0.788, 0.679, 0.483, 0.762, 0.924, 0.242, 0.313, 0.591, 0.887, 0.688, 0.541, 0.014, 0.903, 0.625, 0.428, 0.098, 0.519, 0.987, 0.562, 0.143, 0.037, 0.86, 0.294, 0.474, 0.356, 0.218, 0.811, 0.971, 0.707, 0.185, 0.054, 0.772, 0.851, 0.444, 0.8, 0.966, 0.14, 0.077, 0.372]
global q = [0.61, 0.817, 0.602, 0.577, 0.823, 0.883, 0.902, 0.21, 0.761, 0.545, 0.91, 0.98, 0.952, 0.8, 0.946, 0.551, 0.372, 0.414, 0.566, 0.86, 0.848, 0.982, 0.859, 0.72, 0.653, 0.7, 0.975, 0.923, 0.954, 0.346, 0.913, 0.885, 0.977, 0.815, 0.881, 0.993, 0.979, 0.309, 0.979, 0.769, 0.558, 0.335, 0.868, 0.965, 0.297, 0.832, 0.911, 0.941, 0.174, 0.689, 0.933, 0.96, 0.992, 0.819, 0.449, 0.888, 0.297, 0.755, 0.937, 0.966, 0.96, 0.876, 0.612, 0.923, 0.928, 0.504, 0.65, 0.678, 0.828, 0.517, 0.799, 0.956, 0.799, 0.6, 0.926, 0.628, 0.528, 0.75, 0.695, 0.299, 0.818, 0.866, 0.624, 0.997, 0.906, 0.913, 0.705, 0.947, 0.78, 0.084, 0.776, 0.971, 0.737, 0.749, 0.615, 0.853, 0.843, 0.971, 0.774, 0.513, 0.981, 0.801, 0.224, 0.593, 0.562, 0.231, 0.907, 0.635, 0.999, 0.974, 0.666, 0.627, 0.916, 0.88, 0.724, 0.479, 0.887, 0.687, 0.932, 0.83, 0.898, 0.89, 0.733, 0.802, 0.597, 0.421, 0.781, 0.985, 0.594, 0.738, 0.88, 0.867, 0.086, 0.955, 0.8, 0.975, 0.898, 0.428, 0.303, 0.549, 0.981, 0.912, 0.679, 0.993, 0.989, 0.468, 0.947, 0.995, 0.924, 0.926, 0.731, 0.743, 0.982, 0.698, 0.94, 0.509, 0.908, 0.519, 0.775, 0.903, 0.989, 0.976, 0.559, 0.984, 0.931, 0.76, 0.73, 0.86, 0.928, 0.513, 0.793, 0.824, 0.817, 0.876, 0.397, 0.927, 0.972, 0.591, 0.999, 0.871, 0.963, 0.915, 0.669, 0.99, 0.908, 0.882, 0.994, 0.987, 0.734, 0.873, 0.915, 0.93, 0.413, 0.94, 0.832, 0.966, 0.709, 0.862, 0.122, 0.958, 0.738, 0.768, 0.489, 0.908, 0.997, 0.872, 0.342, 0.463, 0.927, 0.565, 0.514, 0.799, 0.835, 0.998, 0.983, 0.941, 0.385, 0.439, 0.788, 0.885, 0.941, 0.906, 0.976, 0.935, 0.722, 0.993]
global origin = 1
global destination = 50