global arcs = [1 13; 1 26; 1 30; 1 40; 1 41; 1 50; 2 12; 2 18; 2 25; 2 37; 2 38; 2 39; 3 31; 3 33; 3 42; 3 46; 4 11; 4 18; 4 32; 4 36; 4 47; 4 48; 5 10; 5 25; 5 32; 5 35; 5 41; 5 49; 5 50; 6 27; 6 41; 6 47; 6 50; 7 3; 7 10; 7 16; 7 18; 7 19; 7 20; 7 21; 7 25; 7 29; 7 49; 7 50; 8 3; 8 34; 8 42; 9 4; 9 17; 9 22; 9 24; 9 32; 9 39; 10 22; 10 24; 10 38; 11 7; 11 15; 11 27; 11 30; 11 35; 11 39; 11 45; 12 8; 12 9; 12 11; 12 15; 12 32; 12 38; 12 39; 13 10; 13 14; 14 17; 14 26; 14 35; 14 44; 14 48; 15 21; 15 24; 15 25; 15 30; 15 41; 15 48; 16 4; 16 7; 16 27; 16 35; 16 40; 17 11; 17 25; 18 5; 18 10; 18 11; 18 14; 18 15; 18 21; 18 29; 18 40; 19 4; 19 5; 19 11; 19 14; 19 22; 19 26; 19 27; 19 38; 19 42; 20 15; 20 29; 20 37; 21 5; 21 9; 21 12; 22 2; 22 11; 22 14; 22 17; 22 34; 23 24; 23 41; 23 43; 24 5; 24 13; 24 16; 24 19; 24 33; 24 38; 24 39; 25 3; 25 6; 25 12; 25 23; 25 28; 26 11; 26 12; 26 18; 26 21; 26 38; 26 45; 26 49; 26 50; 27 9; 27 14; 27 22; 27 30; 27 44; 28 4; 28 16; 28 20; 28 22; 28 32; 29 4; 29 9; 29 11; 29 15; 29 40; 30 6; 30 22; 31 27; 31 35; 32 12; 32 20; 32 27; 32 37; 32 41; 32 44; 32 49; 33 21; 33 24; 33 25; 33 46; 33 47; 34 11; 34 17; 34 18; 34 24; 34 32; 34 35; 34 36; 34 48; 35 2; 35 3; 35 4; 35 9; 35 37; 35 43; 35 44; 36 13; 36 15; 36 30; 36 40; 36 42; 37 5; 37 9; 37 32; 37 46; 38 3; 38 32; 38 34; 39 7; 39 16; 39 19; 39 37; 40 23; 40 27; 40 29; 40 32; 40 43; 41 11; 41 23; 41 29; 41 30; 41 45; 42 2; 42 41; 43 11; 43 13; 43 37; 44 4; 44 22; 44 37; 44 39; 44 50; 45 5; 45 9; 45 15; 45 29; 45 40; 45 47; 46 6; 46 22; 46 31; 46 39; 47 35; 48 2; 48 6; 48 31; 48 35; 48 44; 49 13; 49 30; 49 35; 49 44]
global d_x = [1.0, 3.0, 5.0, 5.0, 2.0, 6.0, 3.0, 4.0, 9.0, 1.0, 4.0, 4.0, 5.0, 3.0, 4.0, 5.0, 7.0, 5.0, 7.0, 5.0, 6.0, 10.0, 7.0, 6.0, 8.0, 10.0, 7.0, 7.0, 2.0, 4.0, 2.0, 5.0, 10.0, 8.0, 9.0, 2.0, 4.0, 1.0, 8.0, 10.0, 5.0, 1.0, 5.0, 8.0, 4.0, 7.0, 5.0, 7.0, 4.0, 9.0, 9.0, 1.0, 9.0, 10.0, 9.0, 2.0, 8.0, 7.0, 9.0, 3.0, 3.0, 4.0, 6.0, 3.0, 8.0, 7.0, 4.0, 6.0, 9.0, 3.0, 2.0, 5.0, 2.0, 3.0, 7.0, 2.0, 6.0, 9.0, 8.0, 5.0, 1.0, 3.0, 8.0, 2.0, 6.0, 3.0, 10.0, 6.0, 6.0, 9.0, 9.0, 4.0, 4.0, 10.0, 5.0, 10.0, 7.0, 1.0, 4.0, 5.0, 9.0, 6.0, 3.0, 9.0, 2.0, 7.0, 5.0, 10.0, 3.0, 2.0, 5.0, 7.0, 6.0, 8.0, 5.0, 3.0, 7.0, 5.0, 7.0, 5.0, 9.0, 5.0, 6.0, 10.0, 2.0, 8.0, 2.0, 3.0, 7.0, 3.0, 8.0, 7.0, 10.0, 10.0, 4.0, 2.0, 9.0, 7.0, 5.0, 3.0, 3.0, 6.0, 2.0, 9.0, 8.0, 6.0, 5.0, 7.0, 5.0, 6.0, 9.0, 3.0, 9.0, 7.0, 9.0, 3.0, 10.0, 5.0, 3.0, 4.0, 3.0, 8.0, 6.0, 8.0, 5.0, 4.0, 6.0, 10.0, 2.0, 10.0, 8.0, 5.0, 6.0, 4.0, 2.0, 5.0, 1.0, 5.0, 9.0, 6.0, 2.0, 8.0, 1.0, 1.0, 1.0, 10.0, 4.0, 6.0, 2.0, 1.0, 6.0, 7.0, 5.0, 10.0, 3.0, 3.0, 7.0, 6.0, 1.0, 9.0, 5.0, 8.0, 2.0, 2.0, 4.0, 3.0, 6.0, 3.0, 7.0, 2.0, 2.0, 8.0, 1.0, 5.0, 8.0, 10.0, 3.0, 3.0, 6.0, 6.0, 7.0, 2.0, 2.0, 10.0, 9.0, 4.0, 5.0, 5.0, 9.0, 7.0, 9.0, 6.0, 7.0, 4.0, 8.0, 7.0, 6.0, 7.0, 7.0, 9.0, 7.0, 5.0, 10.0]
global b_x = 5
global d_y = [5.0, 8.0, 2.0, 5.0, 3.0, 7.0, 10.0, 2.0, 9.0, 9.0, 4.0, 7.0, 10.0, 4.0, 5.0, 4.0, 1.0, 6.0, 5.0, 2.0, 10.0, 9.0, 1.0, 4.0, 7.0, 7.0, 1.0, 2.0, 10.0, 5.0, 3.0, 5.0, 1.0, 7.0, 8.0, 1.0, 3.0, 6.0, 2.0, 10.0, 9.0, 7.0, 2.0, 9.0, 9.0, 3.0, 2.0, 8.0, 1.0, 10.0, 8.0, 9.0, 2.0, 5.0, 3.0, 1.0, 3.0, 6.0, 4.0, 7.0, 5.0, 5.0, 5.0, 8.0, 2.0, 4.0, 8.0, 6.0, 3.0, 7.0, 4.0, 4.0, 9.0, 5.0, 8.0, 5.0, 6.0, 4.0, 10.0, 4.0, 3.0, 3.0, 9.0, 2.0, 4.0, 6.0, 8.0, 7.0, 9.0, 1.0, 2.0, 4.0, 9.0, 7.0, 3.0, 9.0, 8.0, 8.0, 8.0, 8.0, 10.0, 9.0, 5.0, 2.0, 10.0, 9.0, 8.0, 5.0, 5.0, 10.0, 10.0, 4.0, 4.0, 5.0, 3.0, 10.0, 3.0, 5.0, 3.0, 2.0, 7.0, 7.0, 10.0, 3.0, 10.0, 10.0, 3.0, 3.0, 7.0, 3.0, 1.0, 10.0, 9.0, 10.0, 5.0, 4.0, 2.0, 10.0, 2.0, 8.0, 8.0, 3.0, 7.0, 1.0, 8.0, 6.0, 2.0, 1.0, 3.0, 9.0, 4.0, 8.0, 1.0, 6.0, 2.0, 2.0, 3.0, 8.0, 1.0, 3.0, 8.0, 10.0, 8.0, 10.0, 10.0, 3.0, 5.0, 10.0, 7.0, 9.0, 9.0, 2.0, 5.0, 7.0, 5.0, 1.0, 1.0, 8.0, 10.0, 3.0, 10.0, 9.0, 2.0, 6.0, 10.0, 7.0, 9.0, 3.0, 5.0, 7.0, 8.0, 10.0, 6.0, 10.0, 9.0, 8.0, 6.0, 7.0, 6.0, 4.0, 3.0, 8.0, 7.0, 9.0, 8.0, 1.0, 6.0, 3.0, 5.0, 9.0, 2.0, 5.0, 1.0, 10.0, 3.0, 6.0, 5.0, 5.0, 3.0, 2.0, 3.0, 6.0, 7.0, 8.0, 6.0, 10.0, 5.0, 10.0, 7.0, 10.0, 3.0, 2.0, 6.0, 10.0, 5.0, 4.0, 6.0, 6.0, 3.0, 3.0, 5.0, 9.0, 6.0]
global b_y = 10
global p = [0.032, 0.193, 0.414, 0.597, 0.683, 0.849, 0.333, 0.758, 0.08, 0.518, 0.108, 0.398, 0.34, 0.507, 0.044, 0.508, 0.475, 0.866, 0.258, 0.046, 0.81, 0.426, 0.671, 0.925, 0.806, 0.439, 0.899, 0.875, 0.737, 0.845, 0.133, 0.936, 0.245, 0.883, 0.861, 0.44, 0.079, 0.365, 0.323, 0.66, 0.865, 0.518, 0.011, 0.501, 0.67, 0.208, 0.693, 0.69, 0.822, 0.12, 0.141, 0.134, 0.549, 0.138, 0.731, 0.371, 0.218, 0.981, 0.765, 0.061, 0.013, 0.887, 0.946, 0.83, 0.608, 0.384, 0.814, 0.947, 0.772, 0.143, 0.436, 0.305, 0.678, 0.481, 0.796, 0.148, 0.2, 0.35, 0.256, 0.071, 0.918, 0.362, 0.302, 0.7, 0.123, 0.912, 0.01, 0.655, 0.873, 0.808, 0.381, 0.561, 0.521, 0.457, 0.606, 0.553, 0.577, 0.826, 0.603, 0.737, 0.966, 0.354, 0.199, 0.574, 0.302, 0.004, 0.238, 0.138, 0.9, 0.777, 0.019, 0.215, 0.164, 0.117, 0.306, 0.504, 0.15, 0.11, 0.498, 0.905, 0.349, 0.201, 0.213, 0.435, 0.798, 0.198, 0.987, 0.931, 0.181, 0.88, 0.266, 0.685, 0.593, 0.708, 0.637, 0.66, 0.973, 0.761, 0.758, 0.733, 0.895, 0.872, 0.831, 0.353, 0.713, 0.8, 0.497, 0.478, 0.81, 0.054, 0.281, 0.807, 0.654, 0.366, 0.838, 0.679, 0.37, 0.576, 0.446, 0.587, 0.479, 0.294, 0.221, 0.311, 0.204, 0.041, 0.121, 0.871, 0.187, 0.743, 0.731, 0.371, 0.064, 0.802, 0.248, 0.943, 0.793, 0.429, 0.428, 0.054, 0.136, 0.886, 0.403, 0.193, 0.112, 0.813, 0.48, 0.12, 0.012, 0.407, 0.406, 0.836, 0.936, 0.969, 0.193, 0.039, 0.722, 0.105, 0.874, 0.888, 0.719, 0.245, 0.959, 0.482, 0.109, 0.057, 0.16, 0.778, 0.322, 0.925, 0.915, 0.464, 0.921, 0.608, 0.225, 0.641, 0.025, 0.838, 0.409, 0.888, 0.658, 0.058, 0.053, 0.155, 0.421, 0.027, 0.175, 0.943, 0.761, 0.946, 0.491, 0.795, 0.253, 0.901, 0.984, 0.491, 0.349, 0.39, 0.357, 0.036, 0.201, 0.366, 0.23]
global q = [0.141, 0.921, 0.72, 0.689, 0.725, 0.856, 0.539, 0.931, 0.093, 0.874, 0.716, 0.919, 0.466, 0.926, 0.916, 0.6, 0.912, 0.984, 0.486, 0.724, 0.904, 0.535, 0.916, 0.946, 0.903, 0.519, 0.971, 0.991, 0.916, 0.845, 0.787, 0.997, 0.265, 0.918, 0.882, 0.691, 0.663, 0.653, 0.414, 0.825, 0.895, 0.645, 0.167, 0.903, 0.746, 0.771, 0.953, 0.802, 0.923, 0.458, 0.326, 0.958, 0.707, 0.441, 0.82, 0.608, 0.989, 0.986, 0.976, 0.45, 0.405, 0.961, 0.946, 0.836, 0.785, 0.51, 0.938, 0.952, 0.88, 0.831, 0.867, 0.519, 0.983, 0.658, 0.822, 0.68, 0.41, 0.489, 0.764, 0.963, 0.965, 0.527, 0.639, 0.958, 0.907, 0.996, 0.142, 0.969, 0.966, 0.979, 0.785, 0.658, 0.904, 0.923, 0.914, 0.636, 0.883, 0.972, 0.786, 0.91, 0.981, 0.821, 0.79, 0.776, 0.439, 0.382, 0.313, 0.819, 0.908, 0.944, 0.547, 0.991, 0.444, 0.935, 0.739, 0.986, 0.383, 0.481, 0.79, 0.928, 0.756, 0.927, 0.448, 0.999, 0.89, 0.705, 0.997, 0.941, 0.596, 0.992, 0.615, 0.7, 0.833, 0.978, 0.844, 0.997, 0.995, 0.854, 0.974, 0.824, 0.961, 0.989, 0.863, 0.49, 0.722, 0.812, 0.788, 0.867, 0.854, 0.238, 0.862, 0.854, 0.732, 0.569, 0.99, 0.72, 0.985, 0.873, 0.887, 0.669, 0.801, 0.656, 0.647, 0.413, 0.722, 0.615, 0.549, 0.998, 0.416, 0.874, 0.827, 0.593, 0.297, 0.808, 0.857, 0.943, 0.927, 0.857, 0.64, 0.543, 0.391, 0.923, 0.887, 0.214, 0.9, 0.914, 0.722, 0.84, 0.538, 0.819, 0.779, 0.931, 0.949, 0.994, 0.927, 0.708, 0.968, 0.38, 0.98, 0.895, 0.975, 0.56, 0.988, 0.735, 0.797, 0.292, 0.219, 0.864, 0.573, 0.999, 0.942, 0.498, 0.983, 0.904, 0.437, 0.953, 0.296, 0.965, 0.781, 0.989, 0.831, 0.076, 0.68, 0.711, 0.57, 0.667, 0.841, 0.976, 0.839, 0.983, 0.594, 0.912, 0.867, 0.966, 0.989, 0.615, 0.836, 0.759, 0.61, 0.363, 0.356, 0.729, 0.829]
global origin = 1
global destination = 50