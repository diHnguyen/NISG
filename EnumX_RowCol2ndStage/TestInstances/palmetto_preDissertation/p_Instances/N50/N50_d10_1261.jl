global arcs = [1 10; 1 16; 2 28; 2 30; 2 41; 3 25; 3 39; 3 42; 4 27; 4 45; 4 47; 4 50; 5 4; 5 13; 5 17; 5 39; 6 11; 6 13; 6 24; 6 40; 6 43; 7 23; 7 27; 7 34; 7 47; 8 3; 8 5; 8 18; 8 19; 8 23; 8 27; 8 35; 8 42; 8 49; 9 25; 9 39; 9 44; 9 46; 9 47; 10 7; 10 14; 11 3; 11 8; 11 10; 11 22; 11 24; 11 29; 11 39; 11 46; 12 42; 12 44; 13 3; 13 14; 13 16; 13 25; 13 28; 13 35; 13 38; 13 39; 13 41; 13 49; 14 10; 14 24; 14 27; 14 47; 15 9; 15 19; 15 24; 15 25; 15 37; 15 46; 16 17; 16 18; 16 22; 16 29; 16 37; 17 3; 17 13; 17 18; 17 23; 17 30; 17 32; 17 33; 17 45; 18 3; 18 44; 18 47; 19 3; 19 4; 19 5; 19 20; 19 35; 20 11; 20 12; 20 13; 20 22; 20 33; 20 48; 20 49; 20 50; 21 2; 21 14; 21 28; 21 36; 22 20; 22 35; 22 41; 22 42; 22 45; 23 14; 23 29; 23 32; 23 45; 23 46; 24 6; 24 30; 24 33; 24 44; 24 50; 25 22; 25 29; 25 43; 26 8; 26 12; 26 14; 26 28; 26 29; 26 31; 26 43; 26 49; 27 4; 27 5; 27 7; 27 10; 27 16; 27 18; 28 2; 28 3; 28 8; 28 21; 28 22; 28 30; 28 31; 28 36; 28 44; 28 46; 29 4; 29 13; 29 26; 29 40; 29 46; 30 6; 30 8; 30 26; 30 27; 30 43; 30 45; 30 47; 31 7; 31 15; 31 27; 31 39; 32 4; 32 22; 32 29; 32 43; 32 44; 33 13; 33 23; 33 32; 33 40; 34 3; 34 29; 34 35; 34 42; 34 43; 34 45; 35 27; 35 28; 35 45; 36 4; 36 5; 36 8; 36 25; 36 43; 37 3; 37 10; 37 20; 37 25; 37 31; 38 2; 38 22; 38 24; 38 34; 38 37; 38 40; 39 3; 39 49; 40 12; 40 19; 40 45; 41 3; 41 9; 41 29; 41 47; 42 26; 42 36; 43 3; 43 21; 43 23; 43 25; 43 28; 44 9; 44 18; 44 23; 45 2; 45 33; 45 36; 46 10; 46 22; 46 32; 46 34; 47 2; 47 4; 47 6; 47 11; 47 24; 47 31; 47 32; 47 44; 48 14; 48 43; 49 4; 49 7; 49 9; 49 22; 49 28; 49 34; 49 44; 49 48]
global d_x = [5.0, 3.0, 6.0, 6.0, 9.0, 4.0, 7.0, 6.0, 10.0, 3.0, 10.0, 10.0, 5.0, 1.0, 10.0, 9.0, 5.0, 6.0, 6.0, 5.0, 9.0, 6.0, 5.0, 10.0, 3.0, 7.0, 8.0, 7.0, 5.0, 6.0, 7.0, 7.0, 3.0, 9.0, 1.0, 10.0, 5.0, 6.0, 6.0, 4.0, 5.0, 10.0, 2.0, 8.0, 1.0, 1.0, 6.0, 3.0, 3.0, 8.0, 1.0, 9.0, 9.0, 7.0, 2.0, 8.0, 9.0, 5.0, 8.0, 7.0, 2.0, 4.0, 10.0, 9.0, 10.0, 1.0, 10.0, 10.0, 2.0, 7.0, 4.0, 6.0, 8.0, 4.0, 4.0, 7.0, 7.0, 9.0, 1.0, 1.0, 6.0, 10.0, 10.0, 7.0, 9.0, 8.0, 5.0, 7.0, 10.0, 7.0, 1.0, 8.0, 3.0, 1.0, 6.0, 8.0, 9.0, 8.0, 6.0, 7.0, 1.0, 9.0, 9.0, 10.0, 3.0, 2.0, 3.0, 6.0, 3.0, 7.0, 10.0, 3.0, 3.0, 5.0, 6.0, 9.0, 4.0, 5.0, 3.0, 10.0, 2.0, 1.0, 9.0, 8.0, 3.0, 6.0, 10.0, 8.0, 4.0, 9.0, 1.0, 1.0, 7.0, 4.0, 10.0, 1.0, 6.0, 1.0, 9.0, 1.0, 8.0, 8.0, 5.0, 10.0, 5.0, 7.0, 9.0, 3.0, 9.0, 6.0, 7.0, 9.0, 8.0, 2.0, 4.0, 8.0, 1.0, 7.0, 1.0, 7.0, 9.0, 6.0, 10.0, 9.0, 4.0, 3.0, 9.0, 1.0, 6.0, 2.0, 8.0, 9.0, 7.0, 8.0, 3.0, 10.0, 1.0, 4.0, 1.0, 5.0, 5.0, 6.0, 7.0, 5.0, 6.0, 1.0, 6.0, 4.0, 5.0, 1.0, 7.0, 4.0, 1.0, 2.0, 2.0, 6.0, 6.0, 6.0, 3.0, 9.0, 6.0, 6.0, 3.0, 4.0, 9.0, 4.0, 6.0, 3.0, 4.0, 1.0, 5.0, 3.0, 2.0, 10.0, 9.0, 3.0, 2.0, 8.0, 4.0, 7.0, 8.0, 3.0, 6.0, 6.0, 10.0, 1.0, 5.0, 2.0, 3.0, 5.0, 1.0, 6.0, 9.0, 8.0, 3.0, 9.0, 6.0, 8.0, 4.0, 3.0]
global b_x = 5
global d_y = [1.0, 4.0, 3.0, 9.0, 4.0, 4.0, 8.0, 10.0, 7.0, 3.0, 3.0, 1.0, 9.0, 4.0, 1.0, 6.0, 9.0, 7.0, 5.0, 1.0, 8.0, 10.0, 6.0, 5.0, 3.0, 4.0, 8.0, 1.0, 9.0, 5.0, 3.0, 10.0, 7.0, 3.0, 9.0, 4.0, 8.0, 10.0, 1.0, 1.0, 6.0, 7.0, 1.0, 2.0, 9.0, 8.0, 4.0, 9.0, 10.0, 9.0, 5.0, 7.0, 7.0, 3.0, 10.0, 5.0, 9.0, 2.0, 8.0, 10.0, 6.0, 6.0, 1.0, 7.0, 10.0, 8.0, 8.0, 1.0, 6.0, 2.0, 6.0, 5.0, 9.0, 9.0, 1.0, 9.0, 9.0, 10.0, 1.0, 2.0, 3.0, 2.0, 1.0, 10.0, 3.0, 10.0, 1.0, 3.0, 9.0, 3.0, 2.0, 2.0, 4.0, 9.0, 6.0, 9.0, 10.0, 6.0, 3.0, 4.0, 5.0, 4.0, 4.0, 3.0, 8.0, 6.0, 3.0, 8.0, 7.0, 3.0, 10.0, 10.0, 8.0, 1.0, 10.0, 4.0, 8.0, 6.0, 6.0, 5.0, 1.0, 2.0, 6.0, 4.0, 6.0, 7.0, 7.0, 7.0, 6.0, 5.0, 9.0, 1.0, 8.0, 3.0, 3.0, 8.0, 8.0, 2.0, 7.0, 4.0, 10.0, 4.0, 4.0, 6.0, 10.0, 5.0, 2.0, 1.0, 1.0, 5.0, 2.0, 9.0, 1.0, 7.0, 8.0, 9.0, 8.0, 8.0, 6.0, 7.0, 10.0, 8.0, 3.0, 9.0, 10.0, 10.0, 8.0, 3.0, 7.0, 7.0, 5.0, 7.0, 10.0, 1.0, 8.0, 3.0, 8.0, 5.0, 5.0, 8.0, 8.0, 4.0, 3.0, 7.0, 4.0, 2.0, 2.0, 8.0, 8.0, 2.0, 1.0, 6.0, 5.0, 2.0, 10.0, 6.0, 10.0, 3.0, 3.0, 9.0, 8.0, 9.0, 10.0, 1.0, 8.0, 5.0, 6.0, 9.0, 9.0, 10.0, 5.0, 4.0, 3.0, 3.0, 2.0, 2.0, 1.0, 3.0, 6.0, 3.0, 4.0, 3.0, 6.0, 4.0, 9.0, 6.0, 9.0, 4.0, 7.0, 6.0, 5.0, 5.0, 5.0, 1.0, 9.0, 7.0, 6.0, 4.0, 3.0, 2.0]
global b_y = 10
global p = [0.085, 0.023, 0.58, 0.815, 0.997, 0.849, 0.044, 0.578, 0.908, 0.621, 0.573, 0.778, 0.746, 0.725, 0.437, 0.745, 0.751, 0.594, 0.023, 0.172, 0.972, 0.206, 0.348, 0.142, 0.982, 0.088, 0.233, 0.213, 0.661, 0.758, 0.336, 0.872, 0.611, 0.746, 0.632, 0.529, 0.62, 0.993, 0.539, 0.788, 0.892, 0.735, 0.041, 0.701, 0.267, 0.532, 0.446, 0.706, 0.199, 0.207, 0.078, 0.388, 0.213, 0.803, 0.587, 0.672, 0.431, 0.197, 0.574, 0.947, 0.064, 0.126, 0.984, 0.138, 0.159, 0.85, 0.776, 0.474, 0.341, 0.616, 0.435, 0.403, 0.773, 0.073, 0.473, 0.446, 0.242, 0.906, 0.17, 0.638, 0.754, 0.251, 0.153, 0.077, 0.311, 0.444, 0.164, 0.684, 0.711, 0.907, 0.267, 0.89, 0.293, 0.12, 0.693, 0.825, 0.738, 0.695, 0.521, 0.047, 0.786, 0.979, 0.161, 0.243, 0.483, 0.022, 0.922, 0.163, 0.206, 0.42, 0.128, 0.862, 0.256, 0.86, 0.674, 0.937, 0.236, 0.84, 0.579, 0.173, 0.018, 0.599, 0.03, 0.487, 0.701, 0.787, 0.491, 0.681, 0.953, 0.203, 0.032, 0.906, 0.829, 0.785, 0.88, 0.005, 0.982, 0.06, 0.355, 0.955, 0.342, 0.66, 0.629, 0.292, 0.736, 0.181, 0.662, 0.113, 0.169, 0.102, 0.216, 0.83, 0.373, 0.451, 0.164, 0.676, 0.869, 0.26, 0.685, 0.772, 0.183, 0.605, 0.459, 0.757, 0.393, 0.729, 0.3, 0.018, 0.271, 0.187, 0.832, 0.16, 0.431, 0.404, 0.831, 0.948, 0.572, 0.468, 0.351, 0.641, 0.754, 0.015, 0.915, 0.365, 0.309, 0.378, 0.732, 0.529, 0.876, 0.041, 0.28, 0.732, 0.447, 0.673, 0.658, 0.254, 0.166, 0.261, 0.579, 0.174, 0.502, 0.015, 0.076, 0.353, 0.35, 0.059, 0.943, 0.327, 0.045, 0.011, 0.974, 0.973, 0.825, 0.485, 0.981, 0.186, 0.846, 0.422, 0.253, 0.056, 0.463, 0.422, 0.743, 0.724, 0.756, 0.487, 0.224, 0.51, 0.533, 0.8, 0.401, 0.134, 0.027, 0.15, 0.402, 0.383, 0.087, 0.175, 0.043, 0.495]
global q = [0.911, 0.127, 0.86, 0.823, 0.997, 0.897, 0.702, 0.725, 0.956, 0.787, 0.879, 0.822, 0.946, 0.847, 0.479, 0.912, 0.822, 0.798, 0.524, 0.37, 0.981, 0.659, 0.715, 0.424, 0.983, 0.987, 0.508, 0.715, 0.68, 0.981, 0.499, 0.957, 0.992, 0.959, 0.855, 0.753, 0.95, 0.997, 0.996, 0.8, 0.993, 0.922, 0.972, 0.967, 0.372, 0.947, 0.915, 0.825, 0.596, 0.227, 0.96, 0.991, 0.314, 0.93, 0.83, 0.716, 0.601, 0.773, 0.596, 0.972, 0.522, 0.643, 0.987, 0.517, 0.476, 0.878, 0.884, 0.995, 0.728, 0.621, 0.52, 0.925, 0.832, 0.281, 0.954, 0.826, 0.844, 0.952, 0.212, 0.738, 0.788, 0.769, 0.569, 0.211, 0.795, 0.629, 0.527, 0.857, 0.957, 0.974, 0.923, 0.891, 0.953, 0.89, 0.982, 0.933, 0.751, 0.753, 0.918, 0.646, 0.859, 0.98, 0.739, 0.977, 0.497, 0.642, 0.941, 0.545, 0.478, 0.581, 0.671, 0.884, 0.632, 0.864, 0.871, 0.958, 0.286, 0.855, 0.75, 0.494, 0.969, 0.825, 0.412, 0.862, 0.977, 0.788, 0.961, 0.961, 0.973, 0.981, 0.367, 0.965, 0.856, 0.793, 0.896, 0.335, 0.994, 0.234, 0.696, 0.977, 0.496, 0.993, 0.986, 0.903, 0.841, 0.382, 0.828, 0.897, 0.272, 0.927, 0.966, 0.95, 0.672, 0.607, 0.599, 0.758, 0.971, 0.698, 0.904, 0.998, 0.263, 0.688, 0.542, 0.901, 0.583, 0.968, 0.657, 0.286, 0.78, 0.243, 0.987, 0.708, 0.748, 0.67, 0.865, 0.956, 0.81, 0.659, 0.412, 0.836, 0.843, 0.089, 0.941, 0.927, 0.829, 0.972, 0.924, 0.85, 0.959, 0.9, 0.386, 0.781, 0.765, 0.85, 0.796, 0.278, 0.398, 0.492, 0.771, 0.45, 0.643, 0.091, 0.748, 0.855, 0.829, 0.961, 0.96, 0.996, 0.733, 0.211, 0.989, 0.989, 0.884, 0.909, 0.989, 0.893, 0.913, 0.597, 0.272, 0.157, 0.659, 0.642, 0.941, 0.732, 0.796, 0.51, 0.948, 0.859, 0.836, 0.979, 0.486, 0.541, 0.061, 0.94, 0.682, 0.592, 0.394, 0.236, 0.508, 0.519]
global origin = 1
global destination = 50