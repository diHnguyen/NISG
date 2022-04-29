global arcs = [1 7; 1 16; 1 28; 2 12; 2 13; 2 16; 2 31; 2 32; 3 11; 4 5; 4 14; 4 21; 4 23; 5 2; 5 13; 5 18; 5 23; 5 31; 6 11; 6 20; 6 30; 6 32; 7 8; 7 11; 7 14; 7 21; 7 28; 8 4; 8 13; 8 26; 9 12; 10 4; 10 7; 10 8; 10 9; 10 15; 10 25; 10 28; 10 30; 10 31; 10 34; 11 10; 11 13; 11 35; 12 20; 12 23; 12 24; 12 35; 13 11; 13 12; 13 14; 13 35; 14 4; 14 8; 14 13; 14 26; 14 27; 14 30; 15 3; 15 4; 15 9; 15 17; 16 24; 17 4; 17 24; 17 27; 18 10; 18 13; 18 25; 18 28; 19 7; 19 13; 19 22; 19 28; 20 15; 20 17; 21 6; 21 14; 21 24; 22 9; 22 18; 22 25; 22 29; 23 2; 23 4; 23 17; 23 20; 24 5; 24 7; 24 23; 25 10; 25 23; 26 13; 26 18; 26 23; 26 29; 26 30; 26 35; 27 9; 27 12; 27 14; 27 18; 27 20; 28 12; 28 26; 28 31; 28 32; 29 5; 29 16; 29 35; 30 5; 30 32; 31 9; 31 10; 31 26; 31 33; 32 10; 32 22; 33 19; 33 23; 34 21]
global d_x = [3.0, 6.0, 7.0, 9.0, 8.0, 7.0, 4.0, 9.0, 4.0, 1.0, 9.0, 9.0, 9.0, 7.0, 9.0, 10.0, 2.0, 7.0, 6.0, 2.0, 9.0, 1.0, 5.0, 5.0, 5.0, 5.0, 2.0, 2.0, 2.0, 6.0, 1.0, 6.0, 2.0, 4.0, 10.0, 5.0, 6.0, 4.0, 3.0, 4.0, 9.0, 8.0, 6.0, 2.0, 2.0, 1.0, 5.0, 3.0, 1.0, 8.0, 6.0, 1.0, 4.0, 10.0, 9.0, 9.0, 8.0, 6.0, 2.0, 2.0, 2.0, 9.0, 8.0, 5.0, 3.0, 6.0, 5.0, 4.0, 10.0, 2.0, 4.0, 5.0, 9.0, 1.0, 9.0, 6.0, 6.0, 10.0, 2.0, 8.0, 6.0, 5.0, 4.0, 8.0, 4.0, 8.0, 3.0, 10.0, 9.0, 1.0, 4.0, 8.0, 8.0, 10.0, 6.0, 1.0, 4.0, 2.0, 3.0, 2.0, 8.0, 9.0, 2.0, 2.0, 4.0, 5.0, 10.0, 10.0, 9.0, 10.0, 8.0, 2.0, 5.0, 5.0, 3.0, 1.0, 6.0, 7.0, 2.0, 2.0, 1.0]
global b_x = 5
global d_y = [4.0, 3.0, 2.0, 9.0, 7.0, 8.0, 2.0, 1.0, 10.0, 8.0, 9.0, 7.0, 8.0, 5.0, 3.0, 2.0, 4.0, 10.0, 10.0, 4.0, 4.0, 7.0, 5.0, 5.0, 2.0, 1.0, 8.0, 7.0, 8.0, 5.0, 1.0, 4.0, 9.0, 4.0, 2.0, 3.0, 2.0, 8.0, 2.0, 6.0, 8.0, 8.0, 9.0, 3.0, 6.0, 6.0, 2.0, 9.0, 6.0, 5.0, 2.0, 6.0, 6.0, 6.0, 3.0, 2.0, 2.0, 1.0, 6.0, 3.0, 2.0, 10.0, 6.0, 1.0, 4.0, 7.0, 8.0, 3.0, 2.0, 8.0, 6.0, 4.0, 6.0, 5.0, 10.0, 3.0, 6.0, 7.0, 2.0, 2.0, 6.0, 2.0, 9.0, 6.0, 8.0, 4.0, 6.0, 2.0, 2.0, 3.0, 8.0, 8.0, 9.0, 1.0, 6.0, 4.0, 10.0, 9.0, 1.0, 10.0, 2.0, 8.0, 2.0, 9.0, 10.0, 9.0, 10.0, 1.0, 1.0, 4.0, 4.0, 7.0, 4.0, 5.0, 4.0, 2.0, 5.0, 6.0, 10.0, 7.0, 8.0]
global b_y = 10
global p = [0.169, 0.081, 0.862, 0.151, 0.802, 0.491, 0.169, 0.576, 0.627, 0.015, 0.672, 0.206, 0.948, 0.474, 0.18, 0.652, 0.538, 0.891, 0.869, 0.032, 0.257, 0.885, 0.114, 0.473, 0.785, 0.212, 0.943, 0.48, 0.122, 0.729, 0.979, 0.053, 0.406, 0.608, 0.733, 0.577, 0.674, 0.572, 0.537, 0.41, 0.249, 0.874, 0.378, 0.079, 0.957, 0.96, 0.102, 0.717, 0.462, 0.651, 0.787, 0.812, 0.926, 0.527, 0.013, 0.295, 0.197, 0.057, 0.443, 0.431, 0.667, 0.503, 0.019, 0.298, 0.737, 0.962, 0.012, 0.027, 0.879, 0.42, 0.764, 0.447, 0.286, 0.732, 0.487, 0.403, 0.825, 0.021, 0.162, 0.784, 0.775, 0.154, 0.131, 0.683, 0.435, 0.937, 0.26, 0.731, 0.266, 0.804, 0.287, 0.214, 0.679, 0.843, 0.574, 0.726, 0.67, 0.5, 0.648, 0.936, 0.232, 0.989, 0.666, 0.799, 0.592, 0.819, 0.049, 0.388, 0.008, 0.199, 0.739, 0.46, 0.915, 0.244, 0.715, 0.727, 0.057, 0.887, 0.283, 0.33, 0.615]
global q = [0.732, 0.884, 0.939, 0.356, 0.812, 0.749, 0.219, 0.808, 0.968, 0.726, 0.832, 0.577, 0.954, 0.883, 0.214, 0.816, 0.796, 0.894, 0.983, 0.445, 0.822, 0.998, 0.795, 0.974, 0.834, 0.435, 0.966, 0.898, 0.178, 0.958, 0.986, 0.082, 0.57, 0.973, 0.828, 0.955, 0.868, 0.966, 0.913, 0.957, 0.728, 0.987, 0.408, 0.436, 0.988, 0.968, 0.63, 0.994, 0.816, 0.936, 0.971, 0.979, 0.972, 0.568, 0.714, 0.852, 0.368, 0.576, 0.562, 0.64, 0.827, 0.593, 0.041, 0.452, 0.751, 0.995, 0.23, 0.273, 0.971, 0.796, 0.807, 0.89, 0.309, 0.872, 0.696, 0.655, 0.832, 0.678, 0.663, 0.886, 0.837, 0.682, 0.312, 0.93, 0.848, 0.957, 0.601, 0.8, 0.711, 0.918, 0.715, 0.54, 0.833, 0.885, 0.625, 0.809, 0.672, 0.774, 0.867, 0.994, 0.438, 0.996, 0.859, 0.953, 0.817, 0.921, 0.844, 0.993, 0.173, 0.821, 0.75, 0.81, 0.916, 0.792, 0.78, 0.943, 0.724, 0.948, 0.834, 0.791, 0.649]
global origin = 1
global destination = 35
