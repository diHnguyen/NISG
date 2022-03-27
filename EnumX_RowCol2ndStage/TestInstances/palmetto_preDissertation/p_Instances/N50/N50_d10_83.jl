global arcs = [1 13; 1 16; 1 18; 1 23; 1 30; 1 34; 1 37; 1 42; 1 46; 2 5; 2 23; 2 25; 2 31; 2 41; 2 48; 3 4; 3 16; 3 36; 4 9; 4 33; 5 6; 5 7; 5 20; 5 38; 6 4; 6 16; 6 18; 6 31; 6 42; 6 45; 7 5; 7 21; 7 36; 7 39; 7 47; 8 5; 8 13; 8 16; 9 3; 9 27; 9 42; 9 48; 9 49; 10 7; 10 38; 11 19; 11 30; 12 2; 12 4; 12 15; 12 46; 12 48; 13 9; 13 10; 13 15; 13 28; 13 44; 14 5; 14 8; 14 12; 14 27; 14 32; 14 35; 14 39; 15 3; 15 12; 15 17; 15 19; 15 24; 15 33; 15 46; 15 50; 16 5; 16 36; 17 5; 17 7; 17 30; 17 48; 18 13; 18 26; 18 41; 19 8; 19 11; 19 14; 19 23; 19 33; 19 41; 19 44; 20 5; 20 13; 20 29; 20 47; 20 50; 21 4; 21 11; 21 26; 21 44; 22 15; 22 31; 22 32; 22 49; 23 9; 23 10; 23 15; 23 31; 23 49; 24 2; 24 12; 24 15; 24 46; 25 45; 26 3; 26 15; 26 37; 27 12; 27 14; 27 42; 27 49; 28 4; 28 11; 28 14; 28 35; 28 39; 28 47; 29 17; 29 21; 29 22; 29 23; 29 24; 30 7; 30 23; 30 26; 30 33; 30 48; 31 14; 31 30; 31 34; 31 39; 31 40; 31 42; 31 43; 31 47; 31 48; 32 13; 32 27; 32 45; 32 46; 32 47; 32 49; 33 3; 33 11; 33 15; 33 23; 34 15; 34 31; 34 36; 34 37; 34 38; 34 41; 34 50; 35 8; 35 11; 35 26; 35 47; 36 12; 36 14; 36 47; 37 4; 37 6; 37 14; 37 24; 37 47; 38 7; 38 20; 38 43; 39 15; 39 18; 39 25; 39 34; 39 44; 39 48; 40 11; 40 18; 40 31; 40 39; 40 43; 41 9; 41 12; 41 33; 41 40; 41 45; 42 4; 42 19; 42 48; 43 7; 43 19; 43 32; 44 8; 44 29; 44 38; 44 42; 44 50; 45 2; 45 7; 45 11; 45 13; 45 19; 45 34; 45 38; 45 49; 46 6; 46 18; 46 41; 46 50; 47 5; 47 18; 47 33; 48 9; 48 29; 48 44; 49 4; 49 29; 49 30; 49 48]
global d_x = [10.0, 1.0, 5.0, 5.0, 10.0, 8.0, 6.0, 6.0, 8.0, 2.0, 2.0, 6.0, 9.0, 7.0, 8.0, 9.0, 2.0, 6.0, 5.0, 4.0, 1.0, 6.0, 7.0, 5.0, 10.0, 3.0, 8.0, 7.0, 9.0, 8.0, 7.0, 5.0, 2.0, 5.0, 6.0, 10.0, 4.0, 9.0, 4.0, 2.0, 8.0, 7.0, 2.0, 5.0, 5.0, 6.0, 3.0, 7.0, 9.0, 9.0, 9.0, 7.0, 3.0, 6.0, 2.0, 2.0, 9.0, 1.0, 5.0, 5.0, 6.0, 4.0, 9.0, 5.0, 10.0, 2.0, 4.0, 8.0, 5.0, 1.0, 8.0, 7.0, 6.0, 4.0, 10.0, 8.0, 10.0, 10.0, 5.0, 1.0, 10.0, 3.0, 4.0, 4.0, 9.0, 2.0, 8.0, 2.0, 5.0, 2.0, 10.0, 1.0, 7.0, 6.0, 8.0, 3.0, 8.0, 4.0, 2.0, 7.0, 1.0, 8.0, 2.0, 5.0, 3.0, 7.0, 3.0, 8.0, 2.0, 3.0, 6.0, 3.0, 1.0, 9.0, 2.0, 6.0, 6.0, 7.0, 4.0, 4.0, 8.0, 1.0, 10.0, 7.0, 1.0, 9.0, 2.0, 2.0, 4.0, 5.0, 8.0, 3.0, 2.0, 5.0, 7.0, 2.0, 7.0, 6.0, 4.0, 2.0, 2.0, 6.0, 6.0, 2.0, 7.0, 5.0, 6.0, 9.0, 3.0, 9.0, 4.0, 9.0, 8.0, 6.0, 7.0, 4.0, 7.0, 5.0, 3.0, 3.0, 7.0, 5.0, 8.0, 1.0, 1.0, 10.0, 1.0, 8.0, 3.0, 10.0, 5.0, 8.0, 7.0, 3.0, 2.0, 8.0, 4.0, 5.0, 4.0, 10.0, 9.0, 8.0, 7.0, 3.0, 2.0, 3.0, 1.0, 5.0, 5.0, 4.0, 5.0, 9.0, 7.0, 6.0, 9.0, 4.0, 10.0, 4.0, 4.0, 1.0, 8.0, 1.0, 10.0, 4.0, 7.0, 7.0, 4.0, 3.0, 5.0, 10.0, 1.0, 9.0, 7.0, 8.0, 10.0, 9.0, 4.0, 3.0, 2.0, 5.0, 6.0, 2.0, 9.0, 1.0]
global b_x = 5
global d_y = [8.0, 6.0, 2.0, 5.0, 5.0, 2.0, 8.0, 5.0, 4.0, 10.0, 8.0, 2.0, 1.0, 2.0, 5.0, 2.0, 8.0, 7.0, 10.0, 4.0, 3.0, 4.0, 10.0, 6.0, 5.0, 1.0, 3.0, 6.0, 1.0, 7.0, 1.0, 2.0, 5.0, 8.0, 2.0, 9.0, 3.0, 1.0, 10.0, 5.0, 7.0, 5.0, 6.0, 3.0, 1.0, 8.0, 8.0, 7.0, 7.0, 9.0, 8.0, 6.0, 4.0, 7.0, 10.0, 9.0, 5.0, 3.0, 2.0, 2.0, 6.0, 7.0, 4.0, 5.0, 6.0, 5.0, 4.0, 5.0, 5.0, 9.0, 5.0, 7.0, 1.0, 4.0, 10.0, 4.0, 9.0, 2.0, 8.0, 8.0, 2.0, 5.0, 10.0, 3.0, 6.0, 1.0, 8.0, 7.0, 7.0, 8.0, 10.0, 9.0, 3.0, 10.0, 6.0, 6.0, 8.0, 9.0, 1.0, 4.0, 3.0, 4.0, 4.0, 6.0, 9.0, 7.0, 10.0, 5.0, 2.0, 7.0, 5.0, 10.0, 9.0, 9.0, 4.0, 7.0, 7.0, 5.0, 4.0, 2.0, 9.0, 1.0, 1.0, 9.0, 1.0, 7.0, 6.0, 2.0, 6.0, 6.0, 4.0, 8.0, 8.0, 2.0, 8.0, 8.0, 7.0, 2.0, 8.0, 2.0, 2.0, 2.0, 8.0, 3.0, 8.0, 9.0, 6.0, 4.0, 9.0, 7.0, 1.0, 4.0, 1.0, 7.0, 3.0, 2.0, 7.0, 3.0, 4.0, 9.0, 8.0, 9.0, 4.0, 6.0, 9.0, 5.0, 9.0, 8.0, 3.0, 4.0, 10.0, 6.0, 2.0, 3.0, 1.0, 5.0, 5.0, 3.0, 4.0, 10.0, 6.0, 9.0, 2.0, 8.0, 5.0, 4.0, 3.0, 1.0, 6.0, 3.0, 3.0, 1.0, 3.0, 6.0, 7.0, 4.0, 1.0, 3.0, 5.0, 1.0, 10.0, 7.0, 6.0, 2.0, 3.0, 7.0, 3.0, 5.0, 1.0, 5.0, 2.0, 4.0, 6.0, 10.0, 10.0, 9.0, 7.0, 5.0, 3.0, 3.0, 10.0, 4.0, 5.0, 4.0]
global b_y = 10
global p = [0.066, 0.015, 0.43, 0.877, 0.71, 0.039, 0.688, 0.408, 0.613, 0.514, 0.499, 0.185, 0.82, 0.512, 0.783, 0.008, 0.166, 0.924, 0.066, 0.819, 0.044, 0.394, 0.887, 0.44, 0.108, 0.399, 0.323, 0.363, 0.639, 0.182, 0.39, 0.027, 0.597, 0.636, 0.474, 0.411, 0.871, 0.161, 0.416, 0.184, 0.75, 0.962, 0.179, 0.318, 0.966, 0.674, 0.15, 0.671, 0.396, 0.496, 0.888, 0.631, 0.912, 0.299, 0.892, 0.842, 0.157, 0.836, 0.118, 0.007, 0.968, 0.979, 0.383, 0.09, 0.813, 0.562, 0.001, 0.517, 0.974, 0.079, 0.408, 0.689, 0.204, 0.281, 0.951, 0.705, 0.263, 0.042, 0.303, 0.254, 0.506, 0.489, 0.003, 0.345, 0.441, 0.249, 0.308, 0.211, 0.707, 0.607, 0.363, 0.963, 0.843, 0.017, 0.299, 0.276, 0.708, 0.55, 0.384, 0.457, 0.777, 0.143, 0.22, 0.657, 0.048, 0.567, 0.779, 0.693, 0.632, 0.03, 0.927, 0.428, 0.948, 0.746, 0.865, 0.752, 0.98, 0.39, 0.076, 0.847, 0.805, 0.434, 0.086, 0.327, 0.088, 0.508, 0.973, 0.305, 0.318, 0.442, 0.863, 0.475, 0.484, 0.479, 0.078, 0.735, 0.662, 0.093, 0.721, 0.939, 0.189, 0.591, 0.596, 0.126, 0.952, 0.631, 0.507, 0.227, 0.171, 0.576, 0.695, 0.597, 0.213, 0.651, 0.445, 0.672, 0.865, 0.783, 0.093, 0.464, 0.429, 0.641, 0.501, 0.488, 0.396, 0.834, 0.752, 0.2, 0.409, 0.571, 0.789, 0.863, 0.047, 0.022, 0.633, 0.853, 0.05, 0.028, 0.668, 0.695, 0.265, 0.381, 0.214, 0.016, 0.742, 0.615, 0.263, 0.362, 0.174, 0.939, 0.047, 0.76, 0.576, 0.259, 0.965, 0.065, 0.496, 0.329, 0.521, 0.017, 0.477, 0.022, 0.451, 0.591, 0.802, 0.728, 0.441, 0.662, 0.074, 0.361, 0.652, 0.794, 0.617, 0.068, 0.878, 0.942, 0.938, 0.032, 0.29, 0.982, 0.729, 0.098, 0.134, 0.354]
global q = [0.223, 0.824, 0.449, 0.995, 0.72, 0.719, 0.859, 0.692, 0.695, 0.625, 0.687, 0.64, 0.916, 0.633, 0.963, 0.63, 0.318, 0.961, 0.372, 0.986, 0.489, 0.507, 0.96, 0.541, 0.708, 0.815, 0.71, 0.382, 0.951, 0.605, 0.761, 0.153, 0.725, 0.642, 0.82, 0.761, 0.911, 0.918, 0.841, 0.379, 0.753, 0.974, 0.996, 0.623, 0.992, 0.682, 0.172, 0.781, 0.559, 0.594, 0.925, 0.643, 0.998, 0.916, 0.914, 0.915, 0.545, 0.91, 0.259, 0.679, 0.983, 0.991, 0.457, 0.139, 0.904, 0.962, 0.42, 0.838, 0.992, 0.464, 0.718, 0.978, 0.814, 0.716, 0.993, 0.825, 0.298, 0.388, 0.998, 0.763, 0.564, 0.61, 0.693, 0.802, 0.478, 0.903, 0.943, 0.682, 0.905, 0.651, 0.781, 0.964, 0.877, 0.897, 0.786, 0.985, 0.962, 0.68, 0.94, 0.515, 0.909, 0.461, 0.655, 0.915, 0.39, 0.952, 0.93, 0.719, 0.909, 0.655, 0.987, 0.624, 0.954, 0.912, 0.875, 0.807, 0.994, 0.778, 0.465, 0.88, 0.951, 0.947, 0.517, 0.994, 0.374, 0.776, 0.99, 0.658, 0.916, 0.556, 0.946, 0.558, 0.676, 0.833, 0.622, 0.933, 0.971, 0.29, 0.83, 0.965, 0.417, 0.921, 0.944, 0.409, 0.987, 0.737, 0.752, 0.561, 0.957, 0.879, 0.842, 0.896, 0.328, 0.968, 0.811, 0.799, 0.929, 0.981, 0.452, 0.798, 0.472, 0.763, 0.918, 0.891, 0.53, 0.956, 0.804, 0.571, 0.55, 0.855, 0.811, 0.963, 0.548, 0.401, 0.84, 0.909, 0.41, 0.426, 0.931, 0.852, 0.882, 0.865, 0.313, 0.724, 0.842, 0.979, 0.461, 0.381, 0.78, 0.965, 0.472, 0.801, 0.637, 0.729, 0.982, 0.129, 0.815, 0.378, 0.851, 0.909, 0.949, 0.116, 0.665, 0.886, 0.835, 0.832, 0.665, 0.728, 0.774, 0.91, 0.688, 0.874, 0.666, 0.555, 0.902, 0.974, 0.938, 0.863, 0.577, 0.999, 0.772, 0.841, 0.758, 0.717]
global origin = 1
global destination = 50