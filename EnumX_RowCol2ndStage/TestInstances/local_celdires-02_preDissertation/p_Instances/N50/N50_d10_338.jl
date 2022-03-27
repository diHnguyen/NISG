global arcs = [1 2; 1 9; 1 34; 1 44; 1 49; 2 48; 3 4; 3 13; 3 32; 3 35; 3 42; 4 5; 4 12; 4 17; 4 48; 5 2; 5 10; 5 38; 5 39; 5 40; 6 12; 6 23; 6 39; 7 4; 7 19; 7 24; 7 29; 7 33; 8 6; 8 12; 8 13; 8 14; 8 19; 8 21; 8 37; 8 39; 9 6; 9 8; 9 11; 9 22; 9 24; 9 27; 9 29; 9 32; 9 35; 9 37; 9 43; 9 44; 10 16; 10 37; 10 50; 11 17; 11 35; 11 49; 12 4; 12 6; 12 14; 12 24; 12 28; 12 49; 13 5; 13 16; 14 2; 14 15; 14 19; 14 28; 14 33; 14 38; 14 50; 15 16; 15 48; 16 20; 16 21; 16 50; 17 2; 17 42; 18 7; 18 14; 18 20; 18 22; 18 32; 18 46; 18 49; 19 35; 19 39; 20 2; 20 7; 20 24; 20 30; 20 31; 20 33; 20 39; 21 6; 21 16; 21 31; 21 41; 21 49; 22 2; 22 16; 22 17; 22 36; 22 42; 23 10; 23 28; 23 35; 23 49; 24 27; 24 28; 25 16; 25 28; 25 36; 25 41; 25 46; 26 6; 26 9; 26 11; 26 28; 26 36; 26 38; 27 7; 27 11; 27 15; 27 35; 27 40; 27 41; 27 45; 27 47; 27 50; 28 2; 28 3; 28 4; 28 9; 28 15; 28 37; 28 38; 29 4; 29 48; 30 2; 30 14; 30 50; 31 22; 31 23; 31 24; 31 43; 31 50; 32 16; 32 19; 32 20; 32 36; 33 18; 33 22; 33 32; 33 48; 34 10; 34 19; 34 29; 34 32; 35 4; 35 9; 35 13; 35 25; 35 31; 35 39; 36 8; 36 10; 36 16; 36 22; 36 44; 37 5; 37 6; 37 11; 37 20; 37 24; 37 44; 38 8; 38 40; 38 41; 39 6; 39 44; 40 21; 40 24; 40 27; 40 46; 40 49; 40 50; 41 17; 41 21; 41 37; 42 6; 42 36; 42 39; 43 3; 43 7; 43 16; 43 22; 43 36; 43 45; 43 48; 44 3; 44 23; 44 40; 44 48; 45 6; 45 13; 45 22; 45 24; 45 26; 45 32; 45 43; 45 49; 46 10; 46 14; 46 31; 46 33; 46 35; 46 37; 46 47; 47 3; 47 7; 47 8; 47 13; 47 18; 47 32; 47 43; 48 20; 48 22; 49 8; 49 11; 49 19; 49 38; 49 39]
global d_x = [6.0, 2.0, 6.0, 6.0, 7.0, 9.0, 8.0, 5.0, 5.0, 9.0, 9.0, 4.0, 5.0, 3.0, 6.0, 3.0, 4.0, 3.0, 2.0, 4.0, 10.0, 3.0, 6.0, 8.0, 3.0, 10.0, 6.0, 1.0, 1.0, 10.0, 1.0, 3.0, 9.0, 3.0, 6.0, 3.0, 5.0, 6.0, 2.0, 7.0, 10.0, 3.0, 9.0, 7.0, 6.0, 9.0, 9.0, 6.0, 6.0, 9.0, 8.0, 10.0, 2.0, 6.0, 7.0, 4.0, 7.0, 3.0, 6.0, 4.0, 9.0, 10.0, 7.0, 2.0, 2.0, 9.0, 7.0, 1.0, 6.0, 4.0, 2.0, 10.0, 1.0, 7.0, 3.0, 7.0, 1.0, 10.0, 8.0, 6.0, 1.0, 4.0, 2.0, 7.0, 3.0, 2.0, 2.0, 6.0, 6.0, 6.0, 2.0, 5.0, 8.0, 5.0, 10.0, 2.0, 6.0, 9.0, 9.0, 7.0, 9.0, 1.0, 7.0, 4.0, 5.0, 7.0, 7.0, 2.0, 10.0, 5.0, 10.0, 1.0, 3.0, 1.0, 5.0, 8.0, 7.0, 8.0, 5.0, 1.0, 4.0, 8.0, 1.0, 1.0, 4.0, 4.0, 3.0, 3.0, 10.0, 5.0, 2.0, 10.0, 6.0, 6.0, 5.0, 8.0, 8.0, 4.0, 8.0, 2.0, 3.0, 5.0, 4.0, 9.0, 1.0, 6.0, 9.0, 2.0, 8.0, 4.0, 4.0, 10.0, 4.0, 9.0, 7.0, 9.0, 6.0, 1.0, 9.0, 1.0, 2.0, 2.0, 8.0, 3.0, 2.0, 8.0, 4.0, 5.0, 3.0, 4.0, 7.0, 7.0, 7.0, 6.0, 2.0, 10.0, 5.0, 4.0, 8.0, 9.0, 8.0, 8.0, 10.0, 3.0, 8.0, 2.0, 3.0, 6.0, 10.0, 5.0, 5.0, 9.0, 10.0, 7.0, 1.0, 2.0, 10.0, 7.0, 1.0, 2.0, 8.0, 1.0, 7.0, 9.0, 9.0, 2.0, 5.0, 6.0, 4.0, 9.0, 10.0, 3.0, 10.0, 3.0, 1.0, 2.0, 6.0, 8.0, 2.0, 2.0, 1.0, 1.0, 2.0, 9.0, 10.0, 10.0, 8.0, 5.0, 4.0, 2.0, 7.0]
global b_x = 5
global d_y = [3.0, 8.0, 3.0, 5.0, 10.0, 4.0, 8.0, 3.0, 9.0, 3.0, 1.0, 1.0, 9.0, 2.0, 9.0, 2.0, 7.0, 5.0, 2.0, 9.0, 1.0, 3.0, 3.0, 2.0, 1.0, 2.0, 3.0, 3.0, 4.0, 4.0, 3.0, 10.0, 7.0, 3.0, 5.0, 7.0, 2.0, 10.0, 4.0, 5.0, 5.0, 8.0, 4.0, 4.0, 9.0, 8.0, 7.0, 7.0, 1.0, 2.0, 3.0, 6.0, 6.0, 8.0, 6.0, 5.0, 5.0, 9.0, 2.0, 3.0, 9.0, 3.0, 4.0, 9.0, 9.0, 3.0, 1.0, 9.0, 5.0, 9.0, 9.0, 6.0, 3.0, 5.0, 6.0, 6.0, 1.0, 6.0, 9.0, 2.0, 8.0, 7.0, 5.0, 9.0, 7.0, 6.0, 2.0, 4.0, 4.0, 5.0, 8.0, 3.0, 2.0, 10.0, 1.0, 8.0, 8.0, 8.0, 3.0, 7.0, 2.0, 2.0, 7.0, 7.0, 10.0, 10.0, 4.0, 9.0, 7.0, 6.0, 3.0, 2.0, 9.0, 9.0, 10.0, 5.0, 2.0, 6.0, 9.0, 9.0, 5.0, 8.0, 4.0, 1.0, 1.0, 7.0, 8.0, 7.0, 7.0, 4.0, 1.0, 8.0, 10.0, 5.0, 3.0, 7.0, 9.0, 4.0, 2.0, 8.0, 5.0, 10.0, 7.0, 4.0, 2.0, 5.0, 3.0, 2.0, 4.0, 10.0, 5.0, 7.0, 5.0, 9.0, 1.0, 10.0, 5.0, 5.0, 8.0, 9.0, 9.0, 10.0, 9.0, 6.0, 2.0, 5.0, 6.0, 4.0, 8.0, 9.0, 5.0, 1.0, 9.0, 2.0, 10.0, 7.0, 6.0, 6.0, 1.0, 8.0, 4.0, 6.0, 6.0, 1.0, 2.0, 3.0, 3.0, 10.0, 1.0, 8.0, 9.0, 1.0, 8.0, 8.0, 1.0, 2.0, 2.0, 2.0, 1.0, 4.0, 6.0, 5.0, 7.0, 3.0, 1.0, 8.0, 8.0, 9.0, 9.0, 3.0, 3.0, 9.0, 9.0, 6.0, 4.0, 7.0, 6.0, 1.0, 4.0, 9.0, 7.0, 4.0, 10.0, 10.0, 4.0, 6.0, 8.0, 3.0, 9.0, 10.0, 3.0]
global b_y = 10
global p = [0.949, 0.6, 0.096, 0.933, 0.481, 0.955, 0.535, 0.669, 0.055, 0.128, 0.795, 0.13, 0.821, 0.748, 0.714, 0.326, 0.666, 0.13, 0.618, 0.114, 0.019, 0.768, 0.336, 0.862, 0.871, 0.459, 0.197, 0.491, 0.766, 0.843, 0.509, 0.198, 0.407, 0.783, 0.112, 0.96, 0.699, 0.613, 0.523, 0.249, 0.624, 0.303, 0.878, 0.889, 0.387, 0.157, 0.859, 0.106, 0.61, 0.816, 0.85, 0.196, 0.211, 0.667, 0.395, 0.739, 0.466, 0.111, 0.696, 0.826, 0.009, 0.283, 0.883, 0.984, 0.649, 0.961, 0.818, 0.41, 0.926, 0.447, 0.329, 0.875, 0.783, 0.766, 0.652, 0.498, 0.903, 0.409, 0.807, 0.526, 0.071, 0.58, 0.425, 0.209, 0.621, 0.08, 0.855, 0.601, 0.076, 0.516, 0.293, 0.555, 0.614, 0.539, 0.79, 0.418, 0.475, 0.554, 0.164, 0.165, 0.398, 0.656, 0.295, 0.086, 0.388, 0.558, 0.247, 0.854, 0.383, 0.429, 0.208, 0.434, 0.062, 0.442, 0.539, 0.971, 0.434, 0.148, 0.526, 0.797, 0.591, 0.85, 0.983, 0.43, 0.211, 0.278, 0.078, 0.306, 0.884, 0.416, 0.347, 0.135, 0.979, 0.152, 0.122, 0.275, 0.321, 0.484, 0.909, 0.72, 0.059, 0.873, 0.166, 0.572, 0.039, 0.229, 0.346, 0.53, 0.259, 0.669, 0.078, 0.32, 0.099, 0.111, 0.333, 0.657, 0.782, 0.636, 0.785, 0.211, 0.521, 0.458, 0.651, 0.748, 0.903, 0.473, 0.356, 0.3, 0.687, 0.932, 0.395, 0.836, 0.897, 0.671, 0.327, 0.724, 0.446, 0.38, 0.149, 0.194, 0.962, 0.249, 0.659, 0.111, 0.146, 0.986, 0.758, 0.859, 0.599, 0.074, 0.483, 0.064, 0.624, 0.65, 0.613, 0.054, 0.221, 0.636, 0.775, 0.412, 0.245, 0.256, 0.928, 0.421, 0.439, 0.702, 0.23, 0.071, 0.181, 0.92, 0.784, 0.789, 0.301, 0.273, 0.999, 0.134, 0.643, 0.778, 0.171, 0.413, 0.717, 0.786, 0.153, 0.538, 0.452, 0.278, 0.325, 0.938, 0.223, 0.446, 0.439]
global q = [0.969, 0.67, 0.137, 0.951, 0.615, 0.976, 0.642, 0.922, 0.111, 0.717, 0.894, 0.464, 0.979, 0.818, 0.932, 0.377, 0.981, 0.955, 0.666, 0.951, 0.672, 0.818, 0.859, 0.874, 0.896, 0.764, 0.265, 0.852, 0.998, 0.87, 0.704, 0.784, 0.941, 0.851, 0.367, 0.999, 0.915, 0.714, 0.594, 0.488, 0.948, 0.469, 0.937, 0.957, 0.733, 0.816, 0.889, 0.302, 0.974, 0.907, 0.96, 0.824, 0.603, 0.954, 0.658, 0.787, 0.893, 0.602, 0.703, 0.883, 0.332, 0.479, 0.947, 0.989, 0.898, 0.973, 0.932, 0.507, 0.933, 0.859, 0.713, 0.92, 0.829, 0.781, 0.969, 0.836, 0.971, 0.69, 0.841, 0.625, 0.632, 0.817, 0.947, 0.556, 0.863, 0.354, 0.973, 0.682, 0.906, 0.826, 0.677, 0.701, 0.807, 0.817, 0.85, 0.789, 0.937, 0.591, 0.37, 0.445, 0.704, 0.861, 0.334, 0.715, 0.82, 0.875, 0.419, 0.95, 0.828, 0.571, 0.691, 0.825, 0.384, 0.561, 0.66, 0.985, 0.907, 0.897, 0.562, 0.866, 0.846, 0.883, 0.995, 0.696, 0.974, 0.649, 0.39, 0.451, 0.96, 0.787, 0.526, 0.587, 0.986, 0.883, 0.829, 0.678, 0.817, 0.918, 0.963, 0.843, 0.375, 0.949, 0.266, 0.939, 0.556, 0.901, 0.779, 0.772, 0.874, 0.682, 0.489, 0.996, 0.872, 0.348, 0.786, 0.864, 0.944, 0.998, 0.789, 0.461, 0.689, 0.633, 0.983, 0.868, 0.904, 0.61, 0.838, 0.809, 0.999, 0.967, 0.986, 0.9, 0.938, 0.988, 0.74, 0.966, 0.871, 0.816, 0.447, 0.459, 0.987, 0.837, 0.945, 0.14, 0.267, 0.992, 0.936, 0.961, 0.778, 0.788, 0.708, 0.903, 0.937, 0.954, 0.801, 0.148, 0.585, 0.956, 0.881, 0.46, 0.954, 0.947, 0.983, 0.487, 0.749, 0.974, 0.467, 0.393, 0.772, 0.945, 0.974, 0.883, 0.439, 0.419, 0.999, 0.918, 0.682, 0.972, 0.678, 0.571, 0.725, 0.896, 0.752, 0.659, 0.59, 0.522, 0.608, 0.982, 0.328, 0.599, 0.749]
global origin = 1
global destination = 50