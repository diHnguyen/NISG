global arcs = [1 46; 2 22; 3 21; 4 27; 4 32; 4 41; 5 2; 5 8; 5 38; 5 40; 5 41; 6 5; 6 8; 6 49; 7 2; 7 6; 7 9; 7 21; 7 22; 7 33; 7 47; 8 4; 8 38; 8 39; 8 43; 9 24; 9 29; 9 32; 9 33; 9 44; 10 15; 10 17; 10 19; 10 33; 10 41; 11 12; 11 22; 11 32; 11 41; 11 45; 11 46; 12 15; 12 16; 12 22; 12 36; 12 46; 13 10; 13 11; 13 12; 13 21; 13 22; 13 35; 13 40; 14 8; 14 38; 14 39; 15 3; 15 8; 15 16; 16 3; 16 8; 16 15; 16 25; 16 33; 16 39; 16 44; 16 49; 17 2; 17 10; 17 11; 17 15; 17 29; 17 34; 17 41; 17 47; 18 13; 18 15; 18 22; 18 23; 18 25; 18 46; 19 12; 19 17; 19 18; 19 21; 20 5; 20 21; 20 26; 20 28; 20 35; 20 38; 20 44; 20 46; 21 33; 22 3; 22 9; 22 16; 22 27; 22 41; 23 19; 23 20; 23 24; 23 38; 24 6; 24 15; 24 34; 24 39; 25 2; 25 15; 25 28; 25 42; 25 48; 26 5; 26 28; 26 31; 26 34; 26 40; 27 2; 27 9; 27 17; 27 21; 27 45; 28 7; 28 9; 28 22; 28 42; 29 17; 30 2; 30 9; 30 24; 30 47; 31 19; 31 22; 31 23; 31 25; 31 38; 31 40; 31 44; 31 48; 32 26; 32 31; 32 38; 32 48; 32 49; 33 47; 34 17; 34 28; 34 36; 35 14; 35 32; 35 37; 35 48; 35 49; 36 10; 36 12; 36 17; 36 27; 36 46; 36 48; 37 15; 37 20; 37 30; 37 34; 38 6; 38 17; 38 18; 38 19; 38 23; 38 27; 38 47; 39 17; 39 36; 39 45; 40 3; 40 7; 40 32; 40 37; 40 48; 41 16; 41 30; 41 34; 41 43; 42 4; 42 5; 42 39; 43 19; 43 21; 43 34; 43 50; 44 3; 44 8; 44 12; 44 22; 44 24; 44 46; 45 7; 45 10; 45 14; 45 23; 45 30; 45 31; 45 36; 45 39; 46 12; 46 26; 46 42; 46 49; 47 10; 47 26; 47 35; 47 36; 47 43; 47 49; 48 20; 48 25; 48 33; 48 38; 48 44; 48 50; 49 12; 49 34; 49 43]
global d_x = [8.0, 1.0, 10.0, 8.0, 6.0, 2.0, 3.0, 4.0, 7.0, 7.0, 1.0, 3.0, 1.0, 2.0, 1.0, 5.0, 9.0, 1.0, 6.0, 8.0, 10.0, 6.0, 10.0, 3.0, 1.0, 2.0, 2.0, 5.0, 7.0, 7.0, 4.0, 2.0, 1.0, 1.0, 10.0, 3.0, 9.0, 3.0, 7.0, 5.0, 9.0, 1.0, 3.0, 5.0, 3.0, 8.0, 6.0, 1.0, 10.0, 10.0, 9.0, 3.0, 2.0, 10.0, 3.0, 9.0, 4.0, 4.0, 3.0, 6.0, 1.0, 7.0, 10.0, 5.0, 1.0, 4.0, 2.0, 10.0, 5.0, 5.0, 5.0, 1.0, 8.0, 3.0, 7.0, 3.0, 8.0, 1.0, 7.0, 3.0, 6.0, 1.0, 6.0, 8.0, 5.0, 2.0, 4.0, 9.0, 2.0, 6.0, 1.0, 6.0, 7.0, 9.0, 3.0, 3.0, 9.0, 9.0, 9.0, 6.0, 9.0, 7.0, 6.0, 10.0, 4.0, 8.0, 2.0, 5.0, 7.0, 4.0, 9.0, 5.0, 2.0, 2.0, 9.0, 2.0, 10.0, 5.0, 1.0, 5.0, 2.0, 1.0, 4.0, 3.0, 7.0, 8.0, 4.0, 8.0, 1.0, 7.0, 8.0, 3.0, 10.0, 6.0, 9.0, 6.0, 10.0, 2.0, 10.0, 2.0, 4.0, 10.0, 5.0, 4.0, 10.0, 9.0, 5.0, 10.0, 2.0, 8.0, 9.0, 1.0, 2.0, 5.0, 4.0, 3.0, 2.0, 1.0, 7.0, 2.0, 8.0, 9.0, 9.0, 8.0, 4.0, 7.0, 3.0, 10.0, 8.0, 5.0, 6.0, 9.0, 4.0, 5.0, 2.0, 7.0, 4.0, 1.0, 1.0, 3.0, 9.0, 9.0, 1.0, 8.0, 7.0, 5.0, 2.0, 2.0, 8.0, 2.0, 1.0, 7.0, 5.0, 1.0, 1.0, 3.0, 8.0, 7.0, 4.0, 6.0, 5.0, 4.0, 7.0, 6.0, 7.0, 7.0, 6.0, 9.0, 4.0, 9.0, 4.0, 7.0, 9.0, 4.0, 8.0, 3.0, 6.0, 2.0, 10.0, 4.0, 10.0, 3.0]
global b_x = 5
global d_y = [8.0, 7.0, 4.0, 7.0, 6.0, 9.0, 2.0, 1.0, 6.0, 6.0, 8.0, 5.0, 8.0, 1.0, 4.0, 9.0, 10.0, 7.0, 2.0, 1.0, 6.0, 4.0, 4.0, 4.0, 7.0, 6.0, 1.0, 3.0, 7.0, 6.0, 2.0, 5.0, 10.0, 7.0, 8.0, 4.0, 4.0, 5.0, 6.0, 1.0, 2.0, 6.0, 1.0, 5.0, 3.0, 7.0, 1.0, 7.0, 5.0, 2.0, 2.0, 8.0, 10.0, 5.0, 8.0, 7.0, 3.0, 2.0, 4.0, 2.0, 1.0, 7.0, 2.0, 4.0, 8.0, 1.0, 1.0, 5.0, 4.0, 2.0, 10.0, 7.0, 9.0, 5.0, 6.0, 5.0, 9.0, 6.0, 3.0, 6.0, 10.0, 1.0, 1.0, 2.0, 4.0, 9.0, 6.0, 2.0, 10.0, 8.0, 6.0, 1.0, 10.0, 6.0, 6.0, 8.0, 3.0, 6.0, 8.0, 4.0, 6.0, 4.0, 6.0, 9.0, 1.0, 7.0, 6.0, 6.0, 5.0, 1.0, 7.0, 2.0, 7.0, 2.0, 1.0, 4.0, 6.0, 10.0, 4.0, 7.0, 9.0, 2.0, 1.0, 5.0, 5.0, 4.0, 2.0, 7.0, 4.0, 4.0, 1.0, 5.0, 7.0, 6.0, 8.0, 1.0, 7.0, 5.0, 4.0, 3.0, 5.0, 9.0, 1.0, 8.0, 7.0, 4.0, 4.0, 10.0, 8.0, 8.0, 1.0, 5.0, 3.0, 5.0, 3.0, 9.0, 5.0, 9.0, 10.0, 7.0, 4.0, 6.0, 2.0, 5.0, 6.0, 6.0, 2.0, 6.0, 5.0, 1.0, 5.0, 4.0, 5.0, 5.0, 10.0, 6.0, 4.0, 6.0, 8.0, 1.0, 10.0, 3.0, 8.0, 1.0, 1.0, 10.0, 7.0, 5.0, 1.0, 1.0, 10.0, 9.0, 1.0, 3.0, 9.0, 6.0, 5.0, 3.0, 4.0, 7.0, 7.0, 1.0, 9.0, 10.0, 4.0, 4.0, 3.0, 3.0, 5.0, 6.0, 1.0, 10.0, 9.0, 5.0, 5.0, 10.0, 7.0, 1.0, 1.0, 1.0, 1.0, 7.0]
global b_y = 10
global p = [0.465, 0.62, 0.564, 0.162, 0.559, 0.596, 0.385, 0.333, 0.056, 0.326, 0.351, 0.383, 0.204, 0.782, 0.996, 0.838, 0.095, 0.752, 0.687, 0.289, 0.133, 0.894, 0.444, 0.997, 0.463, 0.61, 0.934, 0.87, 0.804, 0.045, 0.751, 0.741, 0.332, 0.641, 0.122, 0.298, 0.117, 0.04, 0.355, 0.077, 0.442, 0.056, 0.411, 0.966, 0.695, 0.038, 0.448, 0.882, 0.071, 0.108, 0.075, 0.238, 0.54, 0.081, 0.771, 0.107, 0.32, 0.382, 0.235, 0.852, 0.834, 0.407, 0.435, 0.952, 0.089, 0.291, 0.617, 0.877, 0.52, 0.539, 0.974, 0.919, 0.149, 0.594, 0.276, 0.565, 0.428, 0.788, 0.109, 0.217, 0.86, 0.116, 0.671, 0.11, 0.96, 0.653, 0.595, 0.978, 0.619, 0.4, 0.751, 0.685, 0.719, 0.915, 0.506, 0.589, 0.404, 0.589, 0.223, 0.32, 0.552, 0.277, 0.587, 0.825, 0.534, 0.557, 0.825, 0.662, 0.808, 0.039, 0.477, 0.107, 0.984, 0.636, 0.637, 0.193, 0.767, 0.262, 0.78, 0.468, 0.363, 0.487, 0.386, 0.586, 0.017, 0.119, 0.163, 0.489, 0.813, 0.548, 0.155, 0.771, 0.877, 0.87, 0.597, 0.206, 0.492, 0.022, 0.604, 0.124, 0.501, 0.003, 0.116, 0.562, 0.189, 0.961, 0.2, 0.483, 0.023, 0.028, 0.455, 0.196, 0.622, 0.743, 0.977, 0.805, 0.43, 0.503, 0.213, 0.698, 0.504, 0.795, 0.379, 0.51, 0.439, 0.957, 0.448, 0.149, 0.759, 0.41, 0.075, 0.967, 0.243, 0.814, 0.11, 0.946, 0.688, 0.26, 0.264, 0.122, 0.112, 0.595, 0.787, 0.688, 0.808, 0.768, 0.824, 0.875, 0.077, 0.26, 0.608, 0.4, 0.727, 0.941, 0.235, 0.918, 0.109, 0.081, 0.929, 0.218, 0.652, 0.86, 0.291, 0.87, 0.754, 0.705, 0.395, 0.91, 0.624, 0.062, 0.398, 0.187, 0.31, 0.141, 0.225, 0.971, 0.123, 0.053, 0.584, 0.72, 0.655, 0.957]
global q = [0.641, 0.778, 0.932, 0.439, 0.73, 0.676, 0.779, 0.421, 0.46, 0.935, 0.645, 0.536, 0.262, 0.875, 0.998, 0.904, 0.291, 0.938, 0.82, 0.969, 0.352, 0.99, 0.723, 0.997, 0.91, 0.814, 0.951, 0.936, 0.976, 0.861, 0.76, 0.893, 0.358, 0.945, 0.722, 0.912, 0.342, 0.488, 0.92, 0.733, 0.894, 0.717, 0.893, 0.98, 0.912, 0.759, 0.592, 0.921, 0.693, 0.431, 0.813, 0.414, 0.558, 0.192, 0.792, 0.799, 0.443, 0.425, 0.303, 0.949, 0.933, 0.808, 0.524, 0.969, 0.725, 0.49, 0.86, 0.94, 0.54, 0.933, 0.986, 0.993, 0.366, 0.781, 0.507, 0.957, 0.804, 0.997, 0.481, 0.481, 0.927, 0.246, 0.972, 0.864, 0.996, 0.899, 0.826, 0.988, 0.746, 0.512, 0.902, 0.991, 0.991, 0.942, 0.623, 0.673, 0.616, 0.781, 0.482, 0.802, 0.811, 0.581, 0.98, 0.964, 0.992, 0.786, 0.97, 0.725, 0.949, 0.596, 0.655, 0.344, 0.985, 0.653, 0.962, 0.246, 0.87, 0.948, 0.781, 0.617, 0.814, 0.804, 0.394, 0.765, 0.39, 0.902, 0.932, 0.8, 0.947, 0.974, 0.176, 0.822, 0.994, 0.945, 0.69, 0.461, 0.704, 0.903, 0.641, 0.752, 0.825, 0.417, 0.243, 0.918, 0.808, 0.976, 0.877, 0.735, 0.101, 0.652, 0.523, 0.76, 0.713, 0.973, 0.986, 0.812, 0.715, 0.636, 0.968, 0.744, 0.734, 0.833, 0.403, 0.955, 0.517, 0.961, 0.457, 0.525, 0.976, 0.676, 0.445, 0.976, 0.624, 0.917, 0.954, 0.965, 0.792, 0.377, 0.337, 0.689, 0.594, 0.684, 0.88, 0.805, 0.989, 0.769, 0.848, 0.948, 0.809, 0.828, 0.882, 0.56, 0.861, 0.985, 0.382, 0.95, 0.734, 0.462, 0.931, 0.776, 0.714, 0.968, 0.89, 0.874, 0.88, 0.81, 0.448, 0.991, 0.996, 0.356, 0.399, 0.952, 0.95, 0.208, 0.474, 0.994, 0.888, 0.32, 0.751, 0.788, 0.837, 0.958]
global origin = 1
global destination = 50