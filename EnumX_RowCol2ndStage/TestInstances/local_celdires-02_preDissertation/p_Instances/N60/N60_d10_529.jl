global arcs = [1 5; 1 6; 1 16; 1 27; 1 45; 2 16; 2 17; 2 19; 2 26; 2 28; 2 49; 2 57; 3 31; 3 42; 3 44; 3 52; 3 54; 3 60; 4 7; 4 13; 4 15; 4 25; 4 32; 4 41; 5 2; 5 9; 5 31; 5 47; 5 57; 6 4; 6 16; 6 22; 6 34; 6 53; 6 54; 7 14; 7 17; 7 29; 7 45; 7 59; 8 2; 8 28; 8 31; 8 38; 8 40; 8 41; 8 43; 8 48; 9 3; 9 30; 9 42; 9 53; 10 3; 10 7; 10 20; 10 42; 10 56; 11 10; 11 15; 11 24; 11 32; 11 39; 11 58; 12 8; 12 14; 12 15; 12 16; 12 39; 12 46; 12 53; 13 5; 13 6; 13 9; 13 14; 13 16; 13 18; 13 23; 13 28; 13 37; 13 43; 13 52; 14 35; 14 37; 14 46; 15 26; 15 27; 15 38; 15 47; 16 4; 16 6; 16 14; 16 22; 16 27; 16 40; 16 50; 17 35; 17 40; 17 52; 18 7; 18 15; 18 19; 18 24; 18 35; 18 41; 18 47; 18 58; 19 18; 19 39; 19 41; 19 54; 19 55; 20 2; 20 53; 21 12; 21 35; 21 45; 22 14; 22 16; 22 23; 22 26; 22 29; 22 47; 23 10; 23 13; 23 16; 23 17; 23 25; 23 30; 23 34; 23 44; 24 6; 24 11; 24 13; 24 38; 24 56; 25 17; 25 32; 25 45; 26 22; 26 46; 26 54; 26 60; 27 11; 27 12; 27 16; 27 19; 27 47; 27 59; 28 6; 28 12; 28 18; 28 50; 29 12; 29 34; 29 54; 30 4; 30 13; 30 17; 30 25; 30 44; 30 50; 31 15; 31 16; 31 32; 31 36; 31 44; 31 45; 31 46; 31 49; 31 54; 32 5; 32 6; 32 26; 32 37; 32 39; 32 46; 32 48; 32 58; 32 60; 33 5; 33 12; 33 18; 33 22; 33 27; 33 31; 33 47; 33 53; 33 59; 34 2; 34 6; 34 28; 34 31; 34 37; 34 39; 35 3; 35 25; 35 29; 35 33; 35 37; 35 45; 35 53; 35 59; 35 60; 36 12; 36 28; 36 32; 36 55; 37 4; 37 11; 37 26; 37 38; 37 60; 38 29; 38 48; 38 55; 39 38; 40 12; 40 55; 41 7; 41 9; 41 29; 41 31; 41 33; 41 34; 42 2; 42 4; 42 17; 42 47; 42 50; 43 4; 43 12; 43 18; 43 26; 43 28; 43 30; 43 35; 43 42; 43 49; 44 8; 44 10; 44 24; 44 26; 44 33; 44 39; 45 11; 45 42; 46 5; 46 8; 46 9; 46 31; 46 39; 46 43; 46 51; 47 10; 47 11; 47 16; 47 27; 47 37; 47 41; 47 58; 47 60; 48 9; 48 11; 48 12; 48 13; 48 31; 48 32; 48 35; 48 44; 48 51; 48 60; 49 15; 49 16; 49 18; 49 23; 50 21; 50 28; 50 31; 50 48; 51 12; 51 24; 51 30; 51 32; 51 34; 51 41; 51 49; 52 7; 52 15; 52 25; 52 27; 52 35; 53 10; 53 16; 53 26; 53 30; 53 32; 53 34; 53 51; 54 8; 54 18; 54 21; 54 52; 54 57; 55 11; 55 32; 55 35; 55 37; 55 44; 55 60; 56 15; 56 31; 57 22; 57 33; 57 41; 58 19; 58 23; 58 25; 58 32; 58 33; 58 37; 58 38; 58 40; 58 41; 58 43; 59 9; 59 12; 59 14]
global d_x = [10.0, 3.0, 4.0, 8.0, 4.0, 7.0, 6.0, 2.0, 2.0, 8.0, 5.0, 1.0, 6.0, 6.0, 10.0, 8.0, 9.0, 1.0, 6.0, 7.0, 4.0, 4.0, 8.0, 4.0, 1.0, 2.0, 4.0, 10.0, 2.0, 3.0, 6.0, 1.0, 1.0, 7.0, 8.0, 2.0, 8.0, 7.0, 7.0, 10.0, 4.0, 9.0, 2.0, 7.0, 6.0, 3.0, 2.0, 10.0, 6.0, 7.0, 6.0, 7.0, 4.0, 4.0, 10.0, 1.0, 7.0, 1.0, 5.0, 4.0, 6.0, 3.0, 7.0, 7.0, 2.0, 6.0, 3.0, 9.0, 5.0, 6.0, 3.0, 1.0, 5.0, 6.0, 4.0, 3.0, 1.0, 10.0, 7.0, 8.0, 7.0, 10.0, 7.0, 10.0, 9.0, 2.0, 2.0, 9.0, 1.0, 5.0, 1.0, 7.0, 4.0, 6.0, 6.0, 8.0, 8.0, 7.0, 9.0, 6.0, 4.0, 2.0, 1.0, 5.0, 3.0, 8.0, 9.0, 1.0, 10.0, 1.0, 6.0, 8.0, 3.0, 6.0, 2.0, 8.0, 2.0, 5.0, 5.0, 1.0, 4.0, 1.0, 5.0, 9.0, 3.0, 7.0, 4.0, 7.0, 5.0, 5.0, 10.0, 4.0, 1.0, 8.0, 7.0, 4.0, 6.0, 8.0, 8.0, 7.0, 7.0, 3.0, 9.0, 5.0, 5.0, 2.0, 8.0, 4.0, 7.0, 3.0, 7.0, 8.0, 10.0, 2.0, 4.0, 8.0, 9.0, 2.0, 5.0, 2.0, 7.0, 10.0, 7.0, 10.0, 7.0, 9.0, 4.0, 2.0, 4.0, 9.0, 8.0, 4.0, 2.0, 5.0, 7.0, 4.0, 5.0, 5.0, 9.0, 6.0, 4.0, 3.0, 3.0, 6.0, 5.0, 2.0, 3.0, 5.0, 2.0, 1.0, 1.0, 9.0, 7.0, 8.0, 3.0, 3.0, 7.0, 4.0, 9.0, 2.0, 7.0, 7.0, 8.0, 7.0, 4.0, 7.0, 4.0, 7.0, 5.0, 3.0, 10.0, 2.0, 4.0, 6.0, 4.0, 10.0, 10.0, 5.0, 2.0, 4.0, 10.0, 1.0, 7.0, 10.0, 10.0, 5.0, 10.0, 9.0, 7.0, 1.0, 2.0, 2.0, 6.0, 1.0, 4.0, 6.0, 10.0, 4.0, 9.0, 8.0, 8.0, 2.0, 3.0, 2.0, 7.0, 6.0, 4.0, 9.0, 3.0, 6.0, 9.0, 8.0, 6.0, 8.0, 4.0, 3.0, 5.0, 3.0, 9.0, 3.0, 5.0, 7.0, 4.0, 8.0, 9.0, 4.0, 7.0, 3.0, 1.0, 6.0, 4.0, 5.0, 10.0, 5.0, 9.0, 6.0, 6.0, 9.0, 3.0, 9.0, 2.0, 10.0, 3.0, 9.0, 1.0, 2.0, 2.0, 2.0, 3.0, 7.0, 7.0, 2.0, 5.0, 9.0, 7.0, 10.0, 5.0, 9.0, 9.0, 3.0, 3.0, 1.0, 7.0, 9.0, 5.0, 10.0, 7.0, 6.0, 10.0, 5.0, 1.0, 6.0, 3.0, 4.0, 9.0, 6.0, 1.0, 2.0, 10.0, 1.0, 2.0, 2.0, 9.0, 4.0, 4.0, 10.0, 2.0]
global b_x = 5
global d_y = [2.0, 3.0, 1.0, 1.0, 9.0, 7.0, 9.0, 10.0, 4.0, 4.0, 9.0, 2.0, 9.0, 3.0, 5.0, 8.0, 3.0, 2.0, 4.0, 1.0, 4.0, 9.0, 5.0, 3.0, 5.0, 8.0, 5.0, 6.0, 7.0, 1.0, 9.0, 4.0, 7.0, 3.0, 10.0, 7.0, 9.0, 8.0, 8.0, 5.0, 3.0, 3.0, 9.0, 5.0, 8.0, 2.0, 1.0, 9.0, 8.0, 5.0, 10.0, 6.0, 7.0, 3.0, 6.0, 8.0, 3.0, 10.0, 10.0, 3.0, 10.0, 8.0, 4.0, 10.0, 8.0, 2.0, 9.0, 9.0, 2.0, 8.0, 9.0, 4.0, 7.0, 10.0, 6.0, 2.0, 5.0, 5.0, 4.0, 5.0, 2.0, 8.0, 6.0, 10.0, 6.0, 10.0, 6.0, 5.0, 7.0, 6.0, 9.0, 3.0, 10.0, 3.0, 4.0, 7.0, 2.0, 2.0, 4.0, 7.0, 8.0, 5.0, 6.0, 2.0, 6.0, 9.0, 3.0, 1.0, 5.0, 8.0, 1.0, 10.0, 1.0, 2.0, 2.0, 2.0, 2.0, 7.0, 2.0, 2.0, 8.0, 4.0, 5.0, 5.0, 4.0, 7.0, 9.0, 4.0, 5.0, 5.0, 5.0, 8.0, 7.0, 7.0, 4.0, 4.0, 10.0, 10.0, 7.0, 6.0, 7.0, 5.0, 1.0, 9.0, 10.0, 4.0, 1.0, 10.0, 6.0, 2.0, 5.0, 4.0, 8.0, 7.0, 9.0, 8.0, 2.0, 1.0, 4.0, 7.0, 4.0, 3.0, 8.0, 10.0, 9.0, 9.0, 1.0, 4.0, 3.0, 1.0, 6.0, 6.0, 10.0, 8.0, 6.0, 3.0, 2.0, 3.0, 5.0, 2.0, 7.0, 8.0, 10.0, 1.0, 7.0, 8.0, 1.0, 1.0, 3.0, 1.0, 5.0, 6.0, 3.0, 9.0, 6.0, 8.0, 9.0, 10.0, 4.0, 9.0, 2.0, 5.0, 8.0, 10.0, 3.0, 4.0, 2.0, 7.0, 2.0, 1.0, 8.0, 9.0, 10.0, 6.0, 5.0, 4.0, 1.0, 8.0, 7.0, 9.0, 2.0, 10.0, 2.0, 5.0, 10.0, 3.0, 3.0, 2.0, 3.0, 9.0, 2.0, 6.0, 2.0, 8.0, 9.0, 3.0, 2.0, 4.0, 5.0, 5.0, 8.0, 5.0, 2.0, 8.0, 1.0, 9.0, 2.0, 5.0, 1.0, 8.0, 5.0, 4.0, 5.0, 6.0, 4.0, 4.0, 9.0, 8.0, 9.0, 10.0, 3.0, 8.0, 2.0, 7.0, 7.0, 7.0, 8.0, 5.0, 1.0, 3.0, 9.0, 8.0, 10.0, 2.0, 9.0, 3.0, 2.0, 8.0, 5.0, 3.0, 8.0, 8.0, 5.0, 3.0, 2.0, 7.0, 2.0, 3.0, 4.0, 6.0, 10.0, 1.0, 4.0, 3.0, 10.0, 10.0, 4.0, 7.0, 6.0, 1.0, 3.0, 1.0, 9.0, 5.0, 10.0, 7.0, 1.0, 8.0, 5.0, 3.0, 2.0, 5.0, 6.0, 8.0, 10.0, 3.0, 1.0, 4.0, 3.0, 4.0, 6.0, 10.0, 5.0, 7.0, 4.0, 3.0, 8.0]
global b_y = 10
global p = [0.857, 0.736, 0.708, 0.558, 0.995, 0.094, 0.737, 0.15, 0.432, 0.66, 0.428, 0.042, 0.972, 0.24, 0.67, 0.897, 0.868, 0.261, 0.331, 0.656, 0.678, 0.323, 0.924, 0.066, 0.486, 0.998, 0.652, 0.855, 0.595, 0.348, 0.393, 0.191, 0.408, 0.442, 0.016, 0.067, 0.806, 0.274, 0.433, 0.692, 0.913, 0.065, 0.49, 0.548, 0.35, 0.841, 0.152, 0.36, 0.678, 0.683, 0.375, 0.389, 0.403, 0.504, 0.652, 0.967, 0.042, 0.542, 0.306, 0.913, 0.517, 0.951, 0.612, 0.581, 0.713, 0.891, 0.211, 0.535, 0.108, 0.758, 0.289, 0.97, 0.513, 0.43, 0.07, 0.377, 0.93, 0.387, 0.438, 0.15, 0.202, 0.172, 0.009, 0.793, 0.197, 0.876, 0.383, 0.282, 0.9, 0.482, 0.54, 0.759, 0.666, 0.587, 0.495, 0.286, 0.563, 0.036, 0.954, 0.296, 0.653, 0.806, 0.322, 0.58, 0.532, 0.047, 0.635, 0.646, 0.519, 0.658, 0.132, 0.662, 0.659, 0.091, 0.899, 0.53, 0.329, 0.893, 0.45, 0.98, 0.323, 0.096, 0.757, 0.769, 0.225, 0.867, 0.021, 0.897, 0.294, 0.966, 0.101, 0.119, 0.005, 0.697, 0.45, 0.988, 0.142, 0.527, 0.789, 0.904, 0.916, 0.413, 0.038, 0.737, 0.819, 0.104, 0.616, 0.826, 0.737, 0.163, 0.257, 0.471, 0.505, 0.438, 0.229, 0.114, 0.689, 0.116, 0.198, 0.598, 0.875, 0.215, 0.165, 0.159, 0.006, 0.79, 0.063, 0.74, 0.278, 0.604, 0.214, 0.652, 0.551, 0.914, 0.176, 0.779, 0.606, 0.047, 0.956, 0.707, 0.062, 0.253, 0.976, 0.757, 0.77, 0.676, 0.072, 0.394, 0.788, 0.724, 0.166, 0.19, 0.79, 0.643, 0.202, 0.941, 0.599, 0.84, 0.382, 0.589, 0.673, 0.335, 0.792, 0.079, 0.554, 0.268, 0.808, 0.16, 0.681, 0.558, 0.652, 0.417, 0.945, 0.198, 0.16, 0.913, 0.72, 0.753, 0.469, 0.811, 0.11, 0.683, 0.919, 0.488, 0.008, 0.568, 0.027, 0.428, 0.347, 0.029, 0.175, 0.363, 0.844, 0.15, 0.285, 0.324, 0.072, 0.733, 0.276, 0.121, 0.278, 0.773, 0.52, 0.618, 0.822, 0.286, 0.175, 0.164, 0.18, 0.035, 0.362, 0.038, 0.68, 0.035, 0.203, 0.803, 0.092, 0.626, 0.568, 0.812, 0.793, 0.486, 0.268, 0.796, 0.249, 0.001, 0.427, 0.486, 0.187, 0.634, 0.706, 0.732, 0.813, 0.706, 0.622, 0.047, 0.52, 0.262, 0.641, 0.364, 0.13, 0.69, 0.382, 0.514, 0.254, 0.235, 0.318, 0.626, 0.908, 0.581, 0.536, 0.512, 0.931, 0.944, 0.994, 0.683, 0.833, 0.867, 0.054, 0.423, 0.801, 0.999, 0.816, 0.76, 0.965, 0.233, 0.638, 0.268, 0.104, 0.611, 0.501, 0.46, 0.548, 0.924, 0.218, 0.089, 0.876, 0.927, 0.97, 0.371, 0.209, 0.603, 0.116, 0.524, 0.851, 0.481, 0.57]
global q = [0.997, 0.984, 0.909, 0.59, 0.999, 0.104, 0.746, 0.934, 0.63, 0.765, 0.612, 0.26, 0.997, 0.694, 0.759, 0.908, 0.943, 0.285, 0.874, 0.72, 0.733, 0.395, 0.984, 0.124, 0.688, 0.998, 0.983, 0.893, 0.749, 0.484, 0.913, 0.57, 0.979, 0.659, 0.711, 0.212, 0.868, 0.558, 0.473, 0.945, 0.982, 0.266, 0.645, 0.585, 0.545, 0.979, 0.755, 0.793, 0.743, 0.958, 0.457, 0.783, 0.893, 0.557, 0.787, 0.983, 0.2, 0.923, 0.749, 0.932, 0.752, 0.957, 0.741, 0.822, 0.993, 0.997, 0.314, 0.757, 0.238, 0.975, 0.386, 0.981, 0.646, 0.976, 0.693, 0.674, 0.955, 0.755, 0.9, 0.476, 0.633, 0.372, 0.277, 0.955, 0.339, 0.925, 0.407, 0.743, 0.981, 0.568, 0.661, 0.821, 0.936, 0.664, 0.646, 0.333, 0.726, 0.39, 0.99, 0.64, 0.871, 0.835, 0.687, 0.836, 0.999, 0.348, 0.959, 0.86, 0.594, 0.832, 0.981, 0.868, 0.67, 0.969, 0.926, 0.95, 0.885, 0.986, 0.559, 0.981, 0.624, 0.906, 0.905, 0.917, 0.468, 0.992, 0.631, 0.999, 0.708, 0.986, 0.548, 0.649, 0.061, 0.948, 0.78, 0.999, 0.402, 0.82, 0.995, 0.913, 0.975, 0.82, 0.37, 0.802, 0.918, 0.227, 0.859, 0.939, 0.892, 0.532, 0.336, 0.662, 0.787, 0.57, 0.584, 0.301, 0.752, 0.322, 0.406, 0.729, 0.972, 0.282, 0.704, 0.867, 0.459, 0.866, 0.398, 0.8, 0.619, 0.963, 0.338, 0.699, 0.613, 0.915, 0.177, 0.879, 0.641, 0.715, 0.982, 0.812, 0.56, 0.861, 0.98, 0.923, 0.817, 0.746, 0.878, 0.608, 0.889, 0.773, 0.343, 0.318, 0.821, 0.931, 0.872, 0.969, 0.814, 0.992, 0.757, 0.648, 0.832, 0.409, 0.958, 0.778, 0.708, 0.299, 0.937, 0.578, 0.868, 0.794, 0.792, 0.783, 0.958, 0.276, 0.217, 0.98, 0.905, 0.893, 0.899, 0.899, 0.21, 0.949, 0.93, 0.775, 0.202, 0.875, 0.383, 0.858, 0.364, 0.849, 0.647, 0.917, 0.85, 0.359, 0.902, 0.458, 0.085, 0.8, 0.67, 0.758, 0.876, 0.84, 0.805, 0.875, 0.87, 0.59, 0.218, 0.83, 0.267, 0.673, 0.824, 0.784, 0.843, 0.098, 0.728, 0.923, 0.46, 0.875, 0.84, 0.872, 0.864, 0.744, 0.339, 0.989, 0.313, 0.243, 0.687, 0.656, 0.477, 0.926, 0.81, 0.741, 0.837, 0.708, 0.875, 0.428, 0.917, 0.348, 0.989, 0.487, 0.69, 0.73, 0.442, 0.704, 0.888, 0.289, 0.641, 0.718, 0.963, 0.881, 0.654, 0.615, 0.931, 0.987, 0.997, 0.747, 0.939, 0.872, 0.255, 0.603, 0.941, 0.999, 0.999, 0.986, 0.991, 0.877, 0.982, 0.604, 0.353, 0.993, 0.978, 0.539, 0.758, 0.949, 0.507, 0.691, 0.96, 0.977, 0.978, 0.559, 0.437, 0.612, 0.988, 0.899, 0.931, 0.557, 0.86]
global origin = 1
global destination = 60