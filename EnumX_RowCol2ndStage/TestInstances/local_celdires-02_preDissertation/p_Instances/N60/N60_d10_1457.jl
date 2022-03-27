global arcs = [1 15; 1 48; 1 56; 1 58; 2 5; 2 37; 2 43; 2 58; 3 18; 3 42; 3 53; 4 6; 4 52; 4 57; 5 3; 5 9; 5 27; 5 43; 5 48; 5 51; 5 59; 6 9; 6 19; 6 29; 6 30; 6 40; 7 12; 7 14; 7 20; 7 46; 7 47; 8 7; 8 13; 8 20; 8 25; 8 29; 8 30; 8 36; 8 38; 8 42; 8 53; 8 55; 9 6; 9 16; 9 24; 9 26; 9 30; 9 40; 9 46; 10 8; 10 12; 10 29; 10 40; 10 46; 10 58; 11 6; 11 14; 11 22; 11 30; 11 33; 11 34; 11 45; 11 59; 11 60; 12 13; 12 33; 12 42; 12 56; 12 57; 13 15; 13 28; 13 43; 13 48; 13 53; 14 4; 14 11; 14 30; 14 31; 14 34; 15 7; 15 12; 15 16; 15 25; 15 37; 15 50; 15 51; 15 57; 16 3; 16 8; 16 19; 16 53; 16 58; 17 25; 17 35; 17 43; 17 58; 18 2; 18 23; 18 32; 18 58; 19 8; 19 27; 19 37; 19 41; 19 45; 19 53; 20 9; 20 10; 20 13; 20 31; 20 32; 21 5; 21 6; 21 7; 21 28; 21 29; 21 37; 21 41; 22 10; 22 11; 22 24; 22 35; 22 53; 22 55; 23 6; 23 21; 23 35; 23 50; 23 59; 24 3; 24 6; 24 20; 24 44; 25 13; 25 26; 25 46; 26 10; 26 19; 26 39; 26 53; 26 55; 27 2; 27 5; 27 9; 27 15; 27 39; 27 46; 27 50; 27 53; 28 2; 28 4; 28 39; 28 43; 29 9; 29 28; 29 32; 29 47; 29 55; 30 25; 30 27; 30 32; 30 55; 31 11; 31 26; 31 29; 31 40; 31 58; 32 3; 32 5; 32 31; 32 36; 32 40; 32 42; 33 11; 33 16; 33 30; 33 45; 33 47; 33 58; 34 22; 34 36; 34 56; 35 4; 35 5; 35 9; 35 11; 35 15; 35 25; 35 26; 35 31; 35 50; 36 17; 36 18; 36 39; 36 42; 36 44; 37 6; 37 21; 37 23; 38 2; 38 5; 38 18; 38 19; 38 25; 38 26; 38 33; 38 39; 38 45; 38 46; 38 53; 39 53; 39 57; 40 17; 40 23; 40 24; 40 28; 40 41; 40 57; 41 21; 41 23; 41 49; 42 6; 42 19; 42 22; 42 53; 42 55; 42 58; 43 12; 43 18; 43 29; 43 48; 43 54; 44 27; 44 47; 44 51; 45 3; 45 9; 45 13; 45 34; 45 44; 46 2; 46 13; 46 14; 46 19; 46 20; 46 36; 46 41; 46 49; 47 4; 47 17; 47 22; 47 50; 47 60; 48 2; 48 3; 48 9; 48 41; 48 58; 49 22; 49 23; 49 26; 49 32; 49 43; 49 53; 49 54; 49 60; 50 2; 50 10; 50 13; 50 59; 51 5; 51 9; 51 16; 51 26; 51 27; 51 49; 52 3; 52 16; 52 26; 52 39; 52 48; 52 49; 53 11; 53 19; 53 21; 53 36; 53 37; 53 43; 53 51; 53 60; 54 5; 54 22; 54 26; 54 29; 54 30; 54 31; 54 39; 54 43; 54 48; 54 51; 55 43; 55 48; 56 7; 56 15; 56 22; 56 26; 56 31; 56 32; 57 3; 57 9; 57 14; 57 46; 57 48; 57 50; 57 51; 58 2; 58 4; 58 18; 58 31; 58 34; 58 47; 58 49; 59 2; 59 3; 59 6; 59 12; 59 19; 59 25; 59 28; 59 42; 59 45; 59 51; 59 56]
global d_x = [1.0, 8.0, 10.0, 4.0, 6.0, 4.0, 5.0, 3.0, 2.0, 2.0, 10.0, 2.0, 10.0, 2.0, 5.0, 8.0, 4.0, 2.0, 9.0, 10.0, 6.0, 10.0, 2.0, 7.0, 2.0, 10.0, 7.0, 7.0, 4.0, 10.0, 10.0, 7.0, 8.0, 8.0, 3.0, 10.0, 1.0, 9.0, 4.0, 1.0, 8.0, 4.0, 5.0, 1.0, 9.0, 7.0, 2.0, 8.0, 5.0, 10.0, 2.0, 6.0, 8.0, 6.0, 9.0, 1.0, 7.0, 8.0, 5.0, 9.0, 7.0, 4.0, 10.0, 5.0, 3.0, 10.0, 10.0, 5.0, 7.0, 4.0, 7.0, 1.0, 9.0, 2.0, 7.0, 3.0, 5.0, 9.0, 4.0, 9.0, 4.0, 2.0, 8.0, 6.0, 7.0, 6.0, 1.0, 2.0, 3.0, 7.0, 6.0, 6.0, 3.0, 3.0, 2.0, 5.0, 1.0, 9.0, 6.0, 3.0, 7.0, 8.0, 7.0, 5.0, 7.0, 10.0, 5.0, 9.0, 8.0, 5.0, 5.0, 8.0, 4.0, 4.0, 9.0, 6.0, 8.0, 8.0, 9.0, 2.0, 4.0, 7.0, 10.0, 9.0, 7.0, 7.0, 8.0, 6.0, 2.0, 5.0, 7.0, 1.0, 10.0, 5.0, 3.0, 9.0, 2.0, 3.0, 1.0, 7.0, 9.0, 9.0, 5.0, 4.0, 9.0, 2.0, 9.0, 2.0, 7.0, 9.0, 5.0, 9.0, 10.0, 6.0, 10.0, 8.0, 6.0, 4.0, 9.0, 6.0, 9.0, 10.0, 9.0, 5.0, 9.0, 7.0, 5.0, 8.0, 7.0, 5.0, 10.0, 9.0, 6.0, 3.0, 8.0, 7.0, 6.0, 2.0, 6.0, 2.0, 9.0, 1.0, 10.0, 3.0, 4.0, 3.0, 2.0, 5.0, 3.0, 3.0, 10.0, 8.0, 5.0, 2.0, 7.0, 4.0, 5.0, 10.0, 1.0, 1.0, 9.0, 2.0, 6.0, 4.0, 10.0, 3.0, 9.0, 6.0, 5.0, 3.0, 10.0, 7.0, 6.0, 4.0, 6.0, 9.0, 4.0, 10.0, 8.0, 7.0, 1.0, 5.0, 5.0, 7.0, 9.0, 4.0, 6.0, 6.0, 9.0, 7.0, 1.0, 9.0, 1.0, 3.0, 2.0, 9.0, 2.0, 1.0, 4.0, 3.0, 4.0, 1.0, 9.0, 4.0, 9.0, 7.0, 6.0, 3.0, 6.0, 8.0, 9.0, 7.0, 4.0, 3.0, 8.0, 5.0, 7.0, 8.0, 2.0, 9.0, 7.0, 10.0, 5.0, 3.0, 3.0, 1.0, 1.0, 3.0, 10.0, 5.0, 2.0, 3.0, 10.0, 10.0, 9.0, 10.0, 6.0, 7.0, 8.0, 9.0, 8.0, 3.0, 5.0, 2.0, 10.0, 6.0, 3.0, 6.0, 7.0, 8.0, 4.0, 6.0, 8.0, 6.0, 4.0, 4.0, 4.0, 3.0, 5.0, 6.0, 1.0, 3.0, 8.0, 10.0, 5.0, 8.0, 6.0, 2.0, 10.0, 4.0, 10.0, 6.0, 2.0, 5.0, 9.0, 6.0, 8.0, 9.0, 3.0, 8.0, 3.0, 7.0, 5.0, 1.0, 9.0, 4.0, 2.0, 3.0, 3.0, 1.0, 7.0, 8.0, 5.0]
global b_x = 5
global d_y = [5.0, 8.0, 5.0, 10.0, 8.0, 9.0, 8.0, 3.0, 8.0, 3.0, 3.0, 8.0, 7.0, 5.0, 10.0, 8.0, 6.0, 10.0, 7.0, 8.0, 6.0, 1.0, 6.0, 2.0, 2.0, 1.0, 6.0, 2.0, 2.0, 3.0, 3.0, 6.0, 9.0, 6.0, 6.0, 9.0, 1.0, 5.0, 8.0, 7.0, 5.0, 1.0, 5.0, 6.0, 7.0, 2.0, 4.0, 2.0, 6.0, 1.0, 7.0, 5.0, 3.0, 5.0, 9.0, 6.0, 6.0, 3.0, 10.0, 8.0, 3.0, 6.0, 1.0, 4.0, 9.0, 4.0, 8.0, 5.0, 6.0, 2.0, 1.0, 7.0, 4.0, 6.0, 3.0, 3.0, 7.0, 2.0, 4.0, 10.0, 9.0, 8.0, 5.0, 2.0, 4.0, 1.0, 5.0, 7.0, 7.0, 7.0, 2.0, 3.0, 9.0, 8.0, 7.0, 6.0, 4.0, 7.0, 1.0, 4.0, 6.0, 4.0, 2.0, 2.0, 2.0, 5.0, 1.0, 8.0, 10.0, 2.0, 3.0, 9.0, 7.0, 8.0, 7.0, 4.0, 8.0, 1.0, 1.0, 7.0, 4.0, 6.0, 7.0, 3.0, 4.0, 9.0, 9.0, 6.0, 4.0, 7.0, 8.0, 5.0, 2.0, 10.0, 2.0, 9.0, 9.0, 10.0, 7.0, 7.0, 4.0, 10.0, 5.0, 3.0, 10.0, 7.0, 3.0, 2.0, 1.0, 8.0, 8.0, 8.0, 8.0, 5.0, 3.0, 1.0, 1.0, 4.0, 1.0, 10.0, 4.0, 9.0, 4.0, 6.0, 6.0, 5.0, 6.0, 4.0, 5.0, 9.0, 2.0, 2.0, 2.0, 1.0, 2.0, 5.0, 1.0, 5.0, 7.0, 4.0, 5.0, 10.0, 5.0, 10.0, 1.0, 5.0, 3.0, 9.0, 5.0, 7.0, 3.0, 9.0, 9.0, 2.0, 8.0, 4.0, 7.0, 4.0, 8.0, 10.0, 7.0, 3.0, 6.0, 8.0, 2.0, 9.0, 8.0, 4.0, 2.0, 9.0, 7.0, 3.0, 1.0, 3.0, 5.0, 3.0, 1.0, 1.0, 2.0, 3.0, 9.0, 7.0, 3.0, 4.0, 1.0, 4.0, 4.0, 6.0, 5.0, 5.0, 8.0, 9.0, 10.0, 1.0, 4.0, 1.0, 7.0, 6.0, 7.0, 6.0, 5.0, 2.0, 1.0, 2.0, 2.0, 9.0, 5.0, 4.0, 4.0, 3.0, 4.0, 2.0, 8.0, 7.0, 6.0, 3.0, 7.0, 10.0, 10.0, 9.0, 10.0, 6.0, 2.0, 6.0, 9.0, 2.0, 2.0, 6.0, 4.0, 4.0, 10.0, 8.0, 6.0, 5.0, 10.0, 8.0, 5.0, 4.0, 8.0, 9.0, 10.0, 3.0, 5.0, 7.0, 1.0, 10.0, 9.0, 7.0, 10.0, 6.0, 3.0, 6.0, 10.0, 2.0, 10.0, 8.0, 8.0, 1.0, 4.0, 10.0, 5.0, 10.0, 1.0, 6.0, 9.0, 9.0, 10.0, 7.0, 3.0, 4.0, 1.0, 6.0, 7.0, 3.0, 9.0, 8.0, 7.0, 2.0, 5.0, 8.0, 6.0, 2.0, 3.0, 9.0, 2.0, 4.0, 9.0, 3.0, 9.0, 7.0, 7.0, 4.0, 8.0]
global b_y = 10
global p = [0.085, 0.894, 0.158, 0.741, 0.241, 0.947, 0.37, 0.512, 0.176, 0.627, 0.504, 0.026, 0.071, 0.244, 0.504, 0.887, 0.561, 0.843, 0.821, 0.652, 0.286, 0.514, 0.488, 0.932, 0.201, 0.117, 0.548, 0.11, 0.471, 0.926, 0.967, 0.876, 0.596, 0.166, 0.296, 0.791, 0.396, 0.074, 0.662, 0.754, 0.922, 0.817, 0.578, 0.683, 0.437, 0.456, 0.38, 0.041, 0.087, 0.726, 0.806, 0.782, 0.087, 0.914, 0.11, 0.757, 0.575, 0.052, 0.357, 0.775, 0.542, 0.475, 0.502, 0.909, 0.547, 0.147, 0.695, 0.39, 0.65, 0.201, 0.562, 0.557, 0.74, 0.427, 0.977, 0.172, 0.784, 0.165, 0.044, 0.394, 0.124, 0.449, 0.563, 0.014, 0.716, 0.802, 0.958, 0.75, 0.059, 0.962, 0.425, 0.686, 0.393, 0.058, 0.454, 0.664, 0.949, 0.944, 0.141, 0.21, 0.868, 0.734, 0.879, 0.376, 0.07, 0.422, 0.687, 0.092, 0.372, 0.431, 0.222, 0.478, 0.254, 0.714, 0.535, 0.422, 0.437, 0.066, 0.295, 0.049, 0.302, 0.18, 0.506, 0.978, 0.48, 0.875, 0.625, 0.921, 0.939, 0.715, 0.052, 0.913, 0.522, 0.196, 0.008, 0.947, 0.784, 0.847, 0.757, 0.378, 0.695, 0.65, 0.977, 0.13, 0.399, 0.613, 0.903, 0.896, 0.807, 0.068, 0.117, 0.121, 0.817, 0.125, 0.347, 0.481, 0.887, 0.981, 0.314, 0.735, 0.396, 0.848, 0.746, 0.678, 0.487, 0.413, 0.505, 0.347, 0.331, 0.073, 0.442, 0.974, 0.829, 0.983, 0.344, 0.937, 0.148, 0.836, 0.404, 0.023, 0.437, 0.513, 0.07, 0.265, 0.498, 0.802, 0.929, 0.604, 0.505, 0.532, 0.426, 0.847, 0.962, 0.34, 0.836, 0.24, 0.207, 0.18, 0.352, 0.86, 0.621, 0.888, 0.816, 0.018, 0.296, 0.123, 0.273, 0.527, 0.016, 0.04, 0.938, 0.778, 0.527, 0.09, 0.838, 0.961, 0.167, 0.111, 0.109, 0.793, 0.086, 0.819, 0.889, 0.278, 0.622, 0.242, 0.855, 0.984, 0.49, 0.677, 0.513, 0.241, 0.453, 0.268, 0.621, 0.806, 0.102, 0.777, 0.584, 0.355, 0.083, 0.808, 0.915, 0.15, 0.331, 0.403, 0.418, 0.934, 0.391, 0.404, 0.003, 0.537, 0.354, 0.029, 0.251, 0.843, 0.015, 0.897, 0.981, 0.402, 0.577, 0.908, 0.968, 0.836, 0.71, 0.124, 0.747, 0.856, 0.488, 0.696, 0.527, 0.235, 0.322, 0.795, 0.446, 0.257, 0.565, 0.488, 0.311, 0.852, 0.391, 0.325, 0.152, 0.698, 0.812, 0.061, 0.793, 0.67, 0.645, 0.049, 0.573, 0.067, 0.611, 0.732, 0.844, 0.325, 0.881, 0.886, 0.041, 0.314, 0.764, 0.542, 0.918, 0.278, 0.451, 0.769, 0.546, 0.18, 0.542, 0.706, 0.498, 0.221, 0.35, 0.58, 0.262, 0.755, 0.163, 0.911, 0.285, 0.808, 0.929, 0.549, 0.771, 0.451, 0.967, 0.303, 0.591, 0.974, 0.297, 0.849, 0.742, 0.636, 0.034]
global q = [0.126, 0.919, 0.214, 0.9, 0.945, 0.987, 0.536, 0.527, 0.576, 0.703, 0.604, 0.942, 0.198, 0.353, 0.69, 0.909, 0.765, 0.981, 0.964, 0.747, 0.951, 0.723, 0.498, 0.983, 0.467, 0.737, 0.635, 0.754, 0.474, 0.969, 0.993, 0.998, 0.623, 0.799, 0.855, 0.918, 0.74, 0.632, 0.813, 0.84, 0.998, 0.981, 0.83, 0.823, 0.537, 0.969, 0.682, 0.076, 0.284, 0.822, 0.933, 0.856, 0.477, 0.953, 0.598, 0.977, 0.931, 0.74, 0.697, 0.855, 0.586, 0.57, 0.532, 0.991, 0.905, 0.893, 0.759, 0.395, 0.668, 0.905, 0.9, 0.682, 0.84, 0.811, 0.985, 0.54, 0.831, 0.609, 0.355, 0.715, 0.291, 0.691, 0.962, 0.203, 0.892, 0.984, 0.974, 0.936, 0.848, 0.989, 0.913, 0.998, 0.997, 0.405, 0.471, 0.906, 0.979, 0.998, 0.242, 0.428, 0.943, 0.921, 0.883, 0.598, 0.993, 0.665, 0.994, 0.131, 0.41, 0.878, 0.867, 0.54, 0.327, 0.776, 0.555, 0.802, 0.88, 0.919, 0.648, 0.93, 0.757, 0.241, 0.61, 0.978, 0.708, 0.884, 0.925, 0.924, 0.965, 0.854, 0.725, 0.938, 0.898, 0.585, 0.506, 0.97, 0.96, 0.883, 0.79, 0.499, 0.762, 0.65, 0.99, 0.457, 0.718, 0.915, 0.998, 0.968, 0.991, 0.662, 0.971, 0.442, 0.91, 0.484, 0.389, 0.753, 0.927, 0.994, 0.733, 0.883, 0.627, 0.975, 0.78, 0.961, 0.863, 0.554, 0.684, 0.535, 0.439, 0.588, 0.971, 0.99, 0.9, 0.985, 0.835, 0.94, 0.905, 0.848, 0.69, 0.593, 0.538, 0.521, 0.621, 0.41, 0.812, 0.851, 0.987, 0.937, 0.506, 0.661, 0.538, 0.962, 0.978, 0.525, 0.907, 0.713, 0.376, 0.881, 0.611, 0.912, 0.76, 0.991, 0.981, 0.839, 0.928, 0.524, 0.61, 0.816, 0.633, 0.251, 0.969, 0.788, 0.665, 0.211, 0.968, 0.977, 0.931, 0.224, 0.723, 0.947, 0.311, 0.926, 0.959, 0.75, 0.824, 0.862, 0.959, 0.998, 0.921, 0.784, 0.755, 0.652, 0.938, 0.654, 0.665, 0.845, 0.304, 0.953, 0.895, 0.489, 0.758, 0.913, 0.97, 0.641, 0.762, 0.655, 0.51, 0.936, 0.612, 0.658, 0.914, 0.865, 0.558, 0.39, 0.701, 0.847, 0.999, 0.916, 0.987, 0.962, 0.967, 0.969, 0.978, 0.915, 0.922, 0.918, 0.792, 0.892, 0.837, 0.897, 0.793, 0.91, 0.89, 0.806, 0.825, 0.923, 0.772, 0.922, 0.971, 0.866, 0.578, 0.439, 0.471, 0.847, 0.941, 0.504, 0.803, 0.706, 0.879, 0.581, 0.827, 0.099, 0.658, 0.845, 0.868, 0.936, 0.911, 0.889, 0.547, 0.392, 0.764, 0.945, 0.937, 0.67, 0.512, 0.961, 0.595, 0.622, 0.804, 0.827, 0.633, 0.399, 0.985, 0.79, 0.799, 0.848, 0.667, 0.941, 0.381, 0.859, 0.949, 0.919, 0.911, 0.938, 0.975, 0.834, 0.703, 0.975, 0.664, 0.943, 0.989, 0.883, 0.64]
global origin = 1
global destination = 60