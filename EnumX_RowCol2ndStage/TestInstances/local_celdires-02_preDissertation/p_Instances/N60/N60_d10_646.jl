global arcs = [1 32; 1 33; 1 37; 1 43; 1 57; 2 10; 2 18; 2 23; 2 35; 2 36; 2 37; 2 38; 2 43; 3 5; 3 7; 3 19; 3 20; 3 33; 3 48; 4 9; 4 12; 4 21; 4 29; 4 32; 4 53; 4 57; 5 4; 5 9; 5 23; 5 32; 5 53; 5 54; 6 19; 6 23; 6 27; 6 51; 7 9; 7 10; 7 19; 7 24; 7 25; 7 27; 7 28; 7 59; 8 4; 8 24; 8 42; 8 57; 8 58; 8 59; 8 60; 9 6; 9 24; 9 25; 9 27; 9 32; 9 40; 9 41; 9 48; 9 50; 9 56; 9 60; 10 29; 10 57; 11 4; 11 7; 11 25; 11 55; 12 4; 12 11; 12 37; 12 38; 12 54; 12 59; 12 60; 13 3; 13 14; 13 17; 13 32; 13 56; 14 19; 14 22; 14 44; 14 59; 15 6; 15 10; 15 20; 15 33; 15 35; 15 50; 16 6; 16 10; 16 45; 16 49; 16 55; 17 9; 17 18; 17 28; 18 38; 18 48; 18 56; 19 4; 19 5; 19 10; 20 2; 20 3; 20 42; 20 43; 20 50; 20 56; 20 60; 21 12; 21 17; 21 25; 21 36; 21 37; 21 39; 21 51; 22 6; 22 10; 22 25; 22 26; 22 27; 23 14; 23 32; 23 43; 23 46; 23 49; 23 60; 24 11; 24 36; 24 39; 24 54; 25 5; 25 8; 25 10; 25 15; 25 56; 26 6; 26 13; 26 29; 26 38; 26 47; 27 7; 27 20; 27 33; 27 39; 27 59; 28 6; 28 12; 28 35; 28 44; 28 45; 28 57; 29 3; 29 9; 29 10; 29 43; 29 60; 30 39; 31 2; 31 3; 31 6; 31 12; 31 13; 31 16; 31 20; 31 24; 31 34; 31 46; 31 59; 32 8; 32 29; 32 35; 32 37; 32 40; 32 47; 33 6; 33 7; 33 13; 33 15; 33 39; 33 40; 33 44; 34 31; 34 39; 35 33; 35 50; 36 4; 36 8; 36 35; 36 39; 36 40; 36 44; 36 49; 36 50; 37 25; 37 38; 37 43; 38 23; 38 35; 38 36; 38 39; 38 46; 38 55; 39 9; 39 20; 39 37; 39 56; 40 4; 40 12; 40 46; 40 50; 40 57; 41 2; 41 14; 41 34; 41 39; 41 59; 42 19; 42 22; 42 23; 42 34; 42 49; 43 8; 43 11; 43 20; 43 23; 43 24; 43 34; 43 35; 43 41; 43 47; 44 4; 44 45; 45 2; 45 27; 45 35; 45 47; 45 49; 46 5; 46 21; 46 23; 46 24; 46 50; 46 53; 46 57; 47 20; 47 24; 47 29; 47 59; 47 60; 48 20; 48 45; 48 53; 49 15; 49 33; 49 36; 49 43; 49 50; 49 51; 49 58; 50 4; 50 15; 50 25; 50 32; 50 34; 50 41; 50 43; 50 45; 50 59; 51 7; 51 12; 51 14; 51 38; 51 44; 52 11; 52 48; 53 5; 53 6; 53 23; 53 28; 53 30; 53 43; 53 44; 53 49; 53 56; 53 58; 53 60; 54 18; 54 50; 54 60; 55 10; 55 18; 55 22; 55 32; 55 41; 55 42; 55 44; 55 46; 56 6; 56 20; 56 29; 56 33; 56 39; 56 45; 56 55; 56 57; 57 6; 57 17; 57 27; 57 28; 57 31; 57 35; 57 38; 57 50; 57 52; 57 55; 58 16; 58 36; 58 45; 59 7; 59 11; 59 14; 59 17; 59 25; 59 39; 59 52]
global d_x = [7.0, 3.0, 6.0, 4.0, 3.0, 6.0, 5.0, 9.0, 10.0, 4.0, 10.0, 2.0, 5.0, 3.0, 3.0, 6.0, 7.0, 4.0, 4.0, 1.0, 5.0, 10.0, 2.0, 8.0, 5.0, 1.0, 4.0, 9.0, 3.0, 10.0, 6.0, 7.0, 1.0, 3.0, 1.0, 3.0, 9.0, 3.0, 7.0, 2.0, 3.0, 7.0, 9.0, 3.0, 9.0, 4.0, 6.0, 8.0, 6.0, 5.0, 5.0, 1.0, 2.0, 2.0, 7.0, 6.0, 6.0, 9.0, 4.0, 2.0, 4.0, 3.0, 7.0, 10.0, 10.0, 6.0, 2.0, 1.0, 1.0, 9.0, 6.0, 4.0, 7.0, 9.0, 4.0, 7.0, 4.0, 2.0, 8.0, 6.0, 2.0, 7.0, 10.0, 8.0, 3.0, 2.0, 10.0, 10.0, 8.0, 5.0, 3.0, 7.0, 3.0, 9.0, 2.0, 8.0, 3.0, 2.0, 8.0, 10.0, 1.0, 2.0, 7.0, 2.0, 4.0, 4.0, 10.0, 6.0, 2.0, 1.0, 3.0, 8.0, 10.0, 2.0, 7.0, 4.0, 9.0, 1.0, 5.0, 7.0, 9.0, 7.0, 7.0, 10.0, 9.0, 3.0, 2.0, 1.0, 2.0, 8.0, 4.0, 10.0, 2.0, 8.0, 9.0, 7.0, 5.0, 5.0, 2.0, 2.0, 9.0, 10.0, 4.0, 8.0, 2.0, 8.0, 5.0, 6.0, 10.0, 6.0, 10.0, 8.0, 1.0, 3.0, 8.0, 1.0, 1.0, 4.0, 3.0, 6.0, 1.0, 5.0, 7.0, 9.0, 9.0, 5.0, 1.0, 8.0, 3.0, 10.0, 5.0, 6.0, 6.0, 10.0, 5.0, 10.0, 2.0, 7.0, 4.0, 6.0, 2.0, 1.0, 2.0, 8.0, 1.0, 8.0, 6.0, 6.0, 2.0, 10.0, 1.0, 6.0, 2.0, 7.0, 1.0, 5.0, 6.0, 3.0, 1.0, 1.0, 9.0, 4.0, 3.0, 8.0, 3.0, 1.0, 9.0, 4.0, 4.0, 6.0, 2.0, 2.0, 5.0, 10.0, 10.0, 10.0, 7.0, 7.0, 4.0, 9.0, 2.0, 10.0, 5.0, 4.0, 9.0, 1.0, 10.0, 6.0, 2.0, 10.0, 5.0, 10.0, 9.0, 6.0, 9.0, 8.0, 7.0, 3.0, 7.0, 4.0, 3.0, 4.0, 1.0, 10.0, 8.0, 10.0, 8.0, 2.0, 5.0, 9.0, 6.0, 8.0, 3.0, 4.0, 6.0, 8.0, 5.0, 2.0, 9.0, 2.0, 1.0, 1.0, 10.0, 5.0, 5.0, 10.0, 7.0, 5.0, 6.0, 2.0, 10.0, 8.0, 9.0, 3.0, 9.0, 2.0, 4.0, 9.0, 10.0, 7.0, 2.0, 9.0, 9.0, 8.0, 7.0, 6.0, 4.0, 9.0, 6.0, 1.0, 4.0, 8.0, 1.0, 10.0, 3.0, 7.0, 5.0, 8.0, 3.0, 9.0, 10.0, 5.0, 10.0, 3.0, 8.0, 3.0, 3.0, 5.0, 4.0, 5.0, 4.0, 6.0, 3.0, 9.0, 4.0, 5.0, 7.0, 6.0, 2.0, 2.0, 2.0, 1.0, 6.0, 2.0, 5.0, 4.0, 7.0, 2.0]
global b_x = 5
global d_y = [3.0, 1.0, 4.0, 6.0, 7.0, 6.0, 4.0, 8.0, 1.0, 3.0, 9.0, 8.0, 2.0, 6.0, 1.0, 9.0, 5.0, 1.0, 8.0, 1.0, 9.0, 8.0, 7.0, 9.0, 9.0, 1.0, 9.0, 7.0, 8.0, 1.0, 7.0, 5.0, 10.0, 5.0, 3.0, 2.0, 3.0, 7.0, 3.0, 8.0, 4.0, 10.0, 3.0, 5.0, 4.0, 9.0, 10.0, 10.0, 8.0, 5.0, 10.0, 9.0, 8.0, 4.0, 9.0, 6.0, 1.0, 10.0, 6.0, 4.0, 9.0, 8.0, 3.0, 4.0, 10.0, 5.0, 2.0, 10.0, 7.0, 6.0, 7.0, 3.0, 10.0, 7.0, 8.0, 10.0, 8.0, 10.0, 3.0, 1.0, 2.0, 3.0, 5.0, 9.0, 6.0, 2.0, 3.0, 5.0, 10.0, 9.0, 7.0, 4.0, 7.0, 8.0, 3.0, 5.0, 9.0, 6.0, 10.0, 7.0, 8.0, 7.0, 3.0, 5.0, 9.0, 6.0, 10.0, 4.0, 10.0, 10.0, 9.0, 10.0, 4.0, 6.0, 9.0, 4.0, 10.0, 5.0, 5.0, 2.0, 10.0, 8.0, 8.0, 3.0, 4.0, 8.0, 4.0, 2.0, 4.0, 1.0, 5.0, 1.0, 1.0, 8.0, 7.0, 8.0, 8.0, 4.0, 7.0, 3.0, 4.0, 5.0, 6.0, 5.0, 5.0, 9.0, 2.0, 6.0, 4.0, 4.0, 4.0, 2.0, 3.0, 9.0, 9.0, 8.0, 4.0, 2.0, 8.0, 8.0, 4.0, 9.0, 6.0, 6.0, 5.0, 8.0, 3.0, 6.0, 8.0, 9.0, 9.0, 9.0, 7.0, 1.0, 7.0, 3.0, 8.0, 9.0, 5.0, 1.0, 8.0, 3.0, 7.0, 6.0, 1.0, 5.0, 3.0, 8.0, 7.0, 10.0, 2.0, 1.0, 4.0, 10.0, 8.0, 3.0, 8.0, 5.0, 5.0, 8.0, 5.0, 10.0, 10.0, 7.0, 2.0, 8.0, 7.0, 3.0, 8.0, 9.0, 1.0, 1.0, 6.0, 7.0, 3.0, 8.0, 3.0, 1.0, 7.0, 8.0, 3.0, 5.0, 8.0, 2.0, 3.0, 5.0, 10.0, 1.0, 8.0, 9.0, 2.0, 8.0, 5.0, 5.0, 5.0, 1.0, 10.0, 5.0, 7.0, 5.0, 7.0, 4.0, 1.0, 1.0, 9.0, 6.0, 9.0, 2.0, 3.0, 8.0, 2.0, 6.0, 1.0, 9.0, 8.0, 2.0, 5.0, 10.0, 9.0, 5.0, 2.0, 3.0, 3.0, 6.0, 9.0, 8.0, 9.0, 6.0, 2.0, 6.0, 2.0, 8.0, 9.0, 1.0, 9.0, 1.0, 9.0, 10.0, 9.0, 3.0, 8.0, 1.0, 4.0, 8.0, 5.0, 2.0, 8.0, 10.0, 5.0, 7.0, 5.0, 3.0, 2.0, 5.0, 3.0, 2.0, 8.0, 3.0, 5.0, 10.0, 6.0, 5.0, 8.0, 7.0, 6.0, 3.0, 7.0, 9.0, 10.0, 6.0, 3.0, 10.0, 6.0, 2.0, 2.0, 9.0, 8.0, 7.0, 7.0, 1.0, 10.0, 5.0, 2.0, 9.0, 5.0, 2.0, 1.0, 3.0]
global b_y = 10
global p = [0.59, 0.331, 0.172, 0.019, 0.764, 0.064, 0.699, 0.503, 0.252, 0.475, 0.341, 0.242, 0.982, 0.008, 0.091, 0.145, 0.82, 0.91, 0.905, 0.393, 0.268, 0.311, 0.061, 0.378, 0.361, 0.173, 0.372, 0.025, 0.577, 0.037, 0.875, 0.474, 0.55, 0.257, 0.273, 0.095, 0.389, 0.994, 0.317, 0.36, 0.346, 0.917, 0.276, 0.359, 0.139, 0.017, 0.851, 0.315, 0.974, 0.852, 0.483, 0.212, 0.584, 0.076, 0.784, 0.557, 0.217, 0.14, 0.926, 0.969, 0.598, 0.84, 0.778, 0.407, 0.268, 0.287, 0.6, 0.102, 0.314, 0.108, 0.444, 0.993, 0.658, 0.326, 0.438, 0.561, 0.454, 0.507, 0.773, 0.846, 0.865, 0.23, 0.27, 0.299, 0.29, 0.953, 0.195, 0.794, 0.266, 0.786, 0.445, 0.103, 0.266, 0.798, 0.574, 0.919, 0.807, 0.523, 0.284, 0.778, 0.063, 0.687, 0.713, 0.494, 0.868, 0.711, 0.496, 0.84, 0.126, 0.684, 0.627, 0.102, 0.293, 0.045, 0.328, 0.952, 0.36, 0.869, 0.783, 0.338, 0.182, 0.983, 0.85, 0.412, 0.391, 0.597, 0.661, 0.055, 0.431, 0.436, 0.94, 0.711, 0.307, 0.49, 0.918, 0.631, 0.468, 0.869, 0.514, 0.937, 0.459, 0.755, 0.883, 0.258, 0.462, 0.837, 0.065, 0.324, 0.346, 0.596, 0.372, 0.041, 0.073, 0.333, 0.394, 0.716, 0.127, 0.299, 0.811, 0.972, 0.908, 0.141, 0.481, 0.32, 0.994, 0.392, 0.491, 0.151, 0.465, 0.51, 0.454, 0.482, 0.607, 0.649, 0.878, 0.048, 0.633, 0.201, 0.097, 0.438, 0.602, 0.025, 0.481, 0.915, 0.531, 0.308, 0.997, 0.415, 0.799, 0.104, 0.166, 0.974, 0.934, 0.113, 0.505, 0.106, 0.719, 0.974, 0.928, 0.141, 0.943, 0.895, 0.373, 0.102, 0.524, 0.123, 0.127, 0.571, 0.367, 0.059, 0.947, 0.613, 0.789, 0.164, 0.573, 0.145, 0.484, 0.701, 0.569, 0.647, 0.321, 0.465, 0.286, 0.167, 0.036, 0.023, 0.028, 0.034, 0.104, 0.922, 0.176, 0.913, 0.992, 0.942, 0.639, 0.558, 0.091, 0.241, 0.094, 0.545, 0.872, 0.854, 0.703, 0.687, 0.928, 0.056, 0.219, 0.867, 0.165, 0.684, 0.086, 0.428, 0.435, 0.846, 0.091, 0.168, 0.09, 0.05, 0.172, 0.123, 0.413, 0.191, 0.796, 0.027, 0.929, 0.044, 0.324, 0.917, 0.78, 0.39, 0.868, 0.279, 0.372, 0.051, 0.55, 0.548, 0.893, 0.768, 0.467, 0.38, 0.711, 0.019, 0.143, 0.959, 0.746, 0.83, 0.989, 0.436, 0.634, 0.617, 0.546, 0.404, 0.151, 0.02, 0.091, 0.876, 0.477, 0.751, 0.01, 0.734, 0.609, 0.263, 0.247, 0.296, 0.381, 0.26, 0.042, 0.738, 0.647, 0.13, 0.442, 0.444, 0.15, 0.497, 0.675, 0.552, 0.514, 0.8, 0.671, 0.516, 0.151, 0.515, 0.043, 0.866, 0.235, 0.061, 0.798, 0.761]
global q = [0.819, 0.364, 0.545, 0.994, 0.994, 0.585, 0.701, 0.928, 0.615, 0.712, 0.622, 0.721, 0.995, 0.293, 0.916, 0.621, 0.966, 0.935, 0.917, 0.555, 0.551, 0.758, 0.689, 0.793, 0.877, 0.943, 0.988, 0.041, 0.808, 0.756, 0.898, 0.744, 0.897, 0.715, 0.574, 0.714, 0.55, 0.996, 0.708, 0.594, 0.686, 0.984, 0.646, 0.917, 0.314, 0.596, 0.965, 0.483, 0.981, 0.93, 0.827, 0.917, 0.763, 0.775, 0.937, 0.847, 0.726, 0.413, 0.959, 0.999, 0.639, 0.965, 0.872, 0.69, 0.343, 0.91, 0.983, 0.576, 0.983, 0.272, 0.617, 0.993, 0.979, 0.752, 0.553, 0.793, 0.972, 0.984, 0.961, 0.985, 0.954, 0.84, 0.899, 0.301, 0.921, 0.982, 0.41, 0.943, 0.942, 0.902, 0.853, 0.477, 0.547, 0.895, 0.876, 0.975, 0.916, 0.951, 0.84, 0.957, 0.15, 0.95, 0.816, 0.957, 0.982, 0.818, 0.632, 0.898, 0.64, 0.887, 0.838, 0.501, 0.659, 0.367, 0.582, 0.97, 0.719, 0.898, 0.837, 0.392, 0.386, 0.991, 0.95, 0.582, 0.601, 0.861, 0.802, 0.673, 0.529, 0.854, 0.964, 0.834, 0.577, 0.683, 0.955, 0.881, 0.89, 0.895, 0.847, 0.988, 0.526, 0.875, 0.972, 0.396, 0.734, 0.899, 0.587, 0.854, 0.848, 0.785, 0.609, 0.457, 0.081, 0.858, 0.445, 0.85, 0.34, 0.418, 0.839, 0.972, 0.989, 0.47, 0.996, 0.34, 0.999, 0.628, 0.819, 0.982, 0.547, 0.7, 0.974, 0.633, 0.88, 0.794, 0.983, 0.601, 0.696, 0.744, 0.681, 0.948, 0.777, 0.123, 0.681, 0.952, 0.824, 0.391, 0.999, 0.746, 0.842, 0.261, 0.61, 0.981, 0.987, 0.77, 0.77, 0.777, 0.733, 0.989, 0.938, 0.457, 0.972, 0.984, 0.845, 0.366, 0.786, 0.613, 0.42, 0.708, 0.804, 0.551, 0.972, 0.862, 0.928, 0.251, 0.962, 0.261, 0.509, 0.74, 0.836, 0.701, 0.345, 0.906, 0.715, 0.446, 0.598, 0.27, 0.622, 0.113, 0.523, 0.972, 0.614, 0.916, 0.997, 0.988, 0.829, 0.883, 0.881, 0.285, 0.85, 0.85, 0.914, 0.944, 0.803, 0.86, 0.993, 0.172, 0.945, 0.94, 0.351, 0.829, 0.724, 0.78, 0.775, 0.869, 0.122, 0.883, 0.374, 0.602, 0.201, 0.567, 0.751, 0.339, 0.878, 0.838, 0.935, 0.632, 0.411, 0.979, 0.852, 0.652, 0.936, 0.54, 0.634, 0.148, 0.867, 0.662, 0.961, 0.833, 0.739, 0.813, 0.9, 0.755, 0.235, 0.988, 0.851, 0.833, 0.992, 0.942, 0.856, 0.756, 0.703, 0.986, 0.617, 0.559, 0.455, 0.918, 0.708, 0.941, 0.278, 0.794, 0.62, 0.287, 0.45, 0.483, 0.921, 0.59, 0.304, 0.89, 0.709, 0.409, 0.477, 0.726, 0.26, 0.615, 0.675, 0.838, 0.869, 0.858, 0.916, 0.692, 0.875, 0.574, 0.544, 0.928, 0.251, 0.179, 0.99, 0.827]
global origin = 1
global destination = 60