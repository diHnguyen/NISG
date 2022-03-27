global arcs = [1 12; 1 38; 1 40; 1 54; 2 6; 2 10; 2 12; 2 15; 2 36; 2 53; 2 54; 3 28; 3 37; 3 38; 3 47; 3 50; 4 12; 4 14; 4 36; 4 54; 5 3; 5 18; 5 33; 5 40; 5 43; 6 17; 6 42; 6 50; 6 52; 6 56; 7 6; 7 9; 7 16; 7 18; 7 22; 7 39; 7 44; 7 48; 7 50; 8 2; 8 4; 8 15; 8 36; 8 37; 8 51; 9 20; 9 32; 9 35; 9 49; 9 50; 9 54; 10 3; 10 16; 10 35; 10 42; 11 6; 11 15; 11 21; 11 31; 11 32; 11 39; 11 44; 11 57; 12 25; 12 29; 12 33; 12 37; 12 51; 13 3; 13 22; 13 37; 13 41; 13 46; 13 56; 14 20; 14 29; 14 53; 15 8; 15 24; 15 32; 15 35; 15 47; 15 56; 16 2; 16 5; 16 31; 16 34; 16 39; 16 42; 16 52; 16 59; 17 2; 17 3; 17 10; 17 32; 17 33; 17 43; 18 2; 18 55; 19 24; 19 45; 19 54; 20 5; 20 27; 20 29; 20 39; 21 24; 21 25; 21 43; 21 49; 21 52; 22 2; 22 10; 22 15; 22 17; 22 18; 22 23; 22 47; 23 4; 23 20; 23 21; 23 26; 23 33; 23 37; 23 44; 23 56; 24 4; 24 5; 24 22; 24 29; 24 32; 24 35; 24 40; 24 49; 24 52; 24 55; 25 6; 25 20; 25 32; 25 50; 26 6; 26 48; 26 52; 26 57; 26 59; 27 18; 27 37; 27 51; 28 13; 28 26; 28 43; 28 48; 29 21; 29 37; 29 54; 30 8; 30 10; 30 33; 30 45; 30 53; 30 59; 31 2; 31 29; 31 53; 32 5; 32 11; 32 12; 32 16; 32 17; 32 18; 32 25; 32 44; 32 51; 33 2; 33 4; 33 44; 34 4; 34 20; 34 22; 34 37; 34 47; 34 57; 35 6; 35 18; 35 24; 35 29; 35 54; 36 3; 36 7; 36 10; 36 35; 36 47; 36 49; 36 52; 36 58; 37 14; 37 23; 37 54; 38 21; 38 32; 38 33; 38 37; 38 55; 39 15; 39 24; 39 25; 39 41; 40 4; 40 11; 40 29; 40 39; 40 49; 40 53; 41 4; 41 12; 41 16; 41 17; 41 19; 41 20; 41 45; 41 59; 42 3; 42 9; 42 15; 42 21; 42 36; 42 46; 42 47; 43 13; 43 19; 43 46; 44 4; 44 5; 44 42; 44 49; 44 58; 45 24; 45 27; 45 29; 45 52; 46 15; 46 18; 46 43; 47 24; 47 33; 47 36; 47 37; 47 53; 48 2; 48 3; 48 25; 48 39; 48 52; 48 56; 49 10; 49 31; 49 47; 49 48; 50 5; 50 6; 50 14; 50 22; 50 31; 50 39; 50 46; 51 9; 51 41; 51 53; 51 54; 52 8; 52 11; 52 13; 52 21; 52 28; 52 33; 52 46; 52 47; 52 54; 52 57; 52 59; 53 14; 53 21; 53 30; 53 34; 53 52; 53 60; 54 3; 54 5; 54 11; 54 17; 54 26; 55 6; 55 7; 55 8; 55 20; 55 33; 55 41; 55 46; 56 26; 56 27; 56 28; 56 29; 56 45; 56 54; 57 7; 57 8; 57 21; 57 43; 57 49; 57 55; 57 56; 57 59; 58 9; 58 16; 58 27; 58 44; 58 47; 58 56; 59 8; 59 37; 59 51]
global d_x = [6.0, 1.0, 1.0, 6.0, 3.0, 5.0, 4.0, 7.0, 4.0, 2.0, 3.0, 5.0, 10.0, 1.0, 4.0, 6.0, 6.0, 8.0, 1.0, 6.0, 9.0, 10.0, 3.0, 5.0, 10.0, 4.0, 10.0, 2.0, 3.0, 10.0, 7.0, 4.0, 6.0, 1.0, 4.0, 10.0, 9.0, 7.0, 10.0, 2.0, 10.0, 9.0, 6.0, 8.0, 10.0, 9.0, 2.0, 5.0, 9.0, 8.0, 9.0, 3.0, 4.0, 6.0, 6.0, 9.0, 6.0, 3.0, 2.0, 1.0, 10.0, 4.0, 10.0, 7.0, 8.0, 2.0, 1.0, 1.0, 8.0, 7.0, 7.0, 7.0, 3.0, 9.0, 7.0, 8.0, 5.0, 2.0, 8.0, 4.0, 6.0, 8.0, 1.0, 4.0, 3.0, 3.0, 9.0, 1.0, 9.0, 10.0, 3.0, 10.0, 10.0, 3.0, 3.0, 10.0, 2.0, 5.0, 3.0, 3.0, 9.0, 8.0, 10.0, 7.0, 9.0, 10.0, 9.0, 6.0, 10.0, 9.0, 3.0, 4.0, 5.0, 6.0, 3.0, 9.0, 6.0, 2.0, 3.0, 9.0, 2.0, 1.0, 3.0, 6.0, 8.0, 5.0, 1.0, 8.0, 10.0, 6.0, 6.0, 3.0, 7.0, 4.0, 10.0, 10.0, 1.0, 4.0, 8.0, 4.0, 6.0, 6.0, 6.0, 7.0, 4.0, 7.0, 5.0, 6.0, 8.0, 4.0, 5.0, 6.0, 9.0, 8.0, 3.0, 2.0, 5.0, 8.0, 3.0, 10.0, 7.0, 9.0, 10.0, 4.0, 2.0, 1.0, 9.0, 9.0, 3.0, 1.0, 7.0, 4.0, 10.0, 6.0, 5.0, 10.0, 2.0, 6.0, 2.0, 6.0, 1.0, 8.0, 4.0, 7.0, 4.0, 6.0, 2.0, 3.0, 9.0, 9.0, 6.0, 7.0, 9.0, 10.0, 8.0, 2.0, 7.0, 7.0, 3.0, 1.0, 7.0, 2.0, 3.0, 9.0, 9.0, 7.0, 3.0, 10.0, 4.0, 8.0, 1.0, 9.0, 2.0, 3.0, 4.0, 4.0, 7.0, 7.0, 4.0, 7.0, 1.0, 10.0, 3.0, 10.0, 5.0, 3.0, 9.0, 2.0, 1.0, 4.0, 8.0, 5.0, 2.0, 6.0, 10.0, 5.0, 2.0, 9.0, 8.0, 5.0, 5.0, 10.0, 7.0, 9.0, 4.0, 6.0, 2.0, 7.0, 5.0, 9.0, 3.0, 5.0, 2.0, 9.0, 6.0, 2.0, 6.0, 10.0, 6.0, 1.0, 7.0, 6.0, 9.0, 10.0, 6.0, 2.0, 9.0, 10.0, 2.0, 5.0, 5.0, 5.0, 6.0, 8.0, 1.0, 10.0, 6.0, 9.0, 4.0, 8.0, 2.0, 8.0, 10.0, 5.0, 1.0, 6.0, 4.0, 4.0, 8.0, 9.0, 7.0, 10.0, 6.0, 9.0, 7.0, 1.0, 2.0, 5.0, 5.0, 2.0, 1.0, 10.0, 3.0, 9.0, 2.0, 7.0, 5.0, 7.0, 8.0, 8.0, 4.0, 2.0, 9.0, 5.0, 8.0, 3.0, 5.0, 4.0, 7.0, 8.0, 10.0]
global b_x = 5
global d_y = [5.0, 9.0, 1.0, 2.0, 7.0, 4.0, 5.0, 6.0, 9.0, 10.0, 5.0, 9.0, 5.0, 5.0, 4.0, 1.0, 1.0, 2.0, 3.0, 7.0, 8.0, 4.0, 6.0, 9.0, 6.0, 3.0, 1.0, 1.0, 6.0, 3.0, 8.0, 5.0, 10.0, 7.0, 4.0, 8.0, 10.0, 3.0, 3.0, 4.0, 8.0, 7.0, 5.0, 10.0, 3.0, 10.0, 10.0, 1.0, 6.0, 9.0, 6.0, 6.0, 3.0, 3.0, 5.0, 3.0, 3.0, 7.0, 9.0, 9.0, 9.0, 3.0, 1.0, 8.0, 8.0, 1.0, 9.0, 4.0, 6.0, 10.0, 6.0, 1.0, 6.0, 6.0, 6.0, 1.0, 10.0, 9.0, 2.0, 6.0, 9.0, 7.0, 8.0, 7.0, 7.0, 6.0, 3.0, 9.0, 4.0, 10.0, 8.0, 10.0, 7.0, 9.0, 10.0, 2.0, 9.0, 4.0, 7.0, 8.0, 4.0, 3.0, 1.0, 8.0, 6.0, 1.0, 2.0, 1.0, 2.0, 8.0, 1.0, 4.0, 10.0, 6.0, 10.0, 10.0, 3.0, 6.0, 1.0, 9.0, 9.0, 2.0, 10.0, 1.0, 7.0, 6.0, 10.0, 2.0, 8.0, 4.0, 2.0, 7.0, 8.0, 6.0, 2.0, 4.0, 3.0, 7.0, 6.0, 7.0, 3.0, 1.0, 4.0, 6.0, 1.0, 7.0, 8.0, 6.0, 8.0, 2.0, 9.0, 2.0, 4.0, 8.0, 9.0, 7.0, 1.0, 4.0, 10.0, 10.0, 3.0, 9.0, 10.0, 2.0, 7.0, 5.0, 5.0, 4.0, 6.0, 3.0, 6.0, 9.0, 9.0, 9.0, 1.0, 7.0, 8.0, 3.0, 1.0, 4.0, 5.0, 10.0, 8.0, 9.0, 8.0, 6.0, 2.0, 10.0, 1.0, 9.0, 3.0, 1.0, 3.0, 7.0, 5.0, 4.0, 1.0, 7.0, 7.0, 8.0, 6.0, 7.0, 1.0, 10.0, 10.0, 6.0, 4.0, 10.0, 9.0, 2.0, 1.0, 4.0, 2.0, 3.0, 10.0, 1.0, 8.0, 9.0, 10.0, 3.0, 6.0, 5.0, 4.0, 6.0, 5.0, 4.0, 5.0, 3.0, 9.0, 5.0, 4.0, 2.0, 8.0, 2.0, 7.0, 1.0, 3.0, 7.0, 4.0, 10.0, 3.0, 10.0, 9.0, 3.0, 8.0, 6.0, 5.0, 2.0, 4.0, 10.0, 4.0, 6.0, 8.0, 1.0, 2.0, 9.0, 1.0, 1.0, 8.0, 6.0, 6.0, 2.0, 3.0, 8.0, 10.0, 9.0, 2.0, 6.0, 7.0, 7.0, 8.0, 3.0, 10.0, 2.0, 6.0, 4.0, 8.0, 6.0, 10.0, 7.0, 2.0, 3.0, 10.0, 2.0, 5.0, 3.0, 6.0, 7.0, 7.0, 10.0, 3.0, 10.0, 2.0, 10.0, 6.0, 10.0, 9.0, 1.0, 7.0, 3.0, 7.0, 10.0, 2.0, 8.0, 4.0, 7.0, 8.0, 4.0, 8.0, 5.0, 8.0, 7.0, 10.0, 6.0, 3.0, 1.0, 2.0, 5.0, 3.0, 6.0, 7.0]
global b_y = 10
global p = [0.916, 0.342, 0.406, 0.169, 0.156, 0.241, 0.389, 0.087, 0.244, 0.158, 0.611, 0.334, 0.345, 0.625, 0.374, 0.511, 0.761, 0.901, 0.67, 0.405, 0.357, 0.015, 0.316, 0.522, 0.003, 0.737, 0.054, 0.147, 0.242, 0.471, 0.203, 0.246, 0.311, 0.626, 0.992, 0.257, 0.677, 0.195, 0.153, 0.039, 0.117, 0.289, 0.159, 0.442, 0.264, 0.324, 0.591, 0.708, 0.725, 0.688, 0.061, 0.836, 0.228, 0.546, 0.949, 0.964, 0.259, 0.975, 0.048, 0.984, 0.384, 0.546, 0.81, 0.737, 0.282, 0.745, 0.693, 0.647, 0.023, 0.719, 0.27, 0.587, 0.594, 0.862, 0.757, 0.833, 0.384, 0.165, 0.458, 0.4, 0.67, 0.556, 0.861, 0.958, 0.846, 0.5, 0.42, 0.576, 0.353, 0.362, 0.625, 0.367, 0.534, 0.34, 0.744, 0.18, 0.575, 0.343, 0.306, 0.844, 0.117, 0.662, 0.692, 0.401, 0.781, 0.763, 0.595, 0.044, 0.608, 0.923, 0.107, 0.952, 0.599, 0.905, 0.699, 0.891, 0.938, 0.007, 0.173, 0.821, 0.2, 0.799, 0.176, 0.001, 0.337, 0.63, 0.86, 0.89, 0.592, 0.423, 0.97, 0.082, 0.085, 0.817, 0.141, 0.576, 0.079, 0.067, 0.014, 0.972, 0.237, 0.349, 0.59, 0.131, 0.467, 0.901, 0.529, 0.259, 0.633, 0.552, 0.602, 0.59, 0.844, 0.825, 0.514, 0.489, 0.246, 0.419, 0.111, 0.454, 0.162, 0.197, 0.056, 0.226, 0.671, 0.596, 0.415, 0.358, 0.767, 0.664, 0.122, 0.136, 0.101, 0.23, 0.906, 0.418, 0.217, 0.246, 0.143, 0.642, 0.801, 0.554, 0.758, 0.032, 0.006, 0.792, 0.257, 0.516, 0.223, 0.315, 0.895, 0.73, 0.5, 0.587, 0.269, 0.74, 0.639, 0.431, 0.927, 0.399, 0.458, 0.923, 0.39, 0.78, 0.961, 0.519, 0.772, 0.878, 0.308, 0.055, 0.416, 0.365, 0.896, 0.631, 0.963, 0.927, 0.245, 0.858, 0.855, 0.787, 0.108, 0.078, 0.371, 0.309, 0.846, 0.232, 0.229, 0.72, 0.922, 0.068, 0.55, 0.89, 0.146, 0.476, 0.236, 0.269, 0.04, 0.047, 0.683, 0.071, 0.926, 0.218, 0.059, 0.262, 0.977, 0.836, 0.192, 0.252, 0.823, 0.154, 0.809, 0.596, 0.491, 0.637, 0.396, 0.19, 0.374, 0.172, 0.749, 0.656, 0.794, 0.335, 0.793, 0.637, 0.227, 0.562, 0.762, 0.724, 0.039, 0.163, 0.386, 0.998, 0.816, 0.059, 0.79, 0.055, 0.782, 0.577, 0.952, 0.863, 0.378, 0.63, 0.819, 0.65, 0.167, 0.495, 0.088, 0.475, 0.515, 0.822, 0.929, 0.259, 0.496, 0.56, 0.059, 0.103, 0.931, 0.953, 0.917, 0.208, 0.768, 0.749, 0.662, 0.705, 0.336, 0.451, 0.921, 0.004, 0.928, 0.057, 0.208, 0.626, 0.896, 0.131, 0.954, 0.076, 0.438, 0.964, 0.656, 0.242, 0.357]
global q = [0.999, 0.364, 0.875, 0.797, 0.549, 0.427, 0.796, 0.673, 0.417, 0.229, 0.991, 0.616, 0.968, 0.659, 0.951, 0.695, 0.963, 0.979, 0.992, 0.837, 0.935, 0.828, 0.339, 0.881, 0.398, 0.994, 0.086, 0.215, 0.976, 0.746, 0.462, 0.647, 0.808, 0.648, 0.993, 0.891, 0.953, 0.85, 0.422, 0.587, 0.672, 0.884, 0.491, 0.49, 0.356, 0.518, 0.868, 0.767, 0.995, 0.934, 0.654, 0.999, 0.372, 0.553, 0.993, 0.976, 0.752, 0.98, 0.813, 0.992, 0.859, 0.915, 0.884, 0.997, 0.769, 0.762, 0.979, 0.9, 0.355, 0.863, 0.389, 0.642, 0.614, 0.897, 0.808, 0.901, 0.583, 0.359, 0.842, 0.457, 0.786, 0.975, 0.913, 0.979, 0.884, 0.873, 0.984, 0.661, 0.896, 0.407, 0.813, 0.805, 0.81, 0.884, 0.966, 0.244, 0.755, 0.432, 0.913, 0.979, 0.235, 0.861, 0.849, 0.485, 0.993, 0.932, 0.689, 0.553, 0.631, 0.956, 0.494, 0.986, 0.793, 0.967, 0.946, 0.909, 0.953, 0.487, 0.486, 0.859, 0.263, 0.909, 0.966, 0.502, 0.6, 0.671, 0.876, 0.988, 0.679, 0.861, 0.999, 0.116, 0.901, 0.963, 0.313, 0.917, 0.152, 0.619, 0.476, 0.988, 0.654, 0.555, 0.77, 0.292, 0.694, 0.935, 0.956, 0.873, 0.96, 0.732, 0.938, 0.882, 0.978, 0.994, 0.909, 0.542, 0.898, 0.796, 0.615, 0.799, 0.931, 0.392, 0.635, 0.953, 0.754, 0.732, 0.632, 0.84, 0.977, 0.972, 0.856, 0.986, 0.716, 0.335, 0.984, 0.672, 0.662, 0.529, 0.21, 0.679, 0.942, 0.64, 0.913, 0.577, 0.157, 0.801, 0.408, 0.982, 0.53, 0.577, 0.98, 0.842, 0.807, 0.982, 0.936, 0.958, 0.975, 0.679, 0.984, 0.425, 0.631, 0.945, 0.587, 0.821, 0.994, 0.62, 0.875, 0.978, 0.723, 0.143, 0.544, 0.551, 0.952, 0.936, 0.978, 0.97, 0.32, 0.951, 0.963, 0.801, 0.772, 0.997, 0.5, 0.421, 0.9, 0.976, 0.282, 0.98, 0.96, 0.769, 0.631, 0.938, 0.797, 0.945, 0.507, 0.403, 0.096, 0.782, 0.965, 0.325, 0.942, 0.52, 0.68, 0.419, 0.998, 0.876, 0.437, 0.58, 0.854, 0.254, 0.919, 0.869, 0.643, 0.818, 0.933, 0.888, 0.51, 0.498, 0.821, 0.99, 0.96, 0.772, 0.952, 0.784, 0.542, 0.661, 0.826, 0.946, 0.68, 0.874, 0.624, 0.999, 0.918, 0.994, 0.872, 0.863, 0.835, 0.93, 0.988, 0.933, 0.633, 0.697, 0.955, 0.786, 0.429, 0.899, 0.239, 0.896, 0.719, 0.949, 0.996, 0.958, 0.766, 0.654, 0.6, 0.726, 0.958, 0.965, 0.93, 0.454, 0.898, 0.881, 0.775, 0.896, 0.957, 0.914, 0.997, 0.849, 0.961, 0.434, 0.413, 0.986, 0.901, 0.243, 0.955, 0.165, 0.44, 0.975, 0.884, 0.815, 0.391]
global origin = 1
global destination = 60