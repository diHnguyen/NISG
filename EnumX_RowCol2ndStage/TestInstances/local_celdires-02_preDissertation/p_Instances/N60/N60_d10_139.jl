global arcs = [1 4; 1 17; 1 33; 1 38; 1 48; 2 21; 2 23; 2 27; 2 41; 2 44; 2 49; 3 8; 3 21; 3 28; 4 9; 4 14; 4 30; 4 36; 4 38; 4 52; 5 8; 5 34; 5 47; 5 50; 5 53; 6 2; 6 8; 6 10; 6 22; 6 43; 7 18; 7 32; 7 35; 7 36; 7 46; 7 50; 8 15; 8 19; 8 25; 8 29; 8 39; 8 42; 8 45; 9 22; 9 24; 9 30; 9 47; 9 51; 10 7; 10 12; 10 28; 10 31; 10 33; 10 38; 10 40; 10 59; 10 60; 11 2; 11 13; 11 18; 11 22; 11 56; 12 26; 12 31; 12 33; 12 41; 12 50; 13 25; 13 52; 13 55; 14 4; 14 17; 14 37; 14 49; 14 55; 14 58; 14 59; 15 14; 15 22; 15 26; 15 39; 16 15; 16 18; 16 21; 16 25; 16 26; 16 29; 16 39; 17 5; 17 6; 17 15; 17 44; 17 51; 17 53; 17 57; 18 10; 18 16; 18 21; 18 33; 18 51; 19 4; 19 47; 19 54; 20 7; 20 22; 20 24; 20 36; 20 37; 20 41; 20 47; 20 48; 21 12; 21 14; 21 15; 21 35; 21 38; 21 46; 21 53; 22 14; 22 40; 22 47; 22 51; 23 60; 24 6; 24 19; 24 25; 24 26; 24 29; 24 48; 24 60; 25 3; 25 4; 25 8; 25 24; 25 29; 25 39; 25 50; 25 55; 25 56; 26 7; 26 13; 26 15; 26 16; 26 32; 26 44; 26 49; 26 59; 27 6; 27 17; 27 23; 27 29; 27 42; 27 43; 27 52; 28 9; 28 36; 28 37; 28 54; 29 12; 29 20; 29 56; 30 23; 30 28; 30 36; 30 37; 30 43; 30 50; 31 2; 31 12; 31 25; 31 36; 31 58; 32 19; 32 31; 32 33; 32 43; 32 47; 32 53; 33 4; 33 8; 33 15; 33 25; 33 40; 33 48; 34 3; 34 11; 34 12; 34 15; 34 18; 34 35; 34 53; 35 3; 35 40; 35 51; 35 55; 36 10; 36 11; 36 12; 36 16; 36 20; 36 22; 36 27; 37 6; 37 30; 37 42; 37 57; 38 7; 38 24; 38 40; 38 41; 39 18; 39 24; 39 28; 39 29; 39 36; 39 47; 39 53; 40 18; 40 22; 40 28; 40 42; 40 43; 41 3; 41 7; 41 8; 41 16; 41 28; 41 34; 41 36; 42 3; 42 16; 42 29; 42 58; 42 59; 43 15; 43 21; 43 32; 43 52; 43 56; 43 58; 44 12; 44 22; 44 23; 44 31; 44 42; 44 50; 45 12; 45 31; 45 38; 45 52; 46 2; 46 8; 46 9; 46 22; 46 24; 46 39; 46 42; 46 49; 46 55; 47 3; 47 25; 47 28; 47 40; 47 50; 48 12; 48 22; 48 32; 48 43; 48 56; 49 18; 49 19; 49 23; 49 29; 49 43; 49 51; 50 22; 50 59; 51 14; 51 37; 51 60; 52 8; 52 11; 52 13; 52 18; 52 31; 52 48; 52 56; 53 3; 53 8; 53 17; 53 34; 53 47; 54 3; 54 20; 54 28; 54 31; 54 32; 54 55; 55 22; 55 50; 55 53; 56 3; 56 7; 56 9; 56 10; 56 20; 56 23; 56 24; 56 48; 56 52; 56 55; 57 3; 57 4; 57 24; 57 28; 57 33; 57 37; 57 47; 57 54; 58 3; 58 7; 58 25; 58 26; 58 31; 58 38; 58 42; 59 5; 59 6; 59 8; 59 13; 59 26; 59 30; 59 36]
global d_x = [3.0, 2.0, 10.0, 10.0, 1.0, 3.0, 8.0, 4.0, 1.0, 7.0, 7.0, 10.0, 10.0, 2.0, 1.0, 6.0, 5.0, 5.0, 2.0, 6.0, 1.0, 2.0, 2.0, 9.0, 3.0, 6.0, 6.0, 8.0, 6.0, 4.0, 4.0, 2.0, 3.0, 1.0, 8.0, 6.0, 6.0, 4.0, 1.0, 5.0, 8.0, 2.0, 5.0, 4.0, 7.0, 6.0, 3.0, 7.0, 2.0, 8.0, 7.0, 3.0, 6.0, 1.0, 9.0, 1.0, 7.0, 2.0, 9.0, 3.0, 1.0, 10.0, 8.0, 6.0, 8.0, 3.0, 6.0, 9.0, 3.0, 6.0, 4.0, 1.0, 7.0, 8.0, 1.0, 1.0, 5.0, 10.0, 7.0, 8.0, 7.0, 3.0, 7.0, 6.0, 8.0, 5.0, 8.0, 5.0, 3.0, 4.0, 9.0, 9.0, 4.0, 10.0, 8.0, 1.0, 8.0, 2.0, 3.0, 1.0, 6.0, 3.0, 7.0, 3.0, 8.0, 5.0, 1.0, 9.0, 6.0, 2.0, 8.0, 6.0, 8.0, 9.0, 1.0, 2.0, 9.0, 2.0, 4.0, 1.0, 6.0, 3.0, 7.0, 10.0, 3.0, 1.0, 10.0, 1.0, 1.0, 4.0, 1.0, 6.0, 5.0, 3.0, 10.0, 7.0, 10.0, 7.0, 2.0, 8.0, 5.0, 9.0, 3.0, 5.0, 6.0, 8.0, 4.0, 3.0, 9.0, 1.0, 6.0, 9.0, 7.0, 7.0, 2.0, 4.0, 3.0, 4.0, 10.0, 5.0, 6.0, 9.0, 6.0, 3.0, 7.0, 6.0, 9.0, 5.0, 9.0, 4.0, 6.0, 1.0, 2.0, 4.0, 5.0, 10.0, 3.0, 9.0, 4.0, 2.0, 5.0, 3.0, 1.0, 9.0, 10.0, 10.0, 2.0, 1.0, 7.0, 4.0, 4.0, 2.0, 5.0, 9.0, 10.0, 2.0, 3.0, 9.0, 4.0, 6.0, 1.0, 3.0, 5.0, 5.0, 8.0, 8.0, 6.0, 7.0, 4.0, 7.0, 9.0, 9.0, 8.0, 4.0, 4.0, 6.0, 3.0, 3.0, 9.0, 7.0, 1.0, 7.0, 8.0, 3.0, 5.0, 7.0, 4.0, 8.0, 6.0, 6.0, 1.0, 3.0, 3.0, 1.0, 3.0, 8.0, 6.0, 6.0, 4.0, 2.0, 4.0, 10.0, 4.0, 5.0, 9.0, 3.0, 5.0, 9.0, 9.0, 5.0, 4.0, 4.0, 1.0, 3.0, 5.0, 4.0, 6.0, 2.0, 6.0, 1.0, 9.0, 5.0, 5.0, 4.0, 7.0, 8.0, 5.0, 10.0, 10.0, 6.0, 4.0, 2.0, 9.0, 2.0, 6.0, 6.0, 3.0, 8.0, 4.0, 8.0, 2.0, 4.0, 7.0, 2.0, 2.0, 6.0, 7.0, 3.0, 2.0, 9.0, 7.0, 1.0, 10.0, 9.0, 7.0, 9.0, 9.0, 5.0, 9.0, 6.0, 2.0, 5.0, 10.0, 7.0, 8.0, 6.0, 2.0, 4.0, 10.0, 9.0, 7.0, 1.0, 7.0, 5.0, 1.0, 2.0, 2.0, 5.0, 10.0, 6.0, 6.0, 4.0, 3.0, 3.0, 2.0, 2.0, 9.0, 4.0, 8.0, 1.0, 2.0, 1.0, 4.0]
global b_x = 5
global d_y = [3.0, 9.0, 3.0, 1.0, 2.0, 6.0, 7.0, 4.0, 5.0, 1.0, 5.0, 10.0, 1.0, 5.0, 2.0, 2.0, 6.0, 6.0, 1.0, 5.0, 10.0, 3.0, 6.0, 8.0, 1.0, 2.0, 5.0, 9.0, 6.0, 6.0, 10.0, 4.0, 7.0, 10.0, 7.0, 10.0, 4.0, 8.0, 10.0, 5.0, 5.0, 1.0, 5.0, 3.0, 5.0, 10.0, 10.0, 9.0, 9.0, 7.0, 9.0, 6.0, 3.0, 9.0, 5.0, 7.0, 6.0, 6.0, 6.0, 2.0, 7.0, 9.0, 8.0, 7.0, 1.0, 7.0, 5.0, 6.0, 1.0, 6.0, 5.0, 7.0, 6.0, 3.0, 3.0, 9.0, 4.0, 6.0, 3.0, 5.0, 1.0, 3.0, 4.0, 8.0, 4.0, 1.0, 7.0, 1.0, 7.0, 8.0, 2.0, 9.0, 10.0, 10.0, 7.0, 1.0, 4.0, 7.0, 1.0, 7.0, 3.0, 10.0, 9.0, 5.0, 4.0, 5.0, 4.0, 6.0, 3.0, 7.0, 10.0, 6.0, 7.0, 10.0, 4.0, 6.0, 4.0, 4.0, 6.0, 1.0, 5.0, 3.0, 7.0, 7.0, 6.0, 9.0, 8.0, 3.0, 6.0, 5.0, 1.0, 6.0, 10.0, 7.0, 4.0, 4.0, 9.0, 8.0, 1.0, 7.0, 3.0, 6.0, 3.0, 8.0, 9.0, 4.0, 4.0, 1.0, 4.0, 4.0, 4.0, 7.0, 3.0, 1.0, 6.0, 3.0, 7.0, 10.0, 5.0, 10.0, 3.0, 5.0, 6.0, 6.0, 6.0, 5.0, 1.0, 9.0, 3.0, 7.0, 9.0, 2.0, 4.0, 7.0, 6.0, 10.0, 1.0, 2.0, 6.0, 6.0, 10.0, 1.0, 1.0, 6.0, 4.0, 6.0, 2.0, 10.0, 9.0, 2.0, 7.0, 9.0, 1.0, 9.0, 2.0, 7.0, 4.0, 1.0, 3.0, 10.0, 6.0, 6.0, 8.0, 10.0, 4.0, 10.0, 7.0, 6.0, 6.0, 5.0, 9.0, 5.0, 1.0, 4.0, 8.0, 8.0, 1.0, 6.0, 4.0, 4.0, 10.0, 9.0, 9.0, 9.0, 7.0, 1.0, 5.0, 5.0, 2.0, 9.0, 9.0, 8.0, 6.0, 10.0, 6.0, 3.0, 5.0, 6.0, 2.0, 6.0, 5.0, 6.0, 8.0, 5.0, 6.0, 4.0, 3.0, 6.0, 7.0, 1.0, 2.0, 8.0, 9.0, 3.0, 3.0, 6.0, 7.0, 5.0, 5.0, 1.0, 8.0, 2.0, 6.0, 10.0, 9.0, 7.0, 2.0, 3.0, 2.0, 7.0, 6.0, 4.0, 7.0, 10.0, 8.0, 9.0, 1.0, 10.0, 9.0, 5.0, 3.0, 6.0, 8.0, 8.0, 8.0, 3.0, 7.0, 5.0, 10.0, 8.0, 5.0, 4.0, 4.0, 7.0, 10.0, 2.0, 5.0, 2.0, 5.0, 4.0, 5.0, 1.0, 7.0, 10.0, 6.0, 9.0, 6.0, 10.0, 6.0, 7.0, 1.0, 3.0, 1.0, 10.0, 7.0, 1.0, 2.0, 6.0, 2.0, 4.0, 10.0, 7.0, 5.0, 3.0, 3.0, 1.0, 3.0, 3.0, 3.0, 7.0, 4.0, 4.0, 6.0]
global b_y = 10
global p = [0.134, 0.742, 0.341, 0.358, 0.601, 0.967, 0.098, 0.695, 0.627, 0.856, 0.928, 0.321, 0.352, 0.101, 0.763, 0.767, 0.766, 0.938, 0.338, 0.151, 0.274, 0.375, 0.666, 0.21, 0.232, 0.632, 0.691, 0.917, 0.359, 0.994, 0.731, 0.459, 0.887, 0.272, 0.082, 0.426, 0.162, 0.296, 0.747, 0.498, 0.936, 0.092, 0.993, 0.92, 0.031, 0.919, 0.107, 0.763, 0.28, 0.972, 0.605, 0.264, 0.911, 0.728, 0.487, 0.54, 0.489, 0.721, 0.799, 0.043, 0.078, 0.769, 0.029, 0.673, 0.129, 0.641, 0.731, 0.539, 0.336, 0.26, 0.04, 0.527, 0.362, 0.619, 0.462, 0.122, 0.115, 0.052, 0.739, 0.341, 0.384, 0.599, 0.728, 0.885, 0.495, 0.75, 0.343, 0.11, 0.924, 0.738, 0.727, 0.945, 0.513, 0.342, 0.606, 0.223, 0.117, 0.882, 0.716, 0.968, 0.522, 0.894, 0.011, 0.309, 0.607, 0.673, 0.72, 0.941, 0.607, 0.191, 0.638, 0.892, 0.864, 0.057, 0.145, 0.411, 0.087, 0.807, 0.679, 0.357, 0.086, 0.009, 0.574, 0.034, 0.045, 0.981, 0.213, 0.67, 0.735, 0.301, 0.556, 0.92, 0.694, 0.442, 0.149, 0.24, 0.478, 0.024, 0.846, 0.602, 0.456, 0.81, 0.878, 0.94, 0.498, 0.749, 0.772, 0.923, 0.433, 0.714, 0.267, 0.575, 0.898, 0.761, 0.334, 0.884, 0.473, 0.186, 0.049, 0.96, 0.997, 0.44, 0.226, 0.238, 0.398, 0.934, 0.594, 0.139, 0.848, 0.381, 0.686, 0.125, 0.144, 0.773, 0.505, 0.198, 0.625, 0.639, 0.719, 0.11, 0.134, 0.023, 0.398, 0.55, 0.437, 0.615, 0.495, 0.731, 0.164, 0.136, 0.766, 0.369, 0.297, 0.096, 0.746, 0.19, 0.066, 0.802, 0.738, 0.039, 0.684, 0.83, 0.014, 0.362, 0.369, 0.549, 0.969, 0.851, 0.33, 0.923, 0.542, 0.281, 0.049, 0.176, 0.623, 0.245, 0.679, 0.26, 0.049, 0.983, 0.214, 0.092, 0.586, 0.321, 0.653, 0.025, 0.389, 0.038, 0.753, 0.286, 0.812, 0.437, 0.56, 0.625, 0.832, 0.329, 0.235, 0.189, 0.861, 0.383, 0.997, 0.045, 0.643, 0.826, 0.782, 0.305, 0.254, 0.246, 0.886, 0.58, 0.122, 0.547, 0.138, 0.35, 0.725, 0.475, 0.385, 0.317, 0.956, 0.679, 0.876, 0.049, 0.747, 0.815, 0.341, 0.785, 0.557, 0.174, 0.015, 0.619, 0.243, 0.943, 0.84, 0.62, 0.394, 0.861, 0.313, 0.361, 0.967, 0.537, 0.961, 0.786, 0.18, 0.218, 0.391, 0.748, 0.217, 0.721, 0.444, 0.584, 0.177, 0.316, 0.159, 0.505, 0.743, 0.845, 0.124, 0.958, 0.227, 0.04, 0.883, 0.863, 0.805, 0.184, 0.472, 0.659, 0.797, 0.609, 0.387, 0.258, 0.874, 0.422, 0.567, 0.315, 0.996, 0.419, 0.528, 0.144, 0.918, 0.681, 0.556, 0.255, 0.789, 0.637, 0.213, 0.17, 0.072, 0.469, 0.292, 0.489, 0.934, 0.554, 0.889]
global q = [0.525, 0.984, 0.742, 0.949, 0.981, 0.982, 0.725, 0.997, 0.675, 0.928, 0.957, 0.397, 0.751, 0.177, 0.87, 0.894, 0.974, 0.974, 0.411, 0.747, 0.297, 0.997, 0.911, 0.75, 0.989, 0.676, 0.853, 0.995, 0.551, 0.999, 0.845, 0.76, 0.968, 0.487, 0.312, 0.915, 0.722, 0.636, 0.923, 0.703, 0.948, 0.199, 0.994, 0.932, 0.713, 0.992, 0.317, 0.787, 0.911, 0.99, 0.922, 0.52, 0.97, 0.898, 0.597, 0.55, 0.604, 0.899, 0.957, 0.673, 0.399, 0.898, 0.92, 0.676, 0.583, 0.728, 0.969, 0.593, 0.642, 0.362, 0.397, 0.802, 0.41, 0.766, 0.985, 0.525, 0.703, 0.055, 0.962, 0.896, 0.602, 0.921, 0.97, 0.99, 0.794, 0.852, 0.616, 0.246, 0.998, 0.847, 0.985, 0.977, 0.871, 0.625, 0.997, 0.366, 0.864, 0.963, 0.792, 0.979, 0.686, 0.982, 0.512, 0.395, 0.671, 0.852, 0.803, 0.979, 0.921, 0.655, 0.987, 0.966, 0.968, 0.074, 0.518, 0.463, 0.394, 0.941, 0.7, 0.674, 0.591, 0.201, 0.98, 0.748, 0.105, 0.996, 0.504, 0.796, 0.738, 0.425, 0.948, 0.937, 0.761, 0.771, 0.988, 0.856, 0.671, 0.275, 0.993, 0.968, 0.9, 0.829, 0.969, 0.978, 0.762, 0.99, 0.905, 0.992, 0.555, 0.834, 0.586, 0.91, 0.956, 0.764, 0.661, 0.962, 0.898, 0.793, 0.916, 0.996, 0.999, 0.642, 0.963, 0.392, 0.534, 0.935, 0.617, 0.38, 0.884, 0.456, 0.805, 0.61, 0.786, 0.931, 0.583, 0.87, 0.873, 0.665, 0.885, 0.393, 0.415, 0.244, 0.801, 0.68, 0.98, 0.94, 0.976, 0.995, 0.856, 0.283, 0.942, 0.962, 0.56, 0.512, 0.911, 0.533, 0.317, 0.82, 0.928, 0.734, 0.965, 0.898, 0.798, 0.5, 0.494, 0.737, 0.971, 0.951, 0.731, 0.996, 0.933, 0.571, 0.914, 0.403, 0.751, 0.578, 0.767, 0.305, 0.923, 0.998, 0.287, 0.121, 0.922, 0.578, 0.714, 0.875, 0.555, 0.974, 0.898, 0.321, 0.918, 0.49, 0.912, 0.701, 0.88, 0.673, 0.828, 0.634, 0.943, 0.918, 0.999, 0.171, 0.895, 0.942, 0.984, 0.489, 0.583, 0.802, 0.927, 0.602, 0.594, 0.919, 0.232, 0.43, 0.748, 0.647, 0.386, 0.464, 0.998, 0.841, 0.929, 0.11, 0.77, 0.846, 0.988, 0.949, 0.95, 0.778, 0.065, 0.681, 0.278, 0.952, 0.901, 0.927, 0.767, 0.939, 0.64, 0.823, 0.985, 0.984, 0.968, 0.897, 0.683, 0.307, 0.445, 0.94, 0.709, 0.774, 0.614, 0.889, 0.852, 0.738, 0.256, 0.9, 0.914, 0.85, 0.617, 0.971, 0.739, 0.196, 0.928, 0.901, 0.988, 0.457, 0.578, 0.862, 0.874, 0.741, 0.862, 0.373, 0.947, 0.529, 0.979, 0.968, 0.998, 0.654, 0.857, 0.816, 0.953, 0.981, 0.855, 0.853, 0.809, 0.812, 0.956, 0.601, 0.869, 0.703, 0.653, 0.535, 0.956, 0.749, 0.93]
global origin = 1
global destination = 60