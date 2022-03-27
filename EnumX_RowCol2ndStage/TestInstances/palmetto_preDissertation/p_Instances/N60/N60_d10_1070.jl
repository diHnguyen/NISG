global arcs = [1 12; 1 13; 1 15; 1 16; 1 31; 1 32; 1 49; 1 60; 2 7; 2 19; 2 26; 2 34; 2 42; 2 58; 3 7; 3 9; 3 14; 3 56; 4 6; 4 7; 4 22; 4 26; 4 34; 4 36; 4 43; 5 3; 5 39; 5 51; 6 4; 6 5; 6 7; 6 27; 6 33; 6 39; 6 44; 6 46; 7 10; 7 14; 7 23; 7 29; 7 32; 7 38; 8 6; 8 13; 8 17; 8 31; 9 32; 9 33; 9 40; 9 43; 10 2; 10 17; 10 34; 10 39; 10 41; 10 45; 10 46; 10 47; 10 50; 10 53; 11 4; 11 8; 11 48; 12 2; 12 15; 12 17; 12 47; 13 9; 13 20; 13 24; 13 43; 13 58; 14 10; 14 34; 14 38; 14 51; 15 12; 15 26; 15 30; 15 36; 15 47; 15 51; 15 60; 16 29; 16 35; 16 36; 17 11; 17 20; 17 30; 17 48; 17 54; 17 55; 18 21; 18 22; 18 27; 18 28; 18 31; 18 43; 19 7; 19 9; 19 20; 19 24; 19 32; 19 48; 19 51; 19 52; 19 58; 20 7; 20 10; 20 11; 20 34; 20 37; 20 46; 21 15; 21 35; 21 36; 21 42; 21 54; 22 11; 22 12; 22 27; 22 34; 22 40; 23 30; 23 38; 23 41; 23 58; 24 7; 24 16; 24 25; 24 30; 24 57; 25 2; 25 8; 25 19; 25 20; 25 22; 25 45; 25 56; 26 27; 26 37; 26 44; 26 52; 27 7; 27 16; 27 18; 27 29; 27 37; 27 58; 28 42; 28 56; 28 59; 29 2; 29 6; 29 8; 29 33; 29 35; 29 42; 29 45; 30 12; 30 15; 30 23; 30 32; 30 40; 30 56; 30 59; 31 5; 31 17; 31 19; 31 26; 31 30; 31 34; 31 54; 32 11; 32 48; 32 56; 32 58; 33 36; 33 37; 34 20; 34 22; 34 26; 34 44; 35 4; 35 12; 35 13; 35 26; 35 29; 35 43; 35 57; 36 17; 36 21; 36 38; 36 42; 36 50; 36 59; 37 35; 37 36; 37 42; 38 3; 38 11; 38 13; 38 15; 38 21; 38 26; 38 55; 38 57; 38 58; 39 11; 39 13; 39 19; 39 24; 39 48; 40 14; 40 31; 40 52; 41 22; 41 36; 41 45; 42 17; 42 24; 42 35; 43 3; 43 7; 43 19; 43 26; 43 45; 44 3; 44 25; 44 60; 45 7; 45 9; 45 51; 45 56; 46 8; 46 12; 46 14; 46 15; 46 27; 46 31; 46 43; 46 45; 46 48; 46 54; 47 2; 47 17; 47 27; 47 36; 47 55; 47 57; 48 11; 49 26; 49 30; 49 38; 49 58; 49 60; 50 26; 50 38; 50 44; 50 54; 51 27; 51 35; 51 46; 52 19; 52 22; 52 53; 53 6; 53 20; 53 23; 53 30; 53 34; 53 39; 53 49; 54 39; 54 42; 54 47; 54 49; 55 3; 55 5; 55 20; 55 44; 56 3; 56 9; 56 14; 56 22; 56 24; 56 31; 56 49; 56 54; 57 6; 57 16; 57 19; 57 22; 57 28; 57 29; 57 42; 57 51; 57 52; 58 7; 58 11; 58 13; 58 14; 58 17; 58 22; 58 24; 58 27; 58 44; 58 45; 58 53; 58 56; 58 57; 59 8; 59 15; 59 20; 59 45; 59 47; 59 51]
global d_x = [1.0, 10.0, 10.0, 6.0, 3.0, 10.0, 10.0, 5.0, 1.0, 9.0, 3.0, 10.0, 2.0, 2.0, 7.0, 1.0, 7.0, 10.0, 8.0, 2.0, 7.0, 9.0, 7.0, 5.0, 10.0, 9.0, 1.0, 5.0, 5.0, 8.0, 3.0, 7.0, 2.0, 3.0, 1.0, 7.0, 4.0, 6.0, 8.0, 2.0, 8.0, 3.0, 2.0, 1.0, 7.0, 10.0, 7.0, 5.0, 10.0, 9.0, 7.0, 2.0, 7.0, 6.0, 9.0, 5.0, 1.0, 10.0, 4.0, 6.0, 1.0, 8.0, 7.0, 1.0, 8.0, 9.0, 6.0, 1.0, 9.0, 8.0, 9.0, 5.0, 4.0, 1.0, 9.0, 10.0, 7.0, 4.0, 2.0, 6.0, 6.0, 5.0, 2.0, 7.0, 7.0, 6.0, 3.0, 8.0, 5.0, 5.0, 4.0, 1.0, 7.0, 10.0, 6.0, 1.0, 1.0, 3.0, 4.0, 10.0, 10.0, 1.0, 4.0, 8.0, 7.0, 1.0, 5.0, 10.0, 1.0, 7.0, 2.0, 9.0, 3.0, 2.0, 3.0, 10.0, 7.0, 7.0, 8.0, 7.0, 10.0, 5.0, 9.0, 7.0, 3.0, 1.0, 8.0, 7.0, 2.0, 3.0, 7.0, 9.0, 2.0, 6.0, 7.0, 7.0, 2.0, 8.0, 4.0, 7.0, 2.0, 7.0, 10.0, 1.0, 4.0, 1.0, 7.0, 1.0, 10.0, 4.0, 7.0, 3.0, 5.0, 3.0, 2.0, 5.0, 4.0, 1.0, 9.0, 3.0, 2.0, 7.0, 1.0, 9.0, 9.0, 8.0, 8.0, 9.0, 3.0, 7.0, 4.0, 4.0, 1.0, 7.0, 5.0, 9.0, 7.0, 10.0, 2.0, 10.0, 10.0, 10.0, 7.0, 4.0, 7.0, 9.0, 5.0, 10.0, 6.0, 9.0, 6.0, 5.0, 2.0, 1.0, 5.0, 5.0, 1.0, 5.0, 7.0, 8.0, 1.0, 1.0, 7.0, 8.0, 8.0, 6.0, 3.0, 8.0, 10.0, 10.0, 5.0, 6.0, 3.0, 6.0, 3.0, 5.0, 5.0, 6.0, 10.0, 4.0, 3.0, 9.0, 1.0, 2.0, 3.0, 3.0, 9.0, 4.0, 3.0, 4.0, 7.0, 2.0, 3.0, 6.0, 8.0, 1.0, 3.0, 8.0, 7.0, 7.0, 3.0, 8.0, 2.0, 5.0, 2.0, 7.0, 5.0, 6.0, 2.0, 3.0, 4.0, 7.0, 2.0, 8.0, 4.0, 7.0, 10.0, 6.0, 10.0, 3.0, 1.0, 3.0, 5.0, 2.0, 2.0, 7.0, 7.0, 7.0, 7.0, 6.0, 8.0, 4.0, 1.0, 9.0, 4.0, 7.0, 1.0, 5.0, 4.0, 4.0, 5.0, 2.0, 5.0, 8.0, 7.0, 7.0, 6.0, 1.0, 2.0, 8.0, 3.0, 5.0, 9.0, 8.0, 5.0, 4.0, 7.0, 8.0, 10.0, 10.0, 6.0, 5.0, 1.0, 10.0, 6.0, 5.0, 8.0, 4.0, 5.0, 7.0, 7.0, 4.0, 5.0, 3.0, 2.0, 6.0, 7.0]
global b_x = 5
global d_y = [8.0, 5.0, 7.0, 1.0, 9.0, 6.0, 3.0, 2.0, 8.0, 1.0, 2.0, 8.0, 10.0, 1.0, 7.0, 5.0, 4.0, 4.0, 1.0, 7.0, 10.0, 6.0, 2.0, 3.0, 6.0, 3.0, 6.0, 2.0, 8.0, 9.0, 4.0, 1.0, 3.0, 8.0, 10.0, 1.0, 4.0, 2.0, 10.0, 10.0, 5.0, 5.0, 10.0, 5.0, 7.0, 6.0, 2.0, 6.0, 4.0, 10.0, 3.0, 4.0, 10.0, 6.0, 2.0, 4.0, 7.0, 5.0, 5.0, 3.0, 7.0, 4.0, 1.0, 9.0, 6.0, 5.0, 4.0, 2.0, 6.0, 8.0, 9.0, 4.0, 5.0, 7.0, 5.0, 8.0, 8.0, 5.0, 4.0, 1.0, 10.0, 6.0, 10.0, 4.0, 3.0, 1.0, 2.0, 4.0, 10.0, 6.0, 1.0, 3.0, 3.0, 4.0, 9.0, 2.0, 6.0, 2.0, 7.0, 10.0, 9.0, 8.0, 10.0, 1.0, 8.0, 3.0, 6.0, 1.0, 6.0, 3.0, 6.0, 10.0, 4.0, 4.0, 9.0, 7.0, 3.0, 2.0, 6.0, 2.0, 4.0, 3.0, 1.0, 6.0, 7.0, 2.0, 1.0, 9.0, 5.0, 8.0, 10.0, 4.0, 10.0, 6.0, 6.0, 3.0, 4.0, 4.0, 8.0, 3.0, 4.0, 10.0, 8.0, 1.0, 4.0, 7.0, 7.0, 7.0, 4.0, 9.0, 2.0, 4.0, 3.0, 7.0, 3.0, 2.0, 6.0, 6.0, 3.0, 8.0, 3.0, 10.0, 3.0, 3.0, 9.0, 9.0, 10.0, 8.0, 8.0, 1.0, 8.0, 8.0, 3.0, 6.0, 4.0, 10.0, 9.0, 10.0, 5.0, 9.0, 9.0, 1.0, 1.0, 6.0, 2.0, 3.0, 7.0, 1.0, 1.0, 7.0, 5.0, 1.0, 6.0, 3.0, 7.0, 5.0, 2.0, 3.0, 3.0, 9.0, 6.0, 10.0, 9.0, 8.0, 2.0, 5.0, 5.0, 10.0, 6.0, 8.0, 7.0, 2.0, 3.0, 9.0, 1.0, 10.0, 1.0, 9.0, 8.0, 6.0, 6.0, 10.0, 2.0, 4.0, 5.0, 5.0, 1.0, 9.0, 10.0, 5.0, 2.0, 9.0, 4.0, 4.0, 3.0, 9.0, 3.0, 2.0, 3.0, 6.0, 9.0, 5.0, 5.0, 5.0, 10.0, 7.0, 2.0, 5.0, 1.0, 7.0, 10.0, 1.0, 9.0, 10.0, 8.0, 3.0, 9.0, 4.0, 2.0, 6.0, 10.0, 4.0, 4.0, 10.0, 5.0, 1.0, 9.0, 5.0, 2.0, 5.0, 10.0, 6.0, 8.0, 7.0, 4.0, 9.0, 3.0, 10.0, 10.0, 9.0, 8.0, 8.0, 9.0, 9.0, 3.0, 4.0, 2.0, 4.0, 9.0, 3.0, 10.0, 2.0, 6.0, 2.0, 9.0, 7.0, 6.0, 5.0, 9.0, 4.0, 2.0, 10.0, 6.0, 6.0, 9.0, 8.0, 4.0, 10.0, 4.0, 3.0, 5.0, 8.0, 10.0, 6.0, 9.0, 1.0, 7.0]
global b_y = 10
global p = [0.705, 0.827, 0.145, 0.172, 0.085, 0.059, 0.18, 0.049, 0.383, 0.026, 0.6, 0.037, 0.57, 0.331, 0.048, 0.726, 0.765, 0.682, 0.133, 0.595, 0.363, 0.87, 0.446, 0.021, 0.598, 0.346, 0.592, 0.11, 0.745, 0.457, 0.248, 0.419, 0.884, 0.12, 0.627, 0.918, 0.916, 0.484, 0.853, 0.275, 0.255, 0.877, 0.846, 0.421, 0.817, 0.54, 0.449, 0.48, 0.674, 0.86, 0.925, 0.07, 0.696, 0.114, 0.265, 0.694, 0.627, 0.563, 0.931, 0.754, 0.217, 0.635, 0.147, 0.38, 0.536, 0.024, 0.12, 0.15, 0.013, 0.878, 0.581, 0.876, 0.127, 0.786, 0.184, 0.12, 0.693, 0.819, 0.47, 0.239, 0.723, 0.279, 0.856, 0.385, 0.029, 0.783, 0.85, 0.235, 0.102, 0.854, 0.975, 0.169, 0.134, 0.799, 0.15, 0.653, 0.242, 0.723, 0.51, 0.048, 0.888, 0.955, 0.371, 0.606, 0.358, 0.354, 0.555, 0.414, 0.789, 0.229, 0.774, 0.002, 0.177, 0.654, 0.67, 0.775, 0.976, 0.833, 0.742, 0.733, 0.435, 0.621, 0.172, 0.244, 0.898, 0.006, 0.895, 0.81, 0.622, 0.416, 0.202, 0.035, 0.552, 0.299, 0.876, 0.016, 0.98, 0.456, 0.569, 0.762, 0.551, 0.208, 0.013, 0.081, 0.794, 0.664, 0.027, 0.516, 0.333, 0.935, 0.796, 0.814, 0.174, 0.906, 0.88, 0.908, 0.71, 0.558, 0.895, 0.137, 0.908, 0.614, 0.654, 0.339, 0.497, 0.261, 0.464, 0.161, 0.392, 0.117, 0.813, 0.921, 0.045, 0.256, 0.304, 0.467, 0.467, 0.729, 0.745, 0.663, 0.242, 0.847, 0.487, 0.048, 0.486, 0.971, 0.952, 0.32, 0.998, 0.151, 0.749, 0.536, 0.649, 0.64, 0.488, 0.581, 0.003, 0.777, 0.379, 0.022, 0.967, 0.19, 0.356, 0.9, 0.2, 0.031, 0.144, 0.702, 0.013, 0.86, 0.092, 0.122, 0.533, 0.609, 0.737, 0.604, 0.509, 0.886, 0.689, 0.099, 0.211, 0.306, 0.538, 0.977, 0.16, 0.935, 0.205, 0.325, 0.417, 0.85, 0.348, 0.975, 0.557, 0.677, 0.643, 0.33, 0.861, 0.432, 0.971, 0.024, 0.672, 0.381, 0.635, 0.834, 0.819, 0.039, 0.159, 0.658, 0.425, 0.679, 0.077, 0.528, 0.292, 0.974, 0.024, 0.054, 0.218, 0.309, 0.043, 0.262, 0.421, 0.942, 0.186, 0.546, 0.474, 0.194, 0.405, 0.415, 0.558, 0.04, 0.874, 0.375, 0.181, 0.252, 0.605, 0.739, 0.306, 0.433, 0.712, 0.337, 0.185, 0.262, 0.124, 0.523, 0.431, 0.497, 0.908, 0.949, 0.493, 0.492, 0.527, 0.76, 0.964, 0.175, 0.51, 0.598, 0.459, 0.334, 0.654, 0.616, 0.853, 0.712, 0.214, 0.726, 0.711, 0.768, 0.614, 0.274, 0.434, 0.064, 0.138, 0.512, 0.022, 0.647, 0.673, 0.522, 0.482]
global q = [0.728, 0.996, 0.362, 0.253, 0.107, 0.164, 0.984, 0.212, 0.608, 0.636, 0.922, 0.303, 0.618, 0.665, 0.393, 0.782, 0.999, 0.896, 0.735, 0.812, 0.688, 0.983, 0.629, 0.355, 0.697, 0.703, 0.718, 0.586, 0.961, 0.916, 0.284, 0.539, 0.961, 0.813, 0.72, 0.938, 0.942, 0.593, 0.903, 0.574, 0.781, 0.925, 0.942, 0.433, 0.824, 0.592, 0.984, 0.753, 0.776, 0.984, 0.96, 0.194, 0.813, 0.691, 0.703, 0.914, 0.818, 0.763, 0.954, 0.801, 0.613, 0.803, 0.538, 0.549, 0.71, 0.839, 0.567, 0.358, 0.725, 0.933, 0.885, 0.911, 0.201, 0.875, 0.391, 0.483, 0.705, 0.864, 0.876, 0.954, 0.798, 0.391, 0.998, 0.575, 0.859, 0.983, 0.919, 0.298, 0.488, 0.932, 0.997, 0.937, 0.799, 0.893, 0.691, 0.983, 0.519, 0.869, 0.545, 0.239, 0.969, 0.977, 0.805, 0.926, 0.994, 0.84, 0.638, 0.718, 0.869, 0.422, 0.982, 0.723, 0.863, 0.722, 0.843, 0.831, 0.987, 0.99, 0.938, 0.963, 0.558, 0.993, 0.989, 0.433, 0.988, 0.463, 0.906, 0.916, 0.723, 0.795, 0.822, 0.43, 0.608, 0.803, 0.969, 0.157, 0.991, 0.746, 0.598, 0.791, 0.913, 0.432, 0.08, 0.791, 0.804, 0.774, 0.557, 0.954, 0.96, 0.962, 0.804, 0.96, 0.928, 0.953, 0.906, 0.948, 0.737, 0.678, 0.921, 0.945, 0.982, 0.638, 0.716, 0.722, 0.815, 0.402, 0.776, 0.429, 0.856, 0.252, 0.824, 0.945, 0.648, 0.356, 0.631, 0.825, 0.687, 0.856, 0.872, 0.933, 0.26, 0.857, 0.678, 0.416, 0.577, 0.996, 0.966, 0.343, 0.999, 0.214, 0.894, 0.714, 0.65, 0.776, 0.556, 0.606, 0.514, 0.825, 0.51, 0.938, 0.991, 0.874, 0.457, 0.995, 0.857, 0.259, 0.279, 0.744, 0.676, 0.946, 0.118, 0.399, 0.959, 0.706, 0.847, 0.808, 0.561, 0.995, 0.763, 0.37, 0.836, 0.612, 0.94, 0.983, 0.418, 0.975, 0.825, 0.393, 0.965, 0.993, 0.81, 0.996, 0.88, 0.872, 0.866, 0.383, 0.969, 0.733, 0.999, 0.693, 0.787, 0.998, 0.865, 0.844, 0.872, 0.548, 0.657, 0.945, 0.919, 0.979, 0.079, 0.687, 0.811, 0.986, 0.817, 0.664, 0.317, 0.412, 0.801, 0.268, 0.73, 0.993, 0.556, 0.815, 0.683, 0.815, 0.513, 0.874, 0.662, 0.056, 0.938, 0.813, 0.518, 0.492, 0.848, 0.784, 0.906, 0.601, 0.899, 0.538, 0.495, 0.414, 0.459, 0.628, 0.897, 0.664, 0.998, 0.974, 0.547, 0.681, 0.853, 0.845, 0.992, 0.426, 0.869, 0.987, 0.521, 0.368, 0.669, 0.958, 0.886, 0.843, 0.628, 0.789, 0.834, 0.92, 0.654, 0.974, 0.769, 0.688, 0.18, 0.888, 0.341, 0.679, 0.836, 0.94, 0.527]
global origin = 1
global destination = 60