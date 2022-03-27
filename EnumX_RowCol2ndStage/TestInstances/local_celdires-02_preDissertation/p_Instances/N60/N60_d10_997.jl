global arcs = [1 15; 1 20; 1 29; 1 60; 2 16; 2 23; 2 48; 3 10; 3 18; 3 34; 3 49; 3 58; 4 15; 4 42; 4 43; 4 52; 5 28; 5 40; 5 49; 5 57; 6 16; 6 26; 6 37; 6 47; 7 12; 7 28; 7 33; 7 34; 7 38; 7 44; 7 49; 7 57; 8 6; 8 17; 8 19; 8 24; 8 35; 8 44; 8 47; 8 57; 9 2; 9 5; 9 27; 9 36; 9 37; 10 2; 10 3; 10 14; 10 18; 10 29; 10 35; 11 3; 11 17; 11 23; 11 40; 11 45; 12 17; 12 25; 12 26; 12 32; 12 54; 12 59; 12 60; 13 8; 13 40; 13 55; 14 2; 14 11; 14 25; 14 27; 14 36; 14 50; 14 56; 15 10; 15 24; 15 29; 15 31; 15 32; 15 40; 15 47; 15 48; 15 53; 15 60; 16 4; 16 6; 16 37; 16 45; 16 48; 16 49; 16 59; 17 2; 17 13; 17 41; 17 50; 18 11; 18 30; 18 49; 18 52; 18 53; 19 2; 19 4; 19 15; 19 26; 19 41; 19 46; 19 53; 19 55; 19 56; 19 60; 20 22; 20 25; 20 38; 20 46; 20 49; 20 50; 20 56; 20 57; 21 14; 21 16; 21 19; 21 28; 21 41; 22 13; 22 23; 22 27; 22 33; 22 36; 22 51; 22 53; 23 5; 23 14; 23 26; 24 23; 24 38; 24 41; 25 13; 25 19; 25 30; 25 41; 25 53; 26 12; 26 21; 26 28; 26 31; 26 37; 27 20; 27 24; 27 35; 27 41; 27 50; 27 51; 27 58; 28 9; 28 13; 28 14; 28 47; 28 50; 29 31; 29 44; 29 45; 29 54; 29 59; 30 20; 30 38; 30 55; 31 18; 31 19; 31 32; 31 48; 31 49; 32 18; 32 41; 32 43; 32 60; 33 3; 33 28; 33 32; 33 42; 33 48; 33 49; 34 16; 34 24; 34 42; 34 45; 34 49; 35 20; 35 40; 35 43; 35 47; 35 54; 35 56; 36 4; 36 38; 36 40; 36 43; 36 57; 36 59; 36 60; 37 4; 37 9; 37 12; 37 16; 37 21; 37 41; 37 47; 37 51; 37 60; 38 13; 38 29; 38 54; 39 2; 39 21; 39 40; 39 41; 39 58; 40 41; 41 10; 41 34; 41 45; 41 51; 42 3; 42 5; 42 15; 42 23; 42 28; 42 31; 42 35; 42 43; 42 51; 42 52; 42 54; 43 4; 43 21; 43 23; 43 32; 43 44; 43 45; 43 46; 43 54; 43 55; 43 57; 44 40; 44 45; 44 50; 44 53; 45 5; 45 6; 45 16; 45 21; 45 32; 45 33; 45 40; 45 43; 45 44; 45 48; 45 60; 46 40; 46 43; 47 16; 47 17; 47 18; 47 21; 47 28; 47 35; 48 11; 48 33; 48 51; 48 52; 48 53; 49 3; 49 5; 49 6; 49 7; 49 15; 49 18; 49 42; 49 52; 49 57; 50 2; 50 6; 50 14; 50 17; 50 21; 50 51; 51 10; 51 33; 51 39; 51 42; 51 55; 52 16; 52 18; 52 28; 52 32; 53 7; 53 9; 53 31; 53 35; 53 48; 54 13; 54 19; 54 26; 55 13; 55 14; 55 34; 55 53; 56 47; 56 48; 57 2; 57 5; 57 8; 57 11; 57 12; 57 21; 57 32; 57 59; 58 6; 58 11; 58 12; 58 19; 58 41; 59 9; 59 11; 59 30; 59 49]
global d_x = [1.0, 9.0, 2.0, 3.0, 7.0, 4.0, 8.0, 7.0, 2.0, 3.0, 2.0, 8.0, 5.0, 8.0, 8.0, 5.0, 10.0, 8.0, 10.0, 10.0, 4.0, 2.0, 6.0, 2.0, 9.0, 1.0, 5.0, 2.0, 5.0, 1.0, 10.0, 4.0, 5.0, 10.0, 6.0, 9.0, 6.0, 10.0, 5.0, 5.0, 3.0, 10.0, 6.0, 6.0, 3.0, 2.0, 7.0, 3.0, 1.0, 2.0, 4.0, 4.0, 8.0, 9.0, 5.0, 7.0, 7.0, 10.0, 7.0, 6.0, 2.0, 3.0, 6.0, 9.0, 1.0, 6.0, 10.0, 4.0, 10.0, 4.0, 8.0, 5.0, 1.0, 8.0, 2.0, 7.0, 4.0, 6.0, 6.0, 10.0, 10.0, 7.0, 1.0, 10.0, 1.0, 9.0, 6.0, 1.0, 3.0, 7.0, 9.0, 3.0, 2.0, 6.0, 10.0, 3.0, 4.0, 3.0, 2.0, 1.0, 4.0, 3.0, 4.0, 4.0, 4.0, 2.0, 10.0, 4.0, 4.0, 3.0, 9.0, 5.0, 3.0, 10.0, 8.0, 7.0, 4.0, 4.0, 1.0, 7.0, 6.0, 5.0, 8.0, 9.0, 1.0, 10.0, 6.0, 4.0, 6.0, 3.0, 7.0, 7.0, 8.0, 10.0, 5.0, 10.0, 1.0, 2.0, 7.0, 7.0, 7.0, 7.0, 7.0, 3.0, 6.0, 1.0, 4.0, 3.0, 1.0, 7.0, 6.0, 5.0, 10.0, 7.0, 4.0, 10.0, 9.0, 10.0, 6.0, 8.0, 10.0, 6.0, 2.0, 9.0, 5.0, 8.0, 1.0, 1.0, 6.0, 8.0, 5.0, 10.0, 10.0, 6.0, 6.0, 1.0, 6.0, 10.0, 9.0, 4.0, 3.0, 9.0, 4.0, 2.0, 6.0, 3.0, 6.0, 2.0, 6.0, 8.0, 5.0, 9.0, 7.0, 9.0, 1.0, 3.0, 3.0, 2.0, 9.0, 10.0, 3.0, 5.0, 4.0, 9.0, 8.0, 7.0, 2.0, 1.0, 5.0, 9.0, 1.0, 10.0, 10.0, 6.0, 3.0, 5.0, 9.0, 8.0, 6.0, 4.0, 8.0, 1.0, 7.0, 8.0, 6.0, 8.0, 2.0, 9.0, 6.0, 1.0, 8.0, 7.0, 8.0, 10.0, 7.0, 7.0, 7.0, 10.0, 8.0, 9.0, 3.0, 8.0, 10.0, 5.0, 8.0, 8.0, 2.0, 10.0, 7.0, 7.0, 2.0, 3.0, 8.0, 8.0, 2.0, 9.0, 7.0, 3.0, 1.0, 2.0, 6.0, 10.0, 10.0, 4.0, 9.0, 7.0, 1.0, 2.0, 6.0, 5.0, 4.0, 2.0, 6.0, 8.0, 10.0, 2.0, 6.0, 2.0, 1.0, 10.0, 10.0, 5.0, 2.0, 3.0, 6.0, 8.0, 6.0, 2.0, 8.0, 6.0, 10.0, 10.0, 2.0, 8.0, 6.0, 3.0, 9.0, 3.0, 4.0, 6.0, 8.0, 3.0, 10.0, 6.0, 9.0, 2.0, 10.0, 5.0, 3.0, 8.0, 9.0, 10.0, 3.0, 4.0, 8.0, 4.0, 10.0, 7.0, 8.0, 1.0, 7.0, 3.0, 9.0, 5.0]
global b_x = 5
global d_y = [10.0, 6.0, 1.0, 8.0, 7.0, 10.0, 1.0, 3.0, 1.0, 1.0, 7.0, 7.0, 1.0, 2.0, 9.0, 3.0, 3.0, 10.0, 5.0, 8.0, 5.0, 1.0, 10.0, 2.0, 4.0, 7.0, 10.0, 2.0, 7.0, 9.0, 4.0, 6.0, 5.0, 6.0, 2.0, 5.0, 4.0, 1.0, 1.0, 9.0, 4.0, 8.0, 6.0, 6.0, 4.0, 9.0, 8.0, 6.0, 2.0, 2.0, 9.0, 5.0, 6.0, 1.0, 4.0, 8.0, 9.0, 5.0, 4.0, 9.0, 7.0, 8.0, 2.0, 4.0, 6.0, 2.0, 9.0, 2.0, 4.0, 8.0, 9.0, 1.0, 4.0, 1.0, 1.0, 10.0, 4.0, 3.0, 6.0, 6.0, 6.0, 8.0, 4.0, 10.0, 9.0, 2.0, 4.0, 1.0, 5.0, 7.0, 6.0, 5.0, 9.0, 9.0, 1.0, 9.0, 6.0, 8.0, 10.0, 9.0, 3.0, 7.0, 7.0, 5.0, 5.0, 2.0, 3.0, 5.0, 5.0, 6.0, 5.0, 3.0, 8.0, 9.0, 1.0, 6.0, 4.0, 1.0, 2.0, 4.0, 7.0, 4.0, 7.0, 10.0, 8.0, 3.0, 6.0, 5.0, 4.0, 4.0, 4.0, 2.0, 3.0, 1.0, 8.0, 6.0, 8.0, 3.0, 4.0, 7.0, 7.0, 7.0, 6.0, 5.0, 4.0, 1.0, 1.0, 7.0, 3.0, 9.0, 2.0, 2.0, 8.0, 6.0, 1.0, 10.0, 1.0, 4.0, 2.0, 9.0, 10.0, 4.0, 5.0, 4.0, 10.0, 4.0, 1.0, 7.0, 1.0, 7.0, 4.0, 5.0, 9.0, 2.0, 2.0, 9.0, 2.0, 3.0, 2.0, 7.0, 1.0, 6.0, 3.0, 10.0, 3.0, 5.0, 4.0, 8.0, 9.0, 3.0, 8.0, 10.0, 6.0, 3.0, 5.0, 1.0, 6.0, 4.0, 7.0, 8.0, 4.0, 7.0, 3.0, 4.0, 6.0, 3.0, 2.0, 9.0, 1.0, 3.0, 3.0, 6.0, 1.0, 9.0, 1.0, 9.0, 5.0, 3.0, 5.0, 4.0, 9.0, 1.0, 3.0, 5.0, 4.0, 6.0, 1.0, 6.0, 5.0, 1.0, 8.0, 10.0, 3.0, 6.0, 4.0, 5.0, 5.0, 2.0, 2.0, 8.0, 4.0, 1.0, 7.0, 5.0, 3.0, 8.0, 8.0, 3.0, 7.0, 3.0, 9.0, 5.0, 10.0, 3.0, 4.0, 5.0, 7.0, 2.0, 10.0, 4.0, 10.0, 1.0, 9.0, 9.0, 9.0, 3.0, 9.0, 2.0, 9.0, 3.0, 4.0, 10.0, 6.0, 10.0, 3.0, 3.0, 1.0, 2.0, 7.0, 10.0, 4.0, 2.0, 4.0, 4.0, 9.0, 3.0, 5.0, 8.0, 4.0, 10.0, 5.0, 2.0, 7.0, 5.0, 3.0, 8.0, 9.0, 8.0, 5.0, 7.0, 9.0, 5.0, 4.0, 1.0, 1.0, 5.0, 4.0, 5.0, 3.0, 2.0, 8.0, 10.0, 6.0, 8.0, 10.0, 4.0, 9.0, 8.0, 5.0, 2.0, 8.0, 5.0, 7.0, 7.0]
global b_y = 10
global p = [0.037, 0.926, 0.228, 0.712, 0.702, 0.387, 0.522, 0.214, 0.295, 0.375, 0.705, 0.26, 0.092, 0.383, 0.074, 0.629, 0.503, 0.908, 0.147, 0.513, 0.536, 0.302, 0.755, 0.227, 0.895, 0.797, 0.968, 0.877, 0.174, 0.327, 0.229, 0.75, 0.071, 0.209, 0.56, 0.338, 0.121, 0.951, 0.619, 0.481, 0.002, 0.18, 0.034, 0.963, 0.069, 0.112, 0.721, 0.792, 0.22, 0.366, 0.4, 0.818, 0.662, 0.316, 0.028, 0.562, 0.422, 0.178, 0.008, 0.715, 0.552, 0.296, 0.647, 0.753, 0.232, 0.876, 0.619, 0.097, 0.185, 0.808, 0.267, 0.446, 0.749, 0.49, 0.186, 0.036, 0.953, 0.595, 0.75, 0.522, 0.747, 0.572, 0.229, 0.168, 0.047, 0.307, 0.194, 0.548, 0.748, 0.33, 0.807, 0.905, 0.642, 0.515, 0.13, 0.636, 0.776, 0.724, 0.5, 0.397, 0.897, 0.623, 0.821, 0.81, 0.131, 0.362, 0.878, 0.894, 0.831, 0.82, 0.691, 0.032, 0.444, 0.809, 0.522, 0.776, 0.189, 0.637, 0.112, 0.443, 0.293, 0.585, 0.064, 0.2, 0.973, 0.817, 0.784, 0.255, 0.986, 0.064, 0.087, 0.688, 0.662, 0.584, 0.738, 0.699, 0.463, 0.634, 0.375, 0.664, 0.177, 0.899, 0.78, 0.733, 0.047, 0.766, 0.904, 0.279, 0.764, 0.703, 0.042, 0.454, 0.067, 0.095, 0.272, 0.844, 0.321, 0.556, 0.348, 0.056, 0.723, 0.529, 0.579, 0.288, 0.613, 0.594, 0.357, 0.617, 0.977, 0.492, 0.198, 0.351, 0.506, 0.59, 0.17, 0.652, 0.818, 0.157, 0.002, 0.272, 0.256, 0.575, 0.45, 0.715, 0.386, 0.902, 0.954, 0.414, 0.858, 0.322, 0.102, 0.843, 0.102, 0.703, 0.563, 0.489, 0.169, 0.254, 0.787, 0.317, 0.445, 0.218, 0.715, 0.749, 0.371, 0.714, 0.953, 0.571, 0.682, 0.018, 0.018, 0.405, 0.234, 0.192, 0.265, 0.254, 0.778, 0.169, 0.931, 0.845, 0.268, 0.906, 0.849, 0.906, 0.718, 0.979, 0.719, 0.986, 0.975, 0.22, 0.499, 0.616, 0.53, 0.676, 0.969, 0.01, 0.338, 0.166, 0.846, 0.693, 0.816, 0.765, 0.647, 0.74, 0.237, 0.639, 0.964, 0.834, 0.873, 0.468, 0.703, 0.596, 0.4, 0.39, 0.376, 0.774, 0.191, 0.628, 0.765, 0.674, 0.921, 0.52, 0.128, 0.839, 0.114, 0.735, 0.766, 0.748, 0.475, 0.551, 0.384, 0.737, 0.495, 0.282, 0.368, 0.395, 0.609, 0.197, 0.628, 0.11, 0.77, 0.817, 0.432, 0.126, 0.266, 0.17, 0.358, 0.033, 0.333, 0.818, 0.459, 0.132, 0.642, 0.294, 0.157, 0.654, 0.169, 0.278, 0.426, 0.314, 0.067, 0.196, 0.02, 0.828, 0.552, 0.287, 0.209, 0.551, 0.197, 0.147, 0.07, 0.062, 0.204, 0.948, 0.328, 0.885, 0.417, 0.336, 0.175, 0.738, 0.751, 0.408, 0.656, 0.569]
global q = [0.371, 0.934, 0.651, 0.883, 0.923, 0.906, 0.9, 0.809, 0.422, 0.941, 0.94, 0.67, 0.374, 0.423, 0.392, 0.652, 0.63, 0.974, 0.874, 0.595, 0.646, 0.891, 0.84, 0.882, 0.912, 0.808, 0.983, 0.992, 0.435, 0.979, 0.57, 0.794, 0.826, 0.457, 0.679, 0.704, 0.383, 0.958, 0.663, 0.801, 0.333, 0.305, 0.781, 0.963, 0.235, 0.921, 0.939, 0.862, 0.639, 0.721, 0.919, 0.9, 0.781, 0.476, 0.527, 0.738, 0.906, 0.618, 0.137, 0.815, 0.717, 0.414, 0.737, 0.975, 0.468, 0.914, 0.622, 0.857, 0.638, 0.84, 0.32, 0.585, 0.89, 0.607, 0.936, 0.398, 0.973, 0.796, 0.985, 0.849, 0.779, 0.581, 0.304, 0.303, 0.653, 0.625, 0.677, 0.747, 0.932, 0.64, 0.825, 0.934, 0.825, 0.728, 0.372, 0.671, 0.956, 0.772, 0.594, 0.503, 0.927, 0.815, 0.974, 0.884, 0.586, 0.952, 0.882, 0.95, 0.913, 0.859, 0.769, 0.243, 0.524, 0.883, 0.869, 0.976, 0.663, 0.818, 0.426, 0.735, 0.807, 0.724, 0.365, 0.648, 0.982, 0.97, 0.86, 0.994, 0.988, 0.834, 0.1, 0.808, 0.69, 0.597, 0.738, 0.792, 0.517, 0.777, 0.897, 0.847, 0.324, 0.918, 0.925, 0.991, 0.297, 0.881, 0.952, 0.949, 0.995, 0.979, 0.841, 0.528, 0.733, 0.73, 0.479, 0.944, 0.487, 0.611, 0.349, 0.783, 0.986, 0.6, 0.803, 0.645, 0.926, 0.957, 0.582, 0.772, 0.991, 0.744, 0.445, 0.675, 0.635, 0.972, 0.396, 0.903, 0.868, 0.526, 0.635, 0.993, 0.585, 0.967, 0.48, 0.776, 0.506, 0.952, 0.974, 0.624, 0.927, 0.861, 0.395, 0.887, 0.212, 0.901, 0.869, 0.853, 0.436, 0.509, 0.945, 0.62, 0.766, 0.224, 0.949, 0.88, 0.448, 0.769, 0.969, 0.741, 0.791, 0.998, 0.526, 0.792, 0.314, 0.369, 0.598, 0.933, 0.966, 0.199, 0.952, 0.864, 0.584, 0.961, 0.924, 0.912, 0.904, 0.998, 0.903, 0.994, 0.983, 0.746, 0.831, 0.871, 0.87, 0.886, 0.983, 0.156, 0.839, 0.271, 0.94, 0.979, 0.821, 0.835, 0.673, 0.898, 0.293, 0.9, 0.984, 0.848, 0.876, 0.758, 0.789, 0.947, 0.802, 0.486, 0.376, 0.924, 0.889, 0.695, 0.793, 0.945, 0.962, 0.975, 0.747, 0.968, 0.744, 0.878, 0.976, 0.795, 0.983, 0.953, 0.526, 0.866, 0.537, 0.613, 0.615, 0.444, 0.952, 0.645, 0.769, 0.21, 0.792, 0.987, 0.941, 0.921, 0.267, 0.344, 0.468, 0.523, 0.626, 0.945, 0.853, 0.96, 0.81, 0.591, 0.923, 0.815, 0.962, 0.604, 0.474, 0.368, 0.09, 0.396, 0.665, 0.936, 0.701, 0.43, 0.688, 0.58, 0.329, 0.742, 0.938, 0.405, 0.862, 0.985, 0.58, 0.917, 0.44, 0.481, 0.601, 0.806, 0.994, 0.96, 0.788, 0.926]
global origin = 1
global destination = 60