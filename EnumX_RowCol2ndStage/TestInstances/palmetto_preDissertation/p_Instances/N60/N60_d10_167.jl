global arcs = [1 5; 1 25; 1 36; 1 40; 1 43; 1 48; 2 6; 2 16; 2 32; 2 37; 2 41; 2 48; 2 51; 3 6; 3 25; 3 32; 3 43; 3 55; 4 10; 4 20; 4 53; 4 60; 5 8; 5 10; 5 15; 5 20; 5 41; 5 45; 5 59; 6 21; 6 30; 6 34; 6 35; 6 43; 6 44; 6 47; 6 57; 7 2; 7 8; 7 11; 7 16; 7 27; 7 33; 7 37; 7 40; 8 10; 8 40; 8 41; 8 59; 9 4; 9 22; 9 34; 9 38; 9 43; 9 51; 9 54; 9 60; 10 8; 10 15; 10 22; 10 30; 10 42; 10 47; 10 54; 10 56; 11 10; 11 12; 11 22; 12 14; 12 21; 12 26; 12 29; 12 30; 12 33; 12 35; 12 40; 12 47; 12 60; 13 11; 13 33; 13 45; 13 48; 13 60; 14 13; 14 19; 14 21; 14 29; 14 47; 15 28; 15 41; 15 48; 16 43; 16 49; 16 51; 17 19; 17 21; 17 28; 17 30; 17 40; 17 47; 17 49; 18 9; 18 21; 18 23; 18 37; 18 39; 18 49; 18 53; 19 14; 19 16; 19 17; 19 21; 19 22; 19 23; 19 25; 19 33; 19 36; 20 7; 20 25; 20 46; 20 52; 21 16; 21 32; 21 34; 21 46; 21 47; 21 52; 22 4; 22 24; 22 29; 22 58; 23 11; 23 31; 23 37; 23 47; 23 57; 24 3; 24 15; 24 16; 24 19; 24 22; 24 27; 24 52; 24 56; 25 10; 25 51; 26 13; 26 19; 26 34; 26 46; 26 60; 27 19; 27 21; 27 24; 27 33; 27 51; 28 17; 28 36; 28 37; 28 58; 29 22; 29 27; 29 34; 29 46; 30 16; 30 18; 30 27; 30 53; 30 58; 31 42; 31 46; 32 2; 32 8; 32 22; 32 34; 32 46; 32 47; 32 49; 33 9; 33 14; 33 16; 33 39; 33 46; 33 50; 33 60; 34 7; 34 24; 34 31; 34 36; 34 54; 34 56; 35 6; 35 25; 35 29; 35 40; 35 42; 35 44; 35 46; 35 53; 35 57; 36 23; 36 29; 36 42; 36 43; 36 44; 36 46; 36 50; 36 52; 37 15; 37 33; 37 35; 37 51; 37 53; 38 6; 38 16; 38 22; 38 42; 38 45; 38 60; 39 6; 39 8; 39 10; 39 12; 40 26; 41 4; 41 16; 41 43; 42 7; 42 8; 42 14; 42 16; 42 29; 42 37; 42 45; 42 49; 42 53; 43 8; 43 11; 43 40; 43 50; 43 54; 43 56; 43 60; 44 2; 44 5; 44 14; 44 24; 44 25; 44 33; 44 39; 44 49; 45 14; 45 25; 45 52; 45 53; 46 8; 46 9; 46 11; 46 24; 46 40; 46 42; 46 43; 46 53; 47 20; 47 26; 47 43; 47 49; 47 54; 47 57; 48 10; 48 36; 48 55; 48 56; 48 57; 49 13; 49 15; 49 37; 49 43; 49 46; 49 59; 50 7; 50 8; 50 10; 50 20; 50 33; 50 35; 51 21; 51 23; 51 29; 51 37; 51 45; 51 48; 52 10; 52 11; 52 23; 52 42; 52 44; 52 57; 52 59; 53 8; 53 10; 53 13; 53 26; 53 30; 53 52; 53 58; 54 22; 54 26; 54 39; 55 9; 55 34; 55 60; 56 15; 56 18; 56 27; 56 58; 56 60; 57 19; 57 21; 57 29; 57 37; 57 55; 58 14; 58 19; 58 28; 58 37; 58 38; 58 40; 58 59; 59 3; 59 8; 59 12]
global d_x = [10.0, 7.0, 2.0, 3.0, 8.0, 4.0, 10.0, 6.0, 3.0, 10.0, 7.0, 4.0, 8.0, 6.0, 8.0, 4.0, 9.0, 9.0, 1.0, 8.0, 3.0, 5.0, 2.0, 1.0, 5.0, 3.0, 3.0, 3.0, 2.0, 6.0, 3.0, 4.0, 9.0, 10.0, 6.0, 2.0, 7.0, 5.0, 5.0, 6.0, 5.0, 8.0, 1.0, 5.0, 9.0, 1.0, 3.0, 7.0, 7.0, 8.0, 5.0, 5.0, 10.0, 7.0, 8.0, 6.0, 3.0, 2.0, 10.0, 5.0, 3.0, 7.0, 1.0, 8.0, 6.0, 9.0, 1.0, 10.0, 7.0, 4.0, 3.0, 8.0, 5.0, 9.0, 3.0, 1.0, 4.0, 5.0, 2.0, 7.0, 8.0, 8.0, 3.0, 3.0, 2.0, 9.0, 3.0, 9.0, 8.0, 4.0, 1.0, 9.0, 8.0, 1.0, 3.0, 8.0, 2.0, 1.0, 6.0, 4.0, 2.0, 3.0, 1.0, 10.0, 9.0, 2.0, 1.0, 7.0, 8.0, 8.0, 4.0, 8.0, 3.0, 10.0, 5.0, 1.0, 7.0, 9.0, 3.0, 7.0, 6.0, 6.0, 7.0, 3.0, 7.0, 8.0, 5.0, 7.0, 1.0, 10.0, 2.0, 2.0, 7.0, 8.0, 3.0, 7.0, 10.0, 2.0, 2.0, 10.0, 7.0, 6.0, 2.0, 9.0, 1.0, 9.0, 10.0, 1.0, 9.0, 7.0, 8.0, 3.0, 10.0, 7.0, 4.0, 7.0, 7.0, 5.0, 1.0, 2.0, 2.0, 6.0, 5.0, 2.0, 1.0, 6.0, 4.0, 10.0, 2.0, 2.0, 1.0, 3.0, 8.0, 4.0, 3.0, 6.0, 10.0, 3.0, 3.0, 1.0, 8.0, 4.0, 8.0, 7.0, 2.0, 6.0, 2.0, 6.0, 6.0, 8.0, 6.0, 5.0, 10.0, 2.0, 4.0, 8.0, 8.0, 5.0, 10.0, 3.0, 1.0, 1.0, 3.0, 1.0, 5.0, 7.0, 1.0, 9.0, 6.0, 4.0, 2.0, 4.0, 3.0, 10.0, 10.0, 1.0, 5.0, 8.0, 9.0, 8.0, 1.0, 9.0, 3.0, 3.0, 4.0, 8.0, 8.0, 2.0, 9.0, 10.0, 3.0, 7.0, 3.0, 8.0, 1.0, 8.0, 5.0, 2.0, 9.0, 6.0, 3.0, 10.0, 9.0, 8.0, 8.0, 1.0, 2.0, 8.0, 2.0, 10.0, 8.0, 2.0, 4.0, 3.0, 4.0, 10.0, 8.0, 1.0, 8.0, 4.0, 8.0, 7.0, 7.0, 7.0, 10.0, 5.0, 5.0, 10.0, 6.0, 2.0, 8.0, 8.0, 1.0, 4.0, 8.0, 10.0, 7.0, 5.0, 9.0, 6.0, 7.0, 5.0, 1.0, 6.0, 9.0, 7.0, 3.0, 5.0, 2.0, 7.0, 3.0, 1.0, 8.0, 2.0, 3.0, 10.0, 4.0, 8.0, 8.0, 7.0, 8.0, 8.0, 4.0, 3.0, 6.0, 10.0, 10.0, 10.0, 5.0, 2.0, 2.0, 9.0, 1.0, 1.0, 2.0, 10.0, 9.0, 7.0, 5.0, 9.0, 10.0, 3.0, 10.0, 5.0, 1.0, 2.0, 6.0, 6.0, 3.0, 2.0, 7.0, 10.0]
global b_x = 5
global d_y = [2.0, 1.0, 2.0, 6.0, 3.0, 4.0, 10.0, 8.0, 4.0, 9.0, 3.0, 8.0, 1.0, 3.0, 6.0, 6.0, 5.0, 1.0, 2.0, 6.0, 9.0, 2.0, 5.0, 3.0, 5.0, 8.0, 2.0, 8.0, 3.0, 10.0, 7.0, 7.0, 9.0, 8.0, 1.0, 5.0, 4.0, 7.0, 5.0, 9.0, 5.0, 1.0, 10.0, 8.0, 10.0, 2.0, 3.0, 3.0, 9.0, 6.0, 2.0, 1.0, 5.0, 9.0, 4.0, 8.0, 2.0, 8.0, 5.0, 5.0, 7.0, 3.0, 2.0, 8.0, 9.0, 4.0, 1.0, 1.0, 7.0, 7.0, 10.0, 4.0, 4.0, 3.0, 4.0, 7.0, 9.0, 4.0, 9.0, 9.0, 4.0, 1.0, 1.0, 2.0, 5.0, 5.0, 1.0, 6.0, 5.0, 5.0, 1.0, 3.0, 10.0, 3.0, 9.0, 1.0, 2.0, 2.0, 9.0, 10.0, 6.0, 4.0, 4.0, 2.0, 8.0, 9.0, 7.0, 8.0, 5.0, 9.0, 6.0, 10.0, 10.0, 6.0, 8.0, 10.0, 7.0, 4.0, 1.0, 7.0, 1.0, 10.0, 9.0, 1.0, 2.0, 7.0, 8.0, 8.0, 5.0, 9.0, 8.0, 7.0, 10.0, 5.0, 4.0, 8.0, 5.0, 3.0, 7.0, 4.0, 8.0, 1.0, 1.0, 9.0, 10.0, 7.0, 6.0, 2.0, 2.0, 3.0, 1.0, 3.0, 6.0, 7.0, 4.0, 1.0, 5.0, 6.0, 8.0, 5.0, 6.0, 8.0, 5.0, 9.0, 3.0, 4.0, 1.0, 7.0, 8.0, 3.0, 6.0, 8.0, 2.0, 8.0, 2.0, 3.0, 2.0, 7.0, 8.0, 1.0, 1.0, 4.0, 1.0, 8.0, 8.0, 3.0, 7.0, 1.0, 2.0, 7.0, 6.0, 2.0, 7.0, 1.0, 2.0, 5.0, 6.0, 8.0, 5.0, 5.0, 5.0, 6.0, 5.0, 6.0, 1.0, 9.0, 10.0, 6.0, 7.0, 3.0, 7.0, 3.0, 8.0, 5.0, 1.0, 2.0, 5.0, 4.0, 3.0, 6.0, 2.0, 10.0, 4.0, 8.0, 6.0, 8.0, 6.0, 1.0, 9.0, 8.0, 5.0, 8.0, 2.0, 1.0, 7.0, 6.0, 1.0, 7.0, 2.0, 10.0, 7.0, 5.0, 2.0, 1.0, 5.0, 3.0, 4.0, 10.0, 4.0, 2.0, 5.0, 9.0, 9.0, 6.0, 5.0, 1.0, 5.0, 4.0, 1.0, 9.0, 8.0, 1.0, 3.0, 4.0, 9.0, 10.0, 6.0, 9.0, 7.0, 3.0, 7.0, 9.0, 3.0, 10.0, 3.0, 8.0, 3.0, 9.0, 3.0, 1.0, 5.0, 2.0, 5.0, 6.0, 7.0, 9.0, 6.0, 4.0, 4.0, 1.0, 4.0, 4.0, 9.0, 5.0, 3.0, 9.0, 8.0, 5.0, 4.0, 6.0, 4.0, 1.0, 8.0, 5.0, 1.0, 10.0, 5.0, 3.0, 5.0, 5.0, 2.0, 1.0, 5.0, 9.0, 4.0, 1.0, 1.0, 2.0, 6.0, 10.0, 8.0, 10.0, 6.0, 2.0, 10.0, 9.0, 1.0, 10.0, 1.0, 1.0, 3.0, 3.0]
global b_y = 10
global p = [0.971, 0.754, 0.028, 0.45, 0.593, 0.648, 0.242, 0.009, 0.457, 0.996, 0.305, 0.704, 0.463, 0.861, 0.768, 0.997, 0.072, 0.256, 0.33, 0.646, 0.561, 0.078, 0.16, 0.333, 0.133, 0.061, 0.308, 0.491, 0.502, 0.961, 0.352, 0.27, 0.387, 0.195, 0.404, 0.416, 0.617, 0.33, 0.161, 0.931, 0.488, 0.856, 0.806, 0.926, 0.422, 0.358, 0.796, 0.968, 0.393, 0.64, 0.444, 0.506, 0.038, 0.168, 0.486, 0.099, 0.477, 0.766, 0.075, 0.223, 0.559, 0.693, 0.025, 0.825, 0.199, 0.747, 0.584, 0.387, 0.619, 0.684, 0.529, 0.227, 0.498, 0.927, 0.7, 0.589, 0.258, 0.375, 0.97, 0.964, 0.067, 0.63, 0.525, 0.239, 0.671, 0.478, 0.481, 0.521, 0.509, 0.77, 0.416, 0.178, 0.524, 0.448, 0.782, 0.88, 0.34, 0.005, 0.752, 0.148, 0.07, 0.171, 0.765, 0.361, 0.111, 0.789, 0.843, 0.756, 0.29, 0.647, 0.624, 0.019, 0.106, 0.074, 0.039, 0.598, 0.656, 0.786, 0.658, 0.082, 0.815, 0.729, 0.092, 0.997, 0.465, 0.5, 0.429, 0.713, 0.903, 0.896, 0.727, 0.83, 0.864, 0.157, 0.385, 0.331, 0.767, 0.015, 0.325, 0.254, 0.9, 0.448, 0.615, 0.815, 0.629, 0.743, 0.468, 0.1, 0.112, 0.058, 0.264, 0.196, 0.054, 0.326, 0.728, 0.334, 0.407, 0.467, 0.145, 0.99, 0.048, 0.952, 0.659, 0.941, 0.418, 0.777, 0.63, 0.094, 0.357, 0.23, 0.309, 0.819, 0.888, 0.921, 0.21, 0.665, 0.632, 0.101, 0.987, 0.933, 0.013, 0.693, 0.148, 0.038, 0.774, 0.338, 0.465, 0.812, 0.924, 0.353, 0.259, 0.119, 0.573, 0.883, 0.258, 0.043, 0.75, 0.671, 0.315, 0.365, 0.656, 0.534, 0.512, 0.771, 0.489, 0.896, 0.279, 0.169, 0.029, 0.5, 0.32, 0.928, 0.359, 0.13, 0.836, 0.242, 0.749, 0.89, 0.657, 0.633, 0.018, 0.407, 0.674, 0.311, 0.138, 0.444, 0.586, 0.805, 0.858, 0.161, 0.621, 0.839, 0.974, 0.225, 0.407, 0.007, 0.666, 0.207, 0.874, 0.632, 0.147, 0.141, 0.758, 0.253, 0.561, 0.178, 0.481, 0.101, 0.264, 0.975, 0.227, 0.763, 0.39, 0.211, 0.019, 0.34, 0.812, 0.036, 0.811, 0.298, 0.157, 0.403, 0.217, 0.326, 0.852, 0.553, 0.782, 0.306, 0.945, 0.164, 0.52, 0.663, 0.429, 0.132, 0.257, 0.936, 0.866, 0.715, 0.907, 0.54, 0.901, 0.403, 0.41, 0.883, 0.066, 0.818, 0.113, 0.315, 0.244, 0.012, 0.456, 0.399, 0.731, 0.698, 0.57, 0.823, 0.258, 0.261, 0.126, 0.229, 0.311, 0.351, 0.838, 0.462, 0.396, 0.778, 0.985, 0.394, 0.641, 0.083, 0.606, 0.877, 0.205, 0.224, 0.228, 0.853, 0.114, 0.259, 0.882, 0.275, 0.533, 0.274, 0.06, 0.083, 0.448, 0.069, 0.836, 0.046, 0.987, 0.7, 0.149, 0.278]
global q = [0.991, 0.759, 0.23, 0.52, 0.836, 0.98, 0.622, 0.368, 0.902, 0.999, 0.746, 0.799, 0.94, 0.928, 0.868, 0.998, 0.946, 0.61, 0.982, 0.877, 0.874, 0.245, 0.221, 0.531, 0.3, 0.621, 0.321, 0.699, 0.963, 0.983, 0.403, 0.721, 0.758, 0.407, 0.561, 0.545, 0.729, 0.667, 0.631, 0.948, 0.97, 0.973, 0.857, 0.991, 0.676, 0.493, 0.905, 0.981, 0.839, 0.978, 0.641, 0.911, 0.043, 0.827, 0.881, 0.821, 0.942, 0.959, 0.193, 0.68, 0.804, 0.933, 0.679, 0.893, 0.21, 0.998, 0.74, 0.98, 0.875, 0.71, 0.825, 0.99, 0.863, 0.963, 0.705, 0.774, 0.724, 0.466, 0.984, 0.964, 0.212, 0.721, 0.828, 0.435, 0.708, 0.504, 0.975, 0.653, 0.709, 0.868, 0.662, 0.831, 0.897, 0.533, 0.845, 0.988, 0.937, 0.778, 0.816, 0.972, 0.237, 0.743, 0.791, 0.648, 0.61, 0.971, 0.948, 0.866, 0.805, 0.891, 0.794, 0.338, 0.126, 0.994, 0.099, 0.948, 0.968, 0.788, 0.801, 0.347, 0.979, 0.901, 0.196, 0.998, 0.705, 0.641, 0.998, 0.828, 0.96, 0.995, 0.881, 0.877, 0.929, 0.664, 0.421, 0.845, 0.781, 0.823, 0.977, 0.689, 0.945, 0.889, 0.885, 0.918, 0.69, 0.954, 0.902, 0.711, 0.333, 0.145, 0.653, 0.205, 0.117, 0.339, 0.808, 0.817, 0.521, 0.839, 0.396, 0.999, 0.646, 0.98, 0.792, 0.987, 0.55, 0.902, 0.669, 0.287, 0.403, 0.509, 0.348, 0.858, 0.892, 0.927, 0.644, 0.866, 0.773, 0.351, 0.997, 0.993, 0.926, 0.868, 0.939, 0.045, 0.784, 0.878, 0.536, 0.954, 0.981, 0.652, 0.843, 0.202, 0.645, 0.949, 0.641, 0.243, 0.944, 0.674, 0.831, 0.611, 0.881, 0.772, 0.556, 0.989, 0.622, 0.965, 0.611, 0.982, 0.399, 0.535, 0.999, 0.976, 0.765, 0.914, 0.955, 0.618, 0.948, 0.988, 0.885, 0.684, 0.284, 0.444, 0.829, 0.852, 0.754, 0.743, 0.966, 0.93, 0.96, 0.511, 0.93, 0.961, 0.992, 0.952, 0.708, 0.372, 0.74, 0.312, 0.9, 0.919, 0.558, 0.519, 0.773, 0.796, 0.681, 0.99, 0.928, 0.67, 0.622, 0.984, 0.537, 0.963, 0.875, 0.868, 0.316, 0.458, 0.872, 0.307, 0.924, 0.985, 0.505, 0.612, 0.782, 0.363, 0.91, 0.755, 0.858, 0.449, 0.961, 0.514, 0.896, 0.954, 0.592, 0.412, 0.662, 0.954, 0.88, 0.8, 0.953, 0.728, 0.949, 0.806, 0.894, 0.926, 0.801, 0.857, 0.333, 0.525, 0.992, 0.118, 0.995, 0.741, 0.735, 0.709, 0.949, 0.947, 0.815, 0.645, 0.545, 0.633, 0.647, 0.784, 0.846, 0.59, 0.992, 0.911, 0.996, 0.894, 0.883, 0.951, 0.883, 0.955, 0.6, 0.797, 0.263, 0.977, 0.251, 0.651, 0.995, 0.855, 0.563, 0.53, 0.635, 0.188, 0.835, 0.094, 0.963, 0.656, 0.998, 0.843, 0.514, 0.356]
global origin = 1
global destination = 60