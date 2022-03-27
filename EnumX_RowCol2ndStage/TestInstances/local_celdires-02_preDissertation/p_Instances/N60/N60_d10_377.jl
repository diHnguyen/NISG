global arcs = [1 3; 1 4; 1 7; 1 29; 1 30; 1 46; 1 48; 1 50; 2 4; 2 7; 2 10; 2 38; 2 47; 3 11; 3 13; 3 23; 3 28; 3 41; 4 21; 4 29; 5 11; 5 25; 5 37; 5 38; 5 54; 6 7; 6 29; 6 35; 6 38; 6 40; 6 59; 7 5; 7 6; 7 22; 7 24; 7 31; 7 33; 7 51; 8 6; 8 7; 8 15; 8 26; 8 33; 8 42; 8 54; 8 57; 9 49; 9 59; 10 16; 10 37; 10 43; 10 48; 11 7; 11 33; 11 41; 12 9; 12 14; 12 17; 12 21; 12 23; 12 47; 12 48; 12 52; 12 60; 13 5; 13 17; 13 21; 13 28; 13 29; 13 52; 13 53; 13 55; 13 60; 14 2; 14 15; 14 19; 14 36; 14 56; 15 5; 15 7; 15 9; 15 10; 15 11; 15 22; 15 24; 15 30; 15 33; 15 34; 15 55; 15 58; 15 60; 16 14; 16 31; 16 36; 16 51; 16 53; 17 7; 17 19; 17 22; 17 29; 17 33; 17 42; 18 5; 18 17; 18 19; 18 20; 18 30; 18 51; 18 53; 19 8; 19 15; 19 16; 19 21; 19 30; 19 32; 19 33; 19 39; 19 43; 19 47; 19 53; 20 4; 20 12; 20 29; 20 33; 20 36; 20 40; 20 42; 21 9; 21 35; 21 50; 21 52; 22 15; 22 16; 22 31; 22 43; 22 47; 22 52; 22 54; 23 7; 23 13; 23 22; 23 27; 23 39; 23 59; 24 28; 24 29; 24 40; 24 44; 24 49; 25 19; 25 48; 25 51; 25 52; 26 8; 26 10; 26 28; 26 36; 26 41; 26 43; 26 54; 26 56; 26 60; 27 8; 27 11; 27 16; 27 31; 27 43; 27 45; 27 58; 28 12; 28 51; 28 52; 28 58; 29 22; 29 48; 29 55; 29 59; 30 13; 30 38; 30 39; 30 60; 31 6; 31 11; 31 59; 32 14; 32 42; 33 10; 33 21; 33 38; 33 52; 33 53; 34 5; 34 10; 34 18; 34 21; 34 24; 34 30; 35 2; 35 5; 35 12; 35 22; 35 30; 35 42; 36 5; 36 23; 36 44; 36 54; 36 55; 37 11; 37 14; 37 15; 37 33; 37 38; 37 49; 38 6; 38 12; 38 13; 38 35; 38 43; 38 49; 38 56; 39 10; 39 11; 39 12; 39 32; 39 42; 39 47; 40 34; 40 38; 40 53; 40 60; 41 4; 41 8; 41 11; 41 13; 41 20; 41 28; 41 39; 41 40; 41 49; 42 11; 42 14; 42 33; 42 36; 43 9; 43 18; 43 24; 43 27; 43 57; 44 3; 44 9; 44 14; 44 15; 44 17; 44 23; 44 50; 44 56; 45 3; 45 7; 45 8; 45 13; 45 26; 45 50; 46 4; 46 18; 46 22; 46 24; 46 27; 46 43; 46 48; 47 7; 47 29; 47 48; 47 53; 48 5; 48 24; 48 25; 48 39; 48 51; 48 58; 49 3; 49 18; 49 27; 49 44; 50 15; 50 18; 50 22; 50 30; 50 38; 50 56; 51 11; 51 15; 51 25; 51 31; 51 42; 51 46; 51 54; 52 19; 52 27; 52 44; 52 57; 53 26; 53 36; 53 37; 53 57; 54 2; 54 3; 54 4; 54 7; 54 14; 54 23; 54 28; 54 34; 54 46; 54 47; 55 2; 55 8; 55 14; 55 26; 55 30; 55 31; 55 35; 55 38; 55 48; 55 57; 56 11; 56 18; 56 32; 56 37; 56 51; 56 58; 57 37; 57 50; 57 51; 57 53; 57 54; 58 17; 58 47; 58 50; 58 57; 59 47; 59 55]
global d_x = [8.0, 2.0, 8.0, 5.0, 8.0, 7.0, 1.0, 1.0, 6.0, 9.0, 3.0, 2.0, 9.0, 6.0, 5.0, 5.0, 5.0, 8.0, 1.0, 2.0, 9.0, 10.0, 9.0, 6.0, 4.0, 10.0, 2.0, 10.0, 3.0, 6.0, 8.0, 7.0, 4.0, 4.0, 5.0, 10.0, 7.0, 4.0, 10.0, 2.0, 1.0, 5.0, 6.0, 8.0, 8.0, 10.0, 3.0, 6.0, 5.0, 5.0, 3.0, 6.0, 3.0, 5.0, 1.0, 7.0, 3.0, 2.0, 2.0, 9.0, 3.0, 7.0, 7.0, 1.0, 3.0, 5.0, 10.0, 7.0, 1.0, 7.0, 4.0, 10.0, 10.0, 2.0, 2.0, 1.0, 1.0, 5.0, 6.0, 10.0, 5.0, 3.0, 8.0, 7.0, 9.0, 7.0, 9.0, 4.0, 9.0, 8.0, 8.0, 8.0, 4.0, 5.0, 10.0, 8.0, 1.0, 3.0, 6.0, 3.0, 6.0, 10.0, 5.0, 7.0, 1.0, 6.0, 6.0, 3.0, 6.0, 2.0, 8.0, 8.0, 4.0, 8.0, 6.0, 9.0, 4.0, 10.0, 6.0, 3.0, 1.0, 6.0, 6.0, 6.0, 3.0, 2.0, 10.0, 1.0, 8.0, 7.0, 1.0, 2.0, 9.0, 8.0, 6.0, 5.0, 5.0, 9.0, 4.0, 5.0, 8.0, 9.0, 9.0, 6.0, 7.0, 1.0, 2.0, 5.0, 5.0, 7.0, 6.0, 10.0, 7.0, 9.0, 3.0, 1.0, 6.0, 2.0, 9.0, 7.0, 4.0, 5.0, 8.0, 7.0, 5.0, 9.0, 2.0, 10.0, 6.0, 2.0, 9.0, 1.0, 10.0, 8.0, 4.0, 5.0, 5.0, 7.0, 9.0, 2.0, 4.0, 5.0, 1.0, 6.0, 10.0, 1.0, 10.0, 1.0, 3.0, 9.0, 4.0, 9.0, 4.0, 8.0, 4.0, 1.0, 4.0, 9.0, 5.0, 8.0, 9.0, 9.0, 10.0, 1.0, 6.0, 6.0, 3.0, 2.0, 8.0, 1.0, 1.0, 9.0, 9.0, 1.0, 3.0, 4.0, 1.0, 7.0, 9.0, 2.0, 7.0, 1.0, 5.0, 6.0, 2.0, 8.0, 6.0, 6.0, 10.0, 9.0, 2.0, 1.0, 7.0, 4.0, 6.0, 4.0, 8.0, 2.0, 1.0, 4.0, 9.0, 5.0, 1.0, 10.0, 9.0, 6.0, 4.0, 4.0, 3.0, 1.0, 4.0, 7.0, 2.0, 1.0, 10.0, 1.0, 6.0, 2.0, 1.0, 10.0, 1.0, 8.0, 4.0, 8.0, 8.0, 2.0, 4.0, 5.0, 10.0, 10.0, 8.0, 1.0, 10.0, 4.0, 9.0, 10.0, 5.0, 9.0, 2.0, 3.0, 1.0, 7.0, 4.0, 4.0, 5.0, 8.0, 6.0, 8.0, 3.0, 3.0, 4.0, 2.0, 3.0, 2.0, 9.0, 2.0, 7.0, 4.0, 10.0, 6.0, 3.0, 3.0, 4.0, 4.0, 2.0, 4.0, 8.0, 10.0, 9.0, 10.0, 7.0, 6.0, 5.0, 8.0, 7.0, 1.0, 8.0, 10.0, 1.0, 6.0, 5.0, 9.0, 8.0, 10.0, 5.0, 9.0, 6.0, 1.0, 9.0, 3.0, 5.0, 1.0, 10.0, 1.0, 9.0, 7.0, 5.0, 9.0, 8.0, 1.0, 10.0, 2.0]
global b_x = 5
global d_y = [1.0, 3.0, 3.0, 5.0, 2.0, 10.0, 2.0, 1.0, 10.0, 10.0, 7.0, 6.0, 4.0, 4.0, 2.0, 8.0, 5.0, 1.0, 3.0, 8.0, 7.0, 8.0, 7.0, 7.0, 8.0, 2.0, 10.0, 1.0, 4.0, 5.0, 1.0, 10.0, 3.0, 1.0, 7.0, 5.0, 10.0, 4.0, 4.0, 9.0, 10.0, 6.0, 3.0, 7.0, 7.0, 1.0, 9.0, 7.0, 10.0, 2.0, 6.0, 7.0, 9.0, 3.0, 9.0, 2.0, 7.0, 2.0, 9.0, 4.0, 4.0, 1.0, 2.0, 9.0, 2.0, 8.0, 8.0, 7.0, 8.0, 4.0, 4.0, 4.0, 10.0, 9.0, 2.0, 7.0, 10.0, 10.0, 7.0, 9.0, 8.0, 4.0, 6.0, 10.0, 1.0, 1.0, 9.0, 10.0, 1.0, 9.0, 9.0, 7.0, 9.0, 10.0, 4.0, 4.0, 5.0, 10.0, 2.0, 1.0, 4.0, 7.0, 10.0, 7.0, 6.0, 4.0, 4.0, 5.0, 3.0, 2.0, 1.0, 2.0, 5.0, 9.0, 2.0, 7.0, 5.0, 3.0, 10.0, 7.0, 4.0, 3.0, 6.0, 3.0, 1.0, 8.0, 1.0, 7.0, 4.0, 1.0, 2.0, 8.0, 9.0, 6.0, 7.0, 10.0, 2.0, 4.0, 6.0, 2.0, 9.0, 1.0, 4.0, 8.0, 2.0, 5.0, 3.0, 7.0, 9.0, 10.0, 3.0, 9.0, 2.0, 2.0, 2.0, 7.0, 5.0, 10.0, 10.0, 5.0, 1.0, 1.0, 7.0, 3.0, 9.0, 4.0, 4.0, 6.0, 9.0, 8.0, 2.0, 6.0, 6.0, 3.0, 9.0, 6.0, 9.0, 4.0, 4.0, 6.0, 9.0, 4.0, 6.0, 9.0, 1.0, 4.0, 3.0, 2.0, 7.0, 5.0, 5.0, 6.0, 8.0, 2.0, 2.0, 2.0, 7.0, 7.0, 8.0, 6.0, 5.0, 6.0, 1.0, 1.0, 1.0, 8.0, 9.0, 5.0, 3.0, 2.0, 5.0, 5.0, 6.0, 8.0, 7.0, 3.0, 6.0, 2.0, 2.0, 7.0, 1.0, 1.0, 1.0, 10.0, 6.0, 9.0, 1.0, 4.0, 8.0, 10.0, 4.0, 1.0, 6.0, 9.0, 8.0, 5.0, 2.0, 1.0, 3.0, 4.0, 7.0, 10.0, 3.0, 4.0, 9.0, 3.0, 8.0, 7.0, 3.0, 10.0, 8.0, 6.0, 6.0, 4.0, 7.0, 7.0, 2.0, 2.0, 7.0, 8.0, 2.0, 4.0, 1.0, 1.0, 5.0, 4.0, 1.0, 3.0, 5.0, 8.0, 7.0, 2.0, 5.0, 10.0, 2.0, 9.0, 1.0, 10.0, 1.0, 6.0, 5.0, 3.0, 5.0, 6.0, 2.0, 7.0, 8.0, 5.0, 8.0, 2.0, 10.0, 5.0, 9.0, 8.0, 5.0, 2.0, 6.0, 6.0, 1.0, 10.0, 6.0, 10.0, 6.0, 5.0, 8.0, 2.0, 4.0, 7.0, 8.0, 3.0, 10.0, 7.0, 9.0, 1.0, 6.0, 5.0, 6.0, 10.0, 5.0, 10.0, 10.0, 7.0, 9.0, 7.0, 7.0, 10.0, 3.0, 3.0, 5.0, 9.0, 6.0, 1.0, 6.0, 10.0, 5.0, 2.0, 3.0, 7.0, 6.0, 6.0, 7.0, 3.0]
global b_y = 10
global p = [0.432, 0.776, 0.145, 0.342, 0.321, 0.243, 0.803, 0.71, 0.464, 0.232, 0.247, 0.457, 0.849, 0.561, 0.501, 0.316, 0.398, 0.942, 0.982, 0.15, 0.33, 0.271, 0.543, 0.278, 0.801, 0.587, 0.444, 0.526, 0.64, 0.081, 0.22, 0.949, 0.41, 0.456, 0.604, 0.889, 0.214, 0.929, 0.392, 0.448, 0.167, 0.644, 0.529, 0.677, 0.753, 0.775, 0.987, 0.854, 0.64, 0.04, 0.957, 0.653, 0.071, 0.05, 0.375, 0.381, 0.601, 0.963, 0.815, 0.831, 0.806, 0.972, 0.522, 0.415, 0.106, 0.775, 0.131, 0.303, 0.289, 0.11, 0.537, 0.096, 0.919, 0.628, 0.542, 0.942, 0.328, 0.484, 0.048, 0.696, 0.272, 0.183, 0.8, 0.424, 0.594, 0.875, 0.277, 0.4, 0.772, 0.056, 0.797, 0.097, 0.432, 0.542, 0.83, 0.099, 0.238, 0.455, 0.512, 0.839, 0.635, 0.793, 0.348, 0.905, 0.069, 0.574, 0.312, 0.189, 0.741, 0.825, 0.126, 0.101, 0.524, 0.795, 0.174, 0.426, 0.882, 0.799, 0.886, 0.648, 0.411, 0.445, 0.922, 0.288, 0.764, 0.207, 0.923, 0.502, 0.118, 0.679, 0.441, 0.826, 0.249, 0.353, 0.298, 0.552, 0.804, 0.234, 0.966, 0.694, 0.535, 0.29, 0.658, 0.941, 0.941, 0.145, 0.106, 0.363, 0.479, 0.573, 0.987, 0.235, 0.686, 0.054, 0.045, 0.851, 0.257, 0.634, 0.785, 0.13, 0.37, 0.163, 0.974, 0.758, 0.882, 0.103, 0.046, 0.035, 0.063, 0.415, 0.575, 0.751, 0.578, 0.191, 0.79, 0.546, 0.502, 0.625, 0.941, 0.19, 0.047, 0.381, 0.488, 0.987, 0.598, 0.432, 0.862, 0.415, 0.154, 0.699, 0.809, 0.244, 0.173, 0.237, 0.351, 0.367, 0.663, 0.845, 0.193, 0.358, 0.764, 0.816, 0.11, 0.211, 0.111, 0.177, 0.361, 0.499, 0.425, 0.848, 0.046, 0.236, 0.029, 0.15, 0.148, 0.201, 0.188, 0.524, 0.681, 0.542, 0.595, 0.734, 0.85, 0.514, 0.297, 0.041, 0.263, 0.223, 0.279, 0.059, 0.982, 0.125, 0.662, 0.1, 0.445, 0.611, 0.827, 0.615, 0.225, 0.305, 0.901, 0.285, 0.535, 0.455, 0.771, 0.881, 0.294, 0.567, 0.123, 0.59, 0.873, 0.454, 0.718, 0.431, 0.222, 0.873, 0.121, 0.747, 0.541, 0.066, 0.922, 0.858, 0.898, 0.309, 0.503, 0.355, 0.212, 0.843, 0.584, 0.142, 0.021, 0.23, 0.213, 0.216, 0.524, 0.839, 0.468, 0.81, 0.397, 0.009, 0.196, 0.363, 0.028, 0.302, 0.455, 0.665, 0.71, 0.683, 0.694, 0.348, 0.586, 0.953, 0.42, 0.329, 0.347, 0.405, 0.561, 0.043, 0.42, 0.7, 0.783, 0.873, 0.078, 0.045, 0.144, 0.674, 0.785, 0.615, 0.707, 0.035, 0.184, 0.826, 0.616, 0.05, 0.476, 0.889, 0.62, 0.647, 0.872, 0.227, 0.013, 0.319, 0.529, 0.597, 0.605, 0.503, 0.796, 0.633, 0.085, 0.382, 0.089, 0.622, 0.272, 0.804, 0.061, 0.276, 0.728, 0.706, 0.447, 0.415, 0.234, 0.813]
global q = [0.69, 0.935, 0.234, 0.931, 0.939, 0.872, 0.813, 0.902, 0.984, 0.609, 0.29, 0.622, 0.89, 0.742, 0.704, 0.81, 0.94, 0.97, 0.994, 0.992, 0.975, 0.737, 0.76, 0.469, 0.959, 0.821, 0.953, 0.952, 0.945, 0.263, 0.744, 0.955, 0.609, 0.526, 0.821, 0.954, 0.335, 0.942, 0.722, 0.532, 0.61, 0.704, 0.554, 0.88, 0.754, 0.851, 0.999, 0.874, 0.825, 0.941, 0.995, 0.876, 0.93, 0.681, 0.585, 0.94, 0.852, 0.988, 0.831, 0.999, 0.921, 0.993, 0.929, 0.924, 0.573, 0.98, 0.999, 0.329, 0.555, 0.335, 0.836, 0.712, 0.968, 0.943, 0.638, 0.969, 0.682, 0.857, 0.579, 0.914, 0.282, 0.192, 0.951, 0.678, 0.644, 0.955, 0.551, 0.529, 0.783, 0.906, 0.829, 0.302, 0.529, 0.863, 0.842, 0.14, 0.548, 0.801, 0.97, 0.878, 0.817, 0.987, 0.596, 0.948, 0.412, 0.701, 0.57, 0.563, 0.858, 0.962, 0.506, 0.607, 0.656, 0.989, 0.973, 0.87, 0.974, 0.836, 0.911, 0.694, 0.731, 0.854, 0.924, 0.996, 0.871, 0.487, 0.979, 0.792, 0.348, 0.737, 0.623, 0.856, 0.465, 0.886, 0.5, 0.787, 0.813, 0.688, 0.979, 0.825, 0.536, 0.793, 0.696, 0.986, 0.984, 0.72, 0.955, 0.934, 0.731, 0.576, 0.999, 0.417, 0.749, 0.984, 0.537, 0.983, 0.777, 0.777, 0.847, 0.211, 0.542, 0.511, 0.98, 0.844, 0.893, 0.68, 0.569, 0.956, 0.605, 0.426, 0.599, 0.879, 0.791, 0.372, 0.819, 0.838, 0.806, 0.78, 0.962, 0.55, 0.9, 0.67, 0.509, 0.987, 0.777, 0.537, 0.868, 0.986, 0.728, 0.943, 0.961, 0.473, 0.605, 0.45, 0.466, 0.818, 0.787, 0.916, 0.264, 0.756, 0.919, 0.883, 0.627, 0.879, 0.89, 0.405, 0.744, 0.963, 0.962, 0.926, 0.915, 0.974, 0.712, 0.444, 0.351, 0.543, 0.879, 0.624, 0.926, 0.997, 0.694, 0.781, 0.941, 0.521, 0.569, 0.184, 0.782, 0.831, 0.392, 0.694, 0.982, 0.337, 0.761, 0.735, 0.749, 0.826, 0.839, 0.678, 0.496, 0.661, 0.901, 0.32, 0.617, 0.535, 0.945, 0.934, 0.383, 0.835, 0.937, 0.844, 0.996, 0.698, 0.8, 0.89, 0.8, 0.903, 0.778, 0.955, 0.83, 0.687, 0.999, 0.95, 0.985, 0.655, 0.931, 0.821, 0.472, 0.929, 0.766, 0.983, 0.917, 0.651, 0.252, 0.79, 0.626, 0.889, 0.645, 0.996, 0.902, 0.619, 0.524, 0.391, 0.49, 0.616, 0.481, 0.969, 0.996, 0.748, 0.715, 0.544, 0.827, 0.996, 0.814, 0.822, 0.569, 0.436, 0.944, 0.697, 0.574, 0.794, 0.858, 0.901, 0.119, 0.559, 0.821, 0.725, 0.867, 0.756, 0.909, 0.371, 0.348, 0.941, 0.718, 0.418, 0.68, 0.995, 0.681, 0.856, 0.934, 0.82, 0.055, 0.609, 0.631, 0.912, 0.726, 0.506, 0.925, 0.895, 0.137, 0.665, 0.628, 0.979, 0.517, 0.853, 0.565, 0.563, 0.894, 0.961, 0.652, 0.729, 0.504, 0.859]
global origin = 1
global destination = 60