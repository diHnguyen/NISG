global arcs = [1 15; 1 23; 1 32; 1 46; 1 52; 2 40; 2 52; 3 19; 3 45; 4 8; 4 30; 4 32; 4 33; 4 39; 4 52; 5 7; 5 12; 5 35; 5 46; 6 13; 6 22; 6 33; 6 37; 6 57; 7 5; 8 2; 8 3; 8 4; 8 7; 8 11; 8 13; 8 33; 8 35; 8 38; 8 42; 8 52; 8 55; 8 58; 9 25; 10 7; 10 12; 10 19; 10 22; 10 28; 10 33; 10 42; 10 53; 10 60; 11 8; 11 23; 11 27; 11 32; 11 34; 12 5; 12 16; 12 24; 12 33; 12 40; 12 41; 13 11; 13 17; 13 33; 13 41; 13 42; 13 56; 14 19; 14 29; 14 36; 14 43; 14 46; 15 9; 15 21; 15 24; 15 49; 15 50; 15 53; 16 11; 16 48; 16 49; 16 59; 17 18; 17 23; 17 33; 17 37; 17 38; 17 41; 17 44; 17 54; 18 6; 19 6; 19 8; 19 15; 19 16; 19 22; 19 30; 19 36; 19 38; 19 40; 19 41; 19 43; 19 45; 19 46; 19 50; 19 58; 20 14; 20 19; 20 29; 20 34; 20 35; 20 38; 20 42; 20 49; 20 58; 21 4; 21 11; 21 16; 21 18; 21 33; 21 54; 22 13; 22 16; 22 34; 22 47; 22 51; 22 58; 23 21; 23 25; 23 36; 24 6; 24 15; 24 30; 24 44; 25 5; 25 10; 25 17; 25 21; 25 24; 25 37; 25 41; 25 42; 25 51; 25 57; 26 2; 26 9; 26 11; 26 18; 26 56; 27 5; 27 14; 27 16; 27 17; 27 26; 27 30; 27 31; 27 41; 27 43; 27 57; 28 38; 28 60; 29 10; 29 12; 29 25; 29 44; 29 52; 29 58; 30 6; 31 10; 31 23; 31 25; 31 30; 31 42; 31 49; 31 53; 32 2; 32 23; 32 39; 32 56; 32 57; 32 58; 33 11; 33 12; 33 21; 33 24; 34 3; 34 6; 34 14; 34 18; 34 21; 34 31; 34 39; 34 42; 34 46; 34 50; 35 8; 35 19; 35 33; 35 39; 35 45; 36 10; 36 12; 36 21; 36 23; 37 9; 37 13; 37 27; 37 38; 37 50; 37 54; 38 22; 38 23; 38 26; 38 35; 38 46; 38 48; 39 8; 39 12; 39 18; 39 27; 39 40; 40 18; 40 24; 40 36; 40 37; 40 58; 41 6; 41 20; 41 48; 41 49; 41 50; 42 11; 42 14; 42 15; 42 23; 42 27; 42 43; 42 44; 42 47; 42 58; 43 28; 43 51; 43 59; 44 4; 44 7; 44 10; 44 13; 44 31; 44 54; 44 55; 44 56; 44 58; 44 60; 45 4; 45 5; 45 48; 45 54; 46 33; 46 36; 46 44; 46 60; 47 24; 47 27; 47 32; 47 40; 47 56; 48 33; 48 40; 48 55; 49 17; 49 30; 49 56; 50 4; 50 28; 50 41; 50 46; 50 47; 50 49; 51 4; 51 9; 51 13; 51 15; 51 24; 51 33; 51 45; 51 50; 51 52; 52 13; 52 16; 52 27; 52 41; 52 45; 53 2; 53 27; 53 35; 53 38; 53 49; 53 54; 54 2; 54 5; 54 13; 54 16; 54 23; 54 32; 54 49; 54 58; 55 12; 55 20; 55 31; 55 44; 55 47; 56 27; 56 42; 57 8; 57 36; 57 43; 57 56; 58 7; 58 15; 58 16; 58 23; 58 45; 58 48; 58 56; 59 5; 59 11; 59 15; 59 26; 59 46; 59 57; 59 60]
global d_x = [2.0, 5.0, 8.0, 8.0, 4.0, 8.0, 7.0, 4.0, 4.0, 2.0, 1.0, 7.0, 2.0, 10.0, 8.0, 9.0, 1.0, 9.0, 2.0, 7.0, 10.0, 7.0, 2.0, 4.0, 10.0, 1.0, 5.0, 2.0, 3.0, 6.0, 3.0, 1.0, 4.0, 7.0, 3.0, 3.0, 3.0, 4.0, 7.0, 3.0, 3.0, 7.0, 7.0, 9.0, 5.0, 7.0, 4.0, 2.0, 8.0, 8.0, 1.0, 7.0, 10.0, 5.0, 4.0, 5.0, 7.0, 2.0, 7.0, 9.0, 1.0, 10.0, 2.0, 10.0, 5.0, 1.0, 9.0, 4.0, 7.0, 7.0, 2.0, 4.0, 9.0, 2.0, 10.0, 7.0, 9.0, 5.0, 2.0, 2.0, 3.0, 1.0, 3.0, 9.0, 2.0, 10.0, 7.0, 1.0, 8.0, 5.0, 4.0, 9.0, 4.0, 5.0, 4.0, 7.0, 2.0, 2.0, 10.0, 9.0, 6.0, 4.0, 1.0, 1.0, 3.0, 4.0, 5.0, 4.0, 10.0, 9.0, 7.0, 7.0, 5.0, 5.0, 4.0, 6.0, 6.0, 9.0, 7.0, 7.0, 10.0, 2.0, 8.0, 1.0, 5.0, 6.0, 3.0, 1.0, 1.0, 8.0, 5.0, 10.0, 8.0, 1.0, 5.0, 10.0, 8.0, 4.0, 7.0, 6.0, 9.0, 8.0, 3.0, 1.0, 9.0, 6.0, 4.0, 1.0, 9.0, 6.0, 5.0, 4.0, 7.0, 4.0, 9.0, 6.0, 8.0, 3.0, 10.0, 1.0, 10.0, 4.0, 6.0, 4.0, 1.0, 3.0, 6.0, 7.0, 2.0, 6.0, 6.0, 4.0, 8.0, 8.0, 1.0, 4.0, 10.0, 6.0, 8.0, 2.0, 9.0, 5.0, 3.0, 5.0, 3.0, 2.0, 1.0, 7.0, 1.0, 3.0, 1.0, 2.0, 5.0, 3.0, 4.0, 5.0, 4.0, 9.0, 7.0, 4.0, 7.0, 9.0, 4.0, 3.0, 10.0, 5.0, 2.0, 6.0, 2.0, 3.0, 9.0, 5.0, 6.0, 8.0, 5.0, 4.0, 9.0, 6.0, 1.0, 7.0, 7.0, 7.0, 5.0, 7.0, 3.0, 7.0, 9.0, 5.0, 3.0, 2.0, 3.0, 9.0, 10.0, 7.0, 6.0, 3.0, 9.0, 4.0, 6.0, 4.0, 10.0, 7.0, 8.0, 9.0, 10.0, 10.0, 9.0, 4.0, 9.0, 5.0, 6.0, 4.0, 2.0, 6.0, 3.0, 4.0, 9.0, 10.0, 7.0, 10.0, 8.0, 9.0, 1.0, 2.0, 10.0, 7.0, 2.0, 4.0, 1.0, 7.0, 4.0, 2.0, 1.0, 7.0, 4.0, 1.0, 9.0, 4.0, 6.0, 8.0, 5.0, 9.0, 2.0, 4.0, 5.0, 1.0, 9.0, 4.0, 8.0, 3.0, 5.0, 3.0, 7.0, 3.0, 1.0, 5.0, 3.0, 2.0, 9.0, 1.0, 7.0, 4.0, 8.0, 3.0, 6.0, 8.0, 8.0, 10.0, 3.0, 8.0, 1.0, 10.0, 2.0, 6.0, 8.0, 5.0, 10.0, 6.0, 2.0, 10.0, 9.0, 10.0, 3.0, 2.0, 2.0, 7.0, 2.0, 9.0, 1.0]
global b_x = 5
global d_y = [7.0, 9.0, 1.0, 7.0, 3.0, 10.0, 4.0, 5.0, 7.0, 7.0, 1.0, 7.0, 7.0, 7.0, 1.0, 9.0, 6.0, 9.0, 1.0, 2.0, 8.0, 10.0, 6.0, 9.0, 10.0, 10.0, 1.0, 6.0, 10.0, 4.0, 6.0, 6.0, 6.0, 1.0, 5.0, 10.0, 1.0, 3.0, 1.0, 7.0, 9.0, 1.0, 9.0, 1.0, 2.0, 1.0, 10.0, 2.0, 7.0, 6.0, 7.0, 4.0, 5.0, 8.0, 3.0, 7.0, 5.0, 3.0, 6.0, 4.0, 3.0, 5.0, 3.0, 8.0, 2.0, 10.0, 2.0, 5.0, 3.0, 8.0, 3.0, 4.0, 10.0, 2.0, 9.0, 10.0, 5.0, 4.0, 2.0, 8.0, 4.0, 8.0, 4.0, 7.0, 10.0, 9.0, 10.0, 9.0, 2.0, 7.0, 6.0, 1.0, 10.0, 2.0, 8.0, 10.0, 5.0, 10.0, 4.0, 9.0, 2.0, 5.0, 10.0, 3.0, 4.0, 4.0, 7.0, 3.0, 10.0, 6.0, 1.0, 5.0, 6.0, 10.0, 1.0, 1.0, 4.0, 1.0, 8.0, 6.0, 1.0, 1.0, 9.0, 1.0, 5.0, 1.0, 9.0, 1.0, 9.0, 6.0, 10.0, 7.0, 5.0, 3.0, 1.0, 1.0, 6.0, 7.0, 8.0, 9.0, 3.0, 8.0, 6.0, 9.0, 8.0, 7.0, 4.0, 10.0, 6.0, 6.0, 10.0, 6.0, 4.0, 9.0, 7.0, 1.0, 8.0, 5.0, 4.0, 8.0, 1.0, 9.0, 2.0, 3.0, 3.0, 10.0, 7.0, 7.0, 9.0, 10.0, 8.0, 1.0, 6.0, 9.0, 10.0, 1.0, 8.0, 2.0, 5.0, 4.0, 10.0, 1.0, 3.0, 9.0, 9.0, 1.0, 8.0, 4.0, 8.0, 8.0, 6.0, 8.0, 4.0, 4.0, 7.0, 1.0, 1.0, 1.0, 6.0, 8.0, 5.0, 10.0, 2.0, 9.0, 6.0, 5.0, 3.0, 1.0, 2.0, 10.0, 2.0, 2.0, 10.0, 10.0, 2.0, 7.0, 3.0, 5.0, 7.0, 2.0, 8.0, 5.0, 9.0, 6.0, 8.0, 2.0, 8.0, 3.0, 4.0, 9.0, 6.0, 3.0, 3.0, 8.0, 8.0, 9.0, 2.0, 1.0, 5.0, 8.0, 9.0, 9.0, 2.0, 5.0, 1.0, 3.0, 7.0, 3.0, 4.0, 9.0, 5.0, 6.0, 10.0, 3.0, 6.0, 9.0, 4.0, 5.0, 4.0, 10.0, 10.0, 8.0, 4.0, 1.0, 4.0, 8.0, 1.0, 3.0, 4.0, 4.0, 8.0, 3.0, 2.0, 5.0, 6.0, 6.0, 5.0, 8.0, 9.0, 2.0, 5.0, 7.0, 9.0, 3.0, 8.0, 4.0, 7.0, 5.0, 7.0, 4.0, 5.0, 5.0, 5.0, 1.0, 5.0, 4.0, 10.0, 7.0, 6.0, 2.0, 8.0, 6.0, 9.0, 8.0, 2.0, 3.0, 9.0, 5.0, 8.0, 1.0, 7.0, 8.0, 6.0, 3.0, 5.0, 1.0, 2.0, 2.0, 6.0, 5.0, 1.0, 7.0, 5.0, 6.0, 1.0, 8.0, 2.0, 1.0, 1.0]
global b_y = 10
global p = [0.754, 0.468, 0.701, 0.109, 0.321, 0.646, 0.24, 0.272, 0.045, 0.965, 0.22, 0.512, 0.012, 0.268, 0.416, 0.763, 0.542, 0.453, 0.404, 0.587, 0.322, 0.862, 0.001, 0.814, 0.484, 0.142, 0.134, 0.625, 0.178, 0.807, 0.57, 0.03, 0.007, 0.163, 0.428, 0.618, 0.073, 0.759, 0.182, 0.641, 0.752, 0.474, 0.027, 0.953, 0.535, 0.991, 0.541, 0.968, 0.17, 0.354, 0.768, 0.222, 0.948, 0.26, 0.444, 0.702, 0.38, 0.686, 0.073, 0.071, 0.013, 0.195, 0.648, 0.02, 0.848, 0.355, 0.02, 0.172, 0.844, 0.135, 0.392, 0.603, 0.272, 0.416, 0.636, 0.739, 0.223, 0.597, 0.413, 0.109, 0.086, 0.222, 0.507, 0.679, 0.421, 0.002, 0.213, 0.042, 0.141, 0.21, 0.786, 0.844, 0.736, 0.135, 0.273, 0.713, 0.905, 0.311, 0.103, 0.783, 0.391, 0.434, 0.765, 0.506, 0.077, 0.733, 0.032, 0.04, 0.636, 0.292, 0.249, 0.28, 0.208, 0.948, 0.877, 0.08, 0.574, 0.355, 0.41, 0.314, 0.095, 0.112, 0.791, 0.017, 0.809, 0.55, 0.334, 0.813, 0.091, 0.69, 0.454, 0.464, 0.519, 0.195, 0.842, 0.955, 0.565, 0.883, 0.063, 0.707, 0.301, 0.413, 0.994, 0.248, 0.452, 0.463, 0.469, 0.851, 0.479, 0.708, 0.534, 0.09, 0.899, 0.288, 0.003, 0.472, 0.846, 0.752, 0.927, 0.474, 0.619, 0.046, 0.955, 0.49, 0.016, 0.439, 0.616, 0.51, 0.929, 0.929, 0.317, 0.19, 0.732, 0.116, 0.68, 0.448, 0.167, 0.512, 0.626, 0.835, 0.238, 0.759, 0.165, 0.569, 0.844, 0.764, 0.671, 0.536, 0.289, 0.533, 0.335, 0.147, 0.313, 0.876, 0.246, 0.99, 0.997, 0.655, 0.763, 0.687, 0.81, 0.974, 0.669, 0.769, 0.98, 0.317, 0.748, 0.045, 0.694, 0.543, 0.507, 0.219, 0.036, 0.034, 0.296, 0.389, 0.268, 0.447, 0.062, 0.11, 0.524, 0.144, 0.906, 0.273, 0.722, 0.888, 0.519, 0.879, 0.551, 0.771, 0.829, 0.237, 0.331, 0.831, 0.711, 0.954, 0.887, 0.482, 0.079, 0.315, 0.98, 0.067, 0.403, 0.47, 0.941, 0.562, 0.193, 0.22, 0.654, 0.889, 0.02, 0.623, 0.515, 0.305, 0.946, 0.249, 0.585, 0.986, 0.779, 0.545, 0.799, 0.383, 0.92, 0.747, 0.993, 0.101, 0.635, 0.571, 0.74, 0.398, 0.056, 0.042, 0.317, 0.307, 0.518, 0.617, 0.179, 0.808, 0.873, 0.433, 0.876, 0.914, 0.825, 0.241, 0.33, 0.255, 0.032, 0.297, 0.847, 0.888, 0.444, 0.369, 0.987, 0.77, 0.937, 0.291, 0.234, 0.564, 0.035, 0.686, 0.616, 0.677, 0.723, 0.918, 0.957, 0.376, 0.352, 0.891, 0.808, 0.078, 0.612, 0.507, 0.479, 0.019, 0.82, 0.88, 0.345, 0.731, 0.356, 0.129, 0.01, 0.283, 0.038, 0.224, 0.986, 0.447, 0.325, 0.003, 0.44]
global q = [0.929, 0.987, 0.775, 0.668, 0.548, 0.913, 0.246, 0.684, 0.267, 0.969, 0.307, 0.927, 0.722, 0.644, 0.865, 0.888, 0.785, 0.488, 0.546, 0.931, 0.416, 0.937, 0.057, 0.941, 0.621, 0.706, 0.893, 0.65, 0.689, 0.903, 0.957, 0.445, 0.039, 0.557, 0.889, 0.883, 0.345, 0.927, 0.843, 0.698, 0.902, 0.66, 0.568, 0.996, 0.753, 0.999, 0.546, 0.991, 0.568, 0.55, 0.856, 0.36, 0.972, 0.63, 0.839, 0.882, 0.382, 0.855, 0.252, 0.616, 0.513, 0.871, 0.728, 0.646, 0.882, 0.867, 0.74, 0.348, 0.915, 0.308, 0.606, 0.89, 0.678, 0.705, 0.901, 0.787, 0.684, 0.814, 0.553, 0.3, 0.772, 0.99, 0.935, 0.751, 0.775, 0.739, 0.465, 0.17, 0.266, 0.928, 0.899, 0.915, 0.878, 0.297, 0.327, 0.816, 0.915, 0.883, 0.536, 0.883, 0.396, 0.82, 0.921, 0.952, 0.395, 0.891, 0.222, 0.37, 0.853, 0.332, 0.807, 0.298, 0.849, 0.949, 0.88, 0.918, 0.937, 0.575, 0.597, 0.325, 0.529, 0.432, 0.957, 0.047, 0.977, 0.93, 0.543, 0.989, 0.179, 0.84, 0.573, 0.564, 0.954, 0.698, 0.933, 0.976, 0.792, 0.884, 0.204, 0.838, 0.401, 0.952, 0.996, 0.367, 0.972, 0.802, 0.66, 0.858, 0.479, 0.762, 0.849, 0.923, 0.904, 0.385, 0.226, 0.834, 0.911, 0.985, 0.999, 0.745, 0.723, 0.078, 0.974, 0.885, 0.909, 0.731, 0.842, 0.755, 0.969, 0.93, 0.639, 0.198, 0.824, 0.456, 0.731, 0.893, 0.362, 0.961, 0.857, 0.984, 0.7, 0.835, 0.732, 0.944, 0.859, 0.941, 0.797, 0.732, 0.942, 0.641, 0.945, 0.297, 0.95, 0.897, 0.275, 0.997, 0.999, 0.969, 0.846, 0.811, 0.819, 0.98, 0.719, 0.873, 0.997, 0.555, 0.943, 0.47, 0.782, 0.574, 0.785, 0.774, 0.922, 0.055, 0.444, 0.729, 0.524, 0.982, 0.814, 0.347, 0.569, 0.44, 0.953, 0.4, 0.92, 0.976, 0.76, 0.979, 0.975, 0.963, 0.926, 0.392, 0.458, 0.957, 0.915, 0.976, 0.906, 0.819, 0.235, 0.989, 0.989, 0.989, 0.76, 0.541, 0.967, 0.565, 0.729, 0.844, 0.948, 0.945, 0.446, 0.988, 0.711, 0.842, 0.961, 0.848, 0.943, 0.996, 0.795, 0.837, 0.873, 0.448, 0.934, 0.99, 0.994, 0.12, 0.854, 0.683, 0.781, 0.618, 0.808, 0.302, 0.458, 0.988, 0.548, 0.89, 0.714, 0.994, 0.937, 0.944, 0.993, 0.943, 0.907, 0.841, 0.918, 0.441, 0.722, 0.426, 0.853, 0.902, 0.604, 0.664, 0.995, 0.83, 0.964, 0.475, 0.999, 0.659, 0.993, 0.753, 0.681, 0.762, 0.795, 0.928, 0.978, 0.935, 0.646, 0.949, 0.901, 0.808, 0.767, 0.543, 0.809, 0.437, 0.954, 0.907, 0.604, 0.973, 0.984, 0.334, 0.921, 0.956, 0.28, 0.356, 0.997, 0.965, 0.857, 0.814, 0.725]
global origin = 1
global destination = 60