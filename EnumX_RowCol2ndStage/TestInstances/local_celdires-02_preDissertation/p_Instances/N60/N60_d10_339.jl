global arcs = [1 13; 1 42; 1 47; 2 4; 2 6; 2 9; 2 15; 2 26; 2 28; 2 35; 2 46; 2 48; 2 55; 3 14; 3 16; 3 35; 3 45; 4 21; 4 29; 4 34; 4 39; 4 43; 4 44; 5 8; 5 18; 5 24; 5 30; 5 56; 5 58; 5 60; 6 21; 6 23; 6 25; 6 35; 6 37; 6 38; 6 53; 7 12; 7 16; 7 21; 7 22; 7 23; 7 26; 7 32; 7 53; 7 55; 7 58; 8 17; 8 20; 8 22; 8 37; 8 47; 8 55; 8 58; 9 13; 9 25; 9 30; 9 31; 9 35; 9 44; 10 5; 10 13; 10 21; 10 52; 10 54; 10 58; 11 18; 11 25; 11 47; 11 49; 11 54; 11 57; 12 4; 12 5; 12 6; 12 27; 12 30; 12 42; 12 53; 12 56; 13 2; 13 47; 13 56; 13 59; 14 5; 14 12; 14 20; 14 24; 14 27; 14 31; 14 53; 15 3; 15 7; 15 19; 16 7; 16 20; 17 6; 17 10; 17 23; 17 40; 17 41; 17 48; 17 49; 18 12; 18 19; 18 33; 18 34; 18 44; 18 48; 19 28; 19 32; 19 33; 19 52; 19 58; 20 3; 20 9; 20 10; 20 17; 20 23; 21 14; 21 22; 21 23; 21 27; 21 33; 21 35; 21 40; 22 33; 22 43; 23 26; 23 30; 23 48; 23 54; 24 9; 24 17; 24 19; 24 29; 24 46; 24 57; 25 21; 25 33; 25 38; 25 42; 25 45; 25 46; 25 52; 25 55; 26 11; 26 16; 26 40; 27 4; 27 28; 27 40; 28 12; 28 13; 28 21; 28 22; 28 36; 28 37; 28 38; 29 4; 29 9; 29 15; 29 18; 29 25; 29 31; 29 33; 29 45; 30 5; 30 14; 31 6; 31 18; 31 24; 31 32; 31 39; 32 21; 32 29; 33 2; 33 5; 33 9; 33 39; 33 43; 33 44; 33 58; 33 59; 33 60; 34 8; 34 13; 35 3; 35 9; 35 31; 35 36; 35 57; 36 2; 36 4; 36 24; 36 31; 36 33; 36 39; 36 40; 36 46; 36 51; 36 53; 36 57; 37 2; 37 17; 37 23; 37 48; 37 50; 37 54; 37 57; 38 3; 38 4; 38 7; 38 13; 38 14; 38 20; 38 25; 38 37; 38 39; 39 6; 39 10; 39 21; 39 44; 40 37; 41 25; 41 28; 41 38; 41 50; 42 24; 42 35; 42 39; 42 44; 42 46; 43 4; 43 8; 43 28; 43 35; 43 36; 43 37; 43 39; 43 46; 43 55; 44 4; 44 21; 44 25; 44 27; 44 35; 44 41; 44 46; 44 50; 45 3; 45 17; 45 26; 45 29; 45 31; 45 41; 46 3; 46 28; 46 36; 46 37; 47 3; 47 6; 47 9; 47 29; 47 33; 47 36; 48 8; 48 32; 48 44; 48 54; 48 57; 49 33; 49 34; 49 35; 49 44; 49 57; 49 59; 49 60; 50 15; 50 31; 50 48; 50 51; 50 52; 50 54; 51 15; 51 19; 51 35; 51 39; 52 10; 52 15; 52 23; 52 28; 52 30; 52 43; 53 2; 53 3; 53 21; 53 33; 53 45; 53 47; 54 2; 54 3; 54 12; 54 22; 54 29; 54 33; 54 58; 55 17; 55 38; 55 57; 55 60; 56 11; 56 35; 56 52; 56 53; 56 59; 57 2; 57 8; 57 43; 57 59; 58 52; 59 8; 59 26; 59 38; 59 39; 59 43; 59 53; 59 56; 59 58]
global d_x = [3.0, 8.0, 9.0, 8.0, 8.0, 6.0, 1.0, 4.0, 3.0, 9.0, 4.0, 4.0, 9.0, 4.0, 2.0, 10.0, 1.0, 5.0, 8.0, 10.0, 7.0, 5.0, 4.0, 1.0, 6.0, 5.0, 6.0, 9.0, 10.0, 7.0, 7.0, 2.0, 8.0, 1.0, 1.0, 9.0, 6.0, 6.0, 10.0, 1.0, 4.0, 3.0, 9.0, 3.0, 9.0, 1.0, 7.0, 4.0, 3.0, 2.0, 6.0, 9.0, 1.0, 5.0, 5.0, 8.0, 3.0, 5.0, 10.0, 7.0, 4.0, 7.0, 8.0, 1.0, 1.0, 9.0, 6.0, 9.0, 7.0, 6.0, 5.0, 9.0, 8.0, 3.0, 8.0, 4.0, 5.0, 4.0, 6.0, 2.0, 10.0, 6.0, 6.0, 1.0, 3.0, 8.0, 10.0, 2.0, 8.0, 8.0, 6.0, 4.0, 1.0, 3.0, 3.0, 10.0, 8.0, 8.0, 8.0, 7.0, 1.0, 8.0, 7.0, 2.0, 10.0, 5.0, 6.0, 9.0, 10.0, 10.0, 2.0, 6.0, 9.0, 10.0, 1.0, 2.0, 1.0, 9.0, 1.0, 1.0, 2.0, 10.0, 10.0, 1.0, 9.0, 2.0, 4.0, 10.0, 5.0, 3.0, 8.0, 1.0, 9.0, 1.0, 10.0, 4.0, 4.0, 1.0, 8.0, 5.0, 8.0, 1.0, 4.0, 6.0, 9.0, 9.0, 6.0, 1.0, 1.0, 1.0, 4.0, 6.0, 6.0, 7.0, 7.0, 1.0, 2.0, 7.0, 5.0, 7.0, 10.0, 7.0, 9.0, 2.0, 5.0, 8.0, 7.0, 9.0, 7.0, 7.0, 5.0, 4.0, 10.0, 9.0, 6.0, 10.0, 2.0, 8.0, 2.0, 3.0, 9.0, 6.0, 10.0, 4.0, 2.0, 9.0, 5.0, 8.0, 5.0, 3.0, 1.0, 1.0, 10.0, 5.0, 10.0, 6.0, 2.0, 6.0, 7.0, 6.0, 3.0, 2.0, 5.0, 9.0, 2.0, 2.0, 6.0, 8.0, 1.0, 10.0, 8.0, 6.0, 2.0, 3.0, 2.0, 8.0, 9.0, 6.0, 1.0, 3.0, 7.0, 3.0, 6.0, 8.0, 10.0, 1.0, 7.0, 1.0, 6.0, 2.0, 3.0, 3.0, 6.0, 10.0, 6.0, 10.0, 10.0, 10.0, 9.0, 2.0, 6.0, 8.0, 5.0, 7.0, 6.0, 4.0, 4.0, 7.0, 6.0, 10.0, 2.0, 10.0, 4.0, 4.0, 9.0, 2.0, 10.0, 6.0, 5.0, 5.0, 7.0, 4.0, 1.0, 3.0, 2.0, 9.0, 8.0, 1.0, 3.0, 4.0, 6.0, 6.0, 6.0, 6.0, 2.0, 5.0, 6.0, 8.0, 9.0, 8.0, 7.0, 1.0, 3.0, 9.0, 6.0, 5.0, 4.0, 7.0, 10.0, 6.0, 1.0, 3.0, 6.0, 6.0, 3.0, 3.0, 7.0, 6.0, 1.0, 10.0, 5.0, 4.0, 4.0, 6.0, 8.0, 1.0, 9.0, 2.0, 7.0, 6.0, 9.0, 7.0, 8.0, 9.0, 4.0, 1.0, 5.0, 1.0, 6.0, 8.0, 3.0, 10.0, 8.0, 2.0, 5.0, 2.0, 5.0, 1.0, 2.0]
global b_x = 5
global d_y = [5.0, 2.0, 2.0, 8.0, 5.0, 2.0, 5.0, 3.0, 3.0, 6.0, 5.0, 1.0, 5.0, 7.0, 3.0, 5.0, 2.0, 10.0, 3.0, 1.0, 9.0, 5.0, 6.0, 1.0, 4.0, 1.0, 7.0, 8.0, 2.0, 3.0, 2.0, 8.0, 9.0, 6.0, 9.0, 8.0, 6.0, 3.0, 5.0, 1.0, 1.0, 4.0, 1.0, 4.0, 3.0, 2.0, 5.0, 9.0, 2.0, 3.0, 8.0, 9.0, 3.0, 5.0, 10.0, 9.0, 9.0, 6.0, 8.0, 9.0, 1.0, 9.0, 8.0, 8.0, 9.0, 9.0, 8.0, 6.0, 6.0, 3.0, 6.0, 2.0, 7.0, 1.0, 3.0, 3.0, 9.0, 3.0, 3.0, 7.0, 8.0, 3.0, 9.0, 9.0, 10.0, 6.0, 7.0, 6.0, 2.0, 6.0, 1.0, 6.0, 1.0, 8.0, 3.0, 9.0, 3.0, 3.0, 1.0, 9.0, 6.0, 9.0, 8.0, 3.0, 9.0, 2.0, 2.0, 4.0, 2.0, 4.0, 2.0, 6.0, 8.0, 10.0, 5.0, 7.0, 9.0, 6.0, 6.0, 9.0, 3.0, 3.0, 10.0, 8.0, 1.0, 10.0, 8.0, 8.0, 1.0, 4.0, 4.0, 7.0, 4.0, 6.0, 2.0, 1.0, 5.0, 6.0, 7.0, 10.0, 10.0, 8.0, 8.0, 4.0, 1.0, 3.0, 9.0, 8.0, 1.0, 2.0, 1.0, 4.0, 10.0, 3.0, 2.0, 2.0, 7.0, 10.0, 1.0, 2.0, 1.0, 10.0, 10.0, 4.0, 2.0, 4.0, 9.0, 3.0, 9.0, 1.0, 5.0, 9.0, 5.0, 10.0, 9.0, 8.0, 1.0, 10.0, 9.0, 3.0, 7.0, 7.0, 6.0, 8.0, 10.0, 9.0, 2.0, 5.0, 6.0, 8.0, 8.0, 2.0, 1.0, 9.0, 10.0, 5.0, 5.0, 10.0, 2.0, 2.0, 8.0, 3.0, 6.0, 7.0, 7.0, 10.0, 7.0, 6.0, 5.0, 10.0, 3.0, 10.0, 1.0, 9.0, 10.0, 9.0, 1.0, 7.0, 5.0, 9.0, 10.0, 10.0, 1.0, 7.0, 8.0, 9.0, 3.0, 10.0, 10.0, 8.0, 2.0, 9.0, 2.0, 7.0, 10.0, 1.0, 2.0, 10.0, 2.0, 9.0, 3.0, 10.0, 3.0, 9.0, 9.0, 7.0, 3.0, 5.0, 7.0, 9.0, 5.0, 6.0, 3.0, 8.0, 1.0, 6.0, 6.0, 6.0, 5.0, 7.0, 10.0, 6.0, 4.0, 1.0, 8.0, 2.0, 5.0, 9.0, 9.0, 2.0, 10.0, 8.0, 7.0, 1.0, 4.0, 9.0, 6.0, 1.0, 1.0, 6.0, 2.0, 6.0, 2.0, 9.0, 5.0, 2.0, 5.0, 2.0, 4.0, 3.0, 9.0, 7.0, 9.0, 6.0, 1.0, 1.0, 4.0, 9.0, 10.0, 10.0, 9.0, 6.0, 9.0, 10.0, 8.0, 10.0, 9.0, 6.0, 1.0, 6.0, 7.0, 6.0, 5.0, 4.0, 7.0, 6.0, 7.0, 2.0, 10.0, 4.0, 8.0, 2.0, 2.0, 2.0, 6.0, 3.0, 9.0, 9.0, 4.0]
global b_y = 10
global p = [0.956, 0.322, 0.078, 0.489, 0.602, 0.548, 0.145, 0.936, 0.855, 0.54, 0.882, 0.427, 0.468, 0.146, 0.646, 0.12, 0.918, 0.158, 0.641, 0.756, 0.432, 0.834, 0.672, 0.357, 0.484, 0.529, 0.808, 0.688, 0.761, 0.531, 0.577, 0.657, 0.42, 0.99, 0.089, 0.172, 0.131, 0.892, 0.187, 0.01, 0.271, 0.208, 0.691, 0.474, 0.974, 0.867, 0.553, 0.467, 0.212, 0.084, 0.437, 0.027, 0.209, 0.531, 0.684, 0.948, 0.317, 0.942, 0.844, 0.7, 0.767, 0.122, 0.53, 0.643, 0.413, 0.984, 0.284, 0.311, 0.958, 0.14, 0.061, 0.376, 0.88, 0.884, 0.288, 0.459, 0.808, 0.79, 0.519, 0.48, 0.873, 0.602, 0.008, 0.627, 0.421, 0.586, 0.945, 0.264, 0.985, 0.872, 0.774, 0.125, 0.135, 0.54, 0.549, 0.768, 0.758, 0.17, 0.114, 0.89, 0.591, 0.309, 0.757, 0.963, 0.602, 0.738, 0.704, 0.914, 0.996, 0.914, 0.263, 0.48, 0.264, 0.145, 0.974, 0.017, 0.317, 0.841, 0.022, 0.944, 0.095, 0.573, 0.345, 0.369, 0.725, 0.021, 0.945, 0.447, 0.102, 0.8, 0.503, 0.502, 0.127, 0.253, 0.788, 0.339, 0.898, 0.32, 0.329, 0.688, 0.853, 0.944, 0.479, 0.918, 0.669, 0.688, 0.433, 0.792, 0.813, 0.292, 0.636, 0.429, 0.709, 0.52, 0.189, 0.956, 0.625, 0.997, 0.01, 0.119, 0.481, 0.8, 0.39, 0.822, 0.963, 0.501, 0.385, 0.712, 0.285, 0.192, 0.955, 0.239, 0.378, 0.359, 0.109, 0.76, 0.136, 0.606, 0.359, 0.493, 0.456, 0.798, 0.662, 0.168, 0.614, 0.662, 0.39, 0.568, 0.882, 0.755, 0.615, 0.93, 0.611, 0.815, 0.281, 0.921, 0.95, 0.229, 0.748, 0.933, 0.897, 0.135, 0.952, 0.896, 0.553, 0.825, 0.018, 0.717, 0.901, 0.439, 0.362, 0.881, 0.709, 0.436, 0.309, 0.264, 0.777, 0.445, 0.414, 0.091, 0.12, 0.142, 0.739, 0.49, 0.252, 0.422, 0.656, 0.593, 0.531, 0.73, 0.702, 0.846, 0.176, 0.575, 0.499, 0.145, 0.947, 0.66, 0.768, 0.433, 0.19, 0.189, 0.575, 0.754, 0.222, 0.991, 0.079, 0.046, 0.099, 0.775, 0.504, 0.911, 0.734, 0.651, 0.987, 0.693, 0.81, 0.617, 0.936, 0.135, 0.772, 0.82, 0.925, 0.559, 0.741, 0.954, 0.547, 0.189, 0.716, 0.5, 0.583, 0.965, 0.655, 0.757, 0.427, 0.207, 0.981, 0.568, 0.874, 0.636, 0.326, 0.927, 0.986, 0.466, 0.804, 0.679, 0.401, 0.871, 0.305, 0.771, 0.718, 0.656, 0.136, 0.752, 0.635, 0.459, 0.78, 0.194, 0.137, 0.187, 0.817, 0.261, 0.01, 0.858, 0.219, 0.405, 0.278, 0.158, 0.362, 0.033, 0.81, 0.919, 0.966, 0.403, 0.359, 0.68, 0.502, 0.942, 0.334, 0.488, 0.029, 0.53, 0.045, 0.771, 0.063, 0.548, 0.643, 0.42, 0.592]
global q = [0.982, 0.953, 0.674, 0.798, 0.787, 0.837, 0.874, 0.945, 0.889, 0.616, 0.906, 0.744, 0.935, 0.464, 0.747, 0.365, 0.942, 0.409, 0.716, 0.85, 0.953, 0.888, 0.799, 0.913, 0.913, 0.53, 0.871, 0.847, 0.902, 0.952, 0.808, 0.952, 0.898, 0.992, 0.714, 0.666, 0.923, 0.986, 0.37, 0.624, 0.422, 0.945, 0.749, 0.993, 0.99, 0.872, 0.833, 0.922, 0.99, 0.206, 0.614, 0.881, 0.348, 0.848, 0.835, 0.997, 0.364, 0.96, 0.997, 0.773, 0.977, 0.739, 0.64, 0.889, 0.814, 0.984, 0.395, 0.854, 0.976, 0.808, 0.934, 0.695, 0.928, 0.945, 0.501, 0.997, 0.871, 0.919, 0.844, 0.999, 0.967, 0.723, 0.703, 0.855, 0.786, 0.667, 0.974, 0.925, 0.986, 0.917, 0.977, 0.768, 0.653, 0.883, 0.826, 0.912, 0.889, 0.326, 0.848, 0.896, 0.975, 0.325, 0.875, 0.974, 0.704, 0.774, 0.787, 0.995, 0.999, 0.975, 0.895, 0.93, 0.907, 0.152, 0.983, 0.133, 0.783, 0.842, 0.376, 0.954, 0.836, 0.737, 0.448, 0.519, 0.886, 0.074, 0.966, 0.867, 0.594, 0.861, 0.571, 0.515, 0.558, 0.976, 0.945, 0.556, 0.944, 0.765, 0.938, 0.764, 0.892, 0.999, 0.887, 0.992, 0.722, 0.796, 0.932, 0.912, 0.977, 0.973, 0.936, 0.657, 0.978, 0.58, 0.633, 0.959, 0.697, 0.997, 0.798, 0.529, 0.927, 0.822, 0.866, 0.981, 0.971, 0.962, 0.476, 0.812, 0.351, 0.207, 0.963, 0.526, 0.758, 0.442, 0.125, 0.77, 0.249, 0.953, 0.925, 0.758, 0.867, 0.82, 0.867, 0.709, 0.622, 0.819, 0.933, 0.897, 0.983, 0.975, 0.805, 0.949, 0.953, 0.858, 0.744, 0.94, 0.959, 0.238, 0.789, 0.952, 0.967, 0.474, 0.99, 0.909, 0.944, 0.897, 0.472, 0.875, 0.923, 0.541, 0.666, 0.98, 0.791, 0.506, 0.745, 0.429, 0.934, 0.594, 0.607, 0.838, 0.218, 0.163, 0.758, 0.878, 0.997, 0.941, 0.684, 0.601, 0.557, 0.862, 0.814, 0.881, 0.485, 0.746, 0.947, 0.209, 0.952, 0.929, 0.863, 0.815, 0.283, 0.576, 0.7, 0.988, 0.485, 0.997, 0.295, 0.308, 0.485, 0.981, 0.789, 0.936, 0.98, 0.818, 0.993, 0.846, 0.994, 0.674, 0.985, 0.694, 0.906, 0.873, 0.947, 0.747, 0.869, 0.976, 0.888, 0.733, 0.996, 0.918, 0.812, 0.981, 0.697, 0.868, 0.685, 0.525, 0.988, 0.623, 0.925, 0.747, 0.885, 0.952, 0.986, 0.647, 0.876, 0.856, 0.62, 0.986, 0.894, 0.909, 0.885, 0.746, 0.793, 0.855, 0.857, 0.705, 0.965, 0.797, 0.139, 0.534, 0.842, 0.67, 0.254, 0.912, 0.743, 0.802, 0.708, 0.697, 0.858, 0.235, 0.848, 0.929, 0.99, 0.403, 0.549, 0.825, 0.563, 0.959, 0.807, 0.841, 0.715, 0.967, 0.27, 0.833, 0.695, 0.887, 0.941, 0.599, 0.687]
global origin = 1
global destination = 60