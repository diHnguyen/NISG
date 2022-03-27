global arcs = [1 4; 1 7; 1 19; 1 22; 1 23; 1 27; 1 56; 1 57; 2 9; 2 26; 2 48; 3 14; 3 22; 3 34; 3 36; 3 39; 3 49; 3 51; 4 5; 4 21; 4 24; 4 29; 4 48; 5 6; 5 10; 5 22; 5 33; 5 39; 5 46; 6 7; 6 17; 6 18; 6 30; 6 31; 6 48; 7 11; 7 13; 7 54; 8 2; 8 14; 8 22; 8 43; 8 44; 9 36; 9 39; 9 54; 10 2; 10 14; 10 15; 10 37; 10 59; 11 9; 11 27; 12 11; 12 15; 12 36; 12 43; 12 50; 13 7; 13 10; 13 11; 13 22; 13 25; 13 27; 13 28; 13 35; 13 41; 13 42; 13 56; 13 59; 14 13; 14 31; 14 33; 14 36; 14 46; 14 57; 15 3; 15 5; 15 24; 15 35; 15 45; 15 55; 15 58; 16 27; 16 36; 16 39; 16 48; 16 54; 17 8; 17 19; 17 22; 17 30; 17 35; 17 39; 17 56; 18 12; 18 17; 18 25; 18 27; 18 43; 18 44; 18 54; 18 59; 19 4; 19 58; 20 12; 20 25; 20 41; 20 59; 21 13; 21 34; 21 36; 21 40; 21 43; 21 51; 22 11; 22 30; 22 39; 22 41; 22 47; 22 57; 22 59; 23 13; 23 25; 23 35; 23 54; 23 55; 23 57; 24 4; 24 10; 24 35; 24 37; 24 46; 24 48; 24 56; 25 4; 25 6; 25 10; 25 11; 25 17; 25 20; 25 48; 25 56; 25 57; 25 59; 26 3; 26 5; 26 7; 26 15; 26 49; 27 13; 27 17; 27 21; 27 37; 27 38; 27 50; 27 53; 27 57; 28 17; 28 19; 28 31; 28 38; 28 40; 28 41; 28 52; 28 56; 29 10; 29 18; 29 36; 29 56; 30 19; 30 39; 30 41; 30 49; 31 18; 31 22; 31 25; 31 37; 31 51; 32 3; 32 11; 32 22; 32 42; 32 54; 32 56; 33 9; 33 18; 33 19; 33 25; 33 31; 33 43; 33 59; 34 13; 35 16; 35 22; 35 45; 36 17; 36 30; 36 35; 36 37; 36 52; 36 53; 36 59; 37 23; 37 33; 37 43; 37 44; 37 49; 37 52; 38 7; 38 30; 38 34; 39 2; 39 9; 39 60; 40 2; 40 3; 40 4; 40 17; 40 24; 40 27; 40 51; 41 4; 41 23; 41 34; 41 44; 42 5; 42 14; 42 33; 42 39; 42 49; 43 3; 43 24; 43 38; 43 42; 43 58; 44 22; 44 25; 44 27; 44 33; 44 37; 44 49; 45 4; 45 11; 45 14; 45 55; 45 56; 46 2; 46 17; 46 18; 46 29; 46 44; 46 52; 46 55; 47 8; 47 12; 47 17; 47 26; 47 30; 47 34; 47 40; 47 41; 47 49; 47 55; 47 60; 48 12; 48 20; 48 57; 49 19; 49 27; 49 28; 50 8; 50 11; 50 22; 50 28; 50 41; 50 43; 50 49; 51 4; 51 12; 51 43; 51 48; 52 9; 52 16; 52 20; 52 32; 52 38; 53 6; 53 20; 53 25; 53 29; 54 14; 54 21; 54 34; 54 42; 55 5; 55 9; 55 14; 55 31; 56 9; 56 13; 56 24; 56 38; 56 43; 57 13; 57 16; 57 41; 57 53; 58 2; 58 25; 59 2; 59 16; 59 23; 59 48; 59 60]
global d_x = [3.0, 9.0, 1.0, 3.0, 8.0, 10.0, 7.0, 3.0, 9.0, 9.0, 6.0, 8.0, 10.0, 8.0, 4.0, 6.0, 8.0, 6.0, 4.0, 1.0, 7.0, 5.0, 5.0, 6.0, 7.0, 6.0, 2.0, 6.0, 7.0, 8.0, 8.0, 9.0, 4.0, 9.0, 10.0, 10.0, 9.0, 7.0, 9.0, 10.0, 1.0, 4.0, 1.0, 10.0, 1.0, 4.0, 7.0, 2.0, 10.0, 1.0, 8.0, 9.0, 2.0, 7.0, 7.0, 8.0, 3.0, 10.0, 2.0, 9.0, 4.0, 3.0, 3.0, 10.0, 5.0, 4.0, 1.0, 8.0, 5.0, 1.0, 7.0, 7.0, 10.0, 5.0, 6.0, 8.0, 10.0, 1.0, 1.0, 2.0, 9.0, 4.0, 10.0, 9.0, 4.0, 9.0, 1.0, 1.0, 4.0, 5.0, 1.0, 5.0, 5.0, 1.0, 4.0, 4.0, 8.0, 9.0, 1.0, 4.0, 2.0, 6.0, 5.0, 1.0, 8.0, 7.0, 7.0, 2.0, 7.0, 10.0, 9.0, 4.0, 7.0, 4.0, 1.0, 9.0, 9.0, 6.0, 5.0, 1.0, 8.0, 2.0, 2.0, 9.0, 10.0, 7.0, 8.0, 2.0, 7.0, 7.0, 4.0, 6.0, 6.0, 5.0, 10.0, 3.0, 7.0, 5.0, 9.0, 1.0, 10.0, 7.0, 5.0, 6.0, 3.0, 1.0, 5.0, 6.0, 9.0, 4.0, 7.0, 9.0, 8.0, 6.0, 6.0, 9.0, 8.0, 5.0, 9.0, 2.0, 2.0, 5.0, 2.0, 3.0, 6.0, 10.0, 9.0, 1.0, 5.0, 1.0, 7.0, 2.0, 6.0, 8.0, 3.0, 9.0, 8.0, 8.0, 3.0, 1.0, 3.0, 6.0, 9.0, 3.0, 4.0, 8.0, 3.0, 7.0, 5.0, 8.0, 8.0, 9.0, 2.0, 9.0, 7.0, 1.0, 9.0, 10.0, 7.0, 10.0, 6.0, 6.0, 9.0, 4.0, 6.0, 7.0, 6.0, 7.0, 10.0, 1.0, 4.0, 3.0, 4.0, 8.0, 3.0, 7.0, 6.0, 6.0, 3.0, 6.0, 1.0, 2.0, 8.0, 3.0, 7.0, 4.0, 6.0, 10.0, 6.0, 1.0, 1.0, 4.0, 5.0, 8.0, 8.0, 6.0, 7.0, 3.0, 10.0, 5.0, 6.0, 9.0, 10.0, 3.0, 5.0, 8.0, 6.0, 5.0, 5.0, 5.0, 6.0, 7.0, 4.0, 2.0, 9.0, 4.0, 8.0, 2.0, 9.0, 5.0, 6.0, 4.0, 3.0, 6.0, 8.0, 8.0, 5.0, 7.0, 5.0, 2.0, 3.0, 10.0, 2.0, 8.0, 8.0, 7.0, 9.0, 4.0, 2.0, 2.0, 9.0, 2.0, 3.0, 3.0, 10.0, 8.0, 7.0, 3.0, 4.0, 4.0, 4.0, 1.0, 2.0, 4.0, 10.0, 8.0, 6.0, 10.0, 8.0, 3.0, 7.0, 2.0, 10.0, 7.0, 4.0, 6.0, 8.0, 5.0, 9.0, 7.0, 4.0, 3.0, 4.0, 7.0, 4.0]
global b_x = 5
global d_y = [5.0, 10.0, 3.0, 5.0, 4.0, 8.0, 2.0, 2.0, 7.0, 8.0, 3.0, 2.0, 3.0, 5.0, 3.0, 8.0, 6.0, 8.0, 2.0, 7.0, 6.0, 2.0, 4.0, 10.0, 8.0, 9.0, 5.0, 1.0, 4.0, 2.0, 9.0, 5.0, 6.0, 6.0, 1.0, 10.0, 8.0, 4.0, 6.0, 8.0, 9.0, 2.0, 8.0, 7.0, 5.0, 3.0, 4.0, 1.0, 1.0, 5.0, 10.0, 3.0, 4.0, 7.0, 9.0, 5.0, 2.0, 4.0, 6.0, 8.0, 1.0, 3.0, 5.0, 3.0, 1.0, 10.0, 4.0, 3.0, 3.0, 5.0, 3.0, 7.0, 7.0, 2.0, 1.0, 1.0, 7.0, 4.0, 7.0, 8.0, 10.0, 8.0, 4.0, 1.0, 10.0, 4.0, 4.0, 1.0, 9.0, 2.0, 7.0, 10.0, 6.0, 8.0, 3.0, 5.0, 7.0, 10.0, 8.0, 9.0, 7.0, 2.0, 10.0, 6.0, 9.0, 8.0, 6.0, 8.0, 3.0, 6.0, 4.0, 1.0, 8.0, 2.0, 2.0, 4.0, 3.0, 5.0, 10.0, 2.0, 8.0, 7.0, 3.0, 7.0, 6.0, 2.0, 4.0, 5.0, 2.0, 3.0, 6.0, 3.0, 3.0, 3.0, 7.0, 9.0, 7.0, 1.0, 5.0, 6.0, 10.0, 6.0, 4.0, 5.0, 9.0, 7.0, 3.0, 8.0, 4.0, 7.0, 8.0, 4.0, 9.0, 1.0, 1.0, 1.0, 6.0, 6.0, 10.0, 9.0, 7.0, 1.0, 6.0, 8.0, 10.0, 2.0, 4.0, 6.0, 8.0, 5.0, 8.0, 9.0, 10.0, 4.0, 10.0, 3.0, 6.0, 1.0, 6.0, 9.0, 6.0, 5.0, 8.0, 10.0, 8.0, 1.0, 6.0, 2.0, 8.0, 4.0, 6.0, 3.0, 6.0, 3.0, 9.0, 7.0, 6.0, 8.0, 4.0, 3.0, 2.0, 10.0, 9.0, 9.0, 1.0, 10.0, 1.0, 4.0, 9.0, 5.0, 7.0, 10.0, 9.0, 10.0, 7.0, 7.0, 4.0, 1.0, 3.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0, 10.0, 4.0, 2.0, 5.0, 8.0, 3.0, 3.0, 5.0, 4.0, 4.0, 9.0, 4.0, 6.0, 8.0, 9.0, 5.0, 2.0, 10.0, 1.0, 8.0, 2.0, 8.0, 8.0, 6.0, 4.0, 1.0, 8.0, 10.0, 2.0, 4.0, 2.0, 5.0, 5.0, 8.0, 4.0, 1.0, 6.0, 1.0, 7.0, 9.0, 9.0, 8.0, 8.0, 9.0, 6.0, 1.0, 4.0, 8.0, 8.0, 4.0, 10.0, 7.0, 5.0, 1.0, 8.0, 3.0, 8.0, 9.0, 4.0, 10.0, 9.0, 2.0, 9.0, 6.0, 5.0, 3.0, 7.0, 5.0, 2.0, 7.0, 7.0, 9.0, 3.0, 4.0, 10.0, 7.0, 6.0, 4.0, 6.0, 9.0, 2.0, 6.0, 3.0, 2.0, 9.0, 8.0, 1.0, 6.0, 1.0, 4.0]
global b_y = 10
global p = [0.79, 0.827, 0.06, 0.292, 0.613, 0.49, 0.003, 0.605, 0.492, 0.665, 0.245, 0.589, 0.703, 0.33, 0.716, 0.868, 0.359, 0.556, 0.939, 0.118, 0.621, 0.349, 0.039, 0.689, 0.942, 0.701, 0.241, 0.138, 0.665, 0.987, 0.089, 0.746, 0.153, 0.444, 0.073, 0.359, 0.08, 0.794, 0.499, 0.841, 0.198, 0.64, 0.418, 0.495, 0.613, 0.637, 0.091, 0.881, 0.107, 0.437, 0.721, 0.382, 0.264, 0.196, 0.549, 0.454, 0.857, 0.957, 0.405, 0.799, 0.633, 0.749, 0.017, 0.559, 0.891, 0.004, 0.112, 0.908, 0.695, 0.979, 0.67, 0.264, 0.452, 0.769, 0.257, 0.894, 0.353, 0.138, 0.825, 0.704, 0.358, 0.818, 0.531, 0.827, 0.912, 0.235, 0.643, 0.966, 0.743, 0.133, 0.114, 0.499, 0.972, 0.752, 0.64, 0.269, 0.856, 0.211, 0.744, 0.594, 0.3, 0.105, 0.658, 0.054, 0.2, 0.234, 0.788, 0.586, 0.749, 0.93, 0.82, 0.188, 0.628, 0.644, 0.474, 0.524, 0.864, 0.842, 0.058, 0.598, 0.03, 0.239, 0.927, 0.215, 0.66, 0.266, 0.093, 0.492, 0.477, 0.914, 0.102, 0.435, 0.08, 0.808, 0.394, 0.742, 0.948, 0.17, 0.009, 0.161, 0.425, 0.017, 0.35, 0.584, 0.494, 0.113, 0.708, 0.685, 0.523, 0.06, 0.092, 0.117, 0.689, 0.418, 0.308, 0.627, 0.445, 0.895, 0.029, 0.885, 0.359, 0.592, 0.744, 0.026, 0.665, 0.328, 0.606, 0.206, 0.456, 0.282, 0.951, 0.998, 0.663, 0.263, 0.824, 0.785, 0.148, 0.867, 0.427, 0.753, 0.775, 0.169, 0.416, 0.299, 0.575, 0.161, 0.486, 0.591, 0.928, 0.125, 0.507, 0.056, 0.412, 0.784, 0.269, 0.153, 0.011, 0.078, 0.779, 0.036, 0.886, 0.336, 0.981, 0.372, 0.409, 0.374, 0.042, 0.401, 0.274, 0.446, 0.683, 0.157, 0.542, 0.378, 0.73, 0.792, 0.768, 0.775, 0.157, 0.278, 0.041, 0.786, 0.692, 0.494, 0.885, 0.564, 0.633, 0.855, 0.884, 0.318, 0.114, 0.906, 0.424, 0.167, 0.2, 0.826, 0.458, 0.145, 0.018, 0.303, 0.077, 0.798, 0.22, 0.996, 0.532, 0.93, 0.713, 0.855, 0.319, 0.741, 0.548, 0.897, 0.014, 0.136, 0.811, 0.446, 0.181, 0.696, 0.697, 0.612, 0.536, 0.545, 0.227, 0.934, 0.307, 0.436, 0.127, 0.602, 0.303, 0.124, 0.682, 0.246, 0.584, 0.999, 0.413, 0.009, 0.225, 0.083, 0.643, 0.546, 0.835, 0.325, 0.176, 0.106, 0.175, 0.477, 0.901, 0.631, 0.188, 0.51, 0.767, 0.806, 0.437, 0.578, 0.674, 0.267, 0.283, 0.657, 0.165, 0.738, 0.788, 0.206, 0.74, 0.193, 0.866, 0.252, 0.91, 0.731, 0.397, 0.74, 0.243, 0.14, 0.225, 0.641, 0.841]
global q = [0.799, 0.946, 0.883, 0.328, 0.635, 0.888, 0.965, 0.904, 0.511, 0.789, 0.364, 0.692, 0.961, 0.986, 0.765, 0.955, 0.69, 0.755, 0.996, 0.32, 0.88, 0.833, 0.875, 0.971, 0.951, 0.801, 0.352, 0.574, 0.945, 0.997, 0.598, 0.788, 0.499, 0.548, 0.552, 0.986, 0.529, 0.884, 0.688, 0.857, 0.475, 0.643, 0.695, 0.658, 0.743, 0.698, 0.911, 0.9, 0.485, 0.665, 0.996, 0.846, 0.56, 0.773, 0.676, 0.73, 0.904, 0.991, 0.578, 0.994, 0.745, 0.945, 0.896, 0.982, 0.934, 0.518, 0.882, 0.961, 0.999, 0.979, 0.99, 0.399, 0.883, 0.872, 0.425, 0.935, 0.383, 0.537, 0.916, 0.978, 0.815, 0.983, 0.841, 0.933, 0.971, 0.456, 0.965, 0.981, 0.983, 0.157, 0.279, 0.748, 0.997, 0.935, 0.815, 0.919, 0.891, 0.718, 0.892, 0.77, 0.589, 0.589, 0.737, 0.912, 0.802, 0.585, 0.836, 0.944, 0.922, 0.953, 0.979, 0.909, 0.942, 0.935, 0.557, 0.905, 0.978, 0.891, 0.998, 0.681, 0.225, 0.934, 0.94, 0.748, 0.976, 0.602, 0.467, 0.798, 0.791, 0.939, 0.479, 0.925, 0.476, 0.91, 0.756, 0.998, 0.994, 0.466, 0.78, 0.54, 0.773, 0.309, 0.474, 0.951, 0.532, 0.703, 0.726, 0.948, 0.748, 0.976, 0.592, 0.532, 0.837, 0.986, 0.646, 0.697, 0.744, 0.911, 0.408, 0.977, 0.868, 0.696, 0.918, 0.169, 0.722, 0.486, 0.843, 0.73, 0.733, 0.933, 0.953, 0.998, 0.875, 0.913, 0.942, 0.943, 0.912, 0.883, 0.692, 0.757, 0.98, 0.542, 0.473, 0.806, 0.759, 0.398, 0.665, 0.94, 0.935, 0.705, 0.574, 0.959, 0.936, 0.914, 0.932, 0.769, 0.103, 0.633, 0.841, 0.105, 0.975, 0.812, 0.995, 0.942, 0.738, 0.68, 0.874, 0.841, 0.521, 0.867, 0.905, 0.48, 0.638, 0.911, 0.865, 0.87, 0.797, 0.919, 0.779, 0.961, 0.338, 0.922, 0.97, 0.847, 0.91, 0.931, 0.829, 0.939, 0.957, 0.856, 0.737, 0.966, 0.886, 0.671, 0.449, 0.947, 0.824, 0.676, 0.106, 0.799, 0.481, 0.824, 0.473, 0.996, 0.668, 0.958, 0.932, 0.951, 0.957, 0.989, 0.818, 0.95, 0.172, 0.369, 0.957, 0.561, 0.249, 0.911, 0.79, 0.631, 0.651, 0.697, 0.393, 0.995, 0.792, 0.653, 0.376, 0.611, 0.754, 0.639, 0.81, 0.846, 0.635, 0.999, 0.703, 0.954, 0.846, 0.881, 0.998, 0.574, 0.893, 0.938, 0.629, 0.612, 0.26, 0.598, 0.994, 0.993, 0.419, 0.975, 0.935, 0.915, 0.777, 0.638, 0.82, 0.416, 0.478, 0.975, 0.711, 0.904, 0.847, 0.21, 0.974, 0.691, 0.927, 0.778, 0.938, 0.988, 0.956, 0.746, 0.376, 0.56, 0.321, 0.894, 0.974]
global origin = 1
global destination = 60