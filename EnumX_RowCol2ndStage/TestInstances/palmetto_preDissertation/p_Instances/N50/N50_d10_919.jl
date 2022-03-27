global arcs = [1 8; 1 9; 1 27; 1 31; 1 34; 2 5; 2 20; 2 30; 2 36; 2 45; 3 19; 3 25; 3 27; 3 28; 3 35; 3 36; 3 37; 3 41; 3 49; 4 6; 4 30; 4 41; 4 48; 5 4; 5 9; 5 21; 5 36; 5 41; 6 14; 6 15; 6 16; 6 17; 6 22; 6 29; 7 4; 7 31; 7 35; 7 37; 7 38; 7 43; 7 44; 7 48; 7 50; 8 36; 9 5; 9 17; 9 18; 9 21; 9 26; 9 33; 10 6; 10 13; 10 31; 10 38; 10 43; 10 50; 11 6; 11 12; 11 28; 11 34; 11 38; 11 50; 12 16; 12 38; 13 5; 13 11; 13 14; 13 15; 13 22; 14 3; 14 7; 14 19; 14 41; 14 43; 14 48; 15 8; 15 10; 15 11; 15 21; 15 37; 15 39; 15 45; 16 4; 16 5; 16 32; 16 33; 16 36; 16 47; 17 3; 17 4; 17 11; 17 16; 17 25; 17 31; 17 38; 17 44; 18 17; 18 26; 18 29; 18 32; 18 37; 19 39; 19 50; 20 28; 20 38; 20 43; 21 26; 21 29; 21 30; 22 24; 22 32; 22 36; 22 46; 22 49; 23 2; 23 3; 23 6; 23 15; 23 44; 24 5; 24 7; 24 45; 24 46; 24 50; 25 4; 25 5; 25 17; 25 19; 25 23; 25 24; 25 40; 26 10; 26 17; 26 18; 26 22; 26 25; 26 29; 26 48; 27 6; 27 7; 27 14; 27 20; 27 23; 27 42; 27 44; 28 40; 29 4; 29 33; 29 37; 29 40; 30 17; 30 22; 30 25; 30 26; 30 27; 31 23; 31 36; 31 37; 32 18; 32 22; 32 35; 32 38; 33 3; 33 13; 33 34; 33 45; 34 4; 34 5; 34 6; 34 9; 34 11; 34 16; 34 38; 34 40; 34 42; 35 3; 35 26; 35 37; 35 47; 36 7; 36 8; 36 19; 36 25; 36 29; 36 30; 36 33; 36 34; 36 37; 36 38; 36 39; 37 13; 37 16; 37 23; 37 28; 37 32; 37 35; 37 47; 37 48; 37 49; 38 17; 38 21; 39 38; 39 42; 39 47; 40 39; 41 16; 41 26; 42 26; 42 32; 42 37; 42 43; 42 44; 43 3; 43 10; 43 25; 43 34; 43 35; 43 41; 43 46; 44 5; 44 7; 44 10; 44 18; 44 36; 44 37; 44 41; 45 8; 45 11; 45 23; 45 32; 45 40; 45 42; 46 22; 46 24; 46 26; 46 47; 47 6; 47 16; 47 17; 47 23; 47 35; 47 39; 47 43; 47 50; 48 11; 48 16; 48 19; 48 26; 48 34; 48 42; 48 47; 49 2; 49 7; 49 16; 49 17; 49 36; 49 47]
global d_x = [3.0, 7.0, 1.0, 5.0, 4.0, 5.0, 7.0, 1.0, 1.0, 10.0, 1.0, 6.0, 2.0, 5.0, 5.0, 8.0, 10.0, 9.0, 1.0, 2.0, 2.0, 7.0, 9.0, 5.0, 9.0, 3.0, 4.0, 1.0, 1.0, 9.0, 5.0, 3.0, 3.0, 2.0, 7.0, 1.0, 5.0, 2.0, 1.0, 4.0, 10.0, 1.0, 4.0, 8.0, 1.0, 8.0, 2.0, 6.0, 10.0, 6.0, 3.0, 5.0, 3.0, 6.0, 8.0, 9.0, 3.0, 7.0, 1.0, 5.0, 2.0, 8.0, 1.0, 8.0, 3.0, 7.0, 7.0, 10.0, 3.0, 7.0, 5.0, 3.0, 2.0, 3.0, 9.0, 3.0, 4.0, 2.0, 1.0, 6.0, 4.0, 8.0, 4.0, 2.0, 10.0, 9.0, 9.0, 8.0, 7.0, 3.0, 9.0, 2.0, 6.0, 4.0, 10.0, 8.0, 2.0, 3.0, 9.0, 5.0, 5.0, 6.0, 5.0, 7.0, 4.0, 10.0, 3.0, 2.0, 10.0, 10.0, 4.0, 5.0, 1.0, 5.0, 5.0, 2.0, 4.0, 10.0, 7.0, 5.0, 4.0, 1.0, 9.0, 4.0, 4.0, 6.0, 7.0, 4.0, 6.0, 1.0, 6.0, 9.0, 4.0, 7.0, 2.0, 9.0, 2.0, 6.0, 1.0, 6.0, 1.0, 10.0, 10.0, 7.0, 5.0, 10.0, 10.0, 9.0, 6.0, 10.0, 3.0, 10.0, 9.0, 4.0, 7.0, 2.0, 2.0, 9.0, 9.0, 5.0, 6.0, 8.0, 9.0, 5.0, 2.0, 10.0, 4.0, 4.0, 10.0, 9.0, 9.0, 5.0, 3.0, 10.0, 1.0, 1.0, 2.0, 7.0, 6.0, 1.0, 4.0, 2.0, 2.0, 8.0, 1.0, 1.0, 1.0, 3.0, 9.0, 3.0, 5.0, 4.0, 2.0, 9.0, 7.0, 1.0, 1.0, 8.0, 1.0, 2.0, 3.0, 9.0, 5.0, 1.0, 7.0, 7.0, 3.0, 7.0, 1.0, 6.0, 8.0, 7.0, 8.0, 8.0, 8.0, 6.0, 6.0, 2.0, 7.0, 3.0, 8.0, 4.0, 3.0, 5.0, 6.0, 6.0, 2.0, 1.0, 6.0, 1.0, 2.0, 1.0, 6.0, 8.0, 1.0, 1.0, 1.0, 2.0, 8.0, 6.0, 3.0, 3.0, 6.0, 2.0, 10.0, 1.0, 9.0, 1.0, 8.0, 9.0, 7.0, 6.0, 6.0, 4.0, 10.0, 1.0, 5.0]
global b_x = 5
global d_y = [8.0, 4.0, 10.0, 5.0, 1.0, 3.0, 2.0, 5.0, 4.0, 2.0, 6.0, 5.0, 8.0, 3.0, 1.0, 7.0, 4.0, 4.0, 7.0, 4.0, 4.0, 6.0, 9.0, 9.0, 5.0, 8.0, 9.0, 9.0, 7.0, 4.0, 2.0, 3.0, 4.0, 1.0, 3.0, 10.0, 3.0, 6.0, 3.0, 2.0, 6.0, 10.0, 3.0, 8.0, 6.0, 4.0, 2.0, 10.0, 7.0, 7.0, 8.0, 4.0, 5.0, 9.0, 8.0, 8.0, 8.0, 9.0, 3.0, 7.0, 4.0, 2.0, 6.0, 2.0, 10.0, 2.0, 4.0, 6.0, 6.0, 6.0, 8.0, 6.0, 6.0, 9.0, 7.0, 4.0, 7.0, 9.0, 5.0, 5.0, 8.0, 2.0, 10.0, 2.0, 1.0, 9.0, 1.0, 2.0, 7.0, 7.0, 1.0, 8.0, 7.0, 4.0, 1.0, 9.0, 8.0, 5.0, 10.0, 2.0, 1.0, 3.0, 6.0, 2.0, 4.0, 2.0, 6.0, 2.0, 3.0, 3.0, 6.0, 7.0, 3.0, 10.0, 3.0, 6.0, 2.0, 4.0, 3.0, 7.0, 4.0, 5.0, 3.0, 9.0, 9.0, 3.0, 8.0, 8.0, 1.0, 1.0, 1.0, 8.0, 9.0, 1.0, 6.0, 4.0, 7.0, 4.0, 8.0, 7.0, 3.0, 6.0, 2.0, 6.0, 10.0, 2.0, 4.0, 8.0, 7.0, 10.0, 8.0, 2.0, 10.0, 6.0, 3.0, 10.0, 1.0, 10.0, 9.0, 1.0, 1.0, 3.0, 10.0, 9.0, 2.0, 8.0, 6.0, 2.0, 5.0, 1.0, 7.0, 3.0, 1.0, 4.0, 9.0, 7.0, 1.0, 7.0, 6.0, 5.0, 2.0, 3.0, 5.0, 5.0, 1.0, 5.0, 7.0, 4.0, 10.0, 5.0, 1.0, 3.0, 8.0, 4.0, 9.0, 10.0, 6.0, 1.0, 1.0, 10.0, 4.0, 4.0, 1.0, 10.0, 8.0, 10.0, 1.0, 1.0, 7.0, 3.0, 6.0, 2.0, 1.0, 5.0, 9.0, 6.0, 4.0, 6.0, 8.0, 3.0, 2.0, 1.0, 7.0, 8.0, 2.0, 5.0, 8.0, 4.0, 4.0, 10.0, 7.0, 8.0, 3.0, 1.0, 6.0, 10.0, 10.0, 9.0, 4.0, 3.0, 2.0, 10.0, 3.0, 9.0, 5.0, 2.0, 2.0, 3.0, 6.0, 8.0, 1.0, 7.0, 5.0, 1.0, 1.0, 2.0, 9.0]
global b_y = 10
global p = [0.12, 0.526, 0.995, 0.884, 0.823, 0.614, 0.536, 0.205, 0.481, 0.97, 0.773, 0.633, 0.442, 0.017, 0.255, 0.335, 0.867, 0.897, 0.234, 0.714, 0.388, 0.932, 0.177, 0.177, 0.617, 0.713, 0.876, 0.93, 0.496, 0.359, 0.408, 0.977, 0.442, 0.493, 0.422, 0.776, 0.206, 0.294, 0.774, 0.424, 0.629, 0.7, 0.937, 0.531, 0.229, 0.659, 0.073, 0.996, 0.539, 0.293, 0.988, 0.317, 0.643, 0.31, 0.48, 0.574, 0.023, 0.801, 0.819, 0.493, 0.494, 0.237, 0.445, 0.983, 0.68, 0.411, 0.614, 0.579, 0.195, 0.385, 0.76, 0.059, 0.745, 0.752, 0.418, 0.094, 0.31, 0.229, 0.735, 0.527, 0.797, 0.773, 0.293, 0.938, 0.161, 0.059, 0.989, 0.49, 0.557, 0.503, 0.067, 0.663, 0.666, 0.084, 0.189, 0.866, 0.997, 0.365, 0.032, 0.819, 0.287, 0.682, 0.953, 0.172, 0.731, 0.681, 0.121, 0.456, 0.664, 0.461, 0.991, 0.445, 0.751, 0.56, 0.474, 0.546, 0.769, 0.244, 0.302, 0.577, 0.177, 0.406, 0.249, 0.613, 0.99, 0.166, 0.486, 0.622, 0.173, 0.096, 0.622, 0.08, 0.39, 0.528, 0.64, 0.358, 0.967, 0.173, 0.305, 0.337, 0.363, 0.097, 0.366, 0.818, 0.885, 0.978, 0.447, 0.45, 0.586, 0.695, 0.975, 0.172, 0.968, 0.372, 0.966, 0.67, 0.274, 0.971, 0.496, 0.516, 0.406, 0.461, 0.696, 0.964, 0.852, 0.267, 0.891, 0.17, 0.321, 0.029, 0.595, 0.268, 0.244, 0.056, 0.75, 0.932, 0.798, 0.51, 0.669, 0.891, 0.002, 0.064, 0.891, 0.118, 0.051, 0.411, 0.791, 0.806, 0.014, 0.094, 0.758, 0.353, 0.381, 0.065, 0.82, 0.281, 0.844, 0.576, 0.096, 0.174, 0.301, 0.288, 0.594, 0.158, 0.09, 0.06, 0.544, 0.886, 0.526, 0.109, 0.303, 0.872, 0.859, 0.257, 0.749, 0.971, 0.794, 0.481, 0.281, 0.141, 0.043, 0.164, 0.544, 0.311, 0.609, 0.862, 0.296, 0.059, 0.492, 0.44, 0.017, 0.172, 0.884, 0.519, 0.932, 0.417, 0.985, 0.979, 0.906, 0.959, 0.403, 0.125, 0.613, 0.821, 0.302, 0.919, 0.307, 0.6, 0.487, 0.157, 0.512, 0.336, 0.322, 0.903, 0.721, 0.978, 0.372]
global q = [0.751, 0.565, 0.997, 0.899, 0.933, 0.938, 0.59, 0.734, 0.818, 0.982, 0.8, 0.94, 0.498, 0.276, 0.457, 0.545, 0.881, 0.963, 0.427, 0.782, 0.513, 0.988, 0.782, 0.327, 0.963, 0.742, 0.959, 0.965, 0.644, 0.479, 0.72, 0.991, 0.505, 0.518, 0.999, 0.86, 0.735, 0.421, 0.995, 0.718, 0.705, 0.803, 0.98, 0.86, 0.989, 0.894, 0.72, 0.998, 0.837, 0.82, 0.998, 0.542, 0.981, 0.627, 0.769, 0.998, 0.659, 0.986, 0.825, 0.971, 0.999, 0.355, 0.914, 0.992, 0.917, 0.471, 0.998, 0.743, 0.576, 0.559, 0.887, 0.49, 0.99, 0.853, 0.627, 0.305, 0.313, 0.868, 0.792, 0.923, 0.936, 0.831, 0.869, 0.949, 0.467, 0.907, 0.991, 0.996, 0.952, 0.677, 0.249, 0.929, 0.947, 0.549, 0.996, 0.913, 0.999, 0.782, 0.542, 0.951, 0.923, 0.918, 0.977, 0.805, 0.739, 0.836, 0.306, 0.947, 0.838, 0.77, 0.996, 0.856, 0.881, 0.917, 0.795, 0.934, 0.871, 0.683, 0.906, 0.937, 0.281, 0.439, 0.604, 0.695, 0.995, 0.539, 0.79, 0.756, 0.41, 0.294, 0.624, 0.08, 0.872, 0.981, 0.859, 0.417, 0.975, 0.384, 0.315, 0.808, 0.558, 0.578, 0.915, 0.955, 0.989, 0.991, 0.51, 0.702, 0.653, 0.995, 0.989, 0.901, 0.978, 0.764, 0.986, 0.936, 0.772, 0.974, 0.702, 0.555, 0.954, 0.581, 0.826, 0.99, 0.969, 0.753, 0.912, 0.197, 0.967, 0.293, 0.74, 0.939, 0.388, 0.302, 0.754, 0.953, 0.977, 0.822, 0.669, 0.933, 0.59, 0.529, 0.932, 0.292, 0.674, 0.777, 0.943, 0.864, 0.048, 0.944, 0.817, 0.477, 0.889, 0.996, 0.96, 0.402, 0.854, 0.726, 0.144, 0.811, 0.449, 0.442, 0.758, 0.297, 0.565, 0.819, 0.951, 0.927, 0.72, 0.455, 0.609, 0.936, 0.935, 0.363, 0.813, 0.98, 0.816, 0.551, 0.364, 0.948, 0.744, 0.177, 0.663, 0.548, 0.642, 0.982, 0.59, 0.91, 0.996, 0.752, 0.852, 0.825, 0.904, 0.949, 0.933, 0.696, 0.993, 0.996, 0.931, 0.993, 0.823, 0.633, 0.912, 0.84, 0.322, 0.987, 0.765, 0.693, 0.736, 0.308, 0.774, 0.833, 0.956, 0.948, 0.767, 0.998, 0.502]
global origin = 1
global destination = 50