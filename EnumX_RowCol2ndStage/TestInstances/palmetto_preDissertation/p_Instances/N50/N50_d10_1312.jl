global arcs = [1 14; 1 24; 1 35; 1 43; 2 3; 2 14; 2 19; 2 28; 2 33; 2 36; 2 38; 3 6; 3 12; 3 19; 3 38; 4 9; 4 12; 4 17; 4 31; 5 13; 5 17; 5 19; 5 25; 5 27; 5 29; 5 37; 5 46; 6 5; 6 14; 6 22; 6 29; 6 30; 6 46; 7 3; 7 4; 7 18; 7 27; 8 4; 8 10; 8 29; 8 45; 9 5; 9 6; 9 16; 9 25; 9 29; 9 46; 10 2; 10 3; 10 5; 10 24; 11 13; 11 14; 11 20; 11 30; 11 33; 11 47; 11 48; 12 2; 12 13; 12 28; 12 29; 12 32; 13 9; 13 15; 13 23; 13 41; 13 42; 13 43; 14 15; 14 16; 14 21; 14 22; 14 23; 14 24; 14 25; 14 28; 14 32; 14 36; 14 38; 14 50; 15 3; 15 7; 15 16; 15 21; 15 30; 15 35; 16 26; 16 40; 16 43; 16 47; 17 14; 17 23; 17 24; 17 48; 18 23; 18 24; 18 37; 18 41; 18 49; 19 20; 20 8; 20 10; 20 26; 20 27; 20 41; 20 47; 21 4; 21 13; 21 16; 21 22; 21 29; 21 36; 22 6; 22 11; 22 37; 22 40; 22 42; 23 15; 23 36; 23 44; 24 7; 24 10; 24 17; 24 26; 24 34; 24 45; 24 49; 25 5; 25 15; 25 16; 26 19; 26 20; 26 23; 26 27; 26 28; 26 31; 27 4; 27 24; 27 50; 28 8; 28 43; 29 12; 29 18; 29 41; 30 23; 30 27; 31 17; 31 19; 31 27; 31 39; 32 15; 32 26; 32 30; 32 34; 32 39; 32 46; 33 7; 33 12; 33 18; 33 35; 33 42; 34 14; 34 17; 34 22; 34 42; 35 2; 35 3; 35 12; 35 25; 35 27; 35 44; 35 47; 36 10; 36 50; 37 10; 37 16; 37 19; 37 33; 37 45; 37 50; 38 9; 38 29; 38 30; 38 37; 38 45; 38 47; 39 2; 39 16; 39 23; 39 40; 40 13; 40 20; 40 24; 40 36; 40 39; 40 47; 41 8; 41 20; 41 22; 41 31; 41 44; 42 13; 42 16; 42 18; 42 20; 42 21; 43 5; 43 11; 43 12; 43 25; 43 29; 43 32; 43 34; 43 46; 44 6; 44 7; 44 14; 44 23; 44 24; 44 49; 45 15; 45 28; 45 38; 45 41; 46 13; 46 23; 46 25; 46 29; 47 8; 47 9; 47 11; 47 35; 47 39; 48 11; 48 13; 48 47; 48 49; 49 4; 49 11; 49 35]
global d_x = [9.0, 4.0, 2.0, 4.0, 5.0, 1.0, 9.0, 1.0, 4.0, 6.0, 10.0, 2.0, 4.0, 5.0, 1.0, 5.0, 7.0, 6.0, 4.0, 3.0, 1.0, 1.0, 9.0, 3.0, 6.0, 6.0, 4.0, 9.0, 10.0, 7.0, 4.0, 10.0, 1.0, 2.0, 4.0, 9.0, 9.0, 3.0, 3.0, 6.0, 9.0, 7.0, 5.0, 2.0, 6.0, 1.0, 4.0, 7.0, 1.0, 6.0, 1.0, 5.0, 1.0, 4.0, 1.0, 6.0, 1.0, 2.0, 2.0, 10.0, 5.0, 10.0, 7.0, 2.0, 6.0, 9.0, 7.0, 5.0, 1.0, 3.0, 10.0, 9.0, 10.0, 9.0, 4.0, 6.0, 5.0, 4.0, 7.0, 5.0, 9.0, 3.0, 6.0, 10.0, 5.0, 3.0, 4.0, 1.0, 7.0, 2.0, 10.0, 3.0, 4.0, 9.0, 10.0, 10.0, 6.0, 9.0, 6.0, 6.0, 3.0, 7.0, 4.0, 8.0, 9.0, 9.0, 7.0, 7.0, 2.0, 3.0, 1.0, 5.0, 4.0, 4.0, 3.0, 1.0, 3.0, 7.0, 8.0, 5.0, 6.0, 4.0, 1.0, 4.0, 5.0, 4.0, 9.0, 1.0, 3.0, 1.0, 3.0, 7.0, 9.0, 4.0, 1.0, 1.0, 10.0, 2.0, 6.0, 1.0, 1.0, 6.0, 10.0, 5.0, 2.0, 9.0, 8.0, 3.0, 10.0, 8.0, 6.0, 1.0, 3.0, 2.0, 7.0, 10.0, 8.0, 4.0, 7.0, 6.0, 10.0, 5.0, 8.0, 9.0, 7.0, 7.0, 10.0, 3.0, 8.0, 2.0, 10.0, 6.0, 4.0, 5.0, 3.0, 8.0, 1.0, 9.0, 9.0, 7.0, 1.0, 3.0, 7.0, 5.0, 6.0, 5.0, 8.0, 6.0, 4.0, 3.0, 4.0, 9.0, 3.0, 2.0, 7.0, 10.0, 5.0, 3.0, 5.0, 5.0, 5.0, 10.0, 2.0, 7.0, 9.0, 6.0, 8.0, 2.0, 10.0, 10.0, 6.0, 2.0, 7.0, 3.0, 10.0, 1.0, 5.0, 5.0, 2.0, 1.0, 9.0, 8.0, 6.0, 8.0, 4.0, 6.0, 9.0, 8.0, 6.0, 10.0, 8.0, 4.0, 7.0, 10.0, 5.0, 3.0, 10.0, 4.0, 3.0, 1.0, 4.0]
global b_x = 5
global d_y = [8.0, 5.0, 9.0, 7.0, 10.0, 5.0, 2.0, 5.0, 7.0, 7.0, 7.0, 3.0, 2.0, 4.0, 8.0, 5.0, 9.0, 6.0, 9.0, 3.0, 3.0, 2.0, 7.0, 10.0, 3.0, 7.0, 9.0, 7.0, 5.0, 1.0, 3.0, 2.0, 9.0, 8.0, 2.0, 6.0, 7.0, 8.0, 1.0, 9.0, 2.0, 9.0, 2.0, 1.0, 9.0, 4.0, 2.0, 6.0, 8.0, 9.0, 7.0, 8.0, 2.0, 8.0, 10.0, 7.0, 2.0, 9.0, 4.0, 5.0, 1.0, 6.0, 6.0, 8.0, 3.0, 7.0, 8.0, 2.0, 4.0, 10.0, 4.0, 6.0, 9.0, 2.0, 1.0, 10.0, 8.0, 4.0, 8.0, 2.0, 8.0, 8.0, 9.0, 5.0, 5.0, 5.0, 2.0, 10.0, 4.0, 7.0, 5.0, 9.0, 1.0, 10.0, 10.0, 1.0, 2.0, 5.0, 10.0, 9.0, 1.0, 4.0, 2.0, 2.0, 10.0, 1.0, 10.0, 2.0, 3.0, 7.0, 5.0, 8.0, 5.0, 8.0, 4.0, 2.0, 6.0, 7.0, 3.0, 1.0, 9.0, 1.0, 10.0, 10.0, 1.0, 3.0, 10.0, 4.0, 2.0, 10.0, 3.0, 6.0, 6.0, 3.0, 5.0, 3.0, 6.0, 6.0, 5.0, 5.0, 4.0, 1.0, 8.0, 9.0, 4.0, 4.0, 8.0, 7.0, 5.0, 1.0, 8.0, 2.0, 3.0, 9.0, 2.0, 9.0, 3.0, 6.0, 5.0, 7.0, 6.0, 2.0, 10.0, 7.0, 9.0, 10.0, 2.0, 9.0, 5.0, 10.0, 4.0, 8.0, 5.0, 4.0, 8.0, 9.0, 1.0, 4.0, 5.0, 4.0, 10.0, 3.0, 8.0, 8.0, 2.0, 2.0, 5.0, 9.0, 1.0, 5.0, 6.0, 5.0, 9.0, 2.0, 6.0, 8.0, 6.0, 8.0, 5.0, 8.0, 10.0, 8.0, 3.0, 7.0, 9.0, 1.0, 10.0, 6.0, 8.0, 8.0, 7.0, 3.0, 3.0, 8.0, 6.0, 4.0, 8.0, 6.0, 2.0, 8.0, 4.0, 10.0, 1.0, 6.0, 6.0, 7.0, 2.0, 5.0, 10.0, 4.0, 6.0, 9.0, 8.0, 5.0, 4.0, 3.0, 8.0, 9.0, 8.0, 6.0, 10.0]
global b_y = 10
global p = [0.419, 0.856, 0.656, 0.409, 0.059, 0.2, 0.35, 0.051, 0.281, 0.781, 0.166, 0.038, 0.267, 0.229, 0.991, 0.698, 0.108, 0.838, 0.897, 0.519, 0.42, 0.25, 0.898, 0.032, 0.869, 0.285, 0.179, 0.083, 0.914, 0.658, 0.64, 0.079, 0.499, 0.848, 0.125, 0.278, 0.888, 0.044, 0.673, 0.483, 0.466, 0.159, 0.945, 0.865, 0.23, 0.831, 0.781, 0.018, 0.812, 0.154, 0.683, 0.263, 0.288, 0.454, 0.648, 0.717, 0.634, 0.13, 0.386, 0.995, 0.294, 0.222, 0.931, 0.712, 0.765, 0.667, 0.136, 0.673, 0.966, 0.404, 0.791, 0.205, 0.492, 0.849, 0.911, 0.236, 0.431, 0.141, 0.72, 0.582, 0.294, 0.589, 0.325, 0.232, 0.82, 0.58, 0.82, 0.897, 0.479, 0.092, 0.538, 0.539, 0.174, 0.309, 0.624, 0.202, 0.246, 0.695, 0.826, 0.161, 0.171, 0.911, 0.598, 0.964, 0.281, 0.236, 0.804, 0.835, 0.451, 0.259, 0.363, 0.449, 0.191, 0.629, 0.303, 0.705, 0.097, 0.895, 0.35, 0.745, 0.608, 0.388, 0.104, 0.413, 0.571, 0.02, 0.144, 0.228, 0.275, 0.296, 0.629, 0.005, 0.974, 0.993, 0.042, 0.108, 0.196, 0.012, 0.398, 0.359, 0.808, 0.697, 0.246, 0.197, 0.208, 0.139, 0.398, 0.964, 0.13, 0.772, 0.167, 0.654, 0.601, 0.969, 0.149, 0.321, 0.787, 0.457, 0.716, 0.534, 0.717, 0.216, 0.016, 0.059, 0.599, 0.227, 0.581, 0.514, 0.488, 0.348, 0.4, 0.27, 0.193, 0.926, 0.905, 0.491, 0.945, 0.265, 0.982, 0.639, 0.098, 0.189, 0.078, 0.169, 0.23, 0.256, 0.413, 0.004, 0.792, 0.66, 0.581, 0.777, 0.522, 0.587, 0.453, 0.147, 0.218, 0.473, 0.145, 0.719, 0.469, 0.7, 0.353, 0.744, 0.493, 0.983, 0.495, 0.781, 0.9, 0.942, 0.051, 0.051, 0.245, 0.577, 0.837, 0.004, 0.034, 0.068, 0.527, 0.805, 0.532, 0.857, 0.758, 0.439, 0.794, 0.38, 0.081, 0.563, 0.674, 0.385, 0.319, 0.561, 0.051, 0.091, 0.637, 0.15, 0.384, 0.003, 0.648, 0.371, 0.081]
global q = [0.62, 0.959, 0.973, 0.915, 0.314, 0.623, 0.789, 0.859, 0.785, 0.808, 0.642, 0.188, 0.417, 0.625, 0.998, 0.974, 0.875, 0.959, 0.914, 0.529, 0.523, 0.68, 0.974, 0.187, 0.912, 0.656, 0.349, 0.279, 0.971, 0.831, 0.866, 0.608, 0.684, 0.877, 0.714, 0.521, 0.911, 0.658, 0.778, 0.605, 0.633, 0.638, 0.973, 0.968, 0.298, 0.868, 0.959, 0.289, 0.97, 0.999, 0.897, 0.883, 0.616, 0.647, 0.775, 0.741, 0.649, 0.772, 0.977, 0.998, 0.872, 0.844, 0.995, 0.805, 0.858, 0.853, 0.49, 0.733, 0.997, 0.889, 0.916, 0.91, 0.709, 0.973, 0.976, 0.364, 0.524, 0.336, 0.747, 0.589, 0.944, 0.589, 0.705, 0.571, 0.971, 0.589, 0.834, 0.982, 0.599, 0.909, 0.742, 0.828, 0.808, 0.704, 0.681, 0.733, 0.885, 0.782, 0.987, 0.915, 0.724, 0.951, 0.817, 0.991, 0.498, 0.578, 0.939, 0.971, 0.614, 0.354, 0.938, 0.721, 0.208, 0.643, 0.568, 0.994, 0.25, 0.99, 0.646, 0.997, 0.748, 0.848, 0.945, 0.547, 0.743, 0.542, 0.165, 0.258, 0.526, 0.924, 0.773, 0.675, 0.997, 0.997, 0.332, 0.871, 0.648, 0.21, 0.893, 0.667, 0.854, 0.902, 0.514, 0.231, 0.763, 0.145, 0.619, 0.995, 0.172, 0.879, 0.708, 0.813, 0.668, 0.976, 0.422, 0.682, 0.863, 0.558, 0.729, 0.892, 0.961, 0.607, 0.838, 0.39, 0.671, 0.383, 0.734, 0.656, 0.759, 0.863, 0.958, 0.933, 0.997, 0.966, 0.977, 0.621, 0.956, 0.499, 0.996, 0.859, 0.403, 0.855, 0.466, 0.208, 0.74, 0.924, 0.548, 0.718, 0.813, 0.677, 0.723, 0.926, 0.807, 0.652, 0.87, 0.394, 0.736, 0.55, 0.913, 0.924, 0.609, 0.703, 0.439, 0.921, 0.641, 0.989, 0.811, 0.822, 0.93, 0.984, 0.838, 0.699, 0.792, 0.718, 0.952, 0.995, 0.901, 0.886, 0.634, 0.894, 0.882, 0.898, 0.92, 0.582, 0.813, 0.728, 0.113, 0.983, 0.694, 0.562, 0.449, 0.817, 0.443, 0.882, 0.871, 0.419, 0.515, 0.903, 0.946, 0.371, 0.597]
global origin = 1
global destination = 50