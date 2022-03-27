global arcs = [1 5; 1 17; 1 26; 2 5; 2 8; 2 12; 2 16; 2 28; 2 37; 3 20; 3 22; 3 25; 4 2; 4 31; 4 36; 4 42; 4 46; 5 2; 5 11; 5 36; 5 49; 6 3; 6 17; 6 26; 6 43; 6 45; 6 47; 7 6; 7 33; 8 33; 8 45; 8 48; 9 12; 9 33; 9 43; 9 46; 10 6; 10 37; 10 38; 10 41; 10 42; 11 2; 11 19; 11 21; 11 29; 11 33; 11 50; 12 8; 12 24; 12 25; 12 41; 12 44; 12 46; 12 48; 12 49; 13 7; 13 23; 13 36; 13 43; 14 3; 14 21; 14 25; 14 29; 14 47; 15 16; 15 18; 15 31; 15 34; 15 37; 15 44; 16 3; 16 7; 16 8; 16 9; 16 38; 16 39; 16 43; 16 50; 17 2; 17 3; 17 5; 17 10; 17 12; 17 22; 18 5; 19 6; 19 15; 19 38; 20 34; 20 38; 20 49; 21 8; 21 11; 21 13; 21 32; 21 36; 21 45; 22 5; 22 11; 22 13; 22 16; 22 28; 22 30; 22 40; 22 42; 22 44; 22 45; 23 9; 23 13; 23 17; 23 19; 23 22; 23 39; 23 40; 23 41; 23 43; 23 45; 23 48; 24 9; 24 23; 24 46; 24 48; 25 47; 25 48; 26 3; 26 5; 26 15; 26 19; 26 27; 27 4; 27 15; 27 38; 27 50; 28 14; 28 25; 28 35; 28 42; 28 45; 29 2; 29 22; 29 25; 29 30; 29 31; 30 16; 31 8; 31 14; 31 25; 31 33; 31 38; 31 39; 32 8; 32 24; 32 31; 32 37; 33 7; 33 17; 33 26; 34 49; 35 8; 35 13; 35 14; 35 17; 35 21; 35 37; 35 38; 36 4; 36 12; 36 19; 36 21; 36 33; 36 34; 36 39; 36 48; 37 7; 37 23; 37 31; 38 6; 38 13; 38 15; 38 42; 38 49; 39 34; 39 35; 40 8; 40 9; 40 11; 40 17; 40 21; 40 24; 40 31; 40 34; 40 36; 40 42; 40 43; 40 46; 40 47; 40 50; 41 3; 41 4; 41 6; 41 28; 41 48; 42 5; 42 30; 42 31; 42 34; 42 41; 42 43; 43 24; 43 31; 43 35; 44 4; 44 9; 44 19; 44 29; 44 47; 45 2; 46 7; 46 11; 46 15; 46 17; 47 10; 47 14; 47 23; 47 29; 47 48; 48 14; 48 25; 48 34; 48 35; 48 39; 48 49; 49 6; 49 8; 49 27; 49 39; 49 47]
global d_x = [2.0, 7.0, 1.0, 6.0, 8.0, 2.0, 3.0, 10.0, 3.0, 2.0, 7.0, 8.0, 1.0, 7.0, 6.0, 3.0, 2.0, 9.0, 1.0, 4.0, 4.0, 10.0, 6.0, 8.0, 6.0, 9.0, 6.0, 5.0, 2.0, 9.0, 8.0, 7.0, 1.0, 2.0, 5.0, 6.0, 7.0, 6.0, 4.0, 3.0, 6.0, 5.0, 3.0, 6.0, 10.0, 2.0, 1.0, 7.0, 2.0, 1.0, 3.0, 10.0, 10.0, 7.0, 10.0, 8.0, 1.0, 1.0, 2.0, 3.0, 10.0, 7.0, 10.0, 7.0, 1.0, 2.0, 4.0, 10.0, 10.0, 10.0, 9.0, 5.0, 4.0, 7.0, 7.0, 5.0, 10.0, 10.0, 8.0, 6.0, 2.0, 1.0, 3.0, 8.0, 10.0, 9.0, 2.0, 6.0, 6.0, 3.0, 9.0, 4.0, 6.0, 4.0, 10.0, 8.0, 2.0, 4.0, 6.0, 1.0, 8.0, 4.0, 10.0, 9.0, 2.0, 3.0, 10.0, 4.0, 3.0, 2.0, 8.0, 1.0, 2.0, 10.0, 1.0, 1.0, 10.0, 4.0, 6.0, 10.0, 4.0, 3.0, 3.0, 8.0, 3.0, 5.0, 5.0, 6.0, 1.0, 3.0, 7.0, 5.0, 6.0, 1.0, 4.0, 8.0, 2.0, 8.0, 2.0, 7.0, 1.0, 3.0, 1.0, 5.0, 2.0, 7.0, 6.0, 1.0, 7.0, 3.0, 8.0, 1.0, 7.0, 9.0, 5.0, 6.0, 3.0, 5.0, 2.0, 7.0, 1.0, 8.0, 1.0, 2.0, 2.0, 7.0, 7.0, 3.0, 8.0, 10.0, 1.0, 8.0, 5.0, 8.0, 9.0, 2.0, 3.0, 2.0, 1.0, 9.0, 7.0, 7.0, 1.0, 2.0, 2.0, 1.0, 6.0, 5.0, 9.0, 3.0, 9.0, 6.0, 3.0, 8.0, 3.0, 5.0, 7.0, 4.0, 9.0, 5.0, 6.0, 10.0, 1.0, 6.0, 7.0, 10.0, 9.0, 5.0, 8.0, 4.0, 9.0, 4.0, 6.0, 6.0, 7.0, 8.0, 9.0, 6.0, 1.0, 9.0, 1.0, 6.0, 1.0, 5.0, 2.0, 10.0, 1.0, 5.0, 3.0, 3.0, 10.0, 5.0, 1.0, 8.0, 8.0, 4.0, 1.0]
global b_x = 5
global d_y = [7.0, 7.0, 4.0, 2.0, 9.0, 6.0, 9.0, 1.0, 9.0, 6.0, 2.0, 3.0, 7.0, 1.0, 2.0, 5.0, 10.0, 10.0, 6.0, 3.0, 6.0, 6.0, 8.0, 8.0, 9.0, 1.0, 9.0, 6.0, 1.0, 10.0, 4.0, 5.0, 1.0, 3.0, 4.0, 6.0, 2.0, 4.0, 4.0, 10.0, 10.0, 3.0, 6.0, 5.0, 1.0, 8.0, 10.0, 2.0, 10.0, 7.0, 4.0, 6.0, 7.0, 8.0, 10.0, 5.0, 7.0, 5.0, 8.0, 9.0, 1.0, 8.0, 2.0, 4.0, 1.0, 3.0, 8.0, 8.0, 8.0, 5.0, 6.0, 3.0, 7.0, 8.0, 6.0, 4.0, 3.0, 5.0, 7.0, 5.0, 10.0, 1.0, 2.0, 10.0, 9.0, 5.0, 5.0, 6.0, 5.0, 10.0, 1.0, 6.0, 8.0, 9.0, 3.0, 10.0, 6.0, 10.0, 10.0, 1.0, 4.0, 1.0, 4.0, 1.0, 4.0, 6.0, 5.0, 3.0, 9.0, 7.0, 5.0, 1.0, 7.0, 2.0, 9.0, 9.0, 2.0, 3.0, 1.0, 8.0, 1.0, 4.0, 2.0, 4.0, 4.0, 4.0, 5.0, 9.0, 6.0, 6.0, 6.0, 5.0, 4.0, 5.0, 2.0, 6.0, 3.0, 9.0, 3.0, 5.0, 2.0, 5.0, 7.0, 6.0, 8.0, 6.0, 9.0, 1.0, 5.0, 8.0, 5.0, 3.0, 8.0, 7.0, 2.0, 9.0, 8.0, 10.0, 10.0, 7.0, 9.0, 3.0, 5.0, 1.0, 3.0, 7.0, 2.0, 2.0, 5.0, 2.0, 6.0, 5.0, 9.0, 6.0, 3.0, 3.0, 6.0, 4.0, 5.0, 2.0, 3.0, 6.0, 3.0, 1.0, 10.0, 5.0, 2.0, 8.0, 3.0, 2.0, 1.0, 4.0, 4.0, 2.0, 2.0, 9.0, 9.0, 3.0, 6.0, 6.0, 2.0, 1.0, 6.0, 3.0, 2.0, 3.0, 3.0, 1.0, 3.0, 10.0, 6.0, 2.0, 10.0, 2.0, 6.0, 9.0, 8.0, 10.0, 9.0, 8.0, 1.0, 6.0, 3.0, 10.0, 4.0, 9.0, 10.0, 6.0, 6.0, 5.0, 1.0, 3.0, 3.0, 10.0, 4.0, 7.0, 8.0]
global b_y = 10
global p = [0.353, 0.33, 0.945, 0.004, 0.548, 0.459, 0.921, 0.633, 0.411, 0.806, 0.551, 0.171, 0.077, 0.313, 0.342, 0.886, 0.238, 0.453, 0.862, 0.702, 0.75, 0.982, 0.925, 0.558, 0.608, 0.568, 0.495, 0.533, 0.233, 0.558, 0.631, 0.235, 0.311, 0.846, 0.629, 0.664, 0.457, 0.134, 0.262, 0.592, 0.086, 0.304, 0.229, 0.116, 0.272, 0.66, 0.74, 0.267, 0.196, 0.784, 0.448, 0.424, 0.705, 0.258, 0.999, 0.325, 0.442, 0.078, 0.134, 0.51, 0.105, 0.79, 0.455, 0.37, 0.226, 0.578, 0.91, 0.491, 0.753, 0.868, 0.987, 0.486, 0.811, 0.695, 0.679, 0.703, 0.275, 0.204, 0.519, 0.396, 0.317, 0.665, 0.447, 0.625, 0.302, 0.276, 0.729, 0.82, 0.237, 0.407, 0.572, 0.068, 0.678, 0.21, 0.094, 0.539, 0.902, 0.494, 0.638, 0.429, 0.084, 0.138, 0.615, 0.816, 0.31, 0.631, 0.023, 0.638, 0.497, 0.013, 0.313, 0.314, 0.604, 0.087, 0.016, 0.12, 0.905, 0.345, 0.742, 0.784, 0.614, 0.458, 0.986, 0.729, 0.889, 0.288, 0.5, 0.848, 0.583, 0.099, 0.038, 0.876, 0.433, 0.799, 0.805, 0.081, 0.479, 0.043, 0.917, 0.86, 0.689, 0.349, 0.799, 0.436, 0.45, 0.176, 0.914, 0.495, 0.488, 0.409, 0.871, 0.816, 0.984, 0.939, 0.917, 0.527, 0.173, 0.89, 0.029, 0.9, 0.715, 0.842, 0.102, 0.093, 0.083, 0.615, 0.502, 0.809, 0.101, 0.299, 0.871, 0.275, 0.855, 0.445, 0.27, 0.067, 0.373, 0.966, 0.631, 0.806, 0.175, 0.603, 0.547, 0.251, 0.428, 0.403, 0.683, 0.76, 0.732, 0.005, 0.974, 0.858, 0.266, 0.823, 0.232, 0.621, 0.047, 0.319, 0.716, 0.645, 0.801, 0.605, 0.389, 0.498, 0.334, 0.151, 0.272, 0.86, 0.095, 0.323, 0.134, 0.858, 0.848, 0.734, 0.315, 0.682, 0.712, 0.725, 0.738, 0.724, 0.981, 0.955, 0.211, 0.378, 0.277, 0.302, 0.263, 0.175, 0.085, 0.149, 0.091, 0.236, 0.428, 0.483, 0.434, 0.466, 0.082]
global q = [0.935, 0.459, 0.996, 0.268, 0.758, 0.5, 0.97, 0.662, 0.454, 0.947, 0.644, 0.68, 0.766, 0.706, 0.586, 0.997, 0.94, 0.742, 0.943, 0.929, 0.966, 0.991, 0.936, 0.924, 0.888, 0.642, 0.916, 0.625, 0.408, 0.752, 0.996, 0.284, 0.922, 0.941, 0.945, 0.887, 0.735, 0.325, 0.273, 0.643, 0.089, 0.525, 0.945, 0.953, 0.364, 0.833, 0.761, 0.409, 0.386, 0.987, 0.832, 0.578, 0.977, 0.798, 0.999, 0.49, 0.75, 0.679, 0.54, 0.539, 0.752, 0.937, 0.764, 0.501, 0.636, 0.716, 0.941, 0.818, 0.778, 0.941, 0.996, 0.537, 0.865, 0.927, 0.795, 0.965, 0.583, 0.268, 0.856, 0.892, 0.548, 0.916, 0.902, 0.677, 0.909, 0.581, 0.967, 0.846, 0.237, 0.932, 0.834, 0.279, 0.873, 0.882, 0.253, 0.829, 0.936, 0.874, 0.718, 0.489, 0.211, 0.328, 0.637, 0.939, 0.951, 0.933, 0.976, 0.69, 0.712, 0.87, 0.502, 0.705, 0.761, 0.419, 0.653, 0.657, 0.975, 0.869, 0.942, 0.856, 0.744, 0.723, 0.99, 0.942, 0.981, 0.42, 0.778, 0.914, 0.876, 0.92, 0.12, 0.996, 0.894, 0.986, 0.834, 0.919, 0.954, 0.827, 0.956, 0.928, 0.797, 0.549, 0.893, 0.479, 0.657, 0.809, 0.925, 0.661, 0.853, 0.778, 0.964, 0.907, 0.992, 0.983, 0.967, 0.573, 0.176, 0.993, 0.103, 0.934, 0.78, 0.972, 0.155, 0.891, 0.688, 0.986, 0.601, 0.809, 0.939, 0.498, 0.958, 0.766, 0.916, 0.517, 0.875, 0.143, 0.701, 0.988, 0.844, 0.877, 0.588, 0.733, 0.863, 0.267, 0.609, 0.749, 0.745, 0.9, 0.854, 0.92, 0.996, 0.946, 0.327, 0.945, 0.705, 0.627, 0.152, 0.815, 0.893, 0.687, 0.964, 0.949, 0.519, 0.78, 0.558, 0.202, 0.511, 0.975, 0.96, 0.578, 0.749, 0.973, 0.954, 0.834, 0.322, 0.973, 0.755, 0.785, 0.843, 0.991, 0.998, 0.966, 0.408, 0.587, 0.451, 0.518, 0.726, 0.935, 0.907, 0.439, 0.173, 0.307, 0.886, 0.581, 0.471, 0.99, 0.638]
global origin = 1
global destination = 50