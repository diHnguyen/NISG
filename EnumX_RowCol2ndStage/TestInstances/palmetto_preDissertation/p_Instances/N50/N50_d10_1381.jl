global arcs = [1 7; 1 16; 1 17; 1 32; 1 33; 2 4; 2 25; 2 29; 2 42; 3 29; 3 38; 4 8; 4 19; 4 25; 4 39; 5 16; 5 34; 5 36; 6 3; 6 29; 6 39; 6 42; 7 9; 7 12; 7 13; 7 14; 7 22; 7 27; 7 48; 8 6; 8 18; 8 25; 8 38; 8 49; 9 11; 9 18; 9 20; 9 26; 9 37; 10 4; 10 15; 10 25; 10 36; 11 4; 11 7; 11 16; 11 30; 11 38; 11 42; 11 43; 11 44; 11 50; 12 2; 12 14; 12 20; 12 21; 12 26; 12 29; 12 40; 12 41; 12 48; 12 49; 13 3; 13 14; 13 16; 13 30; 13 34; 13 36; 13 43; 14 9; 15 12; 15 18; 15 32; 16 7; 17 7; 17 16; 17 28; 17 45; 18 11; 18 13; 18 16; 18 38; 18 40; 18 42; 19 10; 19 22; 19 24; 19 26; 19 36; 19 39; 20 16; 20 21; 20 27; 20 28; 20 42; 20 50; 21 2; 21 6; 21 8; 21 17; 21 20; 21 24; 22 25; 23 5; 23 6; 23 14; 23 16; 23 24; 23 32; 23 40; 23 50; 24 3; 24 11; 24 22; 24 33; 24 46; 24 49; 25 7; 25 26; 26 7; 26 13; 26 17; 27 12; 27 26; 27 28; 27 39; 27 46; 28 4; 28 5; 28 10; 28 36; 28 37; 29 2; 29 4; 29 5; 29 15; 29 37; 29 42; 30 25; 30 29; 30 44; 31 6; 31 28; 31 33; 31 42; 32 7; 32 13; 33 10; 33 13; 33 18; 33 22; 33 27; 33 30; 33 36; 34 2; 34 42; 34 50; 35 15; 35 24; 35 27; 36 2; 36 5; 36 11; 36 23; 36 24; 36 37; 36 40; 36 42; 36 49; 37 16; 37 35; 37 50; 38 5; 38 25; 38 30; 38 49; 38 50; 39 3; 39 20; 39 43; 40 6; 40 26; 41 17; 41 18; 41 20; 41 33; 41 34; 41 40; 41 47; 42 5; 42 22; 42 23; 42 44; 43 2; 43 9; 43 38; 44 5; 44 16; 44 18; 44 27; 44 36; 44 38; 44 41; 44 43; 44 47; 44 50; 45 19; 45 21; 45 26; 45 30; 45 37; 45 49; 46 12; 46 15; 46 29; 47 4; 47 6; 47 12; 47 15; 47 38; 48 8; 48 31; 48 38; 48 40; 48 42; 48 44; 49 4; 49 5; 49 7; 49 11; 49 16; 49 23; 49 34; 49 43]
global d_x = [6.0, 3.0, 8.0, 1.0, 3.0, 1.0, 5.0, 3.0, 4.0, 1.0, 8.0, 5.0, 6.0, 4.0, 5.0, 3.0, 4.0, 10.0, 6.0, 10.0, 8.0, 2.0, 3.0, 3.0, 5.0, 8.0, 9.0, 1.0, 5.0, 3.0, 4.0, 10.0, 9.0, 2.0, 1.0, 4.0, 10.0, 7.0, 9.0, 8.0, 9.0, 2.0, 1.0, 8.0, 4.0, 9.0, 5.0, 9.0, 4.0, 3.0, 6.0, 9.0, 7.0, 3.0, 2.0, 2.0, 7.0, 10.0, 2.0, 2.0, 4.0, 2.0, 10.0, 3.0, 7.0, 6.0, 8.0, 1.0, 10.0, 9.0, 7.0, 4.0, 4.0, 10.0, 5.0, 8.0, 5.0, 4.0, 7.0, 8.0, 9.0, 4.0, 7.0, 7.0, 10.0, 10.0, 7.0, 2.0, 5.0, 2.0, 2.0, 2.0, 3.0, 10.0, 6.0, 2.0, 8.0, 7.0, 2.0, 3.0, 1.0, 5.0, 4.0, 5.0, 1.0, 9.0, 4.0, 6.0, 5.0, 3.0, 10.0, 2.0, 3.0, 10.0, 5.0, 4.0, 1.0, 9.0, 7.0, 1.0, 2.0, 9.0, 9.0, 9.0, 1.0, 10.0, 8.0, 5.0, 7.0, 9.0, 4.0, 1.0, 8.0, 8.0, 3.0, 7.0, 1.0, 2.0, 6.0, 10.0, 5.0, 10.0, 8.0, 5.0, 8.0, 8.0, 2.0, 6.0, 9.0, 8.0, 4.0, 6.0, 2.0, 7.0, 5.0, 4.0, 10.0, 3.0, 6.0, 6.0, 10.0, 10.0, 9.0, 6.0, 4.0, 8.0, 10.0, 6.0, 9.0, 4.0, 2.0, 9.0, 8.0, 6.0, 4.0, 7.0, 8.0, 9.0, 2.0, 4.0, 10.0, 10.0, 8.0, 5.0, 7.0, 6.0, 10.0, 2.0, 10.0, 6.0, 6.0, 4.0, 5.0, 3.0, 1.0, 10.0, 5.0, 10.0, 4.0, 6.0, 1.0, 2.0, 3.0, 10.0, 4.0, 10.0, 10.0, 3.0, 4.0, 3.0, 8.0, 8.0, 2.0, 8.0, 8.0, 10.0, 9.0, 7.0, 3.0, 10.0, 6.0, 9.0, 3.0, 8.0, 7.0, 3.0, 7.0, 10.0, 10.0, 10.0, 5.0, 10.0, 4.0, 6.0]
global b_x = 5
global d_y = [4.0, 7.0, 8.0, 8.0, 8.0, 9.0, 2.0, 8.0, 6.0, 1.0, 6.0, 9.0, 1.0, 4.0, 3.0, 7.0, 1.0, 6.0, 2.0, 3.0, 8.0, 8.0, 1.0, 5.0, 8.0, 7.0, 1.0, 3.0, 8.0, 2.0, 1.0, 4.0, 9.0, 10.0, 5.0, 3.0, 10.0, 8.0, 10.0, 6.0, 4.0, 4.0, 7.0, 4.0, 1.0, 9.0, 5.0, 9.0, 10.0, 9.0, 3.0, 6.0, 9.0, 5.0, 5.0, 6.0, 2.0, 9.0, 8.0, 2.0, 9.0, 7.0, 10.0, 2.0, 7.0, 10.0, 1.0, 1.0, 6.0, 2.0, 7.0, 3.0, 3.0, 5.0, 3.0, 8.0, 9.0, 1.0, 2.0, 1.0, 5.0, 4.0, 9.0, 9.0, 4.0, 2.0, 1.0, 8.0, 10.0, 6.0, 9.0, 7.0, 4.0, 3.0, 3.0, 8.0, 7.0, 1.0, 8.0, 1.0, 6.0, 4.0, 7.0, 7.0, 10.0, 5.0, 5.0, 6.0, 5.0, 7.0, 5.0, 8.0, 5.0, 9.0, 9.0, 6.0, 6.0, 3.0, 8.0, 10.0, 1.0, 1.0, 9.0, 3.0, 3.0, 2.0, 8.0, 8.0, 5.0, 9.0, 1.0, 8.0, 9.0, 10.0, 10.0, 6.0, 3.0, 2.0, 3.0, 7.0, 2.0, 9.0, 8.0, 3.0, 6.0, 5.0, 5.0, 7.0, 7.0, 3.0, 7.0, 9.0, 9.0, 1.0, 1.0, 5.0, 2.0, 10.0, 3.0, 8.0, 4.0, 3.0, 7.0, 10.0, 6.0, 9.0, 8.0, 5.0, 10.0, 1.0, 10.0, 3.0, 3.0, 8.0, 10.0, 8.0, 8.0, 7.0, 7.0, 7.0, 3.0, 3.0, 10.0, 3.0, 9.0, 9.0, 1.0, 10.0, 5.0, 7.0, 4.0, 4.0, 2.0, 10.0, 6.0, 6.0, 1.0, 8.0, 6.0, 4.0, 7.0, 2.0, 5.0, 9.0, 7.0, 10.0, 10.0, 7.0, 1.0, 7.0, 10.0, 3.0, 3.0, 5.0, 7.0, 3.0, 6.0, 4.0, 3.0, 3.0, 3.0, 7.0, 2.0, 9.0, 6.0, 2.0, 4.0, 8.0, 1.0, 6.0, 3.0, 8.0, 4.0, 7.0]
global b_y = 10
global p = [0.252, 0.258, 0.876, 0.58, 0.105, 0.335, 0.672, 0.362, 0.109, 0.554, 0.474, 0.727, 0.062, 0.699, 0.083, 0.638, 0.237, 0.107, 0.572, 0.415, 0.991, 0.751, 0.94, 0.256, 0.397, 0.144, 0.525, 0.41, 0.142, 0.391, 0.989, 0.851, 0.636, 0.428, 0.245, 0.613, 0.448, 0.801, 0.662, 0.56, 0.053, 0.898, 0.46, 0.87, 0.161, 0.684, 0.402, 0.719, 0.296, 0.682, 0.519, 0.775, 0.191, 0.524, 0.885, 0.618, 0.61, 0.827, 0.805, 0.315, 0.623, 0.599, 0.372, 0.622, 0.787, 0.567, 0.846, 0.437, 0.453, 0.894, 0.51, 0.608, 0.251, 0.358, 0.099, 0.882, 0.944, 0.352, 0.881, 0.943, 0.688, 0.447, 0.518, 0.479, 0.315, 0.594, 0.043, 0.415, 0.183, 0.346, 0.458, 0.94, 0.387, 0.688, 0.502, 0.66, 0.364, 0.062, 0.561, 0.762, 0.105, 0.306, 0.536, 0.542, 0.923, 0.984, 0.986, 0.848, 0.727, 0.912, 0.267, 0.426, 0.834, 0.885, 0.411, 0.241, 0.222, 0.67, 0.654, 0.044, 0.253, 0.023, 0.484, 0.762, 0.785, 0.509, 0.532, 0.905, 0.153, 0.342, 0.487, 0.548, 0.872, 0.375, 0.462, 0.033, 0.739, 0.352, 0.984, 0.639, 0.563, 0.133, 0.55, 0.37, 0.402, 0.057, 0.714, 0.996, 0.045, 0.508, 0.049, 0.895, 0.041, 0.408, 0.965, 0.699, 0.126, 0.453, 0.033, 0.018, 0.279, 0.023, 0.627, 0.967, 0.794, 0.369, 0.691, 0.754, 0.841, 0.924, 0.599, 0.309, 0.469, 0.476, 0.907, 0.383, 0.088, 0.216, 0.047, 0.414, 0.952, 0.586, 0.977, 0.386, 0.868, 0.956, 0.465, 0.46, 0.754, 0.813, 0.618, 0.205, 0.005, 0.349, 0.626, 0.387, 0.044, 0.348, 0.628, 0.325, 0.212, 0.712, 0.981, 0.729, 0.988, 0.182, 0.373, 0.252, 0.012, 0.34, 0.202, 0.571, 0.43, 0.268, 0.846, 0.419, 0.825, 0.099, 0.302, 0.655, 0.706, 0.792, 0.523, 0.322, 0.762, 0.693, 0.163, 0.764, 0.867, 0.666, 0.346, 0.905, 0.084, 0.451]
global q = [0.858, 0.422, 0.94, 0.697, 0.343, 0.914, 0.736, 0.45, 0.274, 0.768, 0.581, 0.785, 0.298, 0.892, 0.945, 0.685, 0.695, 0.901, 0.735, 0.946, 0.999, 0.808, 0.974, 0.583, 0.714, 0.577, 0.97, 0.711, 0.78, 0.998, 0.99, 0.969, 0.989, 0.649, 0.278, 0.64, 0.717, 0.95, 0.885, 0.839, 0.56, 0.904, 0.615, 0.875, 0.923, 0.776, 0.5, 0.985, 0.923, 0.686, 0.93, 0.917, 0.612, 0.565, 0.964, 0.97, 0.629, 0.975, 0.933, 0.616, 0.743, 0.913, 0.732, 0.968, 0.941, 0.698, 0.852, 0.94, 0.609, 0.899, 0.936, 0.692, 0.786, 0.675, 0.525, 0.97, 0.968, 0.767, 0.913, 0.992, 0.821, 0.686, 0.7, 0.785, 0.438, 0.808, 0.074, 0.579, 0.425, 0.37, 0.578, 0.991, 0.57, 0.832, 0.845, 0.782, 0.659, 0.925, 0.766, 0.824, 0.318, 0.668, 0.662, 0.815, 0.934, 0.987, 0.996, 0.951, 0.779, 0.949, 0.305, 0.8, 0.928, 0.941, 0.911, 0.444, 0.968, 0.853, 0.798, 0.056, 0.707, 0.847, 0.997, 0.919, 0.961, 0.609, 0.777, 0.993, 0.177, 0.62, 0.817, 0.552, 0.931, 0.686, 0.828, 0.45, 0.948, 0.669, 0.999, 0.946, 0.576, 0.271, 0.745, 0.765, 0.943, 0.361, 0.991, 0.997, 0.8, 0.555, 0.72, 0.949, 0.582, 0.853, 0.98, 0.783, 0.98, 0.475, 0.512, 0.436, 0.881, 0.229, 0.832, 0.997, 0.996, 0.768, 0.927, 0.767, 0.849, 0.936, 0.939, 0.455, 0.779, 0.651, 0.943, 0.904, 0.796, 0.765, 0.327, 0.776, 0.974, 0.994, 0.984, 0.689, 0.883, 0.979, 0.886, 0.798, 0.764, 0.845, 0.873, 0.64, 0.974, 0.865, 0.99, 0.6, 0.98, 0.801, 0.745, 0.812, 0.8, 0.903, 0.988, 0.767, 0.997, 0.611, 0.719, 0.351, 0.402, 0.425, 0.907, 0.674, 0.696, 0.416, 0.904, 0.654, 0.881, 0.737, 0.852, 0.946, 0.818, 0.814, 0.732, 0.79, 0.895, 0.846, 0.292, 0.791, 0.888, 0.77, 0.713, 0.911, 0.726, 0.854]
global origin = 1
global destination = 50