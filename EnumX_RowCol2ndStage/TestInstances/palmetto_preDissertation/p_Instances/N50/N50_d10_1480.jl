global arcs = [1 13; 1 18; 1 26; 1 35; 1 39; 2 3; 2 7; 2 8; 2 28; 2 31; 2 35; 2 43; 2 45; 3 4; 3 21; 3 35; 3 40; 3 42; 3 46; 4 2; 4 26; 4 35; 4 36; 4 37; 4 41; 5 8; 5 30; 5 34; 5 50; 6 2; 6 5; 6 18; 7 5; 7 6; 7 8; 7 17; 7 23; 7 24; 7 42; 8 6; 8 24; 8 35; 8 36; 8 46; 8 49; 9 7; 9 11; 9 16; 9 29; 9 30; 9 37; 9 38; 9 43; 9 50; 10 12; 10 28; 10 31; 11 3; 11 13; 11 16; 11 25; 11 39; 11 42; 11 43; 11 49; 12 3; 12 7; 12 8; 12 15; 12 17; 12 22; 12 23; 12 24; 12 31; 12 47; 12 49; 13 5; 13 20; 13 25; 13 29; 13 45; 13 49; 14 5; 14 36; 14 40; 14 43; 14 50; 15 20; 15 22; 15 25; 16 5; 16 18; 16 29; 16 35; 16 48; 17 37; 17 49; 18 9; 18 15; 18 20; 18 26; 18 38; 18 41; 18 45; 19 21; 19 29; 19 37; 20 12; 20 41; 20 44; 20 48; 21 10; 21 16; 21 28; 21 29; 21 33; 21 39; 22 3; 22 12; 22 15; 22 19; 22 28; 22 31; 22 34; 22 47; 23 2; 23 6; 23 27; 23 28; 23 47; 23 50; 24 12; 24 18; 24 35; 24 44; 25 15; 25 38; 25 39; 26 4; 26 6; 26 13; 26 15; 26 19; 26 27; 26 28; 26 46; 27 11; 27 32; 27 43; 27 48; 28 11; 28 15; 28 17; 28 20; 28 32; 28 33; 28 36; 28 45; 29 2; 29 22; 29 33; 29 35; 29 42; 29 43; 30 4; 30 16; 30 32; 30 50; 31 6; 31 10; 31 34; 31 38; 31 40; 31 45; 31 46; 31 48; 32 14; 32 35; 32 43; 32 44; 33 3; 33 7; 33 26; 33 27; 33 29; 33 35; 34 10; 34 23; 34 26; 35 33; 35 34; 35 39; 35 42; 35 43; 35 46; 36 8; 36 23; 36 27; 36 29; 36 33; 36 43; 37 11; 37 20; 37 45; 37 47; 37 50; 38 4; 38 12; 38 14; 38 21; 38 23; 38 39; 39 8; 39 24; 39 27; 39 37; 39 46; 40 15; 40 24; 40 35; 41 8; 41 14; 41 35; 41 48; 42 8; 42 11; 42 12; 42 33; 43 2; 43 19; 43 22; 43 23; 43 29; 43 38; 43 45; 44 7; 44 29; 44 32; 44 47; 45 15; 45 27; 45 33; 45 48; 46 8; 46 14; 46 28; 46 31; 46 37; 46 49; 47 14; 47 39; 47 43; 47 44; 48 2; 48 14; 48 27; 48 31; 48 47; 49 17; 49 25; 49 27; 49 31; 49 32; 49 35; 49 36; 49 41; 49 46]
global d_x = [8.0, 5.0, 8.0, 9.0, 2.0, 6.0, 6.0, 8.0, 3.0, 7.0, 4.0, 4.0, 5.0, 8.0, 5.0, 4.0, 6.0, 10.0, 3.0, 2.0, 9.0, 10.0, 7.0, 8.0, 8.0, 2.0, 10.0, 10.0, 1.0, 3.0, 8.0, 10.0, 9.0, 5.0, 7.0, 10.0, 4.0, 2.0, 5.0, 6.0, 5.0, 6.0, 6.0, 8.0, 8.0, 1.0, 7.0, 6.0, 9.0, 4.0, 7.0, 2.0, 8.0, 4.0, 7.0, 10.0, 1.0, 8.0, 2.0, 1.0, 4.0, 9.0, 4.0, 10.0, 8.0, 10.0, 6.0, 10.0, 4.0, 8.0, 2.0, 8.0, 3.0, 10.0, 8.0, 7.0, 10.0, 10.0, 1.0, 9.0, 6.0, 8.0, 5.0, 6.0, 1.0, 8.0, 1.0, 10.0, 7.0, 4.0, 2.0, 2.0, 2.0, 8.0, 2.0, 5.0, 4.0, 5.0, 1.0, 1.0, 9.0, 3.0, 1.0, 10.0, 2.0, 5.0, 4.0, 6.0, 7.0, 6.0, 6.0, 9.0, 5.0, 6.0, 10.0, 8.0, 4.0, 2.0, 6.0, 3.0, 8.0, 2.0, 5.0, 6.0, 7.0, 8.0, 10.0, 6.0, 5.0, 4.0, 2.0, 10.0, 8.0, 8.0, 9.0, 9.0, 3.0, 5.0, 3.0, 7.0, 3.0, 5.0, 1.0, 1.0, 6.0, 5.0, 7.0, 9.0, 7.0, 8.0, 2.0, 9.0, 5.0, 10.0, 4.0, 4.0, 5.0, 3.0, 7.0, 3.0, 4.0, 6.0, 2.0, 3.0, 4.0, 10.0, 8.0, 10.0, 4.0, 2.0, 5.0, 3.0, 6.0, 2.0, 6.0, 1.0, 6.0, 9.0, 9.0, 2.0, 9.0, 8.0, 6.0, 8.0, 5.0, 9.0, 3.0, 8.0, 7.0, 4.0, 9.0, 10.0, 4.0, 7.0, 9.0, 5.0, 7.0, 2.0, 2.0, 1.0, 3.0, 2.0, 10.0, 4.0, 6.0, 4.0, 3.0, 9.0, 4.0, 9.0, 5.0, 3.0, 7.0, 10.0, 6.0, 9.0, 10.0, 6.0, 10.0, 8.0, 6.0, 2.0, 1.0, 4.0, 5.0, 8.0, 1.0, 1.0, 5.0, 2.0, 8.0, 3.0, 4.0, 5.0, 7.0, 9.0, 5.0, 6.0, 2.0, 6.0, 1.0, 9.0, 10.0, 6.0, 2.0, 5.0, 8.0, 2.0, 6.0, 7.0, 4.0, 10.0, 1.0, 1.0, 5.0, 3.0, 5.0, 8.0, 5.0, 9.0, 8.0, 1.0, 9.0, 5.0, 6.0, 7.0, 3.0]
global b_x = 5
global d_y = [5.0, 9.0, 2.0, 5.0, 6.0, 6.0, 6.0, 9.0, 6.0, 3.0, 6.0, 3.0, 1.0, 2.0, 4.0, 8.0, 4.0, 9.0, 1.0, 5.0, 5.0, 8.0, 6.0, 8.0, 6.0, 4.0, 5.0, 1.0, 9.0, 6.0, 10.0, 5.0, 8.0, 2.0, 9.0, 7.0, 7.0, 1.0, 8.0, 10.0, 5.0, 10.0, 9.0, 3.0, 8.0, 1.0, 7.0, 3.0, 6.0, 3.0, 9.0, 4.0, 1.0, 7.0, 8.0, 2.0, 7.0, 9.0, 4.0, 4.0, 5.0, 6.0, 1.0, 8.0, 3.0, 2.0, 4.0, 3.0, 9.0, 9.0, 8.0, 1.0, 4.0, 1.0, 5.0, 6.0, 6.0, 7.0, 8.0, 2.0, 9.0, 10.0, 9.0, 7.0, 10.0, 5.0, 1.0, 4.0, 4.0, 5.0, 5.0, 1.0, 4.0, 3.0, 10.0, 2.0, 1.0, 3.0, 10.0, 2.0, 6.0, 2.0, 5.0, 8.0, 4.0, 10.0, 4.0, 2.0, 5.0, 2.0, 2.0, 7.0, 7.0, 5.0, 7.0, 5.0, 10.0, 6.0, 6.0, 9.0, 2.0, 10.0, 5.0, 5.0, 9.0, 4.0, 4.0, 7.0, 4.0, 4.0, 8.0, 2.0, 7.0, 7.0, 5.0, 9.0, 5.0, 3.0, 1.0, 6.0, 5.0, 9.0, 6.0, 3.0, 2.0, 7.0, 1.0, 10.0, 8.0, 10.0, 4.0, 4.0, 10.0, 3.0, 8.0, 10.0, 1.0, 10.0, 6.0, 7.0, 9.0, 8.0, 6.0, 3.0, 8.0, 6.0, 10.0, 4.0, 4.0, 3.0, 2.0, 1.0, 3.0, 8.0, 3.0, 7.0, 10.0, 3.0, 7.0, 7.0, 6.0, 9.0, 5.0, 9.0, 2.0, 10.0, 7.0, 10.0, 10.0, 2.0, 8.0, 2.0, 1.0, 10.0, 7.0, 8.0, 7.0, 8.0, 7.0, 3.0, 1.0, 5.0, 5.0, 10.0, 9.0, 9.0, 4.0, 2.0, 6.0, 3.0, 10.0, 4.0, 8.0, 8.0, 3.0, 8.0, 1.0, 5.0, 8.0, 7.0, 9.0, 2.0, 7.0, 8.0, 4.0, 1.0, 9.0, 3.0, 4.0, 5.0, 7.0, 1.0, 10.0, 2.0, 9.0, 2.0, 10.0, 9.0, 8.0, 7.0, 3.0, 10.0, 8.0, 9.0, 2.0, 7.0, 3.0, 5.0, 2.0, 8.0, 4.0, 8.0, 6.0, 3.0, 1.0, 7.0, 2.0, 2.0, 8.0, 3.0, 4.0, 5.0, 3.0, 4.0, 1.0, 1.0, 1.0]
global b_y = 10
global p = [0.333, 0.504, 0.514, 0.343, 0.236, 0.299, 0.821, 0.942, 0.265, 0.834, 0.354, 0.959, 0.128, 0.124, 0.715, 0.845, 0.283, 0.816, 0.569, 0.169, 0.217, 0.294, 0.063, 0.655, 0.193, 0.261, 0.988, 0.415, 0.907, 0.337, 0.566, 0.023, 0.858, 0.555, 0.033, 0.992, 0.393, 0.585, 0.461, 0.66, 0.835, 0.152, 0.911, 0.42, 0.651, 0.866, 0.581, 0.755, 0.325, 0.758, 0.942, 0.957, 0.47, 0.942, 0.004, 0.307, 0.792, 0.697, 0.346, 0.252, 0.789, 0.394, 0.333, 0.338, 0.933, 0.808, 0.272, 0.804, 0.131, 0.877, 0.157, 0.256, 0.126, 0.086, 0.994, 0.834, 0.789, 0.227, 0.922, 0.942, 0.32, 0.762, 0.182, 0.174, 0.233, 0.95, 0.487, 0.052, 0.241, 0.152, 0.422, 0.623, 0.084, 0.6, 0.111, 0.019, 0.734, 0.321, 0.002, 0.741, 0.751, 0.586, 0.72, 0.186, 0.381, 0.768, 0.76, 0.019, 0.66, 0.616, 0.311, 0.694, 0.546, 0.317, 0.557, 0.985, 0.041, 0.175, 0.795, 0.9, 0.597, 0.111, 0.677, 0.115, 0.197, 0.186, 0.392, 0.212, 0.056, 0.922, 0.4, 0.714, 0.001, 0.148, 0.023, 0.488, 0.823, 0.991, 0.921, 0.571, 0.911, 0.205, 0.283, 0.256, 0.639, 0.455, 0.118, 0.451, 0.873, 0.344, 0.341, 0.763, 0.175, 0.961, 0.875, 0.955, 0.917, 0.548, 0.239, 0.966, 0.431, 0.768, 0.675, 0.719, 0.955, 0.281, 0.076, 0.794, 0.872, 0.854, 0.848, 0.225, 0.891, 0.417, 0.791, 0.86, 0.564, 0.136, 0.927, 0.876, 0.156, 0.199, 0.577, 0.224, 0.014, 0.992, 0.701, 0.718, 0.225, 0.721, 0.81, 0.835, 0.565, 0.23, 0.145, 0.748, 0.981, 0.118, 0.941, 0.737, 0.793, 0.091, 0.614, 0.835, 0.883, 0.912, 0.193, 0.707, 0.257, 0.507, 0.19, 0.528, 0.858, 0.479, 0.554, 0.624, 0.957, 0.413, 0.459, 0.638, 0.116, 0.817, 0.871, 0.274, 0.044, 0.024, 0.467, 0.64, 0.087, 0.515, 0.932, 0.636, 0.986, 0.286, 0.673, 0.042, 0.984, 0.245, 0.415, 0.687, 0.2, 0.144, 0.334, 0.48, 0.819, 0.286, 0.316, 0.868, 0.12, 0.997, 0.765, 0.403, 0.026, 0.108, 0.07, 0.959, 0.612, 0.702, 0.496, 0.579, 0.847, 0.357, 0.335, 0.308, 0.807, 0.349, 0.173]
global q = [0.357, 0.61, 0.849, 0.612, 0.756, 0.446, 0.972, 0.976, 0.357, 0.863, 0.729, 0.962, 0.305, 0.225, 0.93, 0.894, 0.988, 0.943, 0.8, 0.715, 0.692, 0.835, 0.111, 0.869, 0.41, 0.6, 0.989, 0.573, 0.971, 0.443, 0.997, 0.671, 0.98, 0.62, 0.217, 0.997, 0.564, 0.754, 0.956, 0.971, 0.913, 0.316, 0.989, 0.923, 0.774, 0.909, 0.591, 0.913, 0.988, 0.968, 0.998, 0.976, 0.944, 0.973, 0.451, 0.98, 0.795, 0.708, 0.375, 0.409, 0.875, 0.986, 0.656, 0.823, 0.965, 0.932, 0.361, 0.962, 0.162, 0.98, 0.738, 0.678, 0.822, 0.61, 0.996, 0.852, 0.986, 0.778, 0.944, 0.998, 0.753, 0.778, 0.455, 0.945, 0.796, 0.999, 0.89, 0.71, 0.951, 0.204, 0.842, 0.744, 0.581, 0.66, 0.821, 0.926, 0.937, 0.791, 0.278, 0.877, 0.935, 0.814, 0.796, 0.818, 0.5, 0.822, 0.933, 0.981, 0.854, 0.831, 0.694, 0.827, 0.772, 0.832, 0.911, 0.994, 0.813, 0.541, 0.985, 0.938, 0.715, 0.877, 0.727, 0.119, 0.662, 0.463, 0.455, 0.99, 0.401, 0.989, 0.905, 0.838, 0.353, 0.259, 0.9, 0.508, 0.883, 0.998, 0.925, 0.877, 0.945, 0.649, 0.749, 0.464, 0.703, 0.876, 0.966, 0.736, 0.875, 0.959, 0.638, 0.944, 0.318, 0.964, 0.929, 0.98, 0.963, 0.574, 0.383, 0.999, 0.699, 0.804, 0.766, 0.749, 0.992, 0.958, 0.232, 0.961, 0.912, 0.942, 0.984, 0.475, 0.955, 0.479, 0.824, 0.947, 0.813, 0.729, 0.933, 0.921, 0.183, 0.334, 0.804, 0.271, 0.034, 0.999, 0.943, 0.801, 0.928, 0.915, 0.919, 0.986, 0.724, 0.63, 0.614, 0.985, 0.991, 0.967, 0.942, 0.873, 0.881, 0.15, 0.739, 0.969, 0.972, 0.968, 0.369, 0.993, 0.326, 0.78, 0.613, 0.533, 0.994, 0.975, 0.57, 0.633, 0.958, 0.719, 0.884, 0.876, 0.712, 0.86, 0.905, 0.575, 0.661, 0.613, 0.588, 0.876, 0.621, 0.835, 0.94, 0.997, 0.999, 0.844, 0.944, 0.736, 0.987, 0.542, 0.438, 0.885, 0.598, 0.854, 0.647, 0.95, 0.86, 0.859, 0.573, 0.883, 0.657, 0.997, 0.881, 0.966, 0.136, 0.172, 0.134, 0.964, 0.703, 0.744, 0.782, 0.727, 0.934, 0.729, 0.339, 0.377, 0.925, 0.481, 0.661]
global origin = 1
global destination = 50