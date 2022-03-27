global arcs = [1 12; 1 17; 1 28; 1 41; 1 42; 1 43; 2 10; 2 29; 2 42; 3 28; 3 29; 4 7; 4 9; 4 15; 4 22; 4 29; 4 39; 4 46; 5 18; 5 25; 6 10; 6 16; 6 32; 6 39; 6 46; 6 47; 7 29; 7 30; 8 6; 8 40; 9 2; 9 6; 9 18; 9 28; 9 32; 10 13; 10 14; 10 27; 10 34; 10 47; 11 2; 11 19; 11 26; 11 44; 12 6; 12 10; 12 35; 12 49; 13 3; 13 19; 13 48; 14 2; 14 3; 14 27; 15 10; 15 35; 15 44; 16 8; 16 14; 16 18; 16 21; 16 27; 16 50; 17 8; 17 10; 17 27; 17 28; 17 48; 18 8; 18 20; 18 25; 18 29; 18 37; 19 24; 19 37; 19 46; 20 34; 20 43; 21 7; 21 14; 21 23; 21 24; 21 38; 21 47; 22 4; 22 9; 22 15; 22 50; 23 8; 23 11; 23 14; 23 22; 23 25; 23 28; 23 31; 23 49; 24 13; 24 26; 25 13; 25 16; 25 20; 25 24; 25 46; 26 33; 26 39; 27 2; 27 13; 27 28; 28 20; 28 21; 28 23; 28 40; 28 47; 29 15; 29 18; 29 22; 29 31; 29 35; 29 43; 29 44; 30 3; 30 18; 30 24; 31 35; 31 45; 31 47; 32 31; 32 46; 32 48; 33 16; 33 20; 33 30; 34 2; 34 27; 34 39; 34 44; 34 45; 34 49; 35 24; 35 34; 36 15; 36 25; 36 44; 36 46; 37 5; 37 8; 37 9; 37 15; 37 16; 37 32; 37 38; 37 46; 38 6; 38 15; 38 25; 38 31; 39 2; 39 3; 39 14; 39 23; 39 24; 39 29; 39 31; 39 34; 40 2; 40 7; 40 12; 40 43; 41 7; 41 19; 41 22; 41 29; 41 33; 42 8; 42 9; 42 30; 42 33; 42 38; 42 46; 43 8; 43 17; 43 33; 43 35; 43 37; 43 41; 44 11; 44 32; 44 39; 44 50; 45 12; 45 15; 45 18; 45 23; 45 25; 45 36; 45 49; 46 9; 46 16; 46 25; 46 36; 47 9; 47 17; 47 18; 48 4; 48 6; 48 17; 48 39; 49 18; 49 19; 49 29; 49 43; 49 50]
global d_x = [4.0, 7.0, 7.0, 10.0, 8.0, 6.0, 5.0, 1.0, 1.0, 7.0, 1.0, 4.0, 10.0, 1.0, 6.0, 3.0, 2.0, 7.0, 2.0, 3.0, 8.0, 6.0, 10.0, 4.0, 10.0, 6.0, 7.0, 6.0, 10.0, 9.0, 5.0, 2.0, 8.0, 5.0, 8.0, 3.0, 6.0, 8.0, 5.0, 8.0, 1.0, 6.0, 10.0, 8.0, 2.0, 2.0, 10.0, 9.0, 7.0, 10.0, 4.0, 1.0, 1.0, 6.0, 5.0, 4.0, 9.0, 8.0, 10.0, 2.0, 5.0, 4.0, 9.0, 1.0, 1.0, 4.0, 4.0, 7.0, 10.0, 2.0, 5.0, 9.0, 8.0, 7.0, 10.0, 2.0, 2.0, 4.0, 7.0, 10.0, 4.0, 3.0, 8.0, 9.0, 3.0, 3.0, 6.0, 10.0, 10.0, 4.0, 2.0, 2.0, 8.0, 5.0, 4.0, 7.0, 7.0, 2.0, 7.0, 1.0, 1.0, 10.0, 9.0, 4.0, 8.0, 2.0, 9.0, 10.0, 4.0, 7.0, 9.0, 9.0, 5.0, 4.0, 10.0, 8.0, 8.0, 9.0, 5.0, 5.0, 1.0, 2.0, 3.0, 8.0, 9.0, 4.0, 10.0, 1.0, 10.0, 7.0, 9.0, 8.0, 9.0, 6.0, 8.0, 9.0, 4.0, 6.0, 7.0, 4.0, 6.0, 8.0, 1.0, 2.0, 9.0, 1.0, 7.0, 7.0, 10.0, 8.0, 3.0, 3.0, 5.0, 2.0, 6.0, 7.0, 6.0, 6.0, 2.0, 1.0, 9.0, 7.0, 5.0, 5.0, 4.0, 6.0, 10.0, 1.0, 6.0, 6.0, 10.0, 8.0, 1.0, 2.0, 6.0, 8.0, 5.0, 5.0, 10.0, 8.0, 9.0, 9.0, 2.0, 3.0, 3.0, 9.0, 3.0, 7.0, 2.0, 3.0, 5.0, 7.0, 10.0, 8.0, 8.0, 4.0, 1.0, 7.0, 6.0, 9.0, 8.0, 1.0, 4.0, 5.0, 2.0, 8.0, 8.0, 1.0, 8.0, 8.0, 1.0, 8.0]
global b_x = 5
global d_y = [5.0, 8.0, 2.0, 4.0, 8.0, 8.0, 8.0, 4.0, 9.0, 10.0, 1.0, 3.0, 3.0, 7.0, 2.0, 9.0, 4.0, 9.0, 8.0, 2.0, 4.0, 3.0, 4.0, 3.0, 7.0, 6.0, 9.0, 3.0, 8.0, 1.0, 1.0, 1.0, 3.0, 9.0, 1.0, 9.0, 6.0, 8.0, 3.0, 6.0, 2.0, 9.0, 6.0, 6.0, 6.0, 8.0, 5.0, 10.0, 5.0, 2.0, 4.0, 4.0, 2.0, 4.0, 2.0, 4.0, 5.0, 6.0, 9.0, 8.0, 6.0, 9.0, 7.0, 7.0, 8.0, 10.0, 6.0, 4.0, 7.0, 8.0, 4.0, 6.0, 10.0, 6.0, 10.0, 9.0, 6.0, 10.0, 5.0, 4.0, 8.0, 1.0, 3.0, 10.0, 6.0, 8.0, 7.0, 7.0, 9.0, 9.0, 9.0, 7.0, 7.0, 7.0, 6.0, 2.0, 6.0, 9.0, 9.0, 3.0, 8.0, 8.0, 2.0, 6.0, 5.0, 10.0, 8.0, 4.0, 8.0, 10.0, 2.0, 7.0, 9.0, 1.0, 6.0, 9.0, 6.0, 8.0, 4.0, 6.0, 5.0, 10.0, 7.0, 1.0, 9.0, 9.0, 5.0, 6.0, 8.0, 9.0, 10.0, 1.0, 5.0, 4.0, 4.0, 4.0, 9.0, 7.0, 4.0, 3.0, 6.0, 3.0, 1.0, 9.0, 10.0, 8.0, 10.0, 2.0, 5.0, 6.0, 8.0, 10.0, 5.0, 1.0, 3.0, 4.0, 4.0, 6.0, 3.0, 1.0, 3.0, 6.0, 9.0, 5.0, 9.0, 3.0, 10.0, 1.0, 2.0, 8.0, 6.0, 9.0, 4.0, 4.0, 3.0, 7.0, 5.0, 3.0, 1.0, 6.0, 4.0, 6.0, 1.0, 3.0, 9.0, 3.0, 1.0, 4.0, 5.0, 1.0, 5.0, 10.0, 9.0, 9.0, 8.0, 8.0, 10.0, 9.0, 10.0, 4.0, 8.0, 10.0, 2.0, 10.0, 9.0, 10.0, 6.0, 9.0, 10.0, 10.0, 10.0, 4.0]
global b_y = 10
global p = [0.387, 0.596, 0.213, 0.153, 0.91, 0.92, 0.591, 0.805, 0.165, 0.581, 0.638, 0.156, 0.451, 0.699, 0.187, 0.764, 0.926, 0.573, 0.208, 0.54, 0.691, 0.445, 0.312, 0.863, 0.937, 0.768, 0.145, 0.931, 0.431, 0.191, 0.13, 0.115, 0.49, 0.389, 0.989, 0.345, 0.801, 0.301, 0.634, 0.264, 0.798, 0.026, 0.054, 0.655, 0.302, 0.21, 0.122, 0.813, 0.866, 0.758, 0.911, 0.764, 0.992, 0.612, 0.325, 0.063, 0.064, 0.216, 0.021, 0.783, 0.914, 0.86, 0.182, 0.687, 0.431, 0.748, 0.649, 0.639, 0.978, 0.659, 0.474, 0.623, 0.012, 0.739, 0.839, 0.996, 0.723, 0.761, 0.492, 0.291, 0.13, 0.909, 0.211, 0.01, 0.175, 0.727, 0.102, 0.535, 0.377, 0.314, 0.155, 0.235, 0.072, 0.566, 0.186, 0.937, 0.461, 0.174, 0.112, 0.181, 0.004, 0.021, 0.183, 0.042, 0.648, 0.813, 0.116, 0.223, 0.074, 0.557, 0.517, 0.694, 0.643, 0.99, 0.635, 0.478, 0.175, 0.007, 0.662, 0.633, 0.21, 0.161, 0.733, 0.009, 0.112, 0.559, 0.29, 0.409, 0.4, 0.931, 0.74, 0.332, 0.845, 0.237, 0.734, 0.028, 0.52, 0.828, 0.26, 0.369, 0.041, 0.574, 0.697, 0.99, 0.781, 0.993, 0.959, 0.335, 0.182, 0.404, 0.788, 0.394, 0.809, 0.495, 0.135, 0.702, 0.336, 0.193, 0.548, 0.233, 0.638, 0.125, 0.402, 0.137, 0.223, 0.551, 0.017, 0.76, 0.829, 0.571, 0.705, 0.136, 0.369, 0.486, 0.402, 0.974, 0.526, 0.169, 0.626, 0.762, 0.691, 0.039, 0.698, 0.338, 0.598, 0.808, 0.782, 0.798, 0.609, 0.758, 0.092, 0.013, 0.207, 0.049, 0.205, 0.549, 0.625, 0.429, 0.279, 0.633, 0.588, 0.751, 0.732, 0.999, 0.497, 0.889, 0.448, 0.642, 0.94, 0.42, 0.662, 0.192]
global q = [0.764, 0.625, 0.951, 0.551, 0.947, 0.997, 0.984, 0.81, 0.236, 0.668, 0.919, 0.638, 0.825, 0.918, 0.684, 0.911, 0.964, 0.776, 0.438, 0.729, 0.81, 0.892, 0.332, 0.916, 0.98, 0.781, 0.854, 0.948, 0.831, 0.641, 0.421, 0.41, 0.79, 0.677, 0.992, 0.487, 0.864, 0.511, 0.805, 0.931, 0.805, 0.716, 0.221, 0.926, 0.514, 0.715, 0.53, 0.966, 0.95, 0.807, 0.96, 0.925, 0.994, 0.843, 0.627, 0.843, 0.791, 0.614, 0.037, 0.884, 0.996, 0.877, 0.928, 0.836, 0.501, 0.899, 0.714, 0.882, 0.989, 0.702, 0.477, 0.751, 0.676, 0.987, 0.966, 0.997, 0.954, 0.987, 0.752, 0.663, 0.972, 0.953, 0.221, 0.209, 0.352, 0.865, 0.501, 0.562, 0.423, 0.839, 0.383, 0.89, 0.882, 0.886, 0.502, 0.974, 0.68, 0.42, 0.217, 0.674, 0.589, 0.335, 0.784, 0.175, 0.924, 0.995, 0.94, 0.806, 0.767, 0.824, 0.988, 0.757, 0.662, 0.991, 0.731, 0.701, 0.801, 0.276, 0.837, 0.828, 0.742, 0.516, 0.871, 0.985, 0.599, 0.559, 0.358, 0.79, 0.953, 0.977, 0.848, 0.852, 0.923, 0.334, 0.906, 0.728, 0.533, 0.858, 0.856, 0.844, 0.594, 0.902, 0.78, 0.991, 0.986, 0.993, 0.995, 0.358, 0.683, 0.606, 0.994, 0.411, 0.949, 0.685, 0.144, 0.775, 0.802, 0.696, 0.76, 0.945, 0.741, 0.511, 0.826, 0.212, 0.714, 0.962, 0.319, 0.992, 0.862, 0.928, 0.99, 0.365, 0.902, 0.737, 0.854, 0.992, 0.743, 0.395, 0.886, 0.947, 0.863, 0.927, 0.973, 0.843, 0.998, 0.99, 0.925, 0.914, 0.866, 0.767, 0.928, 0.535, 0.309, 0.499, 0.233, 0.967, 0.82, 0.944, 0.403, 0.758, 0.992, 0.76, 0.749, 0.999, 0.832, 0.955, 0.672, 0.704, 0.952, 0.442, 0.962, 0.655]
global origin = 1
global destination = 50