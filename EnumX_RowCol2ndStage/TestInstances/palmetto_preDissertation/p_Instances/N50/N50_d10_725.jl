global arcs = [1 34; 1 37; 1 45; 1 47; 2 13; 2 26; 2 41; 3 7; 3 24; 3 28; 3 33; 3 38; 3 43; 3 46; 3 47; 4 8; 4 10; 4 26; 4 40; 5 2; 5 11; 5 22; 5 40; 5 42; 5 45; 5 49; 6 2; 6 37; 6 40; 6 42; 7 6; 7 22; 7 25; 7 27; 7 28; 7 32; 7 34; 8 10; 8 19; 8 22; 9 7; 9 20; 9 30; 9 47; 10 9; 10 16; 10 20; 11 2; 11 4; 11 13; 11 19; 11 27; 11 35; 11 50; 12 6; 12 8; 12 14; 12 15; 12 27; 12 33; 12 36; 12 44; 12 48; 13 15; 13 35; 13 44; 14 2; 14 15; 14 18; 14 28; 15 4; 15 22; 15 43; 15 45; 16 21; 16 24; 16 29; 16 37; 16 42; 16 45; 16 49; 17 6; 17 9; 17 22; 17 27; 17 37; 18 20; 18 23; 18 50; 19 9; 19 17; 19 20; 19 27; 19 29; 19 36; 19 39; 20 7; 20 48; 21 33; 21 42; 21 48; 22 8; 22 10; 22 28; 22 39; 22 42; 23 11; 23 19; 23 32; 23 43; 23 45; 23 49; 24 8; 24 10; 24 25; 24 26; 24 38; 24 43; 24 44; 25 3; 25 8; 25 9; 25 16; 25 41; 25 44; 26 2; 26 13; 26 18; 26 23; 26 28; 26 29; 27 15; 27 29; 27 32; 28 7; 28 9; 28 11; 28 35; 29 7; 29 12; 29 30; 29 32; 29 33; 29 42; 30 29; 30 35; 30 39; 30 41; 30 42; 30 45; 31 5; 31 22; 31 34; 31 36; 31 37; 31 49; 32 8; 32 38; 33 3; 33 31; 33 34; 33 35; 33 45; 34 4; 34 6; 34 12; 34 23; 34 25; 34 38; 34 41; 35 10; 35 14; 35 26; 35 39; 35 40; 36 28; 36 46; 37 6; 37 27; 37 33; 37 38; 37 45; 37 46; 38 2; 38 10; 38 14; 38 20; 38 24; 38 30; 38 50; 39 4; 39 21; 39 45; 40 6; 40 21; 40 26; 41 13; 41 27; 41 32; 41 47; 42 6; 42 10; 42 12; 42 18; 42 22; 42 46; 42 50; 43 6; 43 27; 43 29; 43 34; 43 36; 43 39; 43 45; 44 7; 44 11; 44 14; 44 15; 44 42; 45 14; 45 47; 46 2; 46 17; 46 32; 46 35; 46 42; 46 50; 47 12; 47 13; 47 36; 47 46; 48 8; 48 9; 48 14; 48 20; 48 22; 48 36; 48 43; 48 44; 49 8; 49 12; 49 22; 49 24; 49 31]
global d_x = [3.0, 3.0, 1.0, 4.0, 2.0, 10.0, 9.0, 10.0, 2.0, 7.0, 10.0, 8.0, 1.0, 9.0, 6.0, 8.0, 2.0, 3.0, 9.0, 2.0, 6.0, 7.0, 6.0, 3.0, 2.0, 1.0, 7.0, 3.0, 5.0, 5.0, 4.0, 8.0, 5.0, 4.0, 9.0, 10.0, 6.0, 9.0, 5.0, 6.0, 9.0, 2.0, 4.0, 7.0, 3.0, 3.0, 10.0, 4.0, 6.0, 5.0, 1.0, 6.0, 2.0, 3.0, 5.0, 10.0, 9.0, 2.0, 5.0, 7.0, 6.0, 6.0, 7.0, 3.0, 8.0, 10.0, 2.0, 4.0, 6.0, 6.0, 9.0, 10.0, 8.0, 7.0, 4.0, 10.0, 1.0, 5.0, 6.0, 8.0, 3.0, 6.0, 7.0, 5.0, 10.0, 6.0, 9.0, 3.0, 6.0, 5.0, 5.0, 8.0, 6.0, 9.0, 8.0, 4.0, 6.0, 8.0, 2.0, 2.0, 5.0, 10.0, 3.0, 7.0, 9.0, 3.0, 8.0, 10.0, 3.0, 7.0, 1.0, 2.0, 7.0, 7.0, 6.0, 4.0, 2.0, 8.0, 4.0, 6.0, 10.0, 6.0, 3.0, 5.0, 7.0, 10.0, 4.0, 8.0, 8.0, 6.0, 2.0, 2.0, 1.0, 10.0, 10.0, 7.0, 5.0, 2.0, 3.0, 5.0, 5.0, 1.0, 3.0, 8.0, 10.0, 1.0, 8.0, 10.0, 5.0, 4.0, 10.0, 4.0, 10.0, 1.0, 10.0, 5.0, 7.0, 6.0, 10.0, 4.0, 9.0, 3.0, 2.0, 8.0, 10.0, 6.0, 8.0, 8.0, 1.0, 1.0, 5.0, 6.0, 4.0, 9.0, 9.0, 1.0, 7.0, 3.0, 6.0, 7.0, 2.0, 2.0, 8.0, 2.0, 5.0, 8.0, 3.0, 4.0, 9.0, 8.0, 9.0, 2.0, 2.0, 5.0, 1.0, 9.0, 1.0, 6.0, 7.0, 8.0, 8.0, 7.0, 10.0, 5.0, 3.0, 5.0, 6.0, 3.0, 6.0, 4.0, 8.0, 6.0, 6.0, 1.0, 10.0, 2.0, 9.0, 1.0, 3.0, 6.0, 3.0, 9.0, 4.0, 1.0, 6.0, 10.0, 5.0, 5.0, 6.0, 1.0, 1.0, 9.0, 6.0, 8.0, 6.0, 2.0, 5.0, 8.0, 2.0, 10.0, 3.0, 9.0, 8.0, 3.0]
global b_x = 5
global d_y = [3.0, 4.0, 1.0, 10.0, 1.0, 2.0, 1.0, 4.0, 6.0, 6.0, 3.0, 8.0, 9.0, 2.0, 8.0, 3.0, 1.0, 3.0, 4.0, 5.0, 10.0, 8.0, 10.0, 2.0, 7.0, 8.0, 1.0, 5.0, 7.0, 9.0, 4.0, 4.0, 6.0, 5.0, 7.0, 10.0, 5.0, 8.0, 6.0, 3.0, 6.0, 1.0, 1.0, 5.0, 10.0, 6.0, 4.0, 8.0, 1.0, 5.0, 7.0, 3.0, 9.0, 4.0, 7.0, 3.0, 2.0, 2.0, 1.0, 3.0, 4.0, 2.0, 3.0, 6.0, 3.0, 7.0, 2.0, 4.0, 3.0, 1.0, 3.0, 4.0, 4.0, 10.0, 9.0, 3.0, 2.0, 10.0, 10.0, 5.0, 3.0, 1.0, 9.0, 5.0, 5.0, 6.0, 6.0, 4.0, 10.0, 3.0, 5.0, 6.0, 8.0, 8.0, 6.0, 2.0, 4.0, 3.0, 1.0, 4.0, 3.0, 5.0, 1.0, 9.0, 3.0, 6.0, 2.0, 5.0, 4.0, 9.0, 8.0, 6.0, 7.0, 1.0, 2.0, 4.0, 4.0, 1.0, 4.0, 8.0, 6.0, 1.0, 9.0, 3.0, 5.0, 10.0, 10.0, 7.0, 5.0, 4.0, 3.0, 1.0, 1.0, 1.0, 7.0, 5.0, 10.0, 8.0, 7.0, 10.0, 4.0, 3.0, 3.0, 9.0, 7.0, 8.0, 7.0, 9.0, 9.0, 8.0, 5.0, 5.0, 5.0, 1.0, 8.0, 10.0, 10.0, 1.0, 8.0, 9.0, 6.0, 2.0, 3.0, 8.0, 6.0, 3.0, 5.0, 5.0, 10.0, 2.0, 5.0, 2.0, 7.0, 7.0, 2.0, 1.0, 1.0, 4.0, 9.0, 5.0, 9.0, 9.0, 10.0, 1.0, 5.0, 4.0, 2.0, 9.0, 3.0, 1.0, 3.0, 7.0, 7.0, 1.0, 4.0, 5.0, 10.0, 6.0, 6.0, 7.0, 1.0, 8.0, 2.0, 7.0, 6.0, 3.0, 2.0, 8.0, 9.0, 7.0, 1.0, 9.0, 5.0, 7.0, 3.0, 3.0, 10.0, 10.0, 6.0, 4.0, 8.0, 7.0, 9.0, 1.0, 6.0, 5.0, 3.0, 10.0, 3.0, 2.0, 10.0, 7.0, 1.0, 10.0, 3.0, 5.0, 2.0, 6.0, 2.0, 4.0, 1.0, 5.0, 6.0, 3.0]
global b_y = 10
global p = [0.066, 0.097, 0.538, 0.061, 0.442, 0.543, 0.659, 0.826, 0.044, 0.12, 0.748, 0.449, 0.844, 0.475, 0.071, 0.347, 0.314, 0.789, 0.02, 0.164, 0.298, 0.647, 0.591, 0.684, 0.47, 0.853, 0.859, 0.66, 0.532, 0.165, 0.733, 0.461, 0.559, 0.609, 0.333, 0.213, 0.597, 0.752, 0.5, 0.331, 0.408, 0.752, 0.701, 0.903, 0.703, 0.156, 0.43, 0.287, 0.144, 0.121, 0.528, 0.394, 0.783, 0.584, 0.106, 0.294, 0.438, 0.057, 0.263, 0.575, 0.49, 0.571, 0.331, 0.826, 0.946, 0.599, 0.502, 0.1, 0.435, 0.997, 0.036, 0.769, 0.697, 0.777, 0.054, 0.698, 0.323, 0.235, 0.374, 0.784, 0.504, 0.97, 0.016, 0.319, 0.57, 0.916, 0.431, 0.481, 0.098, 0.698, 0.924, 0.894, 0.846, 0.157, 0.005, 0.527, 0.914, 0.612, 0.998, 0.015, 0.014, 0.301, 0.488, 0.033, 0.164, 0.762, 0.58, 0.801, 0.827, 0.98, 0.057, 0.585, 0.156, 0.157, 0.888, 0.218, 0.063, 0.108, 0.037, 0.815, 0.946, 0.519, 0.143, 0.901, 0.709, 0.204, 0.492, 0.979, 0.638, 0.24, 0.079, 0.917, 0.674, 0.167, 0.705, 0.338, 0.145, 0.037, 0.913, 0.051, 0.763, 0.672, 0.951, 0.223, 0.59, 0.426, 0.07, 0.132, 0.094, 0.506, 0.04, 0.081, 0.136, 0.335, 0.423, 0.404, 0.764, 0.09, 0.952, 0.266, 0.295, 0.412, 0.705, 0.677, 0.6, 0.766, 0.278, 0.682, 0.809, 0.531, 0.772, 0.699, 0.652, 0.406, 0.665, 0.621, 0.188, 0.499, 0.853, 0.203, 0.845, 0.221, 0.08, 0.635, 0.187, 0.376, 0.411, 0.574, 0.287, 0.229, 0.789, 0.223, 0.384, 0.006, 0.997, 0.953, 0.115, 0.787, 0.494, 0.492, 0.31, 0.061, 0.668, 0.787, 0.198, 0.965, 0.521, 0.9, 0.203, 0.457, 0.212, 0.9, 0.767, 0.517, 0.252, 0.291, 0.325, 0.993, 0.513, 0.683, 0.955, 0.07, 0.557, 0.193, 0.572, 0.687, 0.794, 0.932, 0.631, 0.959, 0.327, 0.475, 0.076, 0.708, 0.021, 0.867, 0.818, 0.63, 0.666, 0.206, 0.845, 0.615, 0.709, 0.648]
global q = [0.089, 0.646, 0.734, 0.838, 0.693, 0.837, 0.844, 0.997, 0.444, 0.845, 0.905, 0.736, 0.893, 0.578, 0.687, 0.347, 0.67, 0.921, 0.442, 0.425, 0.346, 0.752, 0.857, 0.765, 0.64, 0.954, 0.874, 0.944, 0.734, 0.731, 0.982, 0.717, 0.925, 0.987, 0.792, 0.311, 0.984, 0.883, 0.72, 0.724, 0.933, 0.859, 0.963, 0.944, 0.887, 0.325, 0.89, 0.878, 0.948, 0.404, 0.633, 0.728, 0.803, 0.676, 0.143, 0.734, 0.705, 0.859, 0.858, 0.635, 0.569, 0.932, 0.801, 0.895, 0.97, 0.761, 0.718, 0.302, 0.749, 0.998, 0.184, 0.777, 0.71, 0.932, 0.808, 0.979, 0.544, 0.682, 0.836, 0.812, 0.93, 0.972, 0.683, 0.481, 0.84, 0.992, 0.634, 0.655, 0.437, 0.845, 0.927, 0.947, 0.945, 0.9, 0.951, 0.782, 0.984, 0.918, 0.998, 0.865, 0.245, 0.756, 0.532, 0.587, 0.451, 0.788, 0.677, 0.907, 0.876, 0.982, 0.75, 0.819, 0.981, 0.692, 0.949, 0.883, 0.346, 0.96, 0.676, 0.974, 0.97, 0.526, 0.542, 0.918, 0.839, 0.208, 0.668, 0.98, 0.87, 0.984, 0.854, 0.928, 0.758, 0.339, 0.832, 0.984, 0.189, 0.633, 0.96, 0.429, 0.912, 0.982, 0.951, 0.75, 0.739, 0.747, 0.471, 0.664, 0.097, 0.937, 0.068, 0.794, 0.381, 0.776, 0.894, 0.903, 0.888, 0.687, 0.998, 0.982, 0.891, 0.634, 0.958, 0.936, 0.923, 0.989, 0.944, 0.722, 0.96, 0.684, 0.955, 0.836, 0.777, 0.812, 0.743, 0.717, 0.556, 0.516, 0.894, 0.719, 0.968, 0.85, 0.149, 0.814, 0.304, 0.412, 0.565, 0.59, 0.804, 0.752, 0.89, 0.48, 0.862, 0.437, 0.998, 0.989, 0.244, 0.966, 0.981, 0.823, 0.909, 0.617, 0.854, 0.865, 0.5, 0.989, 0.842, 0.998, 0.42, 0.987, 0.256, 0.979, 0.875, 0.625, 0.958, 0.587, 0.56, 0.995, 0.551, 0.722, 0.965, 0.369, 0.571, 0.237, 0.916, 0.71, 0.842, 0.953, 0.863, 0.963, 0.72, 0.546, 0.272, 0.823, 0.946, 0.899, 0.825, 0.997, 0.913, 0.817, 0.873, 0.695, 0.949, 0.692]
global origin = 1
global destination = 50