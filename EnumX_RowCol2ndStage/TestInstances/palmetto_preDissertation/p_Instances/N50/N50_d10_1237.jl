global arcs = [1 22; 1 25; 1 30; 1 34; 1 35; 2 4; 2 23; 2 24; 2 40; 2 48; 3 11; 3 26; 3 35; 4 10; 4 46; 4 47; 4 49; 5 6; 5 7; 5 10; 5 25; 5 27; 5 42; 5 43; 5 44; 5 49; 6 4; 6 11; 6 20; 6 34; 7 3; 7 40; 8 2; 8 15; 8 19; 8 36; 9 3; 9 15; 9 26; 10 6; 10 9; 10 20; 10 31; 10 42; 11 2; 11 4; 11 23; 11 25; 11 47; 11 49; 12 9; 12 34; 12 37; 13 15; 13 16; 13 38; 13 45; 13 47; 14 11; 14 22; 14 25; 14 27; 14 32; 14 37; 14 39; 14 48; 14 49; 15 11; 15 18; 15 34; 16 31; 16 37; 16 45; 16 46; 17 4; 17 14; 17 23; 17 33; 17 41; 17 45; 18 4; 18 12; 18 15; 18 29; 18 35; 18 44; 19 12; 19 16; 19 39; 19 46; 19 47; 20 26; 20 35; 20 38; 20 39; 21 8; 21 10; 21 28; 21 33; 21 38; 22 13; 22 27; 22 37; 22 49; 23 11; 23 16; 23 26; 23 34; 24 9; 24 16; 24 19; 24 20; 24 32; 24 40; 25 5; 25 10; 25 15; 25 32; 25 34; 25 38; 26 10; 26 29; 26 45; 26 50; 27 26; 27 35; 27 36; 28 4; 28 11; 28 22; 28 23; 29 28; 29 33; 29 41; 30 7; 30 9; 30 13; 30 14; 30 22; 30 25; 30 39; 31 4; 31 7; 31 9; 31 10; 31 12; 31 13; 31 16; 31 29; 32 9; 32 13; 32 19; 32 30; 32 45; 33 11; 33 13; 33 16; 33 36; 34 12; 34 15; 34 21; 34 22; 35 9; 35 24; 35 29; 35 40; 35 47; 36 9; 36 11; 36 14; 36 35; 36 44; 37 7; 37 8; 37 10; 37 13; 37 17; 37 28; 38 6; 38 14; 38 26; 38 28; 38 43; 38 46; 39 17; 39 23; 39 26; 39 28; 39 30; 39 31; 39 43; 40 6; 40 13; 40 18; 40 24; 40 35; 40 36; 40 39; 40 47; 41 2; 41 16; 41 39; 41 43; 41 46; 42 4; 42 32; 42 33; 42 40; 42 45; 43 11; 43 19; 43 27; 44 7; 44 10; 44 24; 44 25; 44 31; 44 37; 44 38; 45 5; 45 11; 45 44; 46 5; 46 12; 46 25; 46 28; 46 45; 47 9; 47 15; 47 26; 47 41; 47 44; 48 10; 48 17; 48 32; 48 34; 48 41; 49 21; 49 28; 49 30]
global d_x = [8.0, 1.0, 10.0, 9.0, 7.0, 6.0, 1.0, 8.0, 7.0, 2.0, 7.0, 3.0, 8.0, 10.0, 9.0, 6.0, 7.0, 7.0, 1.0, 4.0, 9.0, 10.0, 3.0, 1.0, 9.0, 8.0, 3.0, 10.0, 8.0, 10.0, 3.0, 6.0, 1.0, 1.0, 4.0, 3.0, 4.0, 2.0, 6.0, 8.0, 9.0, 10.0, 6.0, 9.0, 8.0, 10.0, 5.0, 2.0, 4.0, 7.0, 10.0, 2.0, 7.0, 8.0, 3.0, 5.0, 2.0, 4.0, 10.0, 3.0, 6.0, 9.0, 10.0, 10.0, 4.0, 8.0, 4.0, 1.0, 2.0, 1.0, 8.0, 8.0, 8.0, 10.0, 5.0, 7.0, 10.0, 2.0, 3.0, 10.0, 1.0, 2.0, 7.0, 3.0, 3.0, 6.0, 1.0, 10.0, 4.0, 10.0, 4.0, 1.0, 2.0, 6.0, 8.0, 7.0, 9.0, 10.0, 2.0, 3.0, 9.0, 3.0, 3.0, 2.0, 5.0, 5.0, 3.0, 7.0, 8.0, 3.0, 6.0, 10.0, 10.0, 7.0, 8.0, 5.0, 2.0, 5.0, 9.0, 3.0, 5.0, 4.0, 2.0, 9.0, 7.0, 3.0, 3.0, 8.0, 6.0, 5.0, 2.0, 7.0, 4.0, 3.0, 4.0, 9.0, 6.0, 2.0, 7.0, 9.0, 3.0, 5.0, 7.0, 9.0, 10.0, 9.0, 8.0, 8.0, 4.0, 1.0, 3.0, 6.0, 6.0, 6.0, 3.0, 3.0, 9.0, 7.0, 3.0, 1.0, 10.0, 1.0, 9.0, 9.0, 7.0, 7.0, 4.0, 5.0, 3.0, 4.0, 4.0, 7.0, 3.0, 10.0, 7.0, 2.0, 8.0, 6.0, 9.0, 9.0, 1.0, 8.0, 9.0, 7.0, 1.0, 1.0, 5.0, 10.0, 10.0, 10.0, 8.0, 3.0, 2.0, 9.0, 10.0, 6.0, 4.0, 10.0, 6.0, 4.0, 2.0, 2.0, 9.0, 1.0, 7.0, 4.0, 2.0, 8.0, 7.0, 6.0, 8.0, 7.0, 6.0, 6.0, 2.0, 10.0, 10.0, 1.0, 6.0, 4.0, 9.0, 8.0, 5.0, 8.0, 3.0, 8.0, 1.0, 9.0, 8.0, 7.0, 5.0, 9.0, 1.0, 7.0, 8.0, 1.0, 3.0, 3.0, 7.0, 5.0]
global b_x = 5
global d_y = [6.0, 6.0, 2.0, 9.0, 7.0, 3.0, 2.0, 3.0, 5.0, 10.0, 1.0, 1.0, 8.0, 6.0, 6.0, 1.0, 9.0, 4.0, 5.0, 1.0, 9.0, 2.0, 7.0, 3.0, 2.0, 6.0, 5.0, 1.0, 9.0, 7.0, 10.0, 6.0, 10.0, 10.0, 5.0, 9.0, 3.0, 10.0, 5.0, 4.0, 7.0, 8.0, 2.0, 4.0, 5.0, 4.0, 5.0, 7.0, 8.0, 10.0, 7.0, 7.0, 6.0, 3.0, 6.0, 3.0, 2.0, 10.0, 6.0, 3.0, 6.0, 10.0, 4.0, 9.0, 10.0, 2.0, 6.0, 4.0, 1.0, 1.0, 1.0, 8.0, 4.0, 2.0, 5.0, 9.0, 10.0, 4.0, 1.0, 9.0, 2.0, 4.0, 7.0, 1.0, 5.0, 7.0, 4.0, 5.0, 8.0, 10.0, 6.0, 9.0, 10.0, 1.0, 8.0, 7.0, 5.0, 6.0, 9.0, 10.0, 3.0, 5.0, 1.0, 4.0, 4.0, 8.0, 2.0, 3.0, 3.0, 4.0, 9.0, 1.0, 7.0, 7.0, 7.0, 1.0, 10.0, 7.0, 4.0, 9.0, 9.0, 7.0, 1.0, 5.0, 10.0, 9.0, 1.0, 5.0, 4.0, 4.0, 1.0, 8.0, 3.0, 6.0, 3.0, 6.0, 1.0, 2.0, 3.0, 1.0, 1.0, 6.0, 9.0, 8.0, 3.0, 8.0, 4.0, 6.0, 8.0, 7.0, 6.0, 3.0, 1.0, 3.0, 3.0, 8.0, 3.0, 10.0, 10.0, 4.0, 10.0, 2.0, 1.0, 2.0, 8.0, 7.0, 7.0, 6.0, 9.0, 7.0, 10.0, 7.0, 2.0, 8.0, 2.0, 3.0, 1.0, 8.0, 3.0, 10.0, 1.0, 8.0, 7.0, 9.0, 7.0, 9.0, 8.0, 4.0, 3.0, 2.0, 9.0, 8.0, 8.0, 4.0, 5.0, 2.0, 8.0, 10.0, 9.0, 7.0, 7.0, 6.0, 1.0, 8.0, 1.0, 3.0, 6.0, 10.0, 2.0, 10.0, 7.0, 3.0, 5.0, 3.0, 10.0, 3.0, 10.0, 5.0, 9.0, 9.0, 1.0, 1.0, 3.0, 3.0, 3.0, 8.0, 10.0, 8.0, 10.0, 8.0, 4.0, 3.0, 8.0, 1.0, 7.0, 4.0, 10.0, 9.0, 4.0, 9.0]
global b_y = 10
global p = [0.1, 0.226, 0.119, 0.964, 0.249, 0.789, 0.573, 0.208, 0.41, 0.563, 0.879, 0.266, 0.808, 0.629, 0.869, 0.981, 0.223, 0.96, 0.527, 0.456, 0.166, 0.18, 0.139, 0.536, 0.434, 0.649, 0.36, 0.018, 0.817, 0.763, 0.559, 0.73, 0.935, 0.45, 0.728, 0.82, 0.062, 0.412, 0.42, 0.059, 0.131, 0.675, 0.808, 0.469, 0.237, 0.282, 0.133, 0.776, 0.298, 0.048, 0.252, 0.873, 0.389, 0.143, 0.979, 0.642, 0.176, 0.146, 0.273, 0.414, 0.92, 0.937, 0.579, 0.502, 0.392, 0.591, 0.974, 0.735, 0.766, 0.146, 0.942, 0.455, 0.527, 0.526, 0.14, 0.64, 0.873, 0.347, 0.21, 0.695, 0.402, 0.965, 0.241, 0.516, 0.166, 0.161, 0.851, 0.4, 0.001, 0.428, 0.924, 0.993, 0.785, 0.722, 0.04, 0.423, 0.699, 0.708, 0.618, 0.041, 0.295, 0.067, 0.188, 0.229, 0.039, 0.655, 0.672, 0.994, 0.126, 0.763, 0.174, 0.322, 0.146, 0.3, 0.71, 0.216, 0.076, 0.481, 0.066, 0.879, 0.303, 0.983, 0.173, 0.729, 0.169, 0.612, 0.303, 0.784, 0.803, 0.071, 0.936, 0.294, 0.821, 0.769, 0.188, 0.352, 0.596, 0.139, 0.911, 0.254, 0.1, 0.611, 0.21, 0.706, 0.295, 0.384, 0.521, 0.569, 0.702, 0.976, 0.886, 0.263, 0.778, 0.136, 0.078, 0.854, 0.426, 0.326, 0.926, 0.815, 0.333, 0.855, 0.451, 0.808, 0.735, 0.614, 0.208, 0.958, 0.56, 0.218, 0.555, 0.508, 0.57, 0.965, 0.664, 0.342, 0.392, 0.157, 0.107, 0.744, 0.151, 0.563, 0.565, 0.086, 0.515, 0.508, 0.555, 0.958, 0.286, 0.96, 0.168, 0.597, 0.313, 0.168, 0.451, 0.432, 0.037, 0.544, 0.498, 0.847, 0.601, 0.811, 0.278, 0.392, 0.736, 0.476, 0.944, 0.407, 0.727, 0.803, 0.893, 0.874, 0.083, 0.014, 0.467, 0.189, 0.382, 0.457, 0.368, 0.818, 0.281, 0.379, 0.235, 0.187, 0.141, 0.62, 0.316, 0.317, 0.515, 0.706, 0.157, 0.354, 0.573, 0.64, 0.932, 0.924, 0.877, 0.357, 0.04, 0.856]
global q = [0.987, 0.686, 0.142, 0.981, 0.549, 0.926, 0.968, 0.693, 0.741, 0.902, 0.924, 0.438, 0.935, 0.959, 0.951, 0.985, 0.84, 0.961, 0.621, 0.648, 0.438, 0.604, 0.755, 0.954, 0.504, 0.743, 0.788, 0.586, 0.979, 0.893, 0.573, 0.838, 0.985, 0.997, 0.876, 0.924, 0.894, 0.747, 0.817, 0.848, 0.473, 0.714, 0.925, 0.704, 0.526, 0.439, 0.689, 0.976, 0.743, 0.125, 0.58, 0.949, 0.483, 0.558, 0.998, 0.781, 0.788, 0.675, 0.701, 0.844, 0.99, 0.956, 0.739, 0.538, 0.823, 0.791, 0.992, 0.898, 0.804, 0.761, 0.97, 0.749, 0.662, 0.72, 0.718, 0.656, 0.909, 0.509, 0.252, 0.972, 0.94, 0.997, 0.628, 0.98, 0.286, 0.55, 0.852, 0.535, 0.586, 0.659, 0.971, 0.998, 0.831, 0.793, 0.176, 0.682, 0.795, 0.777, 0.953, 0.945, 0.531, 0.297, 0.53, 0.714, 0.555, 0.765, 0.954, 0.997, 0.363, 0.779, 0.226, 0.449, 0.783, 0.581, 0.767, 0.55, 0.641, 0.791, 0.291, 0.994, 0.82, 0.991, 0.254, 0.888, 0.898, 0.933, 0.969, 0.94, 0.882, 0.466, 0.966, 0.953, 0.839, 0.978, 0.75, 0.954, 0.839, 0.412, 0.925, 0.397, 0.77, 0.884, 0.814, 0.945, 0.887, 0.715, 0.915, 0.921, 0.814, 0.981, 0.992, 0.579, 0.931, 0.661, 0.49, 0.909, 0.645, 0.539, 0.94, 0.967, 0.873, 0.959, 0.956, 0.883, 0.951, 0.901, 0.398, 0.966, 0.703, 0.76, 0.771, 0.61, 0.929, 0.979, 0.827, 0.772, 0.656, 0.464, 0.309, 0.901, 0.491, 0.849, 0.622, 0.51, 0.879, 0.581, 0.556, 0.994, 0.875, 0.961, 0.472, 0.7, 0.514, 0.422, 0.827, 0.441, 0.832, 0.94, 0.925, 0.994, 0.765, 0.911, 0.847, 0.462, 0.834, 0.868, 0.984, 0.632, 0.727, 0.841, 0.979, 0.883, 0.965, 0.191, 0.523, 0.705, 0.466, 0.89, 0.51, 0.922, 0.626, 0.758, 0.629, 0.768, 0.717, 0.767, 0.526, 0.53, 0.563, 0.824, 0.878, 0.69, 0.925, 0.949, 0.946, 0.954, 0.898, 0.845, 0.773, 0.996]
global origin = 1
global destination = 50