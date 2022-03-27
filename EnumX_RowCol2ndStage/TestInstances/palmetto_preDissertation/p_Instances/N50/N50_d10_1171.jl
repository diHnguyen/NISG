global arcs = [1 3; 1 5; 1 7; 1 19; 1 21; 1 23; 1 29; 1 30; 1 39; 1 41; 1 46; 2 12; 2 24; 2 25; 2 26; 2 27; 2 47; 3 9; 3 10; 3 26; 4 3; 5 7; 5 24; 5 25; 5 27; 5 38; 6 24; 7 3; 7 5; 7 13; 7 18; 7 40; 8 11; 8 12; 8 21; 8 34; 8 39; 9 8; 9 11; 9 13; 9 14; 9 15; 9 24; 9 26; 9 31; 9 36; 10 15; 10 21; 10 29; 10 30; 10 32; 10 38; 10 44; 10 49; 11 17; 11 27; 11 43; 11 46; 12 5; 12 7; 12 13; 12 35; 13 6; 13 24; 13 29; 14 7; 14 20; 14 25; 14 28; 14 35; 14 38; 14 40; 15 3; 15 12; 15 16; 15 23; 15 32; 15 41; 15 42; 15 47; 16 24; 16 35; 17 2; 17 9; 17 21; 17 22; 17 35; 17 43; 18 23; 18 24; 18 33; 18 35; 18 42; 18 44; 18 46; 18 49; 19 3; 19 11; 19 15; 19 28; 19 36; 19 41; 20 9; 20 17; 20 37; 20 42; 20 49; 20 50; 21 2; 21 3; 21 12; 21 33; 22 12; 22 13; 22 19; 22 23; 22 28; 22 34; 22 37; 22 43; 22 46; 22 47; 23 7; 23 13; 23 14; 23 20; 23 24; 23 43; 24 4; 24 20; 24 29; 24 33; 24 38; 24 48; 25 10; 25 23; 25 31; 25 40; 25 45; 26 6; 26 22; 26 24; 26 30; 26 38; 26 41; 26 48; 27 21; 27 33; 27 48; 28 13; 28 20; 28 21; 28 42; 28 49; 28 50; 29 25; 29 37; 29 38; 29 43; 30 6; 30 16; 30 46; 31 7; 31 27; 31 47; 32 3; 32 15; 32 37; 32 39; 32 40; 33 2; 33 9; 33 25; 33 28; 33 32; 33 41; 34 6; 34 21; 34 25; 34 27; 34 31; 34 38; 35 11; 35 17; 35 41; 35 49; 36 9; 36 27; 36 29; 36 33; 36 41; 37 11; 37 24; 37 29; 37 46; 38 5; 38 6; 38 15; 38 23; 39 6; 39 11; 39 15; 40 15; 40 20; 40 23; 40 29; 40 46; 41 6; 41 14; 41 36; 41 39; 42 4; 42 6; 42 11; 42 17; 42 36; 42 40; 43 16; 43 17; 43 24; 43 25; 43 29; 43 32; 44 2; 44 9; 44 31; 44 32; 44 34; 44 40; 44 47; 45 17; 45 34; 45 42; 46 23; 46 24; 46 44; 46 50; 47 5; 47 24; 47 43; 47 46; 47 49; 48 7; 48 20; 48 31; 48 32; 48 45; 48 46; 49 14; 49 16; 49 26; 49 30; 49 32; 49 50]
global d_x = [9.0, 10.0, 10.0, 8.0, 9.0, 3.0, 9.0, 3.0, 6.0, 3.0, 6.0, 5.0, 7.0, 7.0, 1.0, 1.0, 7.0, 1.0, 4.0, 8.0, 8.0, 4.0, 8.0, 1.0, 10.0, 4.0, 2.0, 10.0, 3.0, 7.0, 5.0, 4.0, 7.0, 9.0, 5.0, 5.0, 5.0, 4.0, 1.0, 7.0, 8.0, 4.0, 10.0, 6.0, 5.0, 10.0, 5.0, 3.0, 1.0, 4.0, 6.0, 3.0, 7.0, 10.0, 10.0, 1.0, 9.0, 2.0, 8.0, 5.0, 6.0, 9.0, 7.0, 3.0, 10.0, 1.0, 2.0, 4.0, 2.0, 9.0, 10.0, 8.0, 10.0, 9.0, 10.0, 9.0, 3.0, 8.0, 1.0, 4.0, 10.0, 7.0, 4.0, 10.0, 8.0, 6.0, 7.0, 7.0, 1.0, 4.0, 8.0, 10.0, 5.0, 7.0, 3.0, 4.0, 4.0, 10.0, 6.0, 4.0, 1.0, 2.0, 3.0, 9.0, 8.0, 4.0, 5.0, 10.0, 10.0, 6.0, 5.0, 10.0, 5.0, 9.0, 8.0, 8.0, 5.0, 4.0, 6.0, 9.0, 6.0, 8.0, 2.0, 8.0, 7.0, 3.0, 6.0, 1.0, 6.0, 5.0, 8.0, 9.0, 10.0, 8.0, 2.0, 9.0, 7.0, 8.0, 5.0, 6.0, 2.0, 5.0, 3.0, 6.0, 10.0, 10.0, 10.0, 5.0, 1.0, 9.0, 7.0, 6.0, 6.0, 7.0, 2.0, 10.0, 7.0, 6.0, 2.0, 7.0, 4.0, 5.0, 8.0, 8.0, 4.0, 8.0, 9.0, 10.0, 4.0, 3.0, 6.0, 6.0, 5.0, 4.0, 8.0, 6.0, 1.0, 6.0, 3.0, 4.0, 1.0, 10.0, 10.0, 10.0, 6.0, 8.0, 5.0, 5.0, 9.0, 5.0, 10.0, 7.0, 6.0, 4.0, 10.0, 1.0, 3.0, 1.0, 10.0, 1.0, 9.0, 7.0, 10.0, 2.0, 5.0, 1.0, 6.0, 9.0, 4.0, 6.0, 3.0, 5.0, 4.0, 10.0, 1.0, 1.0, 1.0, 9.0, 10.0, 3.0, 9.0, 3.0, 3.0, 6.0, 8.0, 10.0, 3.0, 3.0, 10.0, 9.0, 7.0, 3.0, 2.0, 3.0, 6.0, 5.0, 3.0, 3.0, 7.0, 4.0, 4.0, 2.0, 8.0, 4.0, 6.0, 8.0, 3.0, 10.0, 4.0, 2.0, 8.0, 2.0, 1.0, 1.0]
global b_x = 5
global d_y = [1.0, 6.0, 7.0, 7.0, 6.0, 3.0, 5.0, 4.0, 10.0, 6.0, 2.0, 5.0, 2.0, 8.0, 9.0, 9.0, 3.0, 2.0, 9.0, 9.0, 3.0, 10.0, 9.0, 2.0, 2.0, 2.0, 3.0, 7.0, 2.0, 5.0, 5.0, 5.0, 10.0, 8.0, 3.0, 4.0, 4.0, 4.0, 3.0, 2.0, 6.0, 8.0, 10.0, 3.0, 8.0, 10.0, 5.0, 3.0, 7.0, 10.0, 7.0, 5.0, 3.0, 5.0, 5.0, 6.0, 6.0, 10.0, 6.0, 4.0, 8.0, 5.0, 2.0, 7.0, 3.0, 2.0, 4.0, 10.0, 7.0, 7.0, 7.0, 7.0, 2.0, 6.0, 6.0, 8.0, 3.0, 7.0, 10.0, 10.0, 2.0, 10.0, 9.0, 7.0, 1.0, 2.0, 2.0, 5.0, 10.0, 5.0, 8.0, 5.0, 8.0, 6.0, 1.0, 6.0, 5.0, 3.0, 8.0, 7.0, 2.0, 4.0, 5.0, 5.0, 7.0, 8.0, 6.0, 10.0, 2.0, 5.0, 2.0, 1.0, 10.0, 7.0, 4.0, 9.0, 6.0, 4.0, 3.0, 3.0, 7.0, 7.0, 5.0, 1.0, 2.0, 9.0, 10.0, 7.0, 8.0, 3.0, 7.0, 5.0, 5.0, 3.0, 6.0, 5.0, 8.0, 2.0, 9.0, 9.0, 3.0, 3.0, 6.0, 8.0, 3.0, 9.0, 2.0, 3.0, 6.0, 4.0, 1.0, 10.0, 8.0, 6.0, 2.0, 3.0, 1.0, 2.0, 10.0, 10.0, 9.0, 3.0, 8.0, 3.0, 4.0, 7.0, 9.0, 1.0, 1.0, 2.0, 7.0, 7.0, 5.0, 5.0, 4.0, 3.0, 10.0, 8.0, 1.0, 1.0, 4.0, 9.0, 3.0, 2.0, 1.0, 1.0, 8.0, 5.0, 4.0, 6.0, 5.0, 8.0, 6.0, 6.0, 7.0, 7.0, 1.0, 8.0, 10.0, 6.0, 7.0, 4.0, 9.0, 4.0, 4.0, 10.0, 6.0, 3.0, 4.0, 9.0, 2.0, 9.0, 7.0, 9.0, 5.0, 8.0, 4.0, 10.0, 4.0, 5.0, 6.0, 10.0, 2.0, 1.0, 5.0, 8.0, 3.0, 7.0, 10.0, 9.0, 8.0, 5.0, 4.0, 8.0, 1.0, 5.0, 4.0, 9.0, 8.0, 9.0, 5.0, 6.0, 1.0, 4.0, 5.0, 7.0, 5.0, 1.0, 1.0, 4.0, 1.0, 4.0, 1.0, 9.0]
global b_y = 10
global p = [0.137, 0.44, 0.317, 0.827, 0.836, 0.064, 0.893, 0.561, 0.822, 0.961, 0.851, 0.044, 0.879, 0.156, 0.329, 0.064, 0.406, 0.658, 0.798, 0.456, 0.395, 0.235, 0.685, 0.519, 0.193, 0.281, 0.405, 0.691, 0.697, 0.251, 0.596, 0.085, 0.802, 0.843, 0.824, 0.269, 0.033, 0.732, 0.112, 0.245, 0.421, 0.433, 0.65, 0.724, 0.067, 0.166, 0.653, 0.942, 0.047, 0.624, 0.408, 0.55, 0.04, 0.89, 0.649, 0.307, 0.812, 0.365, 0.042, 0.493, 0.68, 0.098, 0.74, 0.761, 0.179, 0.852, 0.549, 0.903, 0.556, 0.508, 0.561, 0.201, 0.927, 0.563, 0.549, 0.552, 0.127, 0.944, 0.327, 0.118, 0.716, 0.202, 0.921, 0.287, 0.612, 0.778, 0.415, 0.141, 0.916, 0.664, 0.984, 0.362, 0.161, 0.108, 0.726, 0.524, 0.237, 0.552, 0.299, 0.79, 0.942, 0.992, 0.758, 0.536, 0.356, 0.289, 0.525, 0.374, 0.338, 0.694, 0.062, 0.269, 0.651, 0.126, 0.962, 0.473, 0.612, 0.524, 0.935, 0.756, 0.75, 0.143, 0.792, 0.922, 0.018, 0.238, 0.557, 0.016, 0.302, 0.289, 0.617, 0.523, 0.78, 0.854, 0.692, 0.187, 0.074, 0.86, 0.897, 0.985, 0.514, 0.241, 0.886, 0.573, 0.973, 0.795, 0.693, 0.96, 0.745, 0.145, 0.754, 0.086, 0.334, 0.872, 0.601, 0.745, 0.295, 0.171, 0.182, 0.544, 0.894, 0.807, 0.785, 0.426, 0.342, 0.35, 0.721, 0.346, 0.65, 0.089, 0.538, 0.767, 0.852, 0.274, 0.3, 0.658, 0.38, 0.555, 0.919, 0.648, 0.825, 0.495, 0.985, 0.676, 0.278, 0.207, 0.596, 0.87, 0.83, 0.165, 0.32, 0.004, 0.554, 0.448, 0.17, 0.632, 0.959, 0.69, 0.854, 0.082, 0.749, 0.872, 0.031, 0.364, 0.601, 0.122, 0.167, 0.549, 0.331, 0.179, 0.025, 0.928, 0.387, 0.05, 0.004, 0.247, 0.941, 0.166, 0.305, 0.337, 0.78, 0.22, 0.166, 0.723, 0.417, 0.062, 0.152, 0.744, 0.626, 0.053, 0.511, 0.715, 0.743, 0.504, 0.049, 0.194, 0.064, 0.601, 0.34, 0.441, 0.098, 0.577, 0.568, 0.588, 0.502, 0.434, 0.477, 0.478, 0.474, 0.528, 0.622, 0.348, 0.962, 0.858]
global q = [0.292, 0.951, 0.727, 0.973, 0.852, 0.964, 0.983, 0.57, 0.839, 0.965, 0.889, 0.609, 0.902, 0.461, 0.875, 0.683, 0.427, 0.714, 0.868, 0.821, 0.72, 0.506, 0.748, 0.549, 0.311, 0.47, 0.816, 0.888, 0.949, 0.827, 0.645, 0.622, 0.897, 0.879, 0.921, 0.488, 0.887, 0.908, 0.747, 0.413, 0.531, 0.882, 0.936, 0.88, 0.624, 0.776, 0.957, 0.954, 0.637, 0.738, 0.892, 0.823, 0.514, 0.939, 0.992, 0.575, 0.836, 0.776, 0.101, 0.6, 0.791, 0.516, 0.928, 0.955, 0.467, 0.918, 0.907, 0.978, 0.863, 0.895, 0.591, 0.619, 0.991, 0.719, 0.909, 0.686, 0.731, 0.979, 0.806, 0.918, 0.778, 0.999, 0.967, 0.997, 0.613, 0.826, 0.664, 0.583, 0.963, 0.849, 0.992, 0.591, 0.884, 0.196, 0.971, 0.832, 0.715, 0.561, 0.357, 0.966, 0.986, 0.997, 0.987, 0.792, 0.881, 0.6, 0.595, 0.905, 0.499, 0.789, 0.252, 0.892, 0.948, 0.642, 0.99, 0.923, 0.917, 0.865, 0.996, 0.918, 0.984, 0.544, 0.794, 0.929, 0.602, 0.894, 0.843, 0.772, 0.689, 0.826, 0.933, 0.915, 0.888, 0.985, 0.726, 0.313, 0.82, 0.89, 0.99, 0.997, 0.811, 0.638, 0.911, 0.613, 0.994, 0.899, 0.854, 0.974, 0.961, 0.644, 0.929, 0.343, 0.464, 0.903, 0.664, 0.762, 0.382, 0.23, 0.644, 0.808, 0.948, 0.909, 0.802, 0.535, 0.355, 0.429, 0.838, 0.903, 0.701, 0.319, 0.662, 0.813, 0.976, 0.502, 0.739, 0.7, 0.981, 0.755, 0.928, 0.962, 0.902, 0.642, 0.997, 0.912, 0.834, 0.822, 0.614, 0.972, 0.963, 0.955, 0.801, 0.026, 0.835, 0.983, 0.867, 0.866, 0.964, 0.802, 0.874, 0.861, 0.902, 0.992, 0.912, 0.916, 0.735, 0.713, 0.193, 0.676, 0.556, 0.398, 0.611, 0.948, 0.414, 0.833, 0.027, 0.594, 0.967, 0.906, 0.728, 0.595, 0.99, 0.712, 0.432, 0.837, 0.553, 0.957, 0.584, 0.937, 0.705, 0.119, 0.936, 0.906, 0.921, 0.687, 0.535, 0.22, 0.315, 0.872, 0.757, 0.606, 0.346, 0.637, 0.711, 0.82, 0.915, 0.499, 0.78, 0.668, 0.957, 0.701, 0.72, 0.798, 0.985, 0.992]
global origin = 1
global destination = 50