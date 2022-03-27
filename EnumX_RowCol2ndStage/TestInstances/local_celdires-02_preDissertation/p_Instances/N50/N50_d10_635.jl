global arcs = [1 17; 1 19; 1 25; 2 3; 2 7; 2 26; 2 40; 3 13; 3 16; 3 37; 3 41; 3 43; 3 48; 4 13; 4 16; 4 23; 4 24; 4 41; 4 45; 5 7; 5 29; 6 11; 6 16; 6 30; 6 37; 7 10; 7 30; 7 37; 7 38; 7 46; 8 5; 8 7; 8 9; 8 11; 8 20; 8 30; 8 32; 8 33; 8 39; 9 8; 9 14; 9 18; 9 19; 9 22; 9 25; 9 50; 10 19; 10 20; 10 23; 10 34; 10 35; 10 42; 10 47; 10 49; 11 2; 11 7; 11 33; 11 34; 11 37; 12 13; 12 17; 12 22; 12 23; 12 32; 12 36; 13 12; 13 18; 14 2; 14 5; 14 11; 14 16; 14 19; 14 30; 14 33; 14 34; 14 40; 15 5; 15 11; 15 21; 15 23; 16 6; 16 20; 16 22; 16 25; 16 27; 16 34; 16 35; 16 37; 17 11; 18 6; 18 38; 18 45; 18 48; 18 50; 19 8; 19 11; 19 24; 19 28; 19 46; 19 47; 20 7; 20 8; 20 18; 20 47; 20 48; 21 27; 21 29; 22 13; 22 20; 22 23; 22 25; 23 3; 23 14; 23 19; 23 24; 23 38; 23 42; 24 3; 24 4; 24 13; 24 44; 24 46; 25 2; 25 11; 25 20; 25 34; 25 37; 25 47; 26 17; 26 18; 26 20; 26 25; 26 27; 26 42; 26 44; 26 47; 26 50; 27 8; 27 10; 27 12; 27 22; 28 7; 28 9; 29 11; 29 28; 29 40; 29 45; 30 12; 30 22; 30 26; 30 29; 31 3; 31 9; 31 23; 31 28; 31 40; 31 42; 31 43; 32 4; 32 7; 32 8; 32 11; 32 14; 32 25; 32 26; 32 37; 32 38; 32 39; 32 44; 33 9; 33 14; 33 24; 33 26; 33 27; 33 48; 33 49; 34 31; 34 36; 34 37; 35 6; 35 36; 35 49; 35 50; 36 3; 36 14; 36 33; 37 4; 38 8; 38 20; 38 27; 38 39; 39 9; 39 20; 39 28; 39 35; 39 43; 39 46; 39 49; 39 50; 40 23; 40 31; 40 32; 40 44; 41 20; 41 21; 42 16; 42 18; 42 39; 43 12; 43 14; 43 23; 43 32; 43 35; 44 3; 44 4; 44 8; 44 18; 44 38; 44 42; 44 49; 45 7; 45 8; 45 43; 45 44; 46 8; 46 11; 46 30; 46 34; 46 40; 46 41; 46 49; 47 4; 47 8; 47 12; 47 15; 47 19; 47 20; 47 30; 47 37; 47 42; 47 44; 48 2; 48 8; 48 17; 48 18; 48 29; 48 44; 48 46; 49 2; 49 8; 49 13; 49 18; 49 25; 49 28; 49 31; 49 34; 49 42]
global d_x = [4.0, 5.0, 9.0, 4.0, 6.0, 2.0, 7.0, 7.0, 2.0, 8.0, 6.0, 9.0, 3.0, 2.0, 9.0, 2.0, 10.0, 8.0, 7.0, 6.0, 5.0, 10.0, 3.0, 8.0, 3.0, 9.0, 4.0, 9.0, 3.0, 10.0, 1.0, 10.0, 2.0, 10.0, 4.0, 1.0, 4.0, 2.0, 8.0, 8.0, 4.0, 2.0, 1.0, 1.0, 9.0, 10.0, 2.0, 1.0, 6.0, 3.0, 6.0, 3.0, 8.0, 3.0, 7.0, 1.0, 2.0, 4.0, 1.0, 6.0, 1.0, 3.0, 5.0, 2.0, 5.0, 7.0, 7.0, 8.0, 7.0, 4.0, 3.0, 9.0, 10.0, 9.0, 7.0, 6.0, 1.0, 1.0, 3.0, 4.0, 4.0, 6.0, 7.0, 3.0, 3.0, 2.0, 6.0, 8.0, 3.0, 8.0, 7.0, 6.0, 4.0, 4.0, 8.0, 3.0, 3.0, 6.0, 7.0, 8.0, 1.0, 4.0, 4.0, 1.0, 3.0, 5.0, 1.0, 10.0, 7.0, 3.0, 5.0, 3.0, 10.0, 9.0, 2.0, 5.0, 7.0, 4.0, 4.0, 7.0, 1.0, 5.0, 5.0, 8.0, 7.0, 5.0, 2.0, 6.0, 2.0, 3.0, 1.0, 3.0, 5.0, 10.0, 3.0, 7.0, 2.0, 1.0, 10.0, 4.0, 6.0, 2.0, 9.0, 9.0, 6.0, 10.0, 1.0, 8.0, 3.0, 3.0, 6.0, 10.0, 6.0, 9.0, 3.0, 5.0, 7.0, 6.0, 7.0, 6.0, 7.0, 5.0, 6.0, 6.0, 2.0, 1.0, 8.0, 7.0, 7.0, 6.0, 7.0, 1.0, 6.0, 3.0, 8.0, 6.0, 10.0, 5.0, 10.0, 3.0, 9.0, 4.0, 7.0, 5.0, 5.0, 10.0, 6.0, 6.0, 5.0, 7.0, 2.0, 8.0, 10.0, 7.0, 4.0, 10.0, 6.0, 6.0, 10.0, 1.0, 10.0, 10.0, 1.0, 2.0, 2.0, 2.0, 5.0, 3.0, 5.0, 2.0, 10.0, 10.0, 8.0, 1.0, 4.0, 4.0, 1.0, 3.0, 5.0, 9.0, 8.0, 10.0, 8.0, 9.0, 6.0, 9.0, 10.0, 8.0, 9.0, 3.0, 2.0, 6.0, 5.0, 6.0, 6.0, 4.0, 10.0, 5.0, 2.0, 10.0, 9.0, 4.0, 9.0, 1.0, 8.0, 2.0, 3.0, 9.0, 10.0, 1.0, 10.0, 4.0, 2.0, 4.0, 2.0, 3.0, 3.0]
global b_x = 5
global d_y = [6.0, 8.0, 8.0, 8.0, 2.0, 2.0, 2.0, 9.0, 8.0, 8.0, 6.0, 4.0, 6.0, 6.0, 9.0, 4.0, 6.0, 9.0, 9.0, 7.0, 10.0, 10.0, 8.0, 2.0, 5.0, 2.0, 9.0, 9.0, 10.0, 2.0, 5.0, 4.0, 4.0, 6.0, 1.0, 2.0, 4.0, 1.0, 10.0, 2.0, 3.0, 2.0, 1.0, 6.0, 9.0, 2.0, 6.0, 6.0, 6.0, 4.0, 10.0, 10.0, 10.0, 3.0, 1.0, 2.0, 5.0, 4.0, 5.0, 7.0, 5.0, 1.0, 2.0, 6.0, 1.0, 7.0, 3.0, 9.0, 2.0, 10.0, 9.0, 1.0, 3.0, 5.0, 1.0, 3.0, 6.0, 10.0, 9.0, 10.0, 2.0, 1.0, 3.0, 10.0, 2.0, 8.0, 5.0, 7.0, 7.0, 4.0, 10.0, 5.0, 8.0, 7.0, 8.0, 8.0, 1.0, 8.0, 5.0, 3.0, 1.0, 9.0, 2.0, 9.0, 5.0, 9.0, 10.0, 7.0, 4.0, 8.0, 1.0, 10.0, 7.0, 6.0, 9.0, 9.0, 10.0, 8.0, 6.0, 7.0, 2.0, 9.0, 3.0, 8.0, 5.0, 2.0, 5.0, 2.0, 9.0, 5.0, 10.0, 9.0, 7.0, 5.0, 8.0, 4.0, 7.0, 9.0, 6.0, 7.0, 6.0, 2.0, 8.0, 6.0, 10.0, 2.0, 9.0, 4.0, 2.0, 8.0, 1.0, 4.0, 1.0, 8.0, 6.0, 10.0, 2.0, 1.0, 10.0, 5.0, 3.0, 6.0, 10.0, 4.0, 1.0, 10.0, 5.0, 3.0, 9.0, 5.0, 10.0, 5.0, 7.0, 3.0, 4.0, 3.0, 9.0, 2.0, 4.0, 1.0, 3.0, 5.0, 9.0, 5.0, 10.0, 3.0, 6.0, 7.0, 10.0, 2.0, 5.0, 9.0, 1.0, 8.0, 3.0, 3.0, 3.0, 9.0, 3.0, 1.0, 10.0, 1.0, 5.0, 7.0, 3.0, 8.0, 6.0, 9.0, 8.0, 3.0, 3.0, 10.0, 6.0, 1.0, 6.0, 5.0, 8.0, 9.0, 1.0, 6.0, 6.0, 4.0, 8.0, 6.0, 6.0, 1.0, 2.0, 8.0, 1.0, 7.0, 6.0, 7.0, 2.0, 2.0, 3.0, 1.0, 7.0, 3.0, 5.0, 7.0, 3.0, 1.0, 4.0, 10.0, 7.0, 10.0, 9.0, 4.0, 6.0, 7.0, 6.0, 6.0, 4.0, 6.0, 1.0, 3.0, 1.0]
global b_y = 10
global p = [0.175, 0.02, 0.659, 0.561, 0.276, 0.319, 0.308, 0.265, 0.008, 0.433, 0.095, 0.674, 0.847, 0.239, 0.501, 0.654, 0.117, 0.81, 0.85, 0.426, 0.749, 0.906, 0.862, 0.646, 0.361, 0.524, 0.502, 0.943, 0.789, 0.746, 0.948, 0.859, 0.4, 0.593, 0.474, 0.396, 0.538, 0.131, 0.814, 0.576, 0.194, 0.546, 0.472, 0.079, 0.767, 0.379, 0.691, 0.598, 0.645, 0.566, 0.043, 0.656, 0.099, 0.33, 0.377, 0.663, 0.054, 0.586, 0.36, 0.893, 0.031, 0.368, 0.432, 0.572, 0.794, 0.56, 0.001, 0.464, 0.549, 0.964, 0.222, 0.99, 0.484, 0.393, 0.62, 0.88, 0.755, 0.378, 0.782, 0.877, 0.068, 0.082, 0.502, 0.295, 0.109, 0.389, 0.602, 0.027, 0.018, 0.977, 0.01, 0.194, 0.054, 0.799, 0.841, 0.724, 0.849, 0.504, 0.294, 0.714, 0.768, 0.222, 0.99, 0.034, 0.817, 0.076, 0.725, 0.096, 0.079, 0.03, 0.98, 0.381, 0.576, 0.146, 0.901, 0.442, 0.713, 0.921, 0.427, 0.338, 0.902, 0.024, 0.121, 0.538, 0.833, 0.76, 0.554, 0.215, 0.655, 0.694, 0.168, 0.144, 0.942, 0.702, 0.299, 0.28, 0.124, 0.934, 0.567, 0.109, 0.876, 0.455, 0.38, 0.486, 0.701, 0.031, 0.97, 0.123, 0.045, 0.593, 0.164, 0.987, 0.14, 0.567, 0.894, 0.259, 0.319, 0.926, 0.97, 0.667, 0.287, 0.429, 0.728, 0.554, 0.507, 0.168, 0.109, 0.347, 0.257, 0.05, 0.939, 0.822, 0.792, 0.642, 0.316, 0.14, 0.717, 0.188, 0.596, 0.786, 0.425, 0.354, 0.433, 0.627, 0.332, 0.852, 0.73, 0.047, 0.831, 0.495, 0.912, 0.628, 0.365, 0.258, 0.522, 0.956, 0.64, 0.85, 0.916, 0.247, 0.862, 0.458, 0.681, 0.132, 0.763, 0.069, 0.303, 0.711, 0.729, 0.775, 0.913, 0.623, 0.301, 0.964, 0.01, 0.866, 0.427, 0.486, 0.344, 0.456, 0.642, 0.529, 0.7, 0.638, 0.846, 0.882, 0.466, 0.168, 0.181, 0.478, 0.477, 0.577, 0.171, 0.055, 0.347, 0.934, 0.431, 0.46, 0.277, 0.751, 0.48, 0.237, 0.413, 0.569, 0.556, 0.604, 0.264, 0.514, 0.986, 0.014, 0.947, 0.882, 0.938, 0.375, 0.019, 0.933, 0.18]
global q = [0.982, 0.827, 0.807, 0.893, 0.821, 0.949, 0.344, 0.798, 0.774, 0.533, 0.901, 0.72, 0.903, 0.557, 0.76, 0.801, 0.822, 0.84, 0.944, 0.702, 0.788, 0.927, 0.992, 0.649, 0.844, 0.852, 0.573, 0.995, 0.953, 0.9, 0.949, 0.867, 0.776, 0.965, 0.558, 0.826, 0.553, 0.56, 0.974, 0.957, 0.303, 0.564, 0.814, 0.395, 0.999, 0.912, 0.91, 0.829, 0.756, 0.733, 0.058, 0.798, 0.398, 0.421, 0.614, 0.787, 0.484, 0.852, 0.451, 0.928, 0.397, 0.699, 0.93, 0.94, 0.799, 0.808, 0.779, 0.877, 0.706, 0.983, 0.909, 0.993, 0.689, 0.67, 0.893, 0.88, 0.781, 0.886, 0.897, 0.933, 0.341, 0.678, 0.92, 0.353, 0.69, 0.59, 0.871, 0.769, 0.852, 0.99, 0.683, 0.382, 0.638, 0.978, 0.971, 0.82, 0.872, 0.696, 0.626, 0.983, 0.872, 0.924, 0.99, 0.269, 0.817, 0.929, 0.928, 0.371, 0.261, 0.435, 0.987, 0.843, 0.685, 0.527, 0.99, 0.743, 0.72, 0.962, 0.847, 0.471, 0.956, 0.855, 0.916, 0.641, 0.86, 0.886, 0.63, 0.801, 0.771, 0.985, 0.897, 0.642, 0.956, 0.983, 0.413, 0.641, 0.339, 0.937, 0.884, 0.188, 0.936, 0.745, 0.71, 0.666, 0.909, 0.742, 0.992, 0.793, 0.832, 0.987, 0.706, 0.995, 0.315, 0.956, 0.909, 0.529, 0.711, 0.966, 0.986, 0.866, 0.789, 0.499, 0.809, 0.821, 0.88, 0.676, 0.259, 0.721, 0.315, 0.978, 0.951, 0.835, 0.995, 0.914, 0.79, 0.405, 0.756, 0.209, 0.672, 0.901, 0.446, 0.792, 0.804, 0.67, 0.513, 0.957, 0.96, 0.222, 0.91, 0.543, 0.971, 0.999, 0.877, 0.648, 0.888, 0.984, 0.842, 0.94, 0.961, 0.727, 0.87, 0.66, 0.974, 0.522, 0.85, 0.729, 0.677, 0.799, 0.775, 0.807, 0.965, 0.641, 0.434, 0.969, 0.566, 0.904, 0.859, 0.916, 0.583, 0.755, 0.671, 0.921, 0.96, 0.856, 0.891, 0.956, 0.793, 0.618, 0.792, 0.563, 0.923, 0.653, 0.512, 0.691, 0.744, 0.987, 0.738, 0.66, 0.807, 0.992, 0.991, 0.469, 0.557, 0.716, 0.727, 0.73, 0.496, 0.85, 0.994, 0.278, 0.987, 0.886, 0.995, 0.613, 0.551, 0.963, 0.619]
global origin = 1
global destination = 50