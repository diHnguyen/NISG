global arcs = [1 5; 1 6; 1 18; 1 24; 2 22; 2 24; 2 37; 3 2; 3 4; 3 7; 3 8; 3 22; 3 31; 3 35; 3 37; 4 3; 4 7; 4 14; 4 33; 4 35; 4 36; 5 6; 5 13; 5 14; 6 13; 6 24; 6 26; 7 9; 7 11; 7 15; 7 23; 7 32; 8 5; 8 11; 9 11; 9 12; 9 18; 9 21; 9 28; 9 38; 10 12; 10 17; 10 28; 10 37; 11 19; 11 22; 12 6; 13 2; 13 17; 13 20; 13 37; 13 40; 14 6; 14 23; 14 28; 14 30; 14 31; 15 8; 15 30; 15 32; 16 2; 16 8; 16 39; 16 40; 17 24; 17 25; 18 5; 18 9; 18 16; 18 26; 18 28; 19 11; 19 21; 19 22; 19 30; 19 33; 20 8; 20 12; 20 18; 20 31; 20 37; 21 15; 21 37; 22 3; 22 7; 22 8; 22 14; 22 20; 22 33; 23 3; 23 4; 23 8; 23 10; 23 13; 23 19; 23 27; 23 35; 23 37; 24 5; 24 6; 24 12; 24 16; 24 17; 24 21; 24 22; 25 34; 26 4; 26 15; 26 38; 27 13; 27 22; 27 29; 28 3; 28 4; 28 13; 28 14; 28 27; 28 37; 29 16; 29 22; 29 25; 29 32; 29 36; 29 39; 30 17; 30 20; 31 12; 31 15; 31 17; 31 37; 32 3; 32 4; 32 5; 32 6; 32 7; 32 9; 32 17; 32 33; 33 9; 33 20; 33 32; 33 36; 33 38; 34 18; 34 27; 35 21; 35 26; 35 30; 35 33; 35 37; 36 7; 36 22; 36 30; 36 40; 37 22; 37 35; 38 4; 38 7; 39 19; 39 24]
global d_x = [3.0, 8.0, 3.0, 7.0, 6.0, 8.0, 10.0, 5.0, 10.0, 1.0, 9.0, 4.0, 3.0, 5.0, 8.0, 4.0, 2.0, 9.0, 10.0, 1.0, 1.0, 4.0, 4.0, 10.0, 6.0, 5.0, 3.0, 5.0, 6.0, 8.0, 1.0, 4.0, 4.0, 5.0, 6.0, 4.0, 9.0, 8.0, 5.0, 9.0, 8.0, 6.0, 8.0, 7.0, 1.0, 8.0, 2.0, 6.0, 8.0, 5.0, 4.0, 8.0, 7.0, 5.0, 7.0, 6.0, 4.0, 6.0, 6.0, 10.0, 8.0, 7.0, 7.0, 2.0, 7.0, 3.0, 1.0, 1.0, 3.0, 5.0, 4.0, 4.0, 10.0, 4.0, 5.0, 1.0, 4.0, 5.0, 5.0, 3.0, 1.0, 5.0, 3.0, 3.0, 7.0, 4.0, 10.0, 3.0, 8.0, 9.0, 7.0, 9.0, 9.0, 10.0, 2.0, 8.0, 8.0, 10.0, 6.0, 2.0, 8.0, 8.0, 8.0, 7.0, 2.0, 5.0, 1.0, 4.0, 5.0, 1.0, 5.0, 5.0, 8.0, 5.0, 9.0, 8.0, 3.0, 7.0, 9.0, 6.0, 4.0, 9.0, 10.0, 4.0, 10.0, 9.0, 9.0, 9.0, 9.0, 7.0, 4.0, 7.0, 5.0, 4.0, 8.0, 6.0, 6.0, 7.0, 2.0, 4.0, 8.0, 3.0, 3.0, 6.0, 3.0, 10.0, 8.0, 9.0, 4.0, 1.0, 3.0, 8.0, 10.0, 4.0, 10.0, 9.0, 2.0, 10.0, 1.0, 8.0]
global b_x = 5
global d_y = [4.0, 1.0, 8.0, 2.0, 3.0, 3.0, 5.0, 7.0, 6.0, 7.0, 1.0, 3.0, 6.0, 5.0, 5.0, 9.0, 2.0, 4.0, 2.0, 10.0, 1.0, 6.0, 2.0, 4.0, 8.0, 1.0, 2.0, 10.0, 9.0, 4.0, 3.0, 5.0, 3.0, 9.0, 9.0, 7.0, 5.0, 6.0, 10.0, 7.0, 6.0, 3.0, 3.0, 8.0, 6.0, 1.0, 5.0, 7.0, 10.0, 1.0, 6.0, 2.0, 6.0, 4.0, 10.0, 10.0, 9.0, 4.0, 6.0, 7.0, 7.0, 1.0, 4.0, 10.0, 3.0, 6.0, 10.0, 9.0, 1.0, 9.0, 8.0, 9.0, 5.0, 10.0, 5.0, 9.0, 1.0, 10.0, 8.0, 2.0, 3.0, 4.0, 3.0, 8.0, 9.0, 6.0, 1.0, 1.0, 6.0, 8.0, 9.0, 5.0, 7.0, 3.0, 7.0, 6.0, 8.0, 9.0, 9.0, 7.0, 8.0, 9.0, 3.0, 5.0, 7.0, 1.0, 8.0, 4.0, 3.0, 2.0, 2.0, 9.0, 3.0, 2.0, 2.0, 8.0, 4.0, 1.0, 6.0, 3.0, 2.0, 9.0, 2.0, 2.0, 7.0, 2.0, 7.0, 1.0, 9.0, 10.0, 2.0, 2.0, 4.0, 3.0, 6.0, 3.0, 8.0, 9.0, 10.0, 8.0, 10.0, 2.0, 1.0, 6.0, 10.0, 6.0, 10.0, 3.0, 5.0, 7.0, 7.0, 10.0, 1.0, 2.0, 10.0, 3.0, 4.0, 1.0, 4.0, 3.0]
global b_y = 10
global p = [0.253, 0.995, 0.003, 0.679, 0.051, 0.565, 0.608, 0.473, 0.107, 0.675, 0.109, 0.116, 0.443, 0.708, 0.666, 0.575, 0.487, 0.034, 0.328, 0.665, 0.439, 0.049, 0.701, 0.59, 0.347, 0.376, 0.56, 0.595, 0.597, 0.691, 0.027, 0.386, 0.01, 0.888, 0.284, 0.469, 0.794, 0.211, 0.484, 0.4, 0.893, 0.737, 0.249, 0.705, 0.999, 0.142, 0.206, 0.224, 0.983, 0.854, 0.758, 0.118, 0.903, 0.746, 0.879, 0.987, 0.496, 0.458, 0.308, 0.1, 0.062, 0.952, 0.885, 0.346, 0.205, 0.782, 0.547, 0.939, 0.618, 0.006, 0.043, 0.945, 0.543, 0.067, 0.496, 0.414, 0.84, 0.575, 0.061, 0.382, 0.618, 0.26, 0.144, 0.257, 0.381, 0.841, 0.369, 0.451, 0.576, 0.046, 0.312, 0.122, 0.083, 0.758, 0.942, 0.931, 0.326, 0.477, 0.262, 0.793, 0.379, 0.194, 0.456, 0.857, 0.833, 0.983, 0.812, 0.903, 0.098, 0.756, 0.272, 0.931, 0.101, 0.848, 0.519, 0.829, 0.41, 0.818, 0.616, 0.521, 0.117, 0.497, 0.931, 0.736, 0.022, 0.043, 0.479, 0.422, 0.65, 0.575, 0.656, 0.349, 0.792, 0.042, 0.522, 0.228, 0.872, 0.925, 0.193, 0.67, 0.2, 0.065, 0.062, 0.276, 0.413, 0.012, 0.869, 0.594, 0.098, 0.42, 0.78, 0.692, 0.7, 0.637, 0.478, 0.212, 0.955, 0.466, 0.371, 0.677]
global q = [0.313, 0.999, 0.591, 0.888, 0.203, 0.683, 0.621, 0.718, 0.565, 0.817, 0.353, 0.974, 0.719, 0.76, 0.823, 0.92, 0.714, 0.696, 0.333, 0.793, 0.563, 0.355, 0.994, 0.732, 0.689, 0.542, 0.992, 0.735, 0.669, 0.99, 0.598, 0.426, 0.211, 0.894, 0.391, 0.616, 0.911, 0.43, 0.789, 0.919, 0.958, 0.992, 0.716, 0.809, 0.999, 0.433, 0.492, 0.995, 0.986, 0.872, 0.768, 0.451, 0.93, 0.972, 0.998, 0.987, 0.613, 0.634, 0.803, 0.509, 0.484, 0.954, 0.968, 0.926, 0.893, 0.897, 0.95, 0.965, 0.931, 0.062, 0.176, 0.974, 0.998, 0.911, 0.733, 0.457, 0.953, 0.693, 0.226, 0.761, 0.872, 0.312, 0.768, 0.348, 0.53, 0.986, 0.663, 0.474, 0.601, 0.518, 0.373, 0.499, 0.618, 0.998, 0.952, 0.966, 0.807, 0.58, 0.596, 0.982, 0.983, 0.466, 0.559, 0.867, 0.945, 0.992, 0.926, 0.924, 0.745, 0.765, 0.882, 0.953, 0.702, 0.903, 0.801, 0.876, 0.866, 0.993, 0.773, 0.598, 0.507, 0.975, 0.932, 0.945, 0.95, 0.834, 0.868, 0.978, 0.711, 0.845, 0.814, 0.378, 0.903, 0.68, 0.869, 0.426, 0.954, 0.953, 0.725, 0.959, 0.357, 0.355, 0.293, 0.655, 0.734, 0.931, 0.915, 0.945, 0.685, 0.998, 0.973, 0.894, 0.778, 0.91, 0.761, 0.746, 0.984, 0.601, 0.439, 0.812]
global origin = 1
global destination = 40