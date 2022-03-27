global arcs = [1 2; 1 22; 1 35; 1 40; 1 41; 1 46; 1 50; 2 6; 2 29; 2 33; 2 42; 2 43; 3 13; 3 17; 3 38; 3 43; 4 7; 4 16; 4 38; 4 39; 4 42; 4 47; 5 8; 5 10; 5 12; 5 24; 5 40; 5 43; 6 24; 6 49; 7 2; 7 8; 7 18; 7 21; 7 30; 7 41; 7 44; 8 17; 8 44; 9 14; 9 16; 9 18; 9 29; 9 38; 9 47; 9 50; 10 8; 10 50; 11 2; 11 3; 11 18; 11 29; 11 35; 11 40; 11 46; 12 7; 12 15; 12 18; 12 21; 12 38; 12 42; 12 46; 12 50; 13 12; 13 33; 14 21; 14 26; 14 27; 14 29; 14 34; 14 37; 14 38; 15 4; 15 41; 15 43; 16 5; 16 30; 16 39; 16 43; 16 47; 17 10; 17 12; 17 18; 18 5; 18 8; 18 15; 18 24; 18 31; 18 34; 18 39; 18 47; 18 48; 19 10; 19 11; 19 43; 19 50; 20 7; 20 11; 20 17; 20 25; 20 44; 21 16; 21 19; 21 23; 22 19; 22 37; 23 5; 23 17; 23 18; 23 26; 24 4; 24 11; 24 15; 24 34; 24 35; 24 39; 24 49; 25 11; 25 14; 25 21; 25 23; 26 5; 26 11; 26 15; 26 23; 26 24; 26 42; 26 43; 27 9; 27 12; 27 17; 28 9; 28 11; 28 21; 28 30; 28 41; 28 44; 29 5; 29 7; 29 10; 29 26; 29 27; 29 30; 29 36; 29 46; 29 50; 30 4; 30 7; 30 11; 30 26; 30 31; 30 35; 30 40; 30 47; 31 9; 31 19; 31 22; 31 28; 31 37; 32 6; 32 9; 32 35; 32 43; 33 8; 33 10; 33 27; 34 2; 34 12; 34 17; 35 7; 35 11; 35 13; 35 16; 35 18; 35 22; 35 27; 35 28; 35 41; 36 27; 36 38; 37 5; 37 13; 37 21; 37 35; 37 49; 38 6; 38 12; 38 16; 38 18; 38 23; 38 27; 38 30; 38 42; 38 50; 39 29; 39 38; 40 7; 40 8; 40 21; 40 22; 40 34; 40 39; 41 15; 41 18; 41 20; 41 37; 41 43; 42 11; 42 13; 42 23; 42 36; 43 12; 43 34; 43 37; 43 45; 44 9; 44 26; 44 28; 44 33; 44 41; 45 4; 45 34; 46 7; 46 10; 46 15; 46 33; 46 38; 46 48; 47 3; 47 30; 47 50; 48 16; 48 17; 48 32; 49 29]
global d_x = [4.0, 8.0, 7.0, 4.0, 5.0, 10.0, 10.0, 10.0, 6.0, 3.0, 10.0, 1.0, 6.0, 5.0, 4.0, 6.0, 6.0, 5.0, 9.0, 9.0, 9.0, 4.0, 2.0, 5.0, 1.0, 10.0, 5.0, 10.0, 4.0, 6.0, 8.0, 9.0, 10.0, 9.0, 6.0, 10.0, 1.0, 1.0, 4.0, 7.0, 3.0, 10.0, 2.0, 1.0, 10.0, 9.0, 9.0, 9.0, 6.0, 4.0, 3.0, 4.0, 4.0, 2.0, 9.0, 10.0, 9.0, 6.0, 4.0, 9.0, 7.0, 8.0, 6.0, 2.0, 3.0, 7.0, 4.0, 1.0, 1.0, 7.0, 2.0, 9.0, 4.0, 9.0, 2.0, 4.0, 8.0, 5.0, 2.0, 7.0, 3.0, 7.0, 7.0, 2.0, 1.0, 7.0, 3.0, 7.0, 5.0, 6.0, 7.0, 5.0, 2.0, 6.0, 2.0, 6.0, 10.0, 1.0, 4.0, 10.0, 8.0, 5.0, 7.0, 5.0, 9.0, 3.0, 5.0, 4.0, 5.0, 1.0, 1.0, 9.0, 8.0, 8.0, 2.0, 9.0, 9.0, 1.0, 10.0, 10.0, 6.0, 7.0, 1.0, 2.0, 3.0, 6.0, 9.0, 9.0, 6.0, 3.0, 4.0, 3.0, 5.0, 3.0, 6.0, 6.0, 2.0, 1.0, 10.0, 10.0, 1.0, 1.0, 4.0, 2.0, 10.0, 5.0, 1.0, 10.0, 2.0, 3.0, 5.0, 6.0, 1.0, 8.0, 3.0, 10.0, 1.0, 2.0, 6.0, 1.0, 6.0, 3.0, 8.0, 2.0, 6.0, 1.0, 5.0, 5.0, 8.0, 5.0, 2.0, 2.0, 4.0, 5.0, 1.0, 4.0, 8.0, 3.0, 6.0, 3.0, 6.0, 6.0, 7.0, 4.0, 4.0, 9.0, 4.0, 8.0, 5.0, 8.0, 5.0, 6.0, 8.0, 1.0, 9.0, 7.0, 2.0, 10.0, 10.0, 7.0, 9.0, 7.0, 2.0, 7.0, 7.0, 6.0, 2.0, 7.0, 10.0, 10.0, 1.0, 4.0, 7.0, 8.0, 5.0, 10.0, 2.0, 2.0, 9.0, 4.0, 6.0, 8.0, 7.0, 4.0, 3.0, 2.0, 10.0, 5.0, 7.0, 9.0, 5.0, 5.0, 2.0, 10.0, 1.0]
global b_x = 5
global d_y = [2.0, 10.0, 6.0, 10.0, 2.0, 3.0, 3.0, 8.0, 5.0, 10.0, 10.0, 10.0, 8.0, 5.0, 6.0, 7.0, 10.0, 9.0, 7.0, 1.0, 6.0, 1.0, 2.0, 5.0, 8.0, 5.0, 6.0, 8.0, 3.0, 2.0, 9.0, 6.0, 2.0, 10.0, 2.0, 2.0, 10.0, 3.0, 8.0, 1.0, 5.0, 2.0, 7.0, 3.0, 9.0, 9.0, 9.0, 8.0, 6.0, 5.0, 4.0, 1.0, 9.0, 9.0, 4.0, 3.0, 8.0, 4.0, 8.0, 6.0, 5.0, 7.0, 1.0, 2.0, 3.0, 9.0, 8.0, 2.0, 2.0, 6.0, 2.0, 10.0, 10.0, 1.0, 5.0, 8.0, 4.0, 4.0, 4.0, 8.0, 4.0, 10.0, 7.0, 4.0, 6.0, 6.0, 9.0, 4.0, 1.0, 10.0, 3.0, 8.0, 2.0, 3.0, 8.0, 6.0, 9.0, 9.0, 2.0, 3.0, 8.0, 10.0, 1.0, 2.0, 6.0, 5.0, 3.0, 4.0, 7.0, 2.0, 5.0, 6.0, 9.0, 4.0, 6.0, 7.0, 2.0, 3.0, 8.0, 2.0, 1.0, 1.0, 6.0, 3.0, 8.0, 9.0, 6.0, 9.0, 9.0, 2.0, 10.0, 10.0, 4.0, 2.0, 4.0, 2.0, 7.0, 10.0, 4.0, 8.0, 7.0, 3.0, 10.0, 6.0, 1.0, 7.0, 6.0, 1.0, 5.0, 4.0, 7.0, 6.0, 5.0, 3.0, 9.0, 7.0, 6.0, 10.0, 9.0, 4.0, 4.0, 2.0, 3.0, 10.0, 9.0, 6.0, 1.0, 3.0, 6.0, 3.0, 10.0, 1.0, 8.0, 3.0, 5.0, 6.0, 4.0, 1.0, 4.0, 3.0, 9.0, 6.0, 3.0, 9.0, 1.0, 5.0, 2.0, 6.0, 6.0, 2.0, 6.0, 6.0, 3.0, 8.0, 9.0, 5.0, 5.0, 1.0, 8.0, 5.0, 1.0, 8.0, 10.0, 1.0, 4.0, 4.0, 1.0, 2.0, 5.0, 4.0, 6.0, 9.0, 5.0, 5.0, 4.0, 5.0, 6.0, 10.0, 7.0, 5.0, 2.0, 7.0, 6.0, 7.0, 3.0, 9.0, 5.0, 7.0, 2.0, 4.0, 10.0, 3.0, 10.0, 8.0, 8.0]
global b_y = 10
global p = [0.832, 0.763, 0.601, 0.416, 0.787, 0.813, 0.315, 0.97, 0.72, 0.668, 0.757, 0.196, 0.165, 0.691, 0.111, 0.387, 0.258, 0.105, 0.586, 0.322, 0.885, 0.572, 0.233, 0.292, 0.116, 0.037, 0.713, 0.199, 0.99, 0.07, 0.371, 0.291, 0.075, 0.406, 0.06, 0.106, 0.398, 0.138, 0.599, 0.804, 0.316, 0.453, 0.728, 0.732, 0.154, 0.573, 0.279, 0.764, 0.171, 0.728, 0.439, 0.863, 0.497, 0.756, 0.401, 0.981, 0.603, 0.221, 0.142, 0.731, 0.946, 0.064, 0.351, 0.924, 0.917, 0.4, 0.906, 0.84, 0.249, 0.674, 0.646, 0.179, 0.595, 0.35, 0.24, 0.783, 0.844, 0.284, 0.721, 0.414, 0.733, 0.68, 0.26, 0.663, 0.705, 0.065, 0.612, 0.49, 0.914, 0.774, 0.637, 0.584, 0.415, 0.073, 0.903, 0.395, 0.162, 0.131, 0.594, 0.623, 0.85, 0.579, 0.808, 0.842, 0.177, 0.317, 0.169, 0.637, 0.291, 0.346, 0.237, 0.228, 0.109, 0.867, 0.103, 0.531, 0.563, 0.013, 0.04, 0.257, 0.163, 0.625, 0.437, 0.622, 0.126, 0.141, 0.152, 0.113, 0.974, 0.392, 0.239, 0.71, 0.403, 0.949, 0.35, 0.671, 0.185, 0.538, 0.221, 0.547, 0.766, 0.404, 0.815, 0.515, 0.605, 0.739, 0.755, 0.811, 0.509, 0.343, 0.083, 0.002, 0.453, 0.983, 0.69, 0.599, 0.365, 0.98, 0.417, 0.386, 0.188, 0.421, 0.97, 0.815, 0.782, 0.817, 0.264, 0.669, 0.851, 0.615, 0.564, 0.012, 0.641, 0.026, 0.777, 0.77, 0.444, 0.974, 0.591, 0.279, 0.293, 0.715, 0.259, 0.332, 0.407, 0.178, 0.928, 0.911, 0.925, 0.624, 0.831, 0.902, 0.419, 0.05, 0.124, 0.708, 0.705, 0.8, 0.345, 0.71, 0.806, 0.58, 0.923, 0.093, 0.423, 0.379, 0.096, 0.606, 0.708, 0.767, 0.367, 0.592, 0.881, 0.861, 0.325, 0.875, 0.495, 0.256, 0.447, 0.477, 0.593, 0.418, 0.773, 0.846, 0.227, 0.398, 0.36, 0.633, 0.092, 0.834, 0.223, 0.715, 0.23, 0.314, 0.243]
global q = [0.874, 0.913, 0.792, 0.815, 0.873, 0.83, 0.339, 0.976, 0.768, 0.719, 0.902, 0.738, 0.342, 0.806, 0.907, 0.772, 0.89, 0.12, 0.633, 0.948, 0.933, 0.994, 0.584, 0.85, 0.934, 0.289, 0.897, 0.209, 0.998, 0.741, 0.895, 0.556, 0.355, 0.9, 0.707, 0.654, 0.91, 0.63, 0.87, 0.823, 0.534, 0.814, 0.974, 0.828, 0.358, 0.701, 0.728, 0.95, 0.812, 0.953, 0.794, 0.965, 0.56, 0.939, 0.937, 0.991, 0.776, 0.894, 0.447, 0.978, 0.949, 0.999, 0.357, 0.965, 0.928, 0.608, 0.993, 0.971, 0.401, 0.769, 0.954, 0.187, 0.705, 0.439, 0.345, 0.959, 0.958, 0.598, 0.971, 0.518, 0.824, 0.79, 0.824, 0.737, 0.824, 0.595, 0.899, 0.506, 0.952, 0.984, 0.811, 0.774, 0.629, 0.503, 0.927, 0.902, 0.443, 0.176, 0.686, 0.664, 0.943, 0.669, 0.972, 0.861, 0.425, 0.568, 0.962, 0.771, 0.412, 0.884, 0.253, 0.721, 0.477, 0.89, 0.924, 0.961, 0.627, 0.933, 0.896, 0.664, 0.864, 0.879, 0.588, 0.765, 0.733, 0.176, 0.623, 0.211, 0.982, 0.55, 0.763, 0.977, 0.618, 0.978, 0.374, 0.682, 0.989, 0.602, 0.524, 0.903, 0.89, 0.746, 0.943, 0.582, 0.624, 0.862, 0.811, 0.996, 0.917, 0.552, 0.407, 0.536, 0.628, 0.997, 0.712, 0.605, 0.96, 0.999, 0.519, 0.838, 0.655, 0.565, 0.985, 0.883, 0.893, 0.954, 0.517, 0.802, 0.892, 0.68, 0.635, 0.724, 0.806, 0.213, 0.821, 0.855, 0.919, 0.974, 0.888, 0.721, 0.785, 0.841, 0.731, 0.44, 0.987, 0.409, 0.982, 0.916, 0.981, 0.654, 0.923, 0.996, 0.758, 0.407, 0.217, 0.864, 0.885, 0.84, 0.795, 0.905, 0.922, 0.758, 0.939, 0.829, 0.915, 0.877, 0.792, 0.821, 0.969, 0.994, 0.552, 0.618, 0.948, 0.895, 0.716, 0.99, 0.625, 0.447, 0.993, 0.792, 0.769, 0.689, 0.829, 0.988, 0.959, 0.648, 0.484, 0.879, 0.938, 0.988, 0.878, 0.82, 0.595, 0.647, 0.535]
global origin = 1
global destination = 50