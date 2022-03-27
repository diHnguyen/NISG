global arcs = [1 2; 1 8; 1 9; 1 50; 2 15; 2 22; 2 34; 2 36; 2 38; 2 41; 2 44; 2 47; 3 26; 3 29; 3 30; 3 47; 3 49; 4 15; 4 16; 4 27; 4 35; 4 49; 5 7; 5 23; 5 24; 5 44; 5 48; 6 13; 6 29; 6 38; 6 46; 7 3; 7 14; 7 17; 7 25; 7 37; 8 3; 8 6; 8 27; 8 38; 8 39; 8 42; 8 43; 8 44; 9 12; 9 23; 9 25; 9 28; 9 38; 9 48; 10 31; 10 43; 11 3; 11 25; 11 32; 11 48; 12 5; 12 6; 12 10; 12 13; 12 31; 12 35; 12 38; 12 39; 12 49; 13 6; 13 12; 13 20; 13 35; 14 4; 14 43; 15 4; 15 7; 15 16; 15 17; 15 18; 15 22; 15 30; 15 39; 15 48; 15 50; 16 2; 16 8; 16 10; 16 25; 16 35; 16 45; 17 3; 17 22; 17 23; 17 24; 17 30; 17 44; 18 22; 18 34; 18 49; 18 50; 19 3; 19 7; 19 15; 19 17; 19 36; 19 38; 20 11; 20 17; 20 21; 20 27; 20 29; 20 45; 21 3; 21 17; 21 19; 21 25; 21 40; 21 42; 22 13; 22 20; 22 31; 22 40; 22 42; 22 46; 23 5; 23 8; 23 27; 23 37; 23 41; 24 6; 24 14; 25 23; 25 28; 25 30; 25 45; 25 49; 26 12; 26 17; 26 19; 26 28; 27 3; 27 20; 27 21; 27 29; 27 32; 27 38; 27 39; 28 10; 28 42; 28 43; 28 44; 29 10; 29 19; 29 20; 29 21; 29 23; 29 26; 29 35; 29 41; 30 28; 30 29; 30 36; 30 43; 31 8; 31 12; 31 16; 31 47; 32 17; 32 21; 32 40; 32 44; 33 30; 33 37; 33 44; 34 4; 34 7; 34 21; 34 45; 35 6; 35 8; 35 18; 35 39; 36 2; 36 8; 36 10; 36 12; 36 33; 36 37; 36 46; 36 49; 37 35; 37 40; 37 43; 37 48; 37 49; 38 7; 38 32; 38 45; 39 13; 39 19; 39 24; 39 36; 40 3; 40 16; 40 18; 40 28; 40 29; 41 7; 41 10; 41 20; 41 24; 41 26; 41 45; 41 46; 42 6; 42 10; 42 11; 42 14; 42 23; 42 43; 43 9; 43 37; 43 38; 43 39; 43 44; 43 49; 44 36; 45 16; 45 23; 45 25; 45 43; 45 50; 46 2; 46 3; 46 12; 46 19; 46 20; 46 38; 46 47; 47 3; 47 25; 47 36; 47 44; 48 5; 48 11; 48 17; 48 28; 48 40; 49 8; 49 19; 49 34; 49 35; 49 38]
global d_x = [1.0, 2.0, 8.0, 7.0, 7.0, 6.0, 3.0, 5.0, 1.0, 3.0, 5.0, 5.0, 5.0, 6.0, 6.0, 1.0, 4.0, 6.0, 3.0, 10.0, 1.0, 10.0, 4.0, 8.0, 8.0, 8.0, 4.0, 5.0, 9.0, 10.0, 4.0, 6.0, 5.0, 1.0, 1.0, 3.0, 5.0, 10.0, 4.0, 6.0, 7.0, 9.0, 2.0, 1.0, 7.0, 5.0, 10.0, 1.0, 8.0, 5.0, 10.0, 2.0, 2.0, 8.0, 2.0, 3.0, 3.0, 2.0, 1.0, 9.0, 5.0, 10.0, 8.0, 8.0, 7.0, 7.0, 3.0, 2.0, 7.0, 9.0, 9.0, 4.0, 5.0, 10.0, 5.0, 7.0, 9.0, 7.0, 4.0, 3.0, 2.0, 4.0, 3.0, 7.0, 5.0, 1.0, 7.0, 1.0, 1.0, 3.0, 5.0, 1.0, 9.0, 1.0, 8.0, 7.0, 4.0, 10.0, 5.0, 10.0, 3.0, 4.0, 2.0, 8.0, 1.0, 8.0, 10.0, 8.0, 7.0, 7.0, 10.0, 3.0, 5.0, 6.0, 6.0, 7.0, 10.0, 8.0, 7.0, 9.0, 3.0, 10.0, 7.0, 7.0, 2.0, 6.0, 4.0, 6.0, 2.0, 7.0, 10.0, 4.0, 8.0, 3.0, 10.0, 3.0, 9.0, 3.0, 2.0, 1.0, 5.0, 10.0, 3.0, 5.0, 2.0, 2.0, 2.0, 4.0, 2.0, 1.0, 9.0, 1.0, 2.0, 3.0, 3.0, 9.0, 4.0, 10.0, 10.0, 4.0, 8.0, 8.0, 7.0, 1.0, 8.0, 1.0, 6.0, 9.0, 8.0, 5.0, 5.0, 9.0, 1.0, 3.0, 6.0, 8.0, 6.0, 5.0, 10.0, 8.0, 5.0, 4.0, 4.0, 4.0, 5.0, 9.0, 3.0, 10.0, 2.0, 1.0, 8.0, 7.0, 7.0, 1.0, 10.0, 5.0, 9.0, 10.0, 7.0, 9.0, 1.0, 3.0, 2.0, 10.0, 5.0, 8.0, 4.0, 2.0, 1.0, 8.0, 7.0, 10.0, 1.0, 10.0, 3.0, 10.0, 4.0, 8.0, 3.0, 9.0, 9.0, 6.0, 8.0, 8.0, 6.0, 9.0, 2.0, 5.0, 5.0, 3.0, 1.0, 8.0, 9.0, 10.0, 1.0, 2.0, 6.0, 9.0, 5.0, 2.0, 4.0, 6.0, 5.0, 4.0, 1.0, 10.0, 2.0, 8.0, 9.0, 2.0]
global b_x = 5
global d_y = [8.0, 2.0, 1.0, 10.0, 5.0, 7.0, 5.0, 3.0, 8.0, 5.0, 6.0, 7.0, 10.0, 7.0, 1.0, 8.0, 9.0, 4.0, 1.0, 2.0, 3.0, 1.0, 6.0, 9.0, 10.0, 6.0, 10.0, 9.0, 2.0, 5.0, 4.0, 8.0, 1.0, 2.0, 9.0, 4.0, 4.0, 1.0, 4.0, 5.0, 1.0, 9.0, 9.0, 8.0, 10.0, 9.0, 1.0, 7.0, 10.0, 1.0, 5.0, 2.0, 10.0, 10.0, 4.0, 8.0, 7.0, 8.0, 2.0, 3.0, 2.0, 3.0, 6.0, 5.0, 5.0, 3.0, 1.0, 8.0, 8.0, 10.0, 6.0, 8.0, 7.0, 6.0, 1.0, 8.0, 7.0, 6.0, 2.0, 7.0, 10.0, 9.0, 5.0, 7.0, 2.0, 3.0, 6.0, 8.0, 6.0, 2.0, 9.0, 7.0, 9.0, 9.0, 5.0, 8.0, 2.0, 7.0, 10.0, 4.0, 4.0, 7.0, 7.0, 6.0, 1.0, 10.0, 3.0, 9.0, 4.0, 1.0, 3.0, 8.0, 3.0, 9.0, 10.0, 1.0, 8.0, 8.0, 10.0, 2.0, 4.0, 4.0, 2.0, 10.0, 8.0, 10.0, 7.0, 7.0, 10.0, 2.0, 4.0, 7.0, 1.0, 10.0, 10.0, 1.0, 9.0, 3.0, 7.0, 9.0, 8.0, 5.0, 8.0, 10.0, 1.0, 8.0, 5.0, 8.0, 8.0, 2.0, 7.0, 4.0, 10.0, 7.0, 5.0, 6.0, 5.0, 2.0, 7.0, 6.0, 6.0, 4.0, 3.0, 5.0, 8.0, 8.0, 2.0, 1.0, 3.0, 4.0, 5.0, 3.0, 10.0, 5.0, 9.0, 10.0, 5.0, 3.0, 10.0, 10.0, 1.0, 5.0, 2.0, 10.0, 2.0, 8.0, 7.0, 8.0, 10.0, 4.0, 10.0, 6.0, 7.0, 8.0, 2.0, 6.0, 10.0, 6.0, 2.0, 1.0, 8.0, 10.0, 1.0, 8.0, 2.0, 6.0, 2.0, 6.0, 3.0, 8.0, 10.0, 1.0, 3.0, 5.0, 3.0, 3.0, 4.0, 10.0, 8.0, 3.0, 4.0, 9.0, 10.0, 10.0, 10.0, 4.0, 8.0, 1.0, 2.0, 5.0, 9.0, 10.0, 8.0, 1.0, 7.0, 1.0, 4.0, 7.0, 2.0, 4.0, 7.0, 9.0, 9.0, 7.0, 3.0, 10.0, 10.0, 6.0, 6.0, 7.0]
global b_y = 10
global p = [0.219, 0.719, 0.264, 0.066, 0.124, 0.903, 0.375, 0.442, 0.092, 0.457, 0.143, 0.323, 0.835, 0.529, 0.822, 0.5, 0.862, 0.714, 0.384, 0.857, 0.713, 0.496, 0.576, 0.254, 0.71, 0.379, 0.354, 0.709, 0.9, 0.303, 0.856, 0.223, 0.044, 0.523, 0.442, 0.013, 0.116, 0.375, 0.736, 0.696, 0.1, 0.05, 0.718, 0.661, 0.868, 0.275, 0.03, 0.467, 0.951, 0.068, 0.353, 0.212, 0.329, 0.802, 0.108, 0.029, 0.328, 0.18, 0.808, 0.853, 0.63, 0.041, 0.935, 0.229, 0.799, 0.915, 0.036, 0.851, 0.221, 0.367, 0.243, 0.001, 0.675, 0.119, 0.645, 0.879, 0.144, 0.966, 0.167, 0.768, 0.764, 0.893, 0.084, 0.764, 0.563, 0.265, 0.536, 0.572, 0.012, 0.65, 0.596, 0.758, 0.528, 0.661, 0.507, 0.976, 0.724, 0.107, 0.716, 0.1, 0.252, 0.629, 0.811, 0.563, 0.577, 0.38, 0.405, 0.376, 0.948, 0.923, 0.646, 0.626, 0.883, 0.837, 0.737, 0.235, 0.082, 0.246, 0.502, 0.327, 0.627, 0.945, 0.485, 0.262, 0.781, 0.888, 0.625, 0.777, 0.459, 0.779, 0.523, 0.108, 0.565, 0.109, 0.226, 0.626, 0.01, 0.389, 0.201, 0.46, 0.433, 0.965, 0.269, 0.288, 0.736, 0.551, 0.796, 0.87, 0.944, 0.752, 0.613, 0.431, 0.128, 0.065, 0.702, 0.517, 0.082, 0.377, 0.86, 0.243, 0.954, 0.464, 0.56, 0.776, 0.449, 0.995, 0.566, 0.669, 0.982, 0.593, 0.826, 0.897, 0.706, 0.269, 0.367, 0.026, 0.657, 0.345, 0.604, 0.793, 0.189, 0.132, 0.204, 0.515, 0.834, 0.983, 0.36, 0.793, 0.746, 0.491, 0.485, 0.316, 0.681, 0.343, 0.111, 0.48, 0.193, 0.458, 0.145, 0.421, 0.125, 0.898, 0.679, 0.753, 0.61, 0.751, 0.671, 0.318, 0.336, 0.495, 0.838, 0.317, 0.422, 0.288, 0.603, 0.537, 0.441, 0.41, 0.978, 0.649, 0.192, 0.825, 0.423, 0.829, 0.085, 0.058, 0.703, 0.147, 0.898, 0.256, 0.101, 0.442, 0.079, 0.996, 0.77, 0.485, 0.843, 0.907, 0.868, 0.855, 0.183, 0.558, 0.361, 0.168, 0.653, 0.386, 0.578, 0.874, 0.915, 0.339]
global q = [0.413, 0.785, 0.412, 0.996, 0.937, 0.935, 0.608, 0.505, 0.544, 0.788, 0.659, 0.605, 0.874, 0.598, 0.902, 0.948, 0.902, 0.989, 0.976, 0.981, 0.788, 0.804, 0.815, 0.749, 0.859, 0.713, 0.385, 0.984, 0.969, 0.726, 0.958, 0.662, 0.322, 0.644, 0.914, 0.939, 0.954, 0.721, 0.952, 0.882, 0.174, 0.209, 0.826, 0.798, 0.951, 0.517, 0.29, 0.681, 0.98, 0.271, 0.981, 0.625, 0.429, 0.932, 0.322, 0.619, 0.4, 0.683, 0.928, 0.963, 0.815, 0.547, 0.979, 0.66, 0.999, 0.937, 0.479, 0.88, 0.287, 0.687, 0.453, 0.102, 0.846, 0.578, 0.812, 0.926, 0.505, 0.978, 0.987, 0.776, 0.929, 0.997, 0.586, 0.953, 0.671, 0.949, 0.666, 0.87, 0.029, 0.974, 0.911, 0.884, 0.737, 0.662, 0.779, 0.987, 0.871, 0.834, 0.945, 0.722, 0.762, 0.808, 0.85, 0.969, 0.582, 0.882, 0.725, 0.76, 0.987, 0.951, 0.815, 0.838, 0.948, 0.967, 0.976, 0.268, 0.825, 0.335, 0.866, 0.748, 0.824, 0.978, 0.98, 0.489, 0.966, 0.9, 0.766, 0.965, 0.739, 0.907, 0.534, 0.729, 0.625, 0.861, 0.67, 0.744, 0.131, 0.43, 0.68, 0.932, 0.938, 0.988, 0.34, 0.634, 0.836, 0.684, 0.862, 0.922, 0.99, 0.824, 0.652, 0.972, 0.309, 0.313, 0.934, 0.694, 0.284, 0.979, 0.948, 0.636, 0.97, 0.927, 0.948, 0.998, 0.956, 0.998, 0.681, 0.955, 0.988, 0.945, 0.885, 0.905, 0.935, 0.911, 0.43, 0.318, 0.866, 0.67, 0.99, 0.845, 0.524, 0.421, 0.677, 0.674, 0.996, 0.986, 0.958, 0.808, 0.812, 0.936, 0.826, 0.443, 0.979, 0.84, 0.115, 0.746, 0.786, 0.962, 0.671, 0.836, 0.641, 0.97, 0.985, 0.808, 0.896, 0.792, 0.911, 0.737, 0.489, 0.657, 0.99, 0.441, 0.612, 0.904, 0.666, 0.929, 0.84, 0.696, 0.991, 0.773, 0.891, 0.856, 0.713, 0.989, 0.593, 0.528, 0.748, 0.269, 0.931, 0.517, 0.785, 0.544, 0.623, 0.998, 0.817, 0.677, 0.872, 0.917, 0.935, 0.875, 0.228, 0.855, 0.795, 0.865, 0.724, 0.736, 0.929, 0.89, 0.991, 0.442]
global origin = 1
global destination = 50