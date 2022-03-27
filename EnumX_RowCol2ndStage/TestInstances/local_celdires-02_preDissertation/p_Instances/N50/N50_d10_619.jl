global arcs = [1 9; 1 16; 1 28; 2 4; 2 13; 2 14; 2 36; 2 39; 2 42; 3 33; 3 44; 4 17; 4 32; 4 46; 5 4; 5 22; 5 28; 5 34; 5 50; 6 2; 6 4; 6 10; 6 12; 6 17; 6 37; 6 45; 7 9; 7 11; 8 2; 8 11; 8 17; 8 28; 8 39; 8 45; 8 47; 8 48; 8 49; 9 2; 9 3; 9 8; 9 10; 9 13; 9 21; 9 27; 9 28; 9 38; 10 15; 10 31; 10 32; 10 36; 11 3; 11 31; 11 36; 11 44; 12 17; 12 19; 12 29; 12 30; 12 42; 13 10; 13 16; 13 17; 13 27; 13 37; 13 45; 13 47; 14 19; 14 29; 15 20; 15 49; 16 35; 16 38; 16 41; 16 45; 16 46; 16 48; 17 6; 17 26; 17 32; 18 2; 18 5; 18 24; 18 34; 18 35; 19 4; 19 12; 19 20; 19 26; 19 35; 19 39; 20 13; 20 26; 20 37; 20 39; 20 41; 21 8; 21 19; 21 40; 21 50; 22 2; 22 24; 22 35; 23 7; 23 15; 23 22; 23 32; 23 37; 23 46; 24 5; 24 14; 24 15; 24 27; 24 38; 25 9; 25 27; 25 36; 25 39; 25 40; 25 42; 25 44; 25 45; 25 46; 25 47; 26 25; 26 36; 26 37; 26 50; 27 3; 27 18; 27 20; 27 23; 27 29; 28 13; 28 17; 28 30; 28 35; 29 5; 29 39; 30 12; 30 19; 30 25; 30 44; 31 12; 31 23; 31 27; 31 42; 32 14; 32 19; 32 25; 32 33; 32 38; 32 40; 33 18; 33 22; 33 26; 33 27; 33 41; 33 46; 33 47; 34 38; 34 47; 35 2; 35 17; 35 30; 35 39; 35 42; 35 44; 36 3; 36 10; 36 13; 36 24; 36 25; 36 26; 36 29; 36 37; 36 49; 37 12; 37 33; 37 45; 38 18; 38 35; 38 44; 38 47; 38 49; 39 2; 39 3; 39 4; 39 5; 39 7; 39 8; 39 10; 39 15; 39 19; 39 20; 39 46; 40 5; 40 9; 40 13; 40 24; 40 31; 40 32; 40 38; 40 43; 41 21; 41 39; 42 9; 42 13; 42 25; 42 31; 42 33; 42 45; 43 4; 43 14; 43 27; 43 29; 44 9; 44 13; 45 11; 45 19; 45 22; 45 32; 45 39; 46 3; 46 20; 46 24; 46 40; 47 14; 47 37; 48 4; 48 18; 48 32; 48 42; 48 49; 49 14; 49 16; 49 28; 49 33; 49 38; 49 42]
global d_x = [8.0, 5.0, 1.0, 6.0, 7.0, 7.0, 8.0, 7.0, 2.0, 8.0, 8.0, 7.0, 10.0, 2.0, 5.0, 1.0, 6.0, 2.0, 7.0, 1.0, 10.0, 8.0, 10.0, 1.0, 2.0, 4.0, 4.0, 7.0, 2.0, 1.0, 5.0, 5.0, 7.0, 3.0, 3.0, 6.0, 3.0, 3.0, 8.0, 3.0, 9.0, 8.0, 1.0, 6.0, 9.0, 4.0, 6.0, 6.0, 2.0, 4.0, 4.0, 5.0, 9.0, 8.0, 1.0, 4.0, 5.0, 3.0, 1.0, 1.0, 1.0, 2.0, 4.0, 4.0, 6.0, 10.0, 4.0, 10.0, 9.0, 10.0, 3.0, 7.0, 9.0, 9.0, 9.0, 2.0, 1.0, 5.0, 5.0, 9.0, 5.0, 8.0, 3.0, 6.0, 4.0, 4.0, 10.0, 7.0, 5.0, 9.0, 3.0, 2.0, 7.0, 10.0, 4.0, 4.0, 8.0, 10.0, 8.0, 9.0, 10.0, 8.0, 10.0, 10.0, 8.0, 1.0, 7.0, 8.0, 1.0, 3.0, 8.0, 7.0, 7.0, 8.0, 6.0, 9.0, 1.0, 3.0, 6.0, 9.0, 2.0, 3.0, 6.0, 10.0, 9.0, 6.0, 1.0, 1.0, 3.0, 10.0, 3.0, 2.0, 2.0, 5.0, 7.0, 5.0, 5.0, 2.0, 10.0, 9.0, 8.0, 1.0, 2.0, 2.0, 10.0, 4.0, 8.0, 4.0, 10.0, 2.0, 1.0, 7.0, 7.0, 5.0, 4.0, 6.0, 1.0, 5.0, 4.0, 8.0, 10.0, 2.0, 8.0, 3.0, 1.0, 3.0, 5.0, 8.0, 9.0, 1.0, 7.0, 8.0, 10.0, 10.0, 4.0, 3.0, 3.0, 8.0, 4.0, 6.0, 10.0, 9.0, 5.0, 8.0, 9.0, 8.0, 9.0, 6.0, 5.0, 6.0, 3.0, 7.0, 4.0, 4.0, 9.0, 6.0, 9.0, 6.0, 8.0, 8.0, 3.0, 8.0, 7.0, 7.0, 8.0, 4.0, 5.0, 6.0, 3.0, 8.0, 3.0, 2.0, 3.0, 3.0, 3.0, 3.0, 7.0, 2.0, 4.0, 10.0, 4.0, 10.0, 2.0, 4.0, 6.0, 4.0, 1.0, 5.0, 8.0, 2.0, 3.0, 2.0, 2.0, 5.0, 9.0, 4.0, 6.0, 9.0, 5.0]
global b_x = 5
global d_y = [9.0, 3.0, 9.0, 8.0, 3.0, 3.0, 7.0, 1.0, 4.0, 2.0, 8.0, 6.0, 3.0, 9.0, 9.0, 6.0, 2.0, 4.0, 3.0, 1.0, 10.0, 8.0, 8.0, 7.0, 9.0, 7.0, 9.0, 1.0, 6.0, 9.0, 7.0, 3.0, 2.0, 2.0, 3.0, 8.0, 6.0, 2.0, 1.0, 10.0, 4.0, 2.0, 5.0, 7.0, 9.0, 2.0, 5.0, 6.0, 6.0, 8.0, 6.0, 5.0, 3.0, 8.0, 4.0, 3.0, 2.0, 4.0, 6.0, 6.0, 5.0, 6.0, 9.0, 2.0, 3.0, 8.0, 9.0, 7.0, 5.0, 1.0, 9.0, 7.0, 10.0, 3.0, 9.0, 9.0, 8.0, 1.0, 10.0, 8.0, 4.0, 6.0, 10.0, 1.0, 4.0, 8.0, 8.0, 1.0, 7.0, 3.0, 9.0, 10.0, 2.0, 10.0, 5.0, 9.0, 4.0, 1.0, 9.0, 1.0, 7.0, 3.0, 1.0, 3.0, 7.0, 7.0, 2.0, 7.0, 8.0, 7.0, 8.0, 8.0, 2.0, 9.0, 5.0, 1.0, 4.0, 5.0, 7.0, 1.0, 3.0, 7.0, 3.0, 10.0, 10.0, 5.0, 10.0, 5.0, 7.0, 6.0, 4.0, 10.0, 4.0, 6.0, 9.0, 9.0, 4.0, 10.0, 3.0, 5.0, 1.0, 3.0, 3.0, 5.0, 8.0, 6.0, 5.0, 10.0, 2.0, 2.0, 6.0, 10.0, 2.0, 5.0, 4.0, 6.0, 5.0, 9.0, 7.0, 2.0, 5.0, 2.0, 5.0, 6.0, 4.0, 6.0, 5.0, 4.0, 2.0, 9.0, 7.0, 10.0, 9.0, 7.0, 6.0, 3.0, 3.0, 8.0, 5.0, 10.0, 8.0, 9.0, 9.0, 8.0, 4.0, 4.0, 2.0, 7.0, 9.0, 5.0, 1.0, 5.0, 3.0, 8.0, 3.0, 8.0, 3.0, 1.0, 5.0, 5.0, 6.0, 10.0, 6.0, 3.0, 10.0, 4.0, 9.0, 5.0, 5.0, 6.0, 1.0, 9.0, 1.0, 1.0, 8.0, 2.0, 1.0, 4.0, 6.0, 6.0, 4.0, 3.0, 8.0, 5.0, 6.0, 4.0, 9.0, 7.0, 2.0, 6.0, 3.0, 10.0, 3.0, 1.0, 2.0, 9.0, 1.0, 8.0, 4.0]
global b_y = 10
global p = [0.389, 0.307, 0.869, 0.758, 0.518, 0.574, 0.469, 0.104, 0.772, 0.799, 0.331, 0.791, 0.252, 0.682, 0.02, 0.586, 0.122, 0.038, 0.57, 0.173, 0.098, 0.47, 0.134, 0.295, 0.356, 0.749, 0.999, 0.124, 0.464, 0.181, 0.667, 0.438, 0.887, 0.269, 0.285, 0.259, 0.312, 0.575, 0.725, 0.459, 0.121, 0.754, 0.468, 0.188, 0.857, 0.387, 0.545, 0.633, 0.339, 0.375, 0.311, 0.748, 0.104, 0.376, 0.958, 0.282, 0.016, 0.142, 0.335, 0.26, 0.231, 0.764, 0.685, 0.632, 0.054, 0.18, 0.211, 0.872, 0.968, 0.463, 0.065, 0.543, 0.659, 0.663, 0.156, 0.562, 0.28, 0.647, 0.417, 0.9, 0.917, 0.604, 0.934, 0.108, 0.188, 0.481, 0.617, 0.262, 0.215, 0.939, 0.119, 0.869, 0.538, 0.23, 0.715, 0.427, 0.186, 0.466, 0.323, 0.967, 0.275, 0.258, 0.359, 0.369, 0.969, 0.679, 0.586, 0.384, 0.439, 0.903, 0.463, 0.972, 0.427, 0.823, 0.093, 0.185, 0.221, 0.888, 0.56, 0.142, 0.947, 0.936, 0.81, 0.367, 0.029, 0.426, 0.465, 0.032, 0.824, 0.595, 0.004, 0.455, 0.121, 0.004, 0.23, 0.31, 0.654, 0.51, 0.163, 0.45, 0.273, 0.118, 0.2, 0.334, 0.064, 0.599, 0.024, 0.162, 0.287, 0.773, 0.792, 0.964, 0.459, 0.546, 0.691, 0.844, 0.636, 0.033, 0.15, 0.634, 0.087, 0.63, 0.186, 0.842, 0.837, 0.336, 0.365, 0.584, 0.018, 0.982, 0.857, 0.123, 0.127, 0.214, 0.895, 0.45, 0.903, 0.712, 0.551, 0.991, 0.036, 0.85, 0.904, 0.059, 0.109, 0.467, 0.671, 0.235, 0.698, 0.825, 0.213, 0.973, 0.822, 0.608, 0.35, 0.959, 0.884, 0.938, 0.46, 0.837, 0.963, 0.673, 0.539, 0.99, 0.091, 0.176, 0.724, 0.031, 0.906, 0.937, 0.633, 0.274, 0.095, 0.17, 0.881, 0.002, 0.858, 0.938, 0.134, 0.081, 0.777, 0.438, 0.714, 0.956, 0.348, 0.998, 0.661, 0.364, 0.595, 0.803, 0.132, 0.69, 0.306, 0.984, 0.307, 0.569, 0.862, 0.107, 0.442]
global q = [0.625, 0.47, 0.974, 0.769, 0.564, 0.739, 0.47, 0.458, 0.786, 0.821, 0.464, 0.966, 0.634, 0.818, 0.021, 0.971, 0.451, 0.89, 0.789, 0.78, 0.961, 0.888, 0.424, 0.637, 0.522, 0.998, 0.999, 0.247, 0.586, 0.841, 0.792, 0.697, 0.943, 0.525, 0.317, 0.287, 0.357, 0.831, 0.827, 0.575, 0.965, 0.99, 0.753, 0.815, 0.881, 0.735, 0.905, 0.928, 0.87, 0.426, 0.983, 0.83, 0.276, 0.465, 0.962, 0.441, 0.514, 0.572, 0.879, 0.536, 0.353, 0.931, 0.837, 0.882, 0.622, 0.358, 0.256, 0.968, 0.985, 0.837, 0.782, 0.972, 0.688, 0.81, 0.278, 0.669, 0.506, 0.731, 0.73, 0.923, 0.984, 0.864, 0.967, 0.136, 0.741, 0.884, 0.926, 0.851, 0.474, 0.983, 0.382, 0.986, 0.894, 0.427, 0.819, 0.869, 0.942, 0.496, 0.55, 0.971, 0.289, 0.696, 0.916, 0.925, 0.981, 0.804, 0.812, 0.987, 0.97, 0.991, 0.558, 0.997, 0.466, 0.938, 0.673, 0.25, 0.539, 0.978, 0.928, 0.983, 0.958, 0.981, 0.975, 0.533, 0.358, 0.718, 0.929, 0.139, 0.931, 0.703, 0.474, 0.995, 0.192, 0.322, 0.498, 0.552, 0.784, 0.51, 0.652, 0.48, 0.755, 0.186, 0.851, 0.369, 0.434, 0.707, 0.253, 0.513, 0.353, 0.988, 0.851, 0.98, 0.757, 0.656, 0.971, 0.871, 0.652, 0.185, 0.671, 0.978, 0.189, 0.8, 0.374, 0.914, 0.841, 0.64, 0.728, 0.758, 0.275, 0.987, 0.911, 0.395, 0.814, 0.704, 0.912, 0.506, 0.958, 0.839, 0.787, 0.994, 0.268, 0.852, 0.937, 0.713, 0.982, 0.483, 0.994, 0.537, 0.875, 0.887, 0.529, 0.978, 0.961, 0.947, 0.803, 0.974, 0.93, 0.997, 0.573, 0.862, 0.966, 0.886, 0.818, 0.995, 0.353, 0.425, 0.823, 0.818, 0.994, 0.969, 0.982, 0.275, 0.848, 0.758, 0.927, 0.052, 0.869, 0.946, 0.368, 0.808, 0.827, 0.879, 0.865, 0.959, 0.639, 0.999, 0.903, 0.967, 0.766, 0.815, 0.483, 0.928, 0.331, 0.999, 0.933, 0.802, 0.958, 0.705, 0.871]
global origin = 1
global destination = 50