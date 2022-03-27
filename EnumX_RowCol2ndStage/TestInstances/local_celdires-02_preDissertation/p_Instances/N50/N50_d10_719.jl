global arcs = [1 6; 1 34; 1 39; 1 44; 2 8; 2 37; 2 38; 3 9; 3 18; 3 21; 3 24; 3 28; 3 29; 3 32; 3 34; 3 39; 4 28; 4 39; 5 8; 5 30; 5 32; 5 46; 5 48; 6 2; 6 3; 6 13; 6 40; 7 2; 7 19; 7 25; 7 26; 7 28; 8 2; 8 3; 8 23; 8 33; 8 39; 8 46; 9 3; 9 7; 9 16; 9 21; 9 39; 9 48; 10 7; 10 19; 10 21; 10 28; 10 47; 11 16; 11 24; 11 35; 11 40; 11 50; 12 9; 12 17; 12 28; 12 30; 12 31; 12 36; 12 38; 12 39; 12 40; 13 6; 13 27; 13 30; 13 38; 14 13; 14 17; 14 18; 14 23; 14 37; 14 39; 14 49; 15 12; 15 32; 15 37; 15 43; 16 22; 16 35; 16 37; 16 39; 16 44; 16 45; 16 46; 16 47; 17 13; 17 19; 17 22; 17 24; 17 31; 17 44; 18 3; 18 8; 18 40; 19 2; 19 14; 19 20; 20 3; 20 11; 20 15; 20 21; 20 33; 20 34; 20 41; 21 7; 21 11; 21 22; 21 25; 21 27; 21 29; 22 2; 22 19; 22 26; 22 29; 22 30; 22 33; 22 37; 22 48; 23 3; 23 4; 23 25; 23 40; 23 45; 24 6; 24 14; 24 19; 24 23; 24 31; 24 41; 24 44; 24 48; 25 7; 25 12; 25 15; 25 24; 25 26; 25 28; 25 32; 25 43; 25 45; 25 46; 25 49; 25 50; 26 19; 26 24; 26 43; 26 50; 27 10; 27 15; 27 20; 27 31; 27 48; 28 4; 28 15; 28 19; 28 31; 28 39; 28 47; 29 7; 29 12; 29 14; 29 17; 29 18; 29 23; 29 30; 29 31; 30 8; 30 19; 30 34; 30 38; 31 7; 31 10; 31 18; 31 29; 31 44; 31 47; 32 3; 32 5; 32 20; 32 29; 33 20; 33 26; 33 34; 33 43; 33 46; 34 9; 34 14; 34 28; 34 29; 35 7; 35 11; 35 12; 35 37; 35 38; 36 9; 36 14; 36 25; 36 48; 37 18; 37 20; 37 21; 37 36; 37 38; 38 8; 38 10; 38 13; 38 20; 38 28; 38 33; 39 2; 39 34; 39 36; 39 45; 40 3; 40 16; 40 26; 40 38; 40 41; 41 34; 42 10; 42 11; 42 20; 42 37; 42 38; 42 39; 43 11; 43 14; 43 42; 43 44; 44 3; 44 18; 44 20; 44 25; 44 26; 44 31; 44 45; 45 3; 45 19; 45 29; 45 46; 46 11; 46 27; 46 32; 46 45; 46 47; 47 2; 47 22; 47 23; 47 28; 48 4; 48 9; 48 17; 48 33; 48 34; 48 36; 48 38; 48 41; 49 7; 49 18; 49 22; 49 25; 49 44; 49 50]
global d_x = [7.0, 9.0, 1.0, 9.0, 9.0, 6.0, 7.0, 4.0, 6.0, 3.0, 1.0, 8.0, 2.0, 5.0, 8.0, 1.0, 6.0, 10.0, 5.0, 3.0, 6.0, 7.0, 3.0, 4.0, 10.0, 5.0, 3.0, 8.0, 9.0, 10.0, 1.0, 3.0, 4.0, 9.0, 4.0, 4.0, 1.0, 10.0, 3.0, 10.0, 1.0, 7.0, 7.0, 10.0, 4.0, 7.0, 10.0, 7.0, 1.0, 10.0, 6.0, 8.0, 10.0, 6.0, 6.0, 4.0, 6.0, 4.0, 9.0, 1.0, 7.0, 9.0, 1.0, 1.0, 10.0, 3.0, 5.0, 3.0, 2.0, 3.0, 2.0, 9.0, 7.0, 10.0, 4.0, 1.0, 2.0, 9.0, 1.0, 8.0, 2.0, 2.0, 8.0, 10.0, 10.0, 4.0, 5.0, 2.0, 7.0, 5.0, 1.0, 7.0, 6.0, 2.0, 4.0, 3.0, 7.0, 9.0, 7.0, 3.0, 6.0, 2.0, 5.0, 1.0, 6.0, 9.0, 7.0, 6.0, 7.0, 5.0, 8.0, 4.0, 4.0, 10.0, 1.0, 9.0, 7.0, 6.0, 10.0, 2.0, 7.0, 2.0, 9.0, 4.0, 9.0, 9.0, 6.0, 2.0, 2.0, 8.0, 2.0, 5.0, 7.0, 9.0, 8.0, 10.0, 4.0, 6.0, 8.0, 10.0, 4.0, 7.0, 2.0, 10.0, 5.0, 8.0, 10.0, 7.0, 3.0, 2.0, 5.0, 4.0, 9.0, 2.0, 3.0, 10.0, 7.0, 3.0, 9.0, 8.0, 2.0, 2.0, 7.0, 7.0, 6.0, 5.0, 1.0, 7.0, 8.0, 6.0, 4.0, 9.0, 7.0, 2.0, 3.0, 10.0, 7.0, 3.0, 7.0, 9.0, 1.0, 10.0, 5.0, 9.0, 10.0, 7.0, 3.0, 6.0, 1.0, 5.0, 3.0, 2.0, 7.0, 2.0, 5.0, 7.0, 6.0, 7.0, 5.0, 4.0, 2.0, 6.0, 8.0, 7.0, 7.0, 10.0, 6.0, 6.0, 8.0, 3.0, 9.0, 2.0, 4.0, 4.0, 1.0, 4.0, 8.0, 1.0, 6.0, 3.0, 9.0, 9.0, 6.0, 2.0, 7.0, 7.0, 3.0, 4.0, 5.0, 7.0, 8.0, 3.0, 2.0, 4.0, 1.0, 5.0, 1.0, 5.0, 1.0, 4.0, 4.0, 1.0, 1.0, 6.0, 2.0, 2.0, 8.0, 8.0, 7.0, 5.0, 7.0, 5.0, 2.0, 1.0, 10.0, 2.0, 9.0, 6.0, 6.0, 2.0, 6.0, 9.0, 4.0, 4.0]
global b_x = 5
global d_y = [6.0, 10.0, 6.0, 10.0, 8.0, 6.0, 8.0, 6.0, 7.0, 5.0, 3.0, 9.0, 5.0, 10.0, 7.0, 10.0, 9.0, 5.0, 9.0, 2.0, 5.0, 9.0, 4.0, 5.0, 2.0, 8.0, 4.0, 3.0, 10.0, 4.0, 7.0, 7.0, 10.0, 7.0, 8.0, 1.0, 6.0, 2.0, 3.0, 4.0, 10.0, 1.0, 5.0, 9.0, 5.0, 3.0, 5.0, 7.0, 3.0, 9.0, 6.0, 6.0, 3.0, 4.0, 7.0, 4.0, 5.0, 1.0, 2.0, 4.0, 8.0, 3.0, 2.0, 4.0, 7.0, 4.0, 2.0, 5.0, 9.0, 10.0, 5.0, 7.0, 2.0, 6.0, 5.0, 2.0, 1.0, 2.0, 1.0, 6.0, 10.0, 8.0, 7.0, 1.0, 8.0, 1.0, 9.0, 10.0, 10.0, 10.0, 8.0, 2.0, 6.0, 3.0, 2.0, 4.0, 8.0, 3.0, 6.0, 6.0, 6.0, 4.0, 7.0, 3.0, 9.0, 10.0, 10.0, 10.0, 3.0, 8.0, 7.0, 1.0, 4.0, 5.0, 1.0, 9.0, 8.0, 10.0, 1.0, 2.0, 3.0, 10.0, 10.0, 3.0, 5.0, 1.0, 10.0, 4.0, 9.0, 6.0, 8.0, 4.0, 10.0, 1.0, 1.0, 7.0, 4.0, 1.0, 3.0, 3.0, 8.0, 5.0, 7.0, 5.0, 4.0, 10.0, 9.0, 3.0, 5.0, 8.0, 3.0, 5.0, 4.0, 9.0, 9.0, 3.0, 10.0, 1.0, 2.0, 10.0, 10.0, 7.0, 10.0, 6.0, 8.0, 7.0, 1.0, 1.0, 5.0, 3.0, 9.0, 10.0, 9.0, 1.0, 2.0, 7.0, 7.0, 7.0, 6.0, 7.0, 3.0, 7.0, 1.0, 7.0, 5.0, 7.0, 3.0, 7.0, 7.0, 7.0, 8.0, 7.0, 4.0, 5.0, 6.0, 9.0, 8.0, 4.0, 4.0, 8.0, 3.0, 10.0, 4.0, 2.0, 3.0, 2.0, 5.0, 1.0, 3.0, 5.0, 3.0, 1.0, 4.0, 4.0, 7.0, 9.0, 4.0, 3.0, 3.0, 5.0, 4.0, 8.0, 8.0, 2.0, 4.0, 8.0, 5.0, 2.0, 6.0, 10.0, 6.0, 7.0, 8.0, 7.0, 3.0, 8.0, 8.0, 7.0, 7.0, 8.0, 5.0, 1.0, 5.0, 10.0, 1.0, 10.0, 9.0, 1.0, 9.0, 8.0, 1.0, 3.0, 5.0, 8.0, 4.0, 10.0, 7.0, 10.0, 4.0, 8.0, 9.0, 9.0, 5.0, 7.0]
global b_y = 10
global p = [0.065, 0.236, 0.745, 0.446, 0.325, 0.304, 0.23, 0.186, 0.452, 0.165, 0.859, 0.767, 0.516, 0.852, 0.196, 0.956, 0.775, 0.399, 0.32, 0.159, 0.542, 0.091, 0.554, 0.971, 0.992, 0.244, 0.082, 0.216, 0.922, 0.771, 0.409, 0.939, 0.243, 0.067, 0.223, 0.237, 0.393, 0.817, 0.085, 0.439, 0.565, 0.705, 0.191, 0.262, 0.072, 0.405, 0.606, 0.241, 0.322, 0.758, 0.155, 0.558, 0.388, 0.417, 0.543, 0.658, 0.236, 0.1, 0.127, 0.334, 0.202, 0.743, 0.845, 0.513, 0.111, 0.477, 0.693, 0.514, 0.733, 0.434, 0.866, 0.9, 0.962, 0.227, 0.274, 0.936, 0.741, 0.424, 0.753, 0.61, 0.968, 0.018, 0.825, 0.302, 0.572, 0.965, 0.342, 0.598, 0.015, 0.431, 0.997, 0.932, 0.864, 0.109, 0.151, 0.593, 0.821, 0.916, 0.964, 0.382, 0.912, 0.946, 0.949, 0.402, 0.624, 0.051, 0.314, 0.638, 0.675, 0.077, 0.696, 0.75, 0.018, 0.991, 0.603, 0.111, 0.871, 0.956, 0.904, 0.765, 0.167, 0.697, 0.151, 0.323, 0.067, 0.553, 0.839, 0.088, 0.464, 0.547, 0.436, 0.19, 0.416, 0.866, 0.622, 0.475, 0.886, 0.233, 0.338, 0.739, 0.195, 0.946, 0.03, 0.403, 0.037, 0.253, 0.1, 0.603, 0.025, 0.579, 0.987, 0.309, 0.837, 0.643, 0.579, 0.982, 0.842, 0.417, 0.452, 0.778, 0.606, 0.719, 0.573, 0.773, 0.256, 0.712, 0.61, 0.558, 0.179, 0.233, 0.128, 0.411, 0.137, 0.229, 0.585, 0.312, 0.745, 0.063, 0.283, 0.641, 0.361, 0.234, 0.935, 0.037, 0.888, 0.284, 0.102, 0.435, 0.737, 0.125, 0.332, 0.254, 0.336, 0.526, 0.235, 0.546, 0.377, 0.64, 0.406, 0.589, 0.118, 0.167, 0.242, 0.18, 0.844, 0.094, 0.967, 0.946, 0.835, 0.771, 0.928, 0.261, 0.019, 0.716, 0.078, 0.424, 0.13, 0.615, 0.451, 0.775, 0.446, 0.121, 0.787, 0.854, 0.695, 0.969, 0.77, 0.754, 0.614, 0.021, 0.562, 0.375, 0.261, 0.772, 0.55, 0.822, 0.342, 0.706, 0.521, 0.935, 0.885, 0.438, 0.218, 0.823, 0.753, 0.929, 0.547, 0.226, 0.728, 0.937, 0.422, 0.323, 0.268, 0.286, 0.023, 0.749, 0.913, 0.819, 0.979, 0.31, 0.279, 0.356, 0.458, 0.949]
global q = [0.513, 0.44, 0.9, 0.632, 0.958, 0.801, 0.553, 0.423, 0.654, 0.708, 0.868, 0.999, 0.635, 0.98, 0.692, 0.961, 0.955, 0.923, 0.329, 0.685, 0.896, 0.395, 0.642, 0.979, 0.993, 0.849, 0.505, 0.535, 0.962, 0.857, 0.983, 0.939, 0.328, 0.996, 0.313, 0.803, 0.63, 0.941, 0.491, 0.843, 0.642, 0.945, 0.553, 0.842, 0.378, 0.847, 0.679, 0.581, 0.635, 0.793, 0.263, 0.716, 0.904, 0.695, 0.635, 0.706, 0.552, 0.141, 0.648, 0.747, 0.829, 0.805, 0.952, 0.888, 0.467, 0.536, 0.856, 0.705, 0.95, 0.796, 0.908, 0.904, 0.979, 0.969, 0.847, 0.985, 0.878, 0.611, 0.827, 0.891, 0.969, 0.833, 0.987, 0.575, 0.758, 0.966, 0.965, 0.784, 0.304, 0.698, 0.999, 0.98, 0.968, 0.273, 0.705, 0.819, 0.991, 0.965, 0.986, 0.445, 0.98, 0.985, 0.965, 0.91, 0.799, 0.406, 0.752, 0.762, 0.706, 0.559, 0.759, 0.771, 0.808, 0.993, 0.63, 0.138, 0.997, 0.993, 0.94, 0.809, 0.905, 0.711, 0.228, 0.458, 0.671, 0.931, 0.846, 0.286, 0.489, 0.625, 0.88, 0.227, 0.521, 0.947, 0.641, 0.765, 0.983, 0.788, 0.611, 0.817, 0.655, 0.954, 0.673, 0.968, 0.465, 0.299, 0.879, 0.637, 0.973, 0.599, 0.988, 0.558, 0.846, 0.887, 0.666, 0.983, 0.862, 0.608, 0.619, 0.98, 0.615, 0.997, 0.877, 0.91, 0.773, 0.99, 0.899, 0.639, 0.407, 0.919, 0.189, 0.785, 0.947, 0.485, 0.674, 0.846, 0.797, 0.749, 0.985, 0.74, 0.372, 0.7, 0.974, 0.537, 0.898, 0.602, 0.787, 0.457, 0.913, 0.846, 0.885, 0.436, 0.654, 0.586, 0.635, 0.566, 0.721, 0.935, 0.65, 0.722, 0.203, 0.895, 0.831, 0.699, 0.939, 0.117, 0.993, 0.983, 0.908, 0.881, 0.956, 0.599, 0.425, 0.739, 0.438, 0.787, 0.353, 0.692, 0.902, 0.991, 0.832, 0.668, 0.889, 0.95, 0.81, 0.989, 0.976, 0.791, 0.738, 0.061, 0.665, 0.532, 0.607, 0.886, 0.971, 0.825, 0.969, 0.837, 0.668, 0.944, 0.915, 0.674, 0.973, 0.995, 0.897, 0.975, 0.626, 0.48, 0.892, 0.962, 0.64, 0.323, 0.89, 0.648, 0.357, 0.893, 0.941, 0.972, 0.99, 0.985, 0.964, 0.711, 0.899, 0.993]
global origin = 1
global destination = 50