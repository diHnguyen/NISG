global arcs = [1 6; 1 11; 1 17; 1 19; 1 35; 2 16; 2 32; 2 36; 2 47; 3 14; 3 16; 3 23; 4 21; 4 32; 5 13; 5 18; 5 19; 5 21; 5 39; 6 23; 6 33; 6 49; 7 2; 7 11; 7 16; 7 38; 8 4; 8 21; 8 23; 8 45; 9 13; 9 19; 9 23; 9 37; 9 45; 10 8; 10 13; 10 14; 10 29; 10 37; 11 5; 11 7; 11 23; 11 41; 11 49; 12 4; 12 6; 12 10; 12 13; 12 17; 12 35; 13 3; 13 8; 13 26; 13 35; 13 43; 14 10; 14 33; 15 4; 15 9; 15 12; 15 13; 15 34; 15 41; 15 46; 15 49; 16 12; 16 18; 16 22; 16 33; 17 6; 17 19; 17 37; 17 39; 18 25; 18 33; 18 40; 19 5; 19 11; 19 15; 19 16; 19 27; 19 30; 20 14; 20 15; 20 31; 20 34; 20 37; 20 48; 21 19; 21 25; 21 44; 21 46; 22 3; 22 6; 22 23; 22 30; 22 31; 22 45; 23 18; 23 26; 23 45; 23 49; 24 2; 24 3; 24 12; 24 13; 24 31; 24 38; 24 42; 25 2; 25 20; 25 26; 26 5; 26 6; 26 11; 26 15; 26 22; 26 29; 26 33; 27 9; 27 15; 27 18; 27 25; 27 32; 27 42; 27 43; 28 13; 28 25; 28 29; 28 43; 29 2; 29 13; 29 19; 29 25; 29 30; 29 37; 29 50; 30 5; 30 11; 30 15; 30 31; 30 47; 31 3; 31 10; 31 11; 31 36; 31 38; 32 15; 32 21; 32 23; 32 24; 32 25; 33 5; 33 12; 33 19; 33 23; 33 36; 34 3; 34 17; 34 23; 34 25; 34 36; 34 46; 34 48; 35 41; 35 44; 36 7; 36 12; 36 33; 36 34; 36 46; 37 31; 37 32; 38 2; 38 29; 38 37; 38 39; 39 3; 39 6; 39 8; 39 35; 39 40; 39 41; 39 43; 39 50; 40 3; 40 13; 40 21; 40 31; 40 38; 41 16; 41 22; 41 30; 41 43; 41 48; 42 10; 42 23; 42 35; 42 39; 42 40; 42 47; 43 14; 43 26; 43 28; 43 29; 43 36; 44 3; 44 6; 44 16; 44 26; 44 28; 45 8; 45 11; 45 12; 45 14; 45 22; 45 23; 45 42; 45 47; 46 6; 46 10; 46 19; 46 27; 46 44; 47 13; 47 14; 47 32; 47 33; 47 37; 48 5; 48 26; 48 29; 48 46; 49 4; 49 13; 49 20; 49 26; 49 28; 49 37; 49 42; 49 45; 49 46]
global d_x = [10.0, 8.0, 4.0, 3.0, 4.0, 6.0, 2.0, 8.0, 3.0, 8.0, 9.0, 3.0, 8.0, 10.0, 1.0, 6.0, 9.0, 7.0, 5.0, 10.0, 1.0, 1.0, 8.0, 1.0, 4.0, 10.0, 4.0, 6.0, 5.0, 9.0, 9.0, 8.0, 7.0, 9.0, 5.0, 1.0, 1.0, 7.0, 10.0, 5.0, 2.0, 4.0, 7.0, 1.0, 4.0, 1.0, 5.0, 4.0, 7.0, 3.0, 3.0, 5.0, 1.0, 5.0, 6.0, 7.0, 9.0, 9.0, 7.0, 2.0, 2.0, 4.0, 6.0, 10.0, 4.0, 4.0, 6.0, 4.0, 4.0, 1.0, 5.0, 4.0, 10.0, 7.0, 3.0, 2.0, 2.0, 6.0, 9.0, 5.0, 2.0, 2.0, 7.0, 9.0, 9.0, 10.0, 3.0, 2.0, 10.0, 7.0, 5.0, 2.0, 1.0, 1.0, 5.0, 1.0, 2.0, 9.0, 5.0, 8.0, 6.0, 9.0, 5.0, 4.0, 6.0, 3.0, 7.0, 3.0, 2.0, 2.0, 5.0, 8.0, 4.0, 2.0, 1.0, 1.0, 6.0, 1.0, 7.0, 10.0, 3.0, 8.0, 4.0, 6.0, 10.0, 2.0, 2.0, 9.0, 3.0, 4.0, 4.0, 6.0, 7.0, 10.0, 1.0, 3.0, 3.0, 7.0, 1.0, 6.0, 4.0, 8.0, 6.0, 3.0, 9.0, 4.0, 10.0, 2.0, 1.0, 5.0, 8.0, 9.0, 1.0, 7.0, 4.0, 1.0, 4.0, 6.0, 7.0, 4.0, 3.0, 3.0, 5.0, 6.0, 4.0, 7.0, 9.0, 10.0, 9.0, 1.0, 10.0, 8.0, 1.0, 4.0, 7.0, 10.0, 5.0, 10.0, 9.0, 1.0, 8.0, 3.0, 5.0, 10.0, 6.0, 10.0, 4.0, 2.0, 1.0, 6.0, 1.0, 10.0, 5.0, 4.0, 5.0, 10.0, 9.0, 4.0, 1.0, 5.0, 3.0, 10.0, 5.0, 5.0, 9.0, 6.0, 2.0, 8.0, 4.0, 9.0, 8.0, 8.0, 4.0, 2.0, 4.0, 7.0, 5.0, 2.0, 1.0, 2.0, 5.0, 9.0, 4.0, 9.0, 9.0, 9.0, 3.0, 3.0, 6.0, 9.0, 7.0, 1.0, 6.0, 7.0, 8.0, 7.0, 6.0, 3.0, 6.0, 7.0, 2.0, 7.0, 10.0]
global b_x = 5
global d_y = [3.0, 10.0, 7.0, 3.0, 8.0, 7.0, 1.0, 5.0, 10.0, 8.0, 2.0, 3.0, 5.0, 1.0, 1.0, 5.0, 5.0, 6.0, 3.0, 2.0, 3.0, 5.0, 2.0, 3.0, 5.0, 7.0, 8.0, 5.0, 9.0, 4.0, 3.0, 10.0, 4.0, 5.0, 2.0, 5.0, 7.0, 3.0, 4.0, 6.0, 10.0, 4.0, 4.0, 10.0, 4.0, 1.0, 10.0, 10.0, 1.0, 4.0, 5.0, 8.0, 1.0, 6.0, 1.0, 3.0, 4.0, 2.0, 8.0, 2.0, 1.0, 4.0, 8.0, 5.0, 6.0, 3.0, 2.0, 4.0, 3.0, 1.0, 7.0, 5.0, 1.0, 2.0, 8.0, 9.0, 8.0, 1.0, 1.0, 5.0, 5.0, 9.0, 8.0, 7.0, 5.0, 8.0, 4.0, 5.0, 6.0, 6.0, 5.0, 10.0, 7.0, 8.0, 2.0, 10.0, 2.0, 2.0, 5.0, 10.0, 7.0, 5.0, 4.0, 5.0, 3.0, 9.0, 8.0, 3.0, 3.0, 9.0, 7.0, 5.0, 6.0, 3.0, 2.0, 8.0, 2.0, 4.0, 2.0, 3.0, 3.0, 9.0, 5.0, 6.0, 10.0, 4.0, 7.0, 10.0, 7.0, 1.0, 9.0, 5.0, 5.0, 2.0, 1.0, 5.0, 6.0, 10.0, 5.0, 1.0, 10.0, 10.0, 2.0, 6.0, 5.0, 4.0, 6.0, 10.0, 3.0, 6.0, 9.0, 9.0, 7.0, 7.0, 8.0, 8.0, 4.0, 2.0, 7.0, 6.0, 1.0, 1.0, 3.0, 10.0, 2.0, 7.0, 4.0, 3.0, 1.0, 8.0, 6.0, 4.0, 4.0, 9.0, 6.0, 5.0, 1.0, 2.0, 4.0, 4.0, 2.0, 5.0, 6.0, 7.0, 1.0, 5.0, 9.0, 5.0, 8.0, 6.0, 8.0, 6.0, 10.0, 4.0, 8.0, 10.0, 10.0, 1.0, 4.0, 7.0, 3.0, 5.0, 10.0, 6.0, 1.0, 1.0, 8.0, 7.0, 4.0, 9.0, 10.0, 5.0, 10.0, 9.0, 8.0, 5.0, 4.0, 2.0, 7.0, 1.0, 1.0, 5.0, 6.0, 1.0, 3.0, 7.0, 3.0, 2.0, 8.0, 3.0, 4.0, 10.0, 4.0, 5.0, 7.0, 7.0, 1.0, 5.0, 4.0, 8.0, 3.0, 10.0, 7.0]
global b_y = 10
global p = [0.064, 0.692, 0.993, 0.616, 0.884, 0.615, 0.745, 0.18, 0.289, 0.51, 0.815, 0.514, 0.356, 0.459, 0.041, 0.215, 0.284, 0.155, 0.223, 0.617, 0.232, 0.954, 0.167, 0.486, 0.512, 0.488, 0.674, 0.095, 0.458, 0.98, 0.384, 0.336, 0.805, 0.823, 0.223, 0.309, 0.825, 0.605, 0.529, 0.159, 0.06, 0.337, 0.478, 0.08, 0.79, 0.28, 0.776, 0.006, 0.22, 0.005, 0.719, 0.753, 0.552, 0.895, 0.894, 0.2, 0.675, 0.49, 0.097, 0.327, 0.354, 0.077, 0.559, 0.278, 0.87, 0.147, 0.446, 0.635, 0.273, 0.222, 0.095, 0.716, 0.944, 0.464, 0.148, 0.081, 0.521, 0.695, 0.916, 0.368, 0.15, 0.413, 0.455, 0.172, 0.477, 0.668, 0.839, 0.963, 0.018, 0.971, 0.056, 0.667, 0.158, 0.865, 0.183, 0.005, 0.757, 0.693, 0.78, 0.227, 0.922, 0.863, 0.54, 0.352, 0.051, 0.292, 0.313, 0.557, 0.225, 0.696, 0.522, 0.092, 0.803, 0.774, 0.43, 0.38, 0.384, 0.642, 0.119, 0.324, 0.021, 0.465, 0.991, 0.93, 0.393, 0.149, 0.124, 0.249, 0.67, 0.349, 0.346, 0.695, 0.731, 0.77, 0.429, 0.775, 0.967, 0.405, 0.809, 0.213, 0.383, 0.689, 0.167, 0.116, 0.206, 0.11, 0.75, 0.247, 0.014, 0.962, 0.272, 0.324, 0.359, 0.945, 0.203, 0.987, 0.918, 0.347, 0.683, 0.038, 0.885, 0.85, 0.511, 0.84, 0.268, 0.461, 0.069, 0.15, 0.045, 0.018, 0.895, 0.119, 0.022, 0.497, 0.001, 0.993, 0.745, 0.812, 0.494, 0.632, 0.5, 0.551, 0.917, 0.178, 0.318, 0.901, 0.498, 0.951, 0.447, 0.112, 0.608, 0.831, 0.281, 0.335, 0.153, 0.691, 0.659, 0.251, 0.031, 0.919, 0.463, 0.5, 0.805, 0.22, 0.37, 0.945, 0.125, 0.152, 0.081, 0.525, 0.142, 0.162, 0.375, 0.727, 0.255, 0.012, 0.456, 0.647, 0.624, 0.935, 0.491, 0.86, 0.848, 0.852, 0.075, 0.479, 0.58, 0.083, 0.328, 0.65, 0.733, 0.581, 0.312, 0.677, 0.307, 0.566, 0.372, 0.859, 0.944, 0.4, 0.054, 0.994, 0.965]
global q = [0.414, 0.752, 0.993, 0.847, 0.942, 0.882, 0.865, 0.381, 0.731, 0.654, 0.968, 0.591, 0.588, 0.545, 0.936, 0.261, 0.867, 0.181, 0.688, 0.895, 0.862, 0.981, 0.383, 0.686, 0.711, 0.711, 0.866, 0.706, 0.955, 0.98, 0.533, 0.912, 0.954, 0.854, 0.914, 0.598, 0.843, 0.664, 0.874, 0.506, 0.298, 0.369, 0.569, 0.554, 0.793, 0.672, 0.827, 0.641, 0.554, 0.181, 0.999, 0.779, 0.556, 0.993, 0.987, 0.436, 0.805, 0.737, 0.336, 0.623, 0.96, 0.298, 0.788, 0.595, 0.985, 0.956, 0.886, 0.776, 0.324, 0.385, 0.241, 0.973, 0.997, 0.724, 0.153, 0.093, 0.592, 0.79, 0.975, 0.716, 0.551, 0.693, 0.793, 0.508, 0.56, 0.683, 0.98, 0.979, 0.715, 0.995, 0.829, 0.971, 0.391, 0.866, 0.262, 0.49, 0.939, 0.843, 0.93, 0.846, 0.928, 0.94, 0.787, 0.53, 0.881, 0.862, 0.773, 0.844, 0.957, 0.7, 0.939, 0.161, 0.824, 0.937, 0.857, 0.612, 0.58, 0.809, 0.731, 0.842, 0.838, 0.517, 0.992, 0.991, 0.797, 0.881, 0.412, 0.893, 0.685, 0.773, 0.957, 0.861, 0.736, 0.948, 0.751, 0.994, 0.99, 0.491, 0.951, 0.505, 0.528, 0.867, 0.36, 0.953, 0.844, 0.937, 0.799, 0.279, 0.739, 0.998, 0.801, 0.33, 0.843, 0.951, 0.672, 0.99, 0.939, 0.658, 0.807, 0.946, 0.991, 0.919, 0.686, 0.892, 0.596, 0.504, 0.619, 0.817, 0.381, 0.511, 0.906, 0.638, 0.591, 0.738, 0.101, 0.995, 0.932, 0.873, 0.727, 0.749, 0.549, 0.704, 0.929, 0.443, 0.697, 0.993, 0.877, 0.958, 0.581, 0.709, 0.805, 0.862, 0.508, 0.796, 0.795, 0.939, 0.731, 0.409, 0.103, 0.986, 0.884, 0.609, 0.968, 0.598, 0.781, 0.972, 0.451, 0.96, 0.543, 0.587, 0.341, 0.597, 0.71, 0.748, 0.447, 0.414, 0.525, 0.774, 0.983, 0.984, 0.533, 0.911, 0.998, 0.982, 0.937, 0.774, 0.764, 0.549, 0.775, 0.942, 0.789, 0.696, 0.794, 0.799, 0.665, 0.775, 0.705, 0.904, 0.981, 0.979, 0.807, 0.996, 0.978]
global origin = 1
global destination = 50