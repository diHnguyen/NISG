global arcs = [1 10; 1 16; 1 28; 1 32; 1 44; 1 50; 2 5; 2 15; 2 39; 2 43; 3 9; 3 15; 3 18; 3 27; 3 47; 4 12; 4 22; 4 36; 4 46; 4 50; 5 6; 5 26; 6 8; 6 23; 6 31; 6 33; 6 35; 6 45; 7 2; 7 11; 7 20; 7 27; 7 28; 7 31; 7 33; 7 49; 8 39; 9 10; 9 14; 9 17; 9 42; 10 16; 10 23; 10 24; 10 30; 10 32; 10 37; 11 5; 11 13; 11 18; 11 31; 12 3; 12 7; 12 16; 12 24; 12 43; 12 48; 13 3; 13 4; 13 11; 13 14; 13 15; 13 34; 13 40; 13 41; 13 49; 14 8; 14 32; 14 35; 15 13; 15 16; 15 19; 15 42; 15 46; 16 8; 16 13; 16 34; 16 43; 16 46; 17 16; 17 22; 17 48; 17 50; 18 29; 18 44; 18 46; 19 2; 19 8; 19 9; 19 15; 19 25; 19 26; 19 31; 19 35; 20 18; 20 23; 20 28; 20 32; 20 38; 20 46; 20 48; 21 6; 21 10; 21 13; 21 15; 21 19; 21 23; 21 33; 22 6; 22 13; 22 15; 22 17; 22 21; 22 24; 22 28; 22 32; 22 34; 22 42; 23 9; 23 25; 23 33; 23 47; 24 3; 24 6; 24 8; 24 9; 24 13; 24 18; 24 36; 24 49; 24 50; 25 9; 25 10; 25 14; 25 15; 25 19; 25 36; 25 43; 25 46; 25 50; 26 7; 26 14; 26 27; 26 30; 27 21; 27 30; 27 41; 27 45; 28 4; 28 13; 28 34; 28 35; 28 49; 29 17; 29 19; 29 31; 30 4; 30 16; 30 25; 30 37; 30 43; 31 2; 31 15; 31 29; 31 32; 31 33; 31 42; 32 22; 32 34; 32 37; 33 19; 33 23; 33 42; 34 4; 34 22; 34 27; 34 42; 34 43; 34 48; 35 25; 35 34; 35 42; 35 44; 36 3; 36 8; 36 13; 36 28; 36 31; 37 5; 38 13; 38 23; 39 7; 39 12; 39 27; 39 41; 39 42; 39 46; 40 22; 40 23; 40 36; 40 41; 40 43; 40 48; 41 2; 41 4; 41 30; 41 31; 41 48; 42 5; 42 19; 42 37; 43 4; 43 12; 43 27; 43 41; 44 7; 44 18; 44 40; 45 3; 45 14; 45 21; 45 22; 45 27; 45 28; 45 30; 45 31; 45 39; 46 19; 46 22; 46 27; 46 39; 47 10; 47 12; 47 20; 47 34; 47 42; 48 43; 49 9; 49 15; 49 18; 49 23; 49 33; 49 50]
global d_x = [7.0, 1.0, 7.0, 7.0, 3.0, 2.0, 7.0, 4.0, 6.0, 5.0, 9.0, 4.0, 10.0, 6.0, 2.0, 8.0, 5.0, 9.0, 9.0, 1.0, 7.0, 3.0, 7.0, 5.0, 1.0, 4.0, 8.0, 1.0, 9.0, 7.0, 8.0, 6.0, 7.0, 10.0, 7.0, 7.0, 5.0, 1.0, 6.0, 2.0, 3.0, 2.0, 9.0, 7.0, 5.0, 9.0, 5.0, 5.0, 2.0, 6.0, 3.0, 1.0, 1.0, 1.0, 7.0, 10.0, 2.0, 3.0, 10.0, 8.0, 7.0, 7.0, 1.0, 9.0, 10.0, 6.0, 5.0, 8.0, 2.0, 5.0, 1.0, 1.0, 6.0, 5.0, 1.0, 4.0, 3.0, 10.0, 8.0, 9.0, 5.0, 4.0, 10.0, 8.0, 1.0, 10.0, 7.0, 2.0, 3.0, 1.0, 4.0, 6.0, 9.0, 3.0, 9.0, 5.0, 4.0, 8.0, 6.0, 4.0, 4.0, 2.0, 8.0, 8.0, 5.0, 7.0, 9.0, 9.0, 1.0, 4.0, 9.0, 5.0, 1.0, 5.0, 9.0, 4.0, 9.0, 3.0, 5.0, 7.0, 6.0, 9.0, 6.0, 6.0, 3.0, 5.0, 3.0, 2.0, 10.0, 4.0, 9.0, 7.0, 7.0, 4.0, 8.0, 4.0, 1.0, 2.0, 1.0, 1.0, 5.0, 4.0, 7.0, 6.0, 5.0, 1.0, 3.0, 8.0, 2.0, 3.0, 6.0, 7.0, 3.0, 4.0, 8.0, 3.0, 7.0, 2.0, 8.0, 1.0, 6.0, 7.0, 3.0, 5.0, 3.0, 7.0, 8.0, 1.0, 2.0, 3.0, 2.0, 5.0, 4.0, 5.0, 9.0, 10.0, 2.0, 7.0, 10.0, 4.0, 6.0, 5.0, 8.0, 4.0, 3.0, 10.0, 5.0, 5.0, 6.0, 3.0, 2.0, 9.0, 9.0, 5.0, 8.0, 3.0, 7.0, 3.0, 9.0, 3.0, 1.0, 8.0, 6.0, 3.0, 8.0, 7.0, 5.0, 8.0, 2.0, 10.0, 3.0, 10.0, 4.0, 5.0, 8.0, 4.0, 10.0, 2.0, 7.0, 10.0, 2.0, 5.0, 6.0, 8.0, 6.0, 10.0, 6.0, 5.0, 7.0, 3.0, 1.0, 1.0, 3.0, 7.0, 7.0, 3.0, 9.0, 10.0, 10.0, 8.0, 8.0, 8.0, 10.0]
global b_x = 5
global d_y = [7.0, 5.0, 1.0, 8.0, 3.0, 6.0, 1.0, 1.0, 6.0, 10.0, 1.0, 1.0, 3.0, 10.0, 6.0, 6.0, 7.0, 1.0, 2.0, 8.0, 1.0, 6.0, 4.0, 5.0, 7.0, 2.0, 7.0, 6.0, 5.0, 6.0, 9.0, 4.0, 7.0, 1.0, 7.0, 9.0, 3.0, 3.0, 2.0, 6.0, 9.0, 9.0, 6.0, 9.0, 7.0, 10.0, 3.0, 9.0, 3.0, 1.0, 10.0, 2.0, 1.0, 7.0, 6.0, 9.0, 7.0, 10.0, 8.0, 8.0, 5.0, 9.0, 3.0, 3.0, 7.0, 7.0, 6.0, 5.0, 9.0, 9.0, 10.0, 7.0, 3.0, 9.0, 4.0, 6.0, 6.0, 10.0, 5.0, 3.0, 7.0, 2.0, 2.0, 5.0, 2.0, 3.0, 1.0, 5.0, 4.0, 10.0, 5.0, 9.0, 7.0, 7.0, 2.0, 1.0, 6.0, 5.0, 4.0, 8.0, 6.0, 6.0, 2.0, 4.0, 10.0, 2.0, 9.0, 3.0, 10.0, 8.0, 3.0, 4.0, 6.0, 2.0, 4.0, 1.0, 2.0, 1.0, 8.0, 6.0, 3.0, 9.0, 2.0, 2.0, 6.0, 7.0, 2.0, 10.0, 1.0, 6.0, 4.0, 2.0, 9.0, 4.0, 10.0, 9.0, 1.0, 7.0, 3.0, 7.0, 2.0, 2.0, 8.0, 8.0, 1.0, 9.0, 3.0, 1.0, 8.0, 5.0, 3.0, 7.0, 6.0, 4.0, 4.0, 1.0, 10.0, 5.0, 9.0, 1.0, 8.0, 9.0, 3.0, 9.0, 3.0, 8.0, 4.0, 4.0, 9.0, 9.0, 3.0, 2.0, 3.0, 4.0, 1.0, 5.0, 5.0, 1.0, 4.0, 1.0, 6.0, 2.0, 5.0, 5.0, 1.0, 10.0, 6.0, 8.0, 8.0, 4.0, 6.0, 4.0, 7.0, 4.0, 8.0, 8.0, 3.0, 2.0, 4.0, 3.0, 3.0, 1.0, 4.0, 9.0, 2.0, 8.0, 2.0, 7.0, 5.0, 6.0, 1.0, 7.0, 4.0, 1.0, 10.0, 5.0, 1.0, 6.0, 2.0, 7.0, 9.0, 3.0, 2.0, 5.0, 8.0, 2.0, 1.0, 7.0, 7.0, 9.0, 1.0, 10.0, 3.0, 8.0, 6.0, 3.0, 5.0, 8.0, 9.0, 7.0, 7.0, 1.0, 1.0]
global b_y = 10
global p = [0.943, 0.28, 0.742, 0.129, 0.128, 0.963, 0.972, 0.91, 0.175, 0.565, 0.329, 0.124, 0.646, 0.466, 0.679, 0.928, 0.135, 0.485, 0.062, 0.793, 0.513, 0.412, 0.376, 0.581, 0.668, 0.125, 0.845, 0.971, 0.556, 0.062, 0.767, 0.058, 0.046, 0.257, 0.752, 0.969, 0.347, 0.08, 0.094, 0.013, 0.738, 0.847, 0.405, 0.193, 0.891, 0.575, 0.356, 0.943, 0.283, 0.915, 0.7, 0.923, 0.848, 0.748, 0.71, 0.773, 0.036, 0.917, 0.735, 0.52, 0.935, 0.344, 0.184, 0.635, 0.267, 0.471, 0.827, 0.445, 0.939, 0.538, 0.964, 0.722, 0.396, 0.804, 0.793, 0.607, 0.728, 0.835, 0.735, 0.608, 0.115, 0.61, 0.254, 0.769, 0.647, 0.566, 0.867, 0.973, 0.789, 0.977, 0.217, 0.345, 0.184, 0.184, 0.018, 0.498, 0.96, 0.42, 0.806, 0.832, 0.174, 0.151, 0.969, 0.175, 0.366, 0.709, 0.343, 0.757, 0.18, 0.221, 0.238, 0.72, 0.595, 0.958, 0.501, 0.327, 0.614, 0.754, 0.724, 0.413, 0.972, 0.829, 0.522, 0.804, 0.802, 0.311, 0.336, 0.568, 0.58, 0.362, 0.5, 0.001, 0.558, 0.937, 0.425, 0.793, 0.335, 0.085, 0.094, 0.195, 0.427, 0.223, 0.291, 0.644, 0.583, 0.206, 0.584, 0.73, 0.689, 0.8, 0.869, 0.981, 0.126, 0.597, 0.313, 0.94, 0.261, 0.377, 0.084, 0.77, 0.9, 0.371, 0.77, 0.552, 0.456, 0.734, 0.547, 0.212, 0.702, 0.372, 0.637, 0.334, 0.52, 0.039, 0.998, 0.47, 0.665, 0.517, 0.373, 0.585, 0.743, 0.859, 0.997, 0.2, 0.807, 0.718, 0.569, 0.395, 0.812, 0.81, 0.847, 0.537, 0.264, 0.663, 0.2, 0.172, 0.178, 0.128, 0.222, 0.521, 0.799, 0.747, 0.422, 0.254, 0.944, 0.797, 0.124, 0.939, 0.286, 0.188, 0.089, 0.907, 0.201, 0.397, 0.431, 0.114, 0.318, 0.393, 0.619, 0.712, 0.046, 0.321, 0.277, 0.131, 0.822, 0.979, 0.465, 0.359, 0.443, 0.167, 0.365, 0.026, 0.635, 0.274, 0.238, 0.624, 0.318, 0.523, 0.169, 0.743, 0.894, 0.863, 0.982]
global q = [0.953, 0.364, 0.927, 0.876, 0.656, 0.98, 0.975, 0.929, 0.643, 0.797, 0.912, 0.608, 0.922, 0.748, 0.772, 0.968, 0.447, 0.71, 0.288, 0.921, 0.991, 0.913, 0.92, 0.932, 0.837, 0.795, 0.934, 0.984, 0.875, 0.539, 0.991, 0.119, 0.505, 0.818, 0.992, 0.992, 0.569, 0.921, 0.834, 0.371, 0.898, 0.917, 0.831, 0.816, 0.963, 0.626, 0.781, 0.977, 0.737, 0.928, 0.855, 0.955, 0.919, 0.856, 0.77, 0.973, 0.592, 0.928, 0.988, 0.68, 0.947, 0.835, 0.866, 0.796, 0.473, 0.512, 0.969, 0.794, 0.994, 0.731, 0.979, 0.907, 0.643, 0.974, 0.867, 0.949, 0.921, 0.845, 0.849, 0.628, 0.565, 0.792, 0.631, 0.923, 0.804, 0.775, 0.877, 0.983, 0.995, 0.985, 0.749, 0.936, 0.267, 0.673, 0.109, 0.667, 0.992, 0.569, 0.998, 0.99, 0.559, 0.68, 0.978, 0.412, 0.658, 0.876, 0.794, 0.807, 0.586, 0.352, 0.842, 0.954, 0.752, 0.962, 0.929, 0.533, 0.991, 0.989, 0.916, 0.44, 0.996, 0.874, 0.876, 0.984, 0.955, 0.673, 0.652, 0.919, 0.941, 0.987, 0.999, 0.365, 0.903, 0.994, 0.821, 0.944, 0.984, 0.63, 0.469, 0.503, 0.461, 0.835, 0.629, 0.965, 0.869, 0.641, 0.615, 0.853, 0.835, 0.887, 0.964, 0.992, 0.357, 0.802, 0.52, 0.981, 0.353, 0.397, 0.105, 0.917, 0.992, 0.591, 0.869, 0.95, 0.533, 0.737, 0.656, 0.456, 0.893, 0.939, 0.753, 0.392, 0.581, 0.947, 0.999, 0.484, 0.733, 0.882, 0.835, 0.667, 0.77, 0.965, 0.997, 0.985, 0.987, 0.782, 0.854, 0.904, 0.943, 0.915, 0.988, 0.649, 0.44, 0.672, 0.387, 0.84, 0.855, 0.395, 0.887, 0.942, 0.808, 0.991, 0.914, 0.571, 0.983, 0.897, 0.174, 0.992, 0.841, 0.579, 0.435, 0.973, 0.353, 0.928, 0.802, 0.443, 0.761, 0.607, 0.831, 0.743, 0.148, 0.674, 0.497, 0.373, 0.887, 0.99, 0.539, 0.454, 0.852, 0.944, 0.721, 0.36, 0.84, 0.799, 0.488, 0.876, 0.867, 0.719, 0.191, 0.852, 0.964, 0.961, 0.994]
global origin = 1
global destination = 50