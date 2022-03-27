global arcs = [1 4; 1 10; 1 12; 1 18; 1 20; 1 41; 1 54; 2 10; 2 18; 2 22; 2 24; 2 32; 2 49; 2 53; 2 54; 3 11; 3 16; 3 32; 3 50; 3 51; 4 3; 4 19; 4 23; 4 29; 5 4; 5 23; 5 46; 5 49; 5 56; 6 11; 6 18; 6 25; 6 44; 6 45; 6 47; 6 48; 6 55; 7 13; 7 22; 8 6; 8 48; 9 3; 9 18; 9 22; 9 29; 9 41; 9 43; 9 44; 10 7; 10 26; 10 35; 10 42; 10 59; 11 17; 11 18; 11 26; 11 29; 11 33; 11 37; 11 39; 11 41; 11 42; 11 45; 11 59; 11 60; 12 4; 12 15; 12 17; 12 22; 12 27; 12 28; 12 52; 13 4; 13 22; 13 24; 13 32; 13 52; 14 12; 14 13; 14 15; 14 20; 14 32; 15 2; 15 21; 15 49; 15 56; 15 58; 16 17; 16 31; 16 49; 16 50; 17 2; 17 3; 17 8; 17 11; 17 20; 18 11; 18 20; 18 22; 18 47; 19 37; 19 38; 20 8; 20 19; 20 28; 20 48; 20 56; 21 12; 21 28; 21 36; 21 38; 21 41; 21 43; 21 53; 22 6; 22 12; 22 23; 22 31; 22 36; 22 38; 22 60; 23 2; 23 11; 23 34; 23 39; 23 44; 24 4; 24 11; 24 15; 24 34; 24 41; 24 45; 25 4; 25 9; 25 10; 25 11; 25 13; 25 20; 25 22; 26 13; 26 24; 26 25; 26 27; 26 47; 27 25; 27 26; 27 49; 27 51; 27 53; 27 58; 28 21; 28 38; 28 49; 28 59; 29 16; 29 37; 29 42; 29 59; 30 13; 30 40; 30 41; 30 47; 30 59; 31 16; 31 19; 31 21; 31 23; 31 35; 31 47; 31 53; 31 60; 32 50; 32 51; 33 11; 33 14; 33 20; 33 24; 33 29; 33 34; 33 38; 33 39; 33 41; 33 52; 33 53; 33 57; 33 60; 34 12; 34 23; 34 27; 34 56; 35 4; 35 27; 35 28; 35 30; 35 42; 36 2; 36 3; 36 9; 36 18; 36 26; 36 42; 36 54; 37 9; 37 27; 37 36; 37 40; 37 42; 38 15; 38 17; 38 29; 38 32; 38 41; 39 2; 39 15; 39 21; 39 24; 39 32; 39 38; 40 2; 40 7; 40 42; 40 43; 40 60; 41 26; 41 47; 42 6; 42 15; 42 43; 42 56; 42 58; 43 39; 43 45; 43 50; 43 53; 43 56; 44 3; 44 4; 44 33; 44 38; 44 41; 45 12; 45 18; 45 21; 45 30; 45 49; 45 51; 45 52; 46 10; 46 38; 47 16; 47 18; 47 19; 47 26; 47 34; 48 20; 48 37; 48 46; 49 4; 49 5; 49 6; 49 9; 49 10; 49 57; 49 60; 50 5; 50 12; 50 13; 50 18; 50 30; 50 31; 50 41; 50 44; 50 47; 51 4; 51 17; 51 40; 52 6; 52 7; 52 14; 52 17; 52 19; 52 24; 52 32; 52 37; 52 41; 52 55; 52 59; 53 11; 53 30; 53 32; 53 47; 54 8; 54 10; 54 12; 54 23; 54 28; 54 35; 54 43; 55 9; 55 13; 55 20; 55 22; 55 23; 55 27; 55 52; 56 2; 56 26; 57 8; 57 11; 57 58; 58 24; 58 31; 58 37; 58 52; 59 2; 59 9; 59 11; 59 18; 59 23; 59 40; 59 53; 59 54]
global d_x = [9.0, 6.0, 4.0, 10.0, 8.0, 6.0, 10.0, 2.0, 1.0, 7.0, 6.0, 3.0, 9.0, 3.0, 2.0, 2.0, 3.0, 3.0, 3.0, 10.0, 6.0, 4.0, 5.0, 2.0, 8.0, 1.0, 8.0, 10.0, 5.0, 9.0, 6.0, 10.0, 4.0, 4.0, 1.0, 6.0, 6.0, 5.0, 9.0, 1.0, 7.0, 6.0, 6.0, 2.0, 2.0, 3.0, 4.0, 1.0, 6.0, 7.0, 8.0, 8.0, 4.0, 1.0, 10.0, 9.0, 1.0, 3.0, 1.0, 10.0, 5.0, 2.0, 2.0, 5.0, 8.0, 10.0, 9.0, 4.0, 10.0, 1.0, 5.0, 5.0, 8.0, 2.0, 6.0, 9.0, 3.0, 6.0, 10.0, 5.0, 5.0, 6.0, 5.0, 3.0, 3.0, 9.0, 7.0, 9.0, 6.0, 4.0, 10.0, 6.0, 10.0, 2.0, 7.0, 6.0, 1.0, 8.0, 3.0, 7.0, 8.0, 4.0, 9.0, 7.0, 9.0, 8.0, 4.0, 2.0, 1.0, 1.0, 7.0, 4.0, 1.0, 10.0, 1.0, 7.0, 6.0, 8.0, 7.0, 2.0, 1.0, 6.0, 4.0, 8.0, 7.0, 8.0, 5.0, 8.0, 8.0, 6.0, 6.0, 9.0, 10.0, 3.0, 2.0, 3.0, 10.0, 1.0, 6.0, 4.0, 4.0, 6.0, 9.0, 3.0, 7.0, 6.0, 10.0, 6.0, 3.0, 5.0, 10.0, 1.0, 5.0, 10.0, 6.0, 4.0, 10.0, 10.0, 8.0, 4.0, 9.0, 7.0, 10.0, 9.0, 8.0, 3.0, 10.0, 8.0, 10.0, 2.0, 10.0, 2.0, 7.0, 4.0, 8.0, 10.0, 4.0, 4.0, 1.0, 2.0, 8.0, 6.0, 1.0, 1.0, 3.0, 10.0, 5.0, 10.0, 7.0, 1.0, 3.0, 6.0, 7.0, 1.0, 4.0, 9.0, 7.0, 5.0, 9.0, 7.0, 4.0, 5.0, 7.0, 2.0, 5.0, 1.0, 4.0, 10.0, 7.0, 7.0, 1.0, 3.0, 1.0, 1.0, 6.0, 9.0, 9.0, 7.0, 9.0, 10.0, 2.0, 4.0, 7.0, 8.0, 1.0, 10.0, 1.0, 2.0, 5.0, 4.0, 8.0, 6.0, 2.0, 10.0, 8.0, 8.0, 10.0, 1.0, 9.0, 10.0, 8.0, 6.0, 3.0, 4.0, 8.0, 9.0, 3.0, 5.0, 3.0, 4.0, 10.0, 2.0, 3.0, 4.0, 4.0, 7.0, 3.0, 9.0, 3.0, 6.0, 4.0, 4.0, 10.0, 10.0, 1.0, 2.0, 9.0, 5.0, 7.0, 2.0, 1.0, 10.0, 5.0, 8.0, 8.0, 3.0, 2.0, 9.0, 10.0, 2.0, 1.0, 7.0, 4.0, 9.0, 10.0, 6.0, 4.0, 1.0, 2.0, 9.0, 7.0, 8.0, 3.0, 5.0, 3.0, 4.0, 3.0, 9.0, 9.0, 2.0, 9.0, 10.0, 2.0, 8.0, 5.0, 5.0, 8.0, 7.0, 6.0, 9.0, 3.0, 1.0, 1.0, 4.0, 10.0, 8.0, 5.0, 8.0, 5.0, 9.0, 7.0, 5.0]
global b_x = 5
global d_y = [9.0, 3.0, 4.0, 3.0, 6.0, 4.0, 8.0, 8.0, 3.0, 5.0, 8.0, 2.0, 2.0, 5.0, 5.0, 5.0, 10.0, 4.0, 7.0, 9.0, 10.0, 9.0, 6.0, 8.0, 8.0, 4.0, 6.0, 10.0, 9.0, 4.0, 8.0, 6.0, 3.0, 10.0, 5.0, 2.0, 6.0, 3.0, 7.0, 4.0, 10.0, 9.0, 7.0, 8.0, 6.0, 6.0, 6.0, 4.0, 2.0, 10.0, 3.0, 7.0, 9.0, 4.0, 7.0, 3.0, 9.0, 9.0, 7.0, 6.0, 5.0, 6.0, 8.0, 2.0, 4.0, 9.0, 10.0, 10.0, 10.0, 9.0, 9.0, 5.0, 10.0, 5.0, 10.0, 2.0, 5.0, 9.0, 5.0, 7.0, 1.0, 1.0, 3.0, 3.0, 3.0, 4.0, 9.0, 8.0, 1.0, 1.0, 10.0, 2.0, 1.0, 2.0, 3.0, 7.0, 2.0, 8.0, 8.0, 7.0, 4.0, 6.0, 5.0, 4.0, 1.0, 2.0, 3.0, 8.0, 7.0, 9.0, 4.0, 6.0, 7.0, 2.0, 4.0, 5.0, 10.0, 2.0, 3.0, 10.0, 8.0, 8.0, 5.0, 2.0, 6.0, 6.0, 9.0, 3.0, 7.0, 5.0, 7.0, 1.0, 2.0, 6.0, 8.0, 9.0, 9.0, 7.0, 4.0, 4.0, 7.0, 10.0, 2.0, 10.0, 10.0, 1.0, 7.0, 9.0, 7.0, 8.0, 5.0, 5.0, 8.0, 9.0, 6.0, 7.0, 3.0, 8.0, 5.0, 1.0, 8.0, 8.0, 6.0, 1.0, 8.0, 9.0, 6.0, 2.0, 9.0, 1.0, 4.0, 8.0, 5.0, 2.0, 4.0, 2.0, 1.0, 2.0, 2.0, 6.0, 8.0, 4.0, 8.0, 2.0, 7.0, 7.0, 2.0, 1.0, 1.0, 1.0, 9.0, 1.0, 9.0, 7.0, 5.0, 9.0, 1.0, 9.0, 5.0, 10.0, 3.0, 3.0, 8.0, 4.0, 3.0, 4.0, 4.0, 3.0, 8.0, 1.0, 7.0, 10.0, 7.0, 2.0, 1.0, 2.0, 1.0, 5.0, 10.0, 10.0, 9.0, 7.0, 4.0, 6.0, 10.0, 2.0, 8.0, 3.0, 5.0, 4.0, 4.0, 8.0, 9.0, 2.0, 8.0, 7.0, 3.0, 4.0, 3.0, 8.0, 3.0, 6.0, 5.0, 10.0, 6.0, 4.0, 8.0, 4.0, 7.0, 8.0, 1.0, 10.0, 9.0, 7.0, 9.0, 2.0, 9.0, 3.0, 8.0, 8.0, 8.0, 1.0, 6.0, 10.0, 1.0, 1.0, 1.0, 8.0, 6.0, 9.0, 4.0, 9.0, 4.0, 2.0, 2.0, 6.0, 10.0, 1.0, 2.0, 3.0, 8.0, 9.0, 2.0, 6.0, 9.0, 10.0, 5.0, 6.0, 8.0, 7.0, 1.0, 5.0, 1.0, 3.0, 7.0, 3.0, 10.0, 7.0, 4.0, 7.0, 2.0, 1.0, 2.0, 2.0, 8.0, 9.0, 1.0, 9.0, 1.0, 3.0, 4.0, 5.0, 4.0, 10.0, 5.0, 2.0, 1.0, 5.0, 8.0, 4.0, 4.0, 9.0]
global b_y = 10
global p = [0.439, 0.495, 0.161, 0.751, 0.039, 0.725, 0.306, 0.111, 0.242, 0.556, 0.828, 0.82, 0.728, 0.089, 0.983, 0.497, 0.406, 0.027, 0.249, 0.488, 0.33, 0.973, 0.539, 0.312, 0.791, 0.307, 0.386, 0.878, 0.651, 0.459, 0.772, 0.974, 0.871, 0.172, 0.186, 0.356, 0.427, 0.017, 0.668, 0.858, 0.372, 0.218, 0.298, 0.954, 0.294, 0.912, 0.181, 0.17, 0.674, 0.755, 0.231, 0.891, 0.592, 0.518, 0.279, 0.392, 0.586, 0.329, 0.043, 0.118, 0.147, 0.734, 0.179, 0.578, 0.871, 0.269, 0.431, 0.124, 0.885, 0.276, 0.978, 0.61, 0.16, 0.902, 0.335, 0.111, 0.959, 0.963, 0.581, 0.861, 0.631, 0.675, 0.919, 0.155, 0.492, 0.043, 0.512, 0.475, 0.784, 0.847, 0.431, 0.849, 0.672, 0.711, 0.996, 0.588, 0.301, 0.943, 0.219, 0.72, 0.38, 0.576, 0.867, 0.276, 0.681, 0.915, 0.615, 0.275, 0.181, 0.275, 0.602, 0.772, 0.95, 0.974, 0.763, 0.905, 0.962, 0.237, 0.28, 0.44, 0.345, 0.698, 0.559, 0.373, 0.191, 0.665, 0.861, 0.094, 0.585, 0.353, 0.249, 0.291, 0.429, 0.302, 0.619, 0.709, 0.427, 0.061, 0.438, 0.725, 0.104, 0.924, 0.429, 0.531, 0.137, 0.285, 0.487, 0.972, 0.056, 0.073, 0.623, 0.311, 0.855, 0.939, 0.343, 0.238, 0.64, 0.13, 0.775, 0.342, 0.326, 0.514, 0.292, 0.1, 0.168, 0.247, 0.184, 0.611, 0.447, 0.383, 0.351, 0.971, 0.703, 0.028, 0.522, 0.336, 0.439, 0.965, 0.085, 0.072, 0.083, 0.311, 0.171, 0.485, 0.57, 0.153, 0.11, 0.605, 0.645, 0.776, 0.373, 0.813, 0.755, 0.651, 0.547, 0.666, 0.282, 0.176, 0.831, 0.456, 0.539, 0.077, 0.124, 0.554, 0.359, 0.868, 0.562, 0.157, 0.265, 0.606, 0.865, 0.069, 0.026, 0.927, 0.593, 0.864, 0.42, 0.534, 0.594, 0.816, 0.841, 0.261, 0.554, 0.778, 0.128, 0.673, 0.204, 0.744, 0.096, 0.923, 0.192, 0.122, 0.873, 0.143, 0.06, 0.415, 0.008, 0.382, 0.183, 0.242, 0.263, 0.052, 0.978, 0.302, 0.091, 0.954, 0.819, 0.932, 0.631, 0.995, 0.703, 0.473, 0.585, 0.531, 0.446, 0.808, 0.379, 0.42, 0.765, 0.396, 0.492, 0.021, 0.957, 0.496, 0.481, 0.863, 0.907, 0.955, 0.671, 0.608, 0.899, 0.767, 0.327, 0.595, 0.443, 0.265, 0.156, 0.201, 0.687, 0.238, 0.545, 0.436, 0.768, 0.227, 0.605, 0.429, 0.098, 0.498, 0.719, 0.5, 0.202, 0.332, 0.475, 0.329, 0.622, 0.89, 0.728, 0.824, 0.363, 0.305, 0.807, 0.299, 0.216, 0.67, 0.791, 0.751, 0.693, 0.108, 0.312, 0.7, 0.051, 0.349, 0.93, 0.711, 0.761, 0.114, 0.745, 0.511, 0.695, 0.602, 0.425, 0.901]
global q = [0.459, 0.779, 0.362, 0.865, 0.223, 0.852, 0.466, 0.652, 0.538, 0.992, 0.918, 0.962, 0.915, 0.51, 0.99, 0.749, 0.865, 0.929, 0.403, 0.91, 0.545, 0.989, 0.682, 0.89, 0.792, 0.365, 0.406, 0.884, 0.678, 0.479, 0.952, 0.974, 0.985, 0.988, 0.672, 0.599, 0.79, 0.784, 0.674, 0.981, 0.63, 0.228, 0.943, 0.968, 0.721, 0.935, 0.924, 0.312, 0.947, 0.937, 0.545, 0.992, 0.747, 0.791, 0.698, 0.92, 0.714, 0.823, 0.047, 0.452, 0.484, 0.862, 0.93, 0.9, 0.987, 0.413, 0.616, 0.192, 0.912, 0.952, 0.99, 0.785, 0.959, 0.956, 0.409, 0.785, 0.97, 0.966, 0.927, 0.989, 0.707, 0.676, 0.921, 0.681, 0.758, 0.576, 0.779, 0.868, 0.94, 0.954, 0.538, 0.939, 0.756, 0.74, 0.998, 0.902, 0.698, 0.961, 0.716, 0.909, 0.954, 0.659, 0.893, 0.377, 0.773, 0.999, 0.93, 0.949, 0.224, 0.358, 0.758, 0.933, 0.97, 0.989, 0.77, 0.94, 0.982, 0.374, 0.864, 0.837, 0.672, 0.967, 0.88, 0.468, 0.268, 0.697, 0.888, 0.972, 0.955, 0.977, 0.529, 0.992, 0.944, 0.41, 0.63, 0.863, 0.837, 0.635, 0.55, 0.895, 0.3, 0.985, 0.452, 0.723, 0.76, 0.585, 0.994, 0.985, 0.887, 0.327, 0.742, 0.627, 0.936, 0.955, 0.813, 0.503, 0.732, 0.787, 0.788, 0.61, 0.463, 0.762, 0.926, 0.502, 0.503, 0.413, 0.351, 0.789, 0.984, 0.719, 0.367, 0.973, 0.874, 0.77, 0.756, 0.535, 0.891, 0.988, 0.513, 0.377, 0.601, 0.395, 0.233, 0.903, 0.955, 0.701, 0.758, 0.962, 0.866, 0.913, 0.564, 0.949, 0.872, 0.665, 0.862, 0.859, 0.44, 0.548, 0.871, 0.711, 0.657, 0.507, 0.706, 0.622, 0.832, 0.876, 0.775, 0.441, 0.288, 0.766, 0.97, 0.738, 0.82, 0.969, 0.707, 0.877, 0.845, 0.93, 0.752, 0.974, 0.966, 0.561, 0.887, 0.929, 0.262, 0.68, 0.329, 0.97, 0.517, 0.99, 0.943, 0.51, 0.907, 0.176, 0.653, 0.707, 0.575, 0.637, 0.701, 0.754, 0.286, 0.44, 0.985, 0.833, 0.365, 0.965, 0.991, 0.992, 0.711, 0.999, 0.835, 0.757, 0.664, 0.76, 0.572, 0.859, 0.761, 0.486, 0.84, 0.523, 0.78, 0.157, 0.993, 0.529, 0.605, 0.893, 0.95, 0.99, 0.902, 0.982, 0.922, 0.842, 0.916, 0.699, 0.869, 0.992, 0.68, 0.696, 0.987, 0.701, 0.941, 0.717, 0.892, 0.67, 0.984, 0.993, 0.898, 0.926, 0.882, 0.946, 0.886, 0.422, 0.885, 0.976, 0.705, 0.952, 0.828, 0.884, 0.545, 0.896, 0.963, 0.304, 0.9, 0.991, 0.98, 0.757, 0.724, 0.293, 0.593, 0.844, 0.317, 0.924, 0.942, 0.769, 0.991, 0.63, 0.994, 0.897, 0.818, 0.642, 0.687, 0.903]
global origin = 1
global destination = 60