global arcs = [1 10; 1 11; 1 18; 1 30; 1 32; 1 36; 1 39; 1 42; 1 50; 2 18; 2 24; 2 38; 2 40; 2 47; 3 6; 3 11; 3 18; 3 27; 3 44; 4 5; 4 7; 4 12; 4 14; 4 30; 4 38; 5 2; 5 16; 5 21; 5 29; 5 35; 5 45; 6 17; 6 31; 6 34; 7 18; 7 23; 8 4; 8 30; 8 46; 9 2; 9 13; 9 32; 9 41; 9 45; 9 49; 10 2; 10 4; 10 12; 10 13; 10 15; 10 23; 10 26; 11 4; 11 20; 11 28; 11 29; 11 36; 11 49; 12 9; 12 24; 12 30; 12 39; 12 47; 13 11; 13 12; 13 24; 13 32; 13 41; 14 2; 14 7; 14 28; 14 29; 14 41; 14 42; 14 48; 15 4; 15 44; 16 3; 16 8; 16 14; 16 23; 16 28; 16 30; 17 4; 17 35; 18 22; 18 27; 18 35; 18 37; 18 45; 18 47; 18 49; 19 5; 19 34; 19 38; 19 41; 19 44; 20 7; 20 11; 20 12; 20 13; 20 16; 20 17; 20 21; 20 22; 20 23; 20 26; 20 38; 20 39; 20 41; 20 43; 20 44; 21 4; 21 20; 21 27; 22 3; 22 18; 22 19; 22 36; 23 4; 23 15; 23 39; 23 49; 24 8; 24 10; 24 27; 24 36; 25 19; 25 20; 25 32; 25 49; 26 5; 26 6; 26 19; 26 38; 26 40; 27 11; 27 20; 27 30; 27 31; 28 8; 28 23; 29 7; 29 9; 29 18; 29 20; 30 23; 30 27; 30 40; 30 46; 31 20; 31 23; 31 32; 31 47; 32 6; 32 33; 33 8; 33 17; 33 21; 33 23; 33 32; 33 37; 33 38; 33 43; 34 2; 34 3; 34 8; 34 11; 34 23; 34 39; 35 11; 35 16; 35 17; 35 20; 35 39; 35 41; 35 46; 36 7; 36 15; 36 22; 36 42; 36 50; 37 16; 37 17; 37 22; 37 29; 37 31; 37 48; 38 2; 38 9; 38 11; 38 21; 38 27; 38 34; 39 8; 39 31; 39 36; 39 46; 40 11; 40 15; 40 33; 40 36; 40 41; 41 7; 41 22; 41 24; 41 47; 41 48; 41 50; 42 4; 42 22; 42 31; 42 34; 42 43; 42 45; 43 34; 43 39; 43 44; 44 2; 44 36; 45 4; 45 11; 45 23; 45 33; 45 36; 45 43; 46 3; 46 5; 46 6; 46 43; 47 13; 47 25; 47 26; 47 33; 47 40; 48 2; 48 5; 48 24; 48 29; 49 3; 49 6; 49 20]
global d_x = [1.0, 4.0, 10.0, 4.0, 1.0, 2.0, 7.0, 6.0, 6.0, 7.0, 8.0, 1.0, 4.0, 1.0, 1.0, 6.0, 6.0, 7.0, 9.0, 6.0, 1.0, 6.0, 8.0, 2.0, 5.0, 3.0, 3.0, 9.0, 9.0, 7.0, 4.0, 4.0, 2.0, 3.0, 5.0, 9.0, 9.0, 8.0, 2.0, 4.0, 6.0, 1.0, 8.0, 2.0, 1.0, 7.0, 2.0, 4.0, 9.0, 1.0, 5.0, 8.0, 4.0, 4.0, 1.0, 7.0, 8.0, 5.0, 5.0, 2.0, 9.0, 4.0, 1.0, 6.0, 7.0, 5.0, 5.0, 3.0, 5.0, 5.0, 1.0, 6.0, 6.0, 10.0, 1.0, 9.0, 7.0, 2.0, 4.0, 3.0, 4.0, 5.0, 6.0, 10.0, 9.0, 2.0, 4.0, 5.0, 9.0, 6.0, 1.0, 4.0, 2.0, 3.0, 8.0, 1.0, 3.0, 3.0, 2.0, 7.0, 10.0, 7.0, 9.0, 6.0, 6.0, 1.0, 1.0, 9.0, 3.0, 2.0, 6.0, 10.0, 10.0, 10.0, 5.0, 1.0, 10.0, 7.0, 1.0, 5.0, 1.0, 1.0, 7.0, 3.0, 8.0, 9.0, 2.0, 4.0, 3.0, 4.0, 1.0, 9.0, 10.0, 7.0, 10.0, 7.0, 5.0, 1.0, 1.0, 9.0, 7.0, 8.0, 8.0, 6.0, 8.0, 3.0, 2.0, 6.0, 4.0, 2.0, 1.0, 7.0, 6.0, 6.0, 10.0, 2.0, 6.0, 9.0, 4.0, 1.0, 8.0, 10.0, 5.0, 4.0, 8.0, 6.0, 10.0, 6.0, 9.0, 4.0, 4.0, 10.0, 7.0, 1.0, 6.0, 3.0, 5.0, 6.0, 7.0, 8.0, 2.0, 7.0, 10.0, 8.0, 1.0, 6.0, 2.0, 4.0, 7.0, 3.0, 4.0, 7.0, 7.0, 1.0, 10.0, 10.0, 8.0, 3.0, 2.0, 8.0, 8.0, 10.0, 1.0, 7.0, 9.0, 8.0, 1.0, 2.0, 9.0, 1.0, 2.0, 9.0, 6.0, 7.0, 5.0, 2.0, 10.0, 7.0, 1.0, 3.0, 8.0, 9.0, 7.0, 2.0, 3.0, 3.0, 4.0, 9.0, 2.0, 2.0, 4.0, 2.0, 1.0, 4.0, 9.0, 1.0, 4.0, 5.0, 2.0, 10.0, 3.0, 5.0]
global b_x = 5
global d_y = [8.0, 1.0, 8.0, 7.0, 8.0, 3.0, 10.0, 6.0, 9.0, 1.0, 4.0, 1.0, 4.0, 9.0, 9.0, 8.0, 7.0, 5.0, 3.0, 7.0, 5.0, 10.0, 3.0, 5.0, 7.0, 5.0, 10.0, 8.0, 2.0, 5.0, 1.0, 5.0, 3.0, 6.0, 1.0, 1.0, 2.0, 5.0, 2.0, 9.0, 6.0, 9.0, 8.0, 10.0, 1.0, 2.0, 1.0, 2.0, 10.0, 10.0, 2.0, 8.0, 7.0, 6.0, 4.0, 5.0, 6.0, 7.0, 3.0, 2.0, 4.0, 8.0, 9.0, 4.0, 1.0, 7.0, 1.0, 5.0, 2.0, 6.0, 6.0, 4.0, 7.0, 7.0, 7.0, 10.0, 3.0, 9.0, 1.0, 4.0, 5.0, 5.0, 6.0, 6.0, 2.0, 3.0, 3.0, 3.0, 7.0, 6.0, 5.0, 9.0, 4.0, 3.0, 4.0, 5.0, 5.0, 10.0, 2.0, 2.0, 7.0, 9.0, 2.0, 8.0, 4.0, 1.0, 6.0, 10.0, 5.0, 9.0, 4.0, 7.0, 3.0, 9.0, 5.0, 4.0, 1.0, 2.0, 10.0, 8.0, 9.0, 3.0, 2.0, 3.0, 4.0, 8.0, 8.0, 4.0, 1.0, 10.0, 5.0, 1.0, 8.0, 3.0, 1.0, 7.0, 4.0, 1.0, 10.0, 8.0, 6.0, 9.0, 4.0, 8.0, 6.0, 7.0, 7.0, 7.0, 7.0, 5.0, 4.0, 6.0, 6.0, 3.0, 4.0, 1.0, 4.0, 2.0, 6.0, 6.0, 6.0, 5.0, 8.0, 8.0, 6.0, 2.0, 10.0, 7.0, 4.0, 8.0, 6.0, 8.0, 3.0, 4.0, 1.0, 6.0, 8.0, 2.0, 6.0, 1.0, 5.0, 3.0, 8.0, 3.0, 10.0, 1.0, 2.0, 1.0, 6.0, 2.0, 5.0, 9.0, 1.0, 10.0, 9.0, 3.0, 8.0, 1.0, 2.0, 3.0, 6.0, 1.0, 3.0, 10.0, 10.0, 6.0, 3.0, 6.0, 3.0, 2.0, 1.0, 3.0, 2.0, 2.0, 3.0, 4.0, 6.0, 10.0, 10.0, 2.0, 7.0, 8.0, 6.0, 6.0, 3.0, 3.0, 7.0, 1.0, 9.0, 10.0, 5.0, 1.0, 2.0, 7.0, 9.0, 6.0, 8.0, 4.0, 7.0, 4.0, 4.0, 5.0]
global b_y = 10
global p = [0.785, 0.692, 0.81, 0.198, 0.519, 0.959, 0.482, 0.071, 0.215, 0.741, 0.637, 0.055, 0.066, 0.101, 0.019, 0.47, 0.546, 0.696, 0.568, 0.16, 0.268, 0.755, 0.897, 0.462, 0.478, 0.055, 0.84, 0.473, 0.085, 0.96, 0.993, 0.304, 0.972, 0.048, 0.855, 0.539, 0.668, 0.158, 0.505, 0.711, 0.084, 0.694, 0.577, 0.906, 0.863, 0.704, 0.732, 0.145, 0.927, 0.566, 0.552, 0.582, 0.415, 0.842, 0.447, 0.578, 0.809, 0.805, 0.647, 0.959, 0.418, 0.733, 0.89, 0.311, 0.73, 0.972, 0.95, 0.968, 0.816, 0.717, 0.349, 0.976, 0.286, 0.98, 0.005, 0.005, 0.412, 0.793, 0.617, 0.261, 0.209, 0.846, 0.458, 0.022, 0.71, 0.825, 0.784, 0.654, 0.475, 0.838, 0.933, 0.158, 0.075, 0.266, 0.669, 0.355, 0.876, 0.039, 0.456, 0.221, 0.447, 0.058, 0.852, 0.577, 0.299, 0.238, 0.438, 0.181, 0.437, 0.783, 0.046, 0.327, 0.548, 0.631, 0.757, 0.972, 0.852, 0.085, 0.992, 0.184, 0.308, 0.255, 0.533, 0.816, 0.41, 0.548, 0.968, 0.191, 0.925, 0.514, 0.64, 0.7, 0.446, 0.045, 0.742, 0.909, 0.23, 0.739, 0.103, 0.234, 0.121, 0.354, 0.879, 0.514, 0.43, 0.197, 0.871, 0.581, 0.952, 0.417, 0.967, 0.625, 0.85, 0.339, 0.482, 0.359, 0.553, 0.754, 0.741, 0.55, 0.943, 0.641, 0.908, 0.808, 0.209, 0.467, 0.753, 0.139, 0.459, 0.908, 0.794, 0.765, 0.796, 0.466, 0.99, 0.799, 0.5, 0.776, 0.331, 0.631, 0.171, 0.478, 0.828, 0.648, 0.939, 0.607, 0.953, 0.341, 0.5, 0.524, 0.449, 0.266, 0.797, 0.112, 0.084, 0.293, 0.741, 0.32, 0.813, 0.961, 0.293, 0.499, 0.938, 0.481, 0.986, 0.289, 0.179, 0.919, 0.034, 0.223, 0.636, 0.446, 0.065, 0.503, 0.394, 0.752, 0.93, 0.999, 0.116, 0.09, 0.427, 0.505, 0.5, 0.284, 0.014, 0.651, 0.749, 0.422, 0.319, 0.87, 0.504, 0.373, 0.352, 0.759, 0.378, 0.481, 0.848, 0.711, 0.056, 0.842, 0.367, 0.301]
global q = [0.839, 0.988, 0.901, 0.733, 0.932, 0.979, 0.529, 0.752, 0.575, 0.772, 0.642, 0.147, 0.787, 0.82, 0.623, 0.822, 0.925, 0.786, 0.672, 0.268, 0.661, 0.853, 0.949, 0.929, 0.743, 0.086, 0.901, 0.833, 0.25, 0.971, 0.995, 0.463, 0.981, 0.59, 0.882, 0.731, 0.816, 0.528, 0.543, 0.731, 0.427, 0.816, 0.948, 0.992, 0.932, 0.943, 0.756, 0.902, 0.978, 0.624, 0.586, 0.959, 0.976, 0.995, 0.457, 0.76, 0.817, 0.997, 0.849, 0.996, 0.77, 0.926, 0.952, 0.667, 0.841, 0.997, 0.994, 0.995, 0.839, 0.974, 0.719, 0.984, 0.753, 0.99, 0.703, 0.985, 0.72, 0.88, 0.88, 0.357, 0.872, 0.923, 0.671, 0.501, 0.733, 0.84, 0.895, 0.842, 0.477, 0.866, 0.946, 0.413, 0.526, 0.382, 0.95, 0.571, 0.969, 0.478, 0.756, 0.226, 0.733, 0.228, 0.885, 0.647, 0.404, 0.496, 0.939, 0.334, 0.49, 0.996, 0.275, 0.482, 0.904, 0.934, 0.833, 0.978, 0.856, 0.993, 0.993, 0.748, 0.652, 0.825, 0.595, 0.96, 0.57, 0.929, 0.987, 0.507, 0.933, 0.902, 0.993, 0.825, 0.691, 0.913, 0.743, 0.915, 0.389, 0.906, 0.254, 0.287, 0.572, 0.538, 0.932, 0.953, 0.721, 0.905, 0.933, 0.778, 0.97, 0.618, 0.973, 0.67, 0.999, 0.598, 0.733, 0.561, 0.896, 0.981, 0.884, 0.715, 0.952, 0.715, 0.928, 0.833, 0.315, 0.494, 0.826, 0.496, 0.866, 0.997, 0.927, 0.822, 0.982, 0.78, 0.996, 0.828, 0.609, 0.836, 0.873, 0.722, 0.199, 0.799, 0.924, 0.836, 0.965, 0.977, 0.986, 0.532, 0.549, 0.649, 0.849, 0.54, 0.985, 0.223, 0.872, 0.651, 0.934, 0.817, 0.934, 0.998, 0.45, 0.564, 0.98, 0.947, 0.986, 0.383, 0.93, 0.998, 0.945, 0.635, 0.805, 0.485, 0.737, 0.906, 0.991, 0.906, 0.999, 0.999, 0.437, 0.639, 0.697, 0.811, 0.663, 0.348, 0.994, 0.7, 0.972, 0.925, 0.607, 0.95, 0.942, 0.46, 0.901, 0.969, 0.529, 0.769, 0.85, 0.945, 0.614, 0.971, 0.804, 0.617]
global origin = 1
global destination = 50