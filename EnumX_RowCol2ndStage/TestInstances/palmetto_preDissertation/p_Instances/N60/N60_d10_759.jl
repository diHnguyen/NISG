global arcs = [1 6; 1 7; 1 10; 1 12; 1 24; 1 57; 2 4; 2 5; 2 24; 2 34; 2 48; 3 4; 3 30; 3 38; 3 41; 3 47; 3 58; 4 12; 4 16; 4 24; 4 25; 4 50; 4 53; 5 8; 5 12; 5 18; 5 20; 5 44; 5 45; 5 54; 6 3; 6 21; 6 36; 7 3; 7 40; 7 56; 7 60; 8 34; 8 49; 9 7; 9 8; 9 21; 9 27; 9 57; 10 11; 10 15; 10 23; 10 31; 10 41; 11 8; 11 16; 11 22; 11 26; 11 30; 11 38; 11 39; 12 17; 12 20; 12 39; 12 43; 12 54; 12 56; 13 23; 13 25; 13 27; 13 33; 13 42; 13 45; 13 51; 13 53; 14 27; 14 35; 14 54; 14 56; 14 58; 15 3; 15 18; 15 26; 15 28; 15 32; 15 47; 16 15; 16 41; 16 42; 16 52; 17 2; 17 5; 17 9; 17 11; 17 19; 17 37; 17 44; 17 48; 18 17; 18 19; 18 35; 18 37; 18 46; 18 51; 18 56; 19 7; 19 11; 19 24; 19 31; 19 33; 19 34; 19 48; 19 49; 20 9; 20 13; 20 14; 20 34; 20 59; 20 60; 21 20; 21 28; 21 38; 21 47; 21 56; 21 57; 22 5; 22 25; 22 31; 22 42; 22 48; 22 51; 22 55; 22 57; 23 2; 23 15; 23 21; 23 24; 23 35; 23 36; 23 52; 24 14; 24 30; 24 35; 24 38; 24 54; 25 24; 25 35; 25 53; 26 10; 26 28; 26 33; 26 40; 26 42; 26 43; 26 49; 26 55; 26 59; 27 2; 27 33; 27 44; 27 45; 27 52; 28 3; 28 14; 28 27; 28 30; 28 31; 28 45; 28 47; 29 9; 29 18; 29 33; 29 51; 29 53; 30 60; 31 20; 31 29; 31 41; 31 45; 31 47; 31 48; 32 27; 32 54; 32 59; 33 3; 33 42; 33 47; 34 6; 34 30; 34 55; 35 7; 35 11; 35 17; 35 23; 35 47; 36 8; 36 11; 36 15; 36 24; 36 37; 36 38; 36 58; 37 5; 37 7; 37 12; 37 20; 37 34; 37 36; 37 47; 38 6; 38 21; 38 41; 38 43; 38 50; 39 7; 39 28; 39 32; 39 42; 39 58; 40 14; 40 22; 40 23; 40 32; 40 53; 40 55; 41 9; 41 48; 41 56; 42 25; 42 31; 42 36; 42 60; 43 5; 43 16; 43 25; 43 27; 43 30; 43 41; 43 46; 43 48; 43 49; 43 51; 43 55; 43 58; 44 2; 44 5; 44 6; 44 14; 44 17; 44 37; 45 12; 45 16; 45 25; 45 31; 45 40; 46 6; 46 16; 46 22; 46 28; 46 30; 46 32; 46 44; 46 50; 47 16; 47 29; 47 31; 47 35; 47 49; 48 15; 48 17; 48 22; 48 28; 48 42; 48 54; 48 57; 49 8; 49 12; 49 24; 49 40; 50 8; 50 9; 50 17; 50 20; 50 23; 50 30; 50 46; 50 52; 50 54; 51 20; 51 34; 51 37; 51 42; 51 44; 51 53; 52 7; 52 11; 52 37; 52 45; 52 47; 52 59; 53 6; 53 33; 54 16; 54 23; 54 52; 54 53; 54 59; 55 19; 55 28; 55 41; 56 2; 56 5; 56 23; 56 26; 56 30; 56 52; 57 32; 57 35; 57 54; 58 31; 58 44; 58 46; 58 53; 59 14; 59 25; 59 53; 59 57]
global d_x = [8.0, 7.0, 2.0, 10.0, 8.0, 5.0, 6.0, 7.0, 10.0, 8.0, 7.0, 6.0, 9.0, 8.0, 1.0, 8.0, 8.0, 9.0, 4.0, 8.0, 9.0, 4.0, 2.0, 10.0, 3.0, 3.0, 1.0, 10.0, 9.0, 2.0, 3.0, 4.0, 10.0, 4.0, 4.0, 10.0, 9.0, 6.0, 3.0, 5.0, 5.0, 9.0, 5.0, 8.0, 3.0, 1.0, 7.0, 4.0, 4.0, 1.0, 1.0, 5.0, 4.0, 4.0, 10.0, 6.0, 9.0, 8.0, 7.0, 10.0, 6.0, 7.0, 1.0, 8.0, 6.0, 2.0, 4.0, 9.0, 1.0, 6.0, 8.0, 3.0, 7.0, 8.0, 7.0, 5.0, 1.0, 6.0, 6.0, 4.0, 5.0, 6.0, 1.0, 10.0, 6.0, 7.0, 6.0, 6.0, 10.0, 8.0, 9.0, 7.0, 3.0, 5.0, 1.0, 5.0, 7.0, 2.0, 5.0, 2.0, 5.0, 1.0, 4.0, 6.0, 5.0, 8.0, 2.0, 5.0, 8.0, 9.0, 8.0, 1.0, 5.0, 10.0, 9.0, 6.0, 4.0, 8.0, 8.0, 7.0, 7.0, 1.0, 8.0, 2.0, 6.0, 8.0, 10.0, 7.0, 5.0, 9.0, 10.0, 5.0, 5.0, 4.0, 5.0, 10.0, 10.0, 2.0, 8.0, 1.0, 9.0, 2.0, 10.0, 6.0, 8.0, 5.0, 6.0, 3.0, 9.0, 8.0, 3.0, 7.0, 6.0, 1.0, 3.0, 6.0, 3.0, 10.0, 4.0, 2.0, 10.0, 3.0, 6.0, 8.0, 1.0, 1.0, 5.0, 10.0, 8.0, 8.0, 4.0, 9.0, 4.0, 9.0, 7.0, 8.0, 2.0, 10.0, 4.0, 1.0, 5.0, 10.0, 10.0, 9.0, 7.0, 2.0, 6.0, 8.0, 4.0, 8.0, 6.0, 8.0, 10.0, 6.0, 2.0, 6.0, 7.0, 4.0, 1.0, 10.0, 3.0, 1.0, 7.0, 4.0, 8.0, 7.0, 4.0, 3.0, 10.0, 2.0, 3.0, 3.0, 6.0, 7.0, 2.0, 6.0, 8.0, 7.0, 10.0, 4.0, 8.0, 5.0, 7.0, 8.0, 5.0, 2.0, 1.0, 5.0, 2.0, 7.0, 2.0, 2.0, 6.0, 1.0, 10.0, 1.0, 8.0, 4.0, 10.0, 4.0, 8.0, 6.0, 4.0, 2.0, 3.0, 3.0, 10.0, 2.0, 2.0, 2.0, 5.0, 4.0, 10.0, 8.0, 7.0, 10.0, 8.0, 3.0, 1.0, 6.0, 4.0, 2.0, 1.0, 5.0, 10.0, 4.0, 5.0, 5.0, 4.0, 5.0, 4.0, 3.0, 3.0, 10.0, 3.0, 3.0, 6.0, 2.0, 9.0, 3.0, 4.0, 7.0, 3.0, 7.0, 7.0, 5.0, 1.0, 10.0, 5.0, 7.0, 7.0, 3.0, 4.0, 7.0, 9.0, 10.0, 5.0, 5.0, 8.0, 7.0, 4.0, 8.0, 7.0, 5.0, 7.0, 4.0, 8.0, 5.0, 3.0, 7.0, 1.0, 2.0, 1.0, 3.0, 8.0, 7.0, 9.0, 10.0, 4.0, 3.0, 8.0, 9.0]
global b_x = 5
global d_y = [8.0, 7.0, 1.0, 5.0, 9.0, 8.0, 2.0, 2.0, 2.0, 1.0, 3.0, 6.0, 2.0, 3.0, 9.0, 7.0, 6.0, 2.0, 6.0, 1.0, 8.0, 6.0, 5.0, 9.0, 9.0, 4.0, 7.0, 5.0, 3.0, 2.0, 6.0, 1.0, 9.0, 5.0, 10.0, 10.0, 4.0, 7.0, 4.0, 4.0, 5.0, 5.0, 3.0, 4.0, 3.0, 3.0, 3.0, 9.0, 8.0, 7.0, 9.0, 1.0, 9.0, 5.0, 3.0, 9.0, 1.0, 4.0, 1.0, 2.0, 5.0, 7.0, 6.0, 3.0, 7.0, 4.0, 4.0, 3.0, 10.0, 4.0, 10.0, 10.0, 1.0, 2.0, 10.0, 7.0, 4.0, 5.0, 7.0, 5.0, 10.0, 8.0, 9.0, 1.0, 7.0, 2.0, 5.0, 9.0, 6.0, 8.0, 1.0, 3.0, 6.0, 10.0, 4.0, 10.0, 3.0, 4.0, 8.0, 9.0, 10.0, 8.0, 6.0, 2.0, 8.0, 5.0, 9.0, 10.0, 5.0, 5.0, 1.0, 7.0, 5.0, 6.0, 1.0, 10.0, 3.0, 2.0, 10.0, 2.0, 5.0, 2.0, 3.0, 1.0, 8.0, 2.0, 6.0, 4.0, 4.0, 4.0, 10.0, 10.0, 10.0, 2.0, 1.0, 5.0, 9.0, 5.0, 1.0, 10.0, 3.0, 4.0, 5.0, 3.0, 7.0, 4.0, 1.0, 10.0, 7.0, 6.0, 8.0, 9.0, 7.0, 8.0, 9.0, 1.0, 4.0, 7.0, 3.0, 2.0, 6.0, 10.0, 10.0, 8.0, 5.0, 8.0, 7.0, 2.0, 3.0, 10.0, 5.0, 10.0, 2.0, 6.0, 2.0, 3.0, 7.0, 9.0, 5.0, 4.0, 4.0, 1.0, 5.0, 9.0, 4.0, 4.0, 10.0, 9.0, 5.0, 10.0, 3.0, 3.0, 6.0, 10.0, 5.0, 2.0, 1.0, 7.0, 4.0, 4.0, 2.0, 1.0, 6.0, 2.0, 4.0, 8.0, 2.0, 9.0, 1.0, 10.0, 2.0, 9.0, 7.0, 8.0, 8.0, 5.0, 2.0, 8.0, 7.0, 5.0, 10.0, 5.0, 4.0, 4.0, 8.0, 9.0, 7.0, 7.0, 2.0, 7.0, 4.0, 2.0, 1.0, 4.0, 7.0, 3.0, 3.0, 9.0, 8.0, 2.0, 1.0, 3.0, 2.0, 5.0, 1.0, 6.0, 9.0, 3.0, 3.0, 4.0, 6.0, 1.0, 5.0, 9.0, 2.0, 4.0, 2.0, 5.0, 4.0, 7.0, 10.0, 5.0, 3.0, 1.0, 5.0, 1.0, 4.0, 5.0, 2.0, 8.0, 6.0, 5.0, 10.0, 7.0, 8.0, 2.0, 3.0, 4.0, 1.0, 1.0, 2.0, 6.0, 4.0, 8.0, 7.0, 2.0, 10.0, 1.0, 5.0, 7.0, 2.0, 10.0, 4.0, 8.0, 5.0, 7.0, 6.0, 4.0, 2.0, 7.0, 10.0, 9.0, 6.0, 5.0, 2.0, 3.0, 3.0, 10.0, 2.0, 10.0, 10.0, 7.0, 10.0, 6.0, 6.0, 4.0, 6.0, 1.0, 10.0, 8.0, 9.0, 2.0]
global b_y = 10
global p = [0.65, 0.604, 0.629, 0.717, 0.642, 0.283, 0.154, 0.743, 0.717, 0.839, 0.575, 0.495, 0.003, 0.147, 0.156, 0.705, 0.555, 0.643, 0.419, 0.972, 0.97, 0.967, 0.14, 0.682, 0.506, 0.613, 0.266, 0.081, 0.513, 0.244, 0.364, 0.201, 0.217, 0.601, 0.158, 0.494, 0.915, 0.227, 0.93, 0.554, 0.895, 0.913, 0.729, 0.871, 0.36, 0.633, 0.288, 0.353, 0.581, 0.969, 0.89, 0.19, 0.311, 0.621, 0.807, 0.237, 0.142, 0.024, 0.381, 0.203, 0.183, 0.489, 0.318, 0.079, 0.542, 0.802, 0.85, 0.755, 0.087, 0.327, 0.172, 0.378, 0.934, 0.452, 0.477, 0.221, 0.649, 0.449, 0.208, 0.914, 0.802, 0.847, 0.722, 0.507, 0.995, 0.43, 0.279, 0.184, 0.608, 0.95, 0.05, 0.764, 0.785, 0.363, 0.58, 0.337, 0.763, 0.376, 0.203, 0.175, 0.436, 0.96, 0.691, 0.709, 0.065, 0.026, 0.128, 0.481, 0.754, 0.677, 0.084, 0.256, 0.112, 0.944, 0.525, 0.689, 0.62, 0.154, 0.001, 0.737, 0.963, 0.261, 0.575, 0.377, 0.644, 0.29, 0.391, 0.061, 0.477, 0.966, 0.772, 0.576, 0.132, 0.924, 0.288, 0.631, 0.916, 0.213, 0.351, 0.44, 0.671, 0.929, 0.782, 0.316, 0.719, 0.222, 0.509, 0.401, 0.741, 0.899, 0.185, 0.005, 0.318, 0.565, 0.26, 0.904, 0.949, 0.41, 0.237, 0.558, 0.016, 0.158, 0.908, 0.333, 0.893, 0.988, 0.586, 0.964, 0.168, 0.215, 0.693, 0.063, 0.493, 0.964, 0.943, 0.129, 0.499, 0.787, 0.552, 0.345, 0.573, 0.586, 0.983, 0.057, 0.314, 0.595, 0.429, 0.52, 0.266, 0.273, 0.425, 0.587, 0.058, 0.688, 0.692, 0.008, 0.717, 0.182, 0.075, 0.286, 0.838, 0.311, 0.451, 0.617, 0.899, 0.041, 0.777, 0.979, 0.348, 0.498, 0.313, 0.252, 0.321, 0.29, 0.176, 0.181, 0.305, 0.86, 0.048, 0.499, 0.617, 0.928, 0.387, 0.305, 0.594, 0.013, 0.24, 0.281, 0.274, 0.026, 0.271, 0.27, 0.267, 0.797, 0.444, 0.588, 0.318, 0.003, 0.001, 0.412, 0.949, 0.778, 0.731, 0.087, 0.931, 0.943, 0.155, 0.303, 0.397, 0.013, 0.019, 0.092, 0.792, 0.357, 0.062, 0.168, 0.476, 0.942, 0.663, 0.669, 0.169, 0.628, 0.107, 0.351, 0.68, 0.218, 0.226, 0.141, 0.067, 0.576, 0.795, 0.346, 0.149, 0.757, 0.191, 0.499, 0.169, 0.188, 0.248, 0.401, 0.722, 0.995, 0.836, 0.687, 0.97, 0.92, 0.599, 0.174, 0.839, 0.105, 0.13, 0.95, 0.205, 0.673, 0.958, 0.82, 0.888, 0.647, 0.321, 0.295, 0.661, 0.247, 0.825, 0.12, 0.345, 0.834, 0.835, 0.313, 0.909, 0.519, 0.772, 0.035, 0.97, 0.494, 0.663, 0.673, 0.844, 0.949, 0.085, 0.434, 0.272, 0.906]
global q = [0.786, 0.853, 0.961, 0.799, 0.894, 0.721, 0.597, 0.808, 0.979, 0.899, 0.971, 0.508, 0.288, 0.245, 0.3, 0.978, 0.623, 0.988, 0.971, 0.994, 0.986, 0.972, 0.507, 0.74, 0.767, 0.696, 0.498, 0.803, 0.937, 0.369, 0.703, 0.433, 0.354, 0.786, 0.408, 0.651, 0.933, 0.728, 0.931, 0.623, 0.897, 0.956, 0.969, 0.887, 0.423, 0.819, 0.319, 0.601, 0.777, 0.984, 0.916, 0.621, 0.947, 0.881, 0.965, 0.336, 0.448, 0.541, 0.876, 0.985, 0.652, 0.968, 0.395, 0.727, 0.929, 0.891, 0.938, 0.903, 0.183, 0.932, 0.601, 0.5, 0.972, 0.77, 0.57, 0.902, 0.753, 0.609, 0.608, 0.977, 0.905, 0.868, 0.916, 0.715, 0.999, 0.775, 0.564, 0.814, 0.829, 0.992, 0.658, 0.859, 0.804, 0.798, 0.661, 0.536, 0.863, 0.776, 0.334, 0.649, 0.447, 0.968, 0.977, 0.891, 0.169, 0.212, 0.742, 0.731, 0.902, 0.902, 0.541, 0.772, 0.173, 0.956, 0.943, 0.823, 0.868, 0.858, 0.842, 0.944, 0.983, 0.696, 0.732, 0.666, 0.865, 0.621, 0.485, 0.145, 0.711, 0.994, 0.898, 0.595, 0.353, 0.938, 0.663, 0.866, 0.964, 0.263, 0.53, 0.567, 0.803, 0.93, 0.906, 0.482, 0.883, 0.935, 0.805, 0.431, 0.761, 0.949, 0.526, 0.838, 0.983, 0.989, 0.555, 0.979, 0.951, 0.872, 0.974, 0.992, 0.955, 0.22, 0.933, 0.978, 0.947, 0.998, 0.897, 0.977, 0.39, 0.826, 0.965, 0.934, 0.763, 0.978, 0.987, 0.273, 0.964, 0.855, 0.914, 0.386, 0.864, 0.916, 0.984, 0.105, 0.732, 0.698, 0.714, 0.94, 0.474, 0.307, 0.829, 0.931, 0.404, 0.699, 0.944, 0.677, 0.855, 0.324, 0.546, 0.342, 0.99, 0.356, 0.862, 0.971, 0.988, 0.086, 0.911, 0.982, 0.659, 0.589, 0.647, 0.438, 0.344, 0.719, 0.509, 0.222, 0.996, 0.92, 0.833, 0.991, 0.631, 0.955, 0.98, 0.578, 0.92, 0.681, 0.386, 0.511, 0.295, 0.345, 0.506, 0.742, 0.613, 0.929, 0.57, 0.924, 0.667, 0.966, 0.45, 0.58, 0.978, 0.986, 0.77, 0.196, 0.987, 0.943, 0.233, 0.497, 0.473, 0.482, 0.215, 0.416, 0.834, 0.854, 0.827, 0.666, 0.538, 0.998, 0.763, 0.803, 0.831, 0.757, 0.178, 0.946, 0.725, 0.829, 0.96, 0.619, 0.298, 0.89, 0.979, 0.721, 0.701, 0.894, 0.957, 0.845, 0.488, 0.27, 0.466, 0.819, 0.737, 0.998, 0.857, 0.978, 0.982, 0.922, 0.876, 0.626, 0.85, 0.68, 0.629, 0.989, 0.832, 0.945, 0.966, 0.883, 0.894, 0.888, 0.406, 0.362, 0.777, 0.643, 0.909, 0.236, 0.723, 0.996, 0.999, 0.458, 0.926, 0.636, 0.898, 0.194, 0.983, 0.958, 0.715, 0.894, 0.975, 0.95, 0.692, 0.536, 0.879, 0.959]
global origin = 1
global destination = 60