global arcs = [1 6; 1 16; 1 18; 1 19; 1 39; 1 48; 2 9; 2 36; 2 38; 3 50; 4 3; 4 25; 4 37; 4 40; 5 2; 5 3; 5 4; 5 20; 5 26; 5 36; 5 45; 6 16; 6 20; 6 28; 6 32; 6 34; 6 50; 7 5; 7 6; 7 42; 8 20; 8 29; 8 37; 9 6; 9 21; 9 24; 9 42; 9 46; 9 49; 10 5; 10 12; 10 21; 10 38; 10 39; 11 28; 11 37; 11 50; 12 17; 12 26; 12 27; 12 36; 12 43; 13 8; 13 32; 13 41; 13 50; 14 4; 14 5; 14 9; 14 37; 14 38; 14 39; 14 40; 15 41; 15 46; 16 17; 16 30; 16 41; 16 50; 17 2; 17 12; 17 16; 17 26; 17 38; 18 7; 18 10; 18 16; 18 25; 18 31; 18 33; 19 5; 19 31; 19 36; 19 42; 19 47; 19 50; 20 14; 20 17; 20 21; 20 31; 20 35; 20 39; 20 45; 20 48; 21 10; 21 13; 21 17; 21 19; 21 31; 21 39; 21 46; 21 48; 22 16; 22 21; 22 35; 23 3; 23 7; 23 9; 23 20; 23 31; 23 42; 23 46; 24 3; 24 8; 24 13; 24 14; 24 18; 24 21; 24 25; 24 28; 25 26; 25 42; 26 2; 26 7; 26 11; 26 19; 26 22; 26 28; 26 35; 26 36; 26 40; 26 47; 27 4; 27 22; 28 8; 28 21; 28 25; 28 31; 29 15; 29 23; 29 40; 30 4; 30 8; 30 14; 30 28; 30 37; 30 45; 31 2; 31 16; 31 23; 31 34; 32 5; 32 11; 32 13; 32 16; 32 26; 32 36; 32 40; 33 9; 33 16; 33 35; 33 37; 33 45; 34 19; 34 24; 34 43; 35 12; 35 18; 35 24; 35 40; 36 7; 36 8; 36 9; 36 23; 36 38; 36 50; 37 10; 37 22; 37 45; 37 46; 37 50; 38 10; 38 14; 38 21; 38 28; 38 36; 38 48; 39 17; 39 43; 40 6; 40 11; 40 12; 40 13; 40 42; 40 44; 40 48; 41 9; 41 10; 41 17; 41 29; 41 36; 41 46; 42 7; 42 26; 42 28; 42 50; 43 9; 43 12; 43 13; 43 19; 43 21; 43 27; 43 32; 43 35; 43 38; 43 39; 43 40; 43 49; 44 20; 44 31; 44 35; 44 41; 45 11; 45 13; 45 18; 45 23; 45 41; 45 47; 46 29; 46 37; 46 38; 46 44; 46 47; 46 48; 46 49; 47 3; 47 20; 47 25; 47 28; 47 45; 48 12; 48 16; 48 22; 48 25; 48 29; 48 32; 48 43; 49 16; 49 24; 49 38; 49 46]
global d_x = [4.0, 2.0, 8.0, 4.0, 8.0, 3.0, 7.0, 4.0, 1.0, 9.0, 10.0, 5.0, 7.0, 9.0, 9.0, 6.0, 6.0, 7.0, 6.0, 2.0, 10.0, 6.0, 1.0, 3.0, 10.0, 7.0, 6.0, 7.0, 3.0, 10.0, 2.0, 2.0, 8.0, 8.0, 8.0, 10.0, 9.0, 8.0, 9.0, 4.0, 7.0, 10.0, 10.0, 4.0, 5.0, 9.0, 2.0, 2.0, 7.0, 9.0, 6.0, 9.0, 8.0, 4.0, 4.0, 2.0, 5.0, 5.0, 7.0, 6.0, 10.0, 2.0, 1.0, 10.0, 10.0, 8.0, 5.0, 6.0, 4.0, 1.0, 7.0, 1.0, 3.0, 2.0, 2.0, 1.0, 4.0, 9.0, 7.0, 10.0, 7.0, 8.0, 8.0, 9.0, 1.0, 10.0, 8.0, 1.0, 2.0, 6.0, 2.0, 1.0, 8.0, 9.0, 5.0, 2.0, 9.0, 10.0, 3.0, 8.0, 3.0, 10.0, 9.0, 7.0, 1.0, 8.0, 1.0, 10.0, 4.0, 6.0, 6.0, 3.0, 7.0, 10.0, 7.0, 4.0, 7.0, 6.0, 3.0, 6.0, 5.0, 4.0, 9.0, 6.0, 6.0, 2.0, 3.0, 7.0, 9.0, 8.0, 9.0, 5.0, 8.0, 7.0, 5.0, 10.0, 7.0, 10.0, 8.0, 10.0, 1.0, 7.0, 8.0, 4.0, 8.0, 1.0, 1.0, 5.0, 3.0, 4.0, 6.0, 9.0, 4.0, 9.0, 5.0, 5.0, 5.0, 8.0, 5.0, 6.0, 2.0, 8.0, 3.0, 8.0, 5.0, 5.0, 7.0, 10.0, 4.0, 5.0, 10.0, 4.0, 7.0, 7.0, 4.0, 10.0, 4.0, 4.0, 4.0, 5.0, 9.0, 7.0, 1.0, 6.0, 3.0, 6.0, 7.0, 5.0, 8.0, 4.0, 6.0, 2.0, 10.0, 3.0, 7.0, 1.0, 3.0, 2.0, 6.0, 10.0, 1.0, 7.0, 4.0, 1.0, 10.0, 8.0, 6.0, 1.0, 2.0, 3.0, 4.0, 10.0, 8.0, 6.0, 5.0, 1.0, 5.0, 1.0, 2.0, 8.0, 5.0, 1.0, 8.0, 10.0, 5.0, 7.0, 5.0, 8.0, 8.0, 1.0, 5.0, 2.0, 4.0, 4.0, 6.0, 5.0, 6.0, 6.0, 5.0, 4.0, 3.0, 9.0, 5.0, 10.0, 1.0, 8.0, 9.0, 1.0, 6.0, 3.0, 10.0]
global b_x = 5
global d_y = [2.0, 9.0, 5.0, 6.0, 7.0, 10.0, 9.0, 1.0, 6.0, 9.0, 6.0, 4.0, 8.0, 1.0, 1.0, 10.0, 6.0, 9.0, 6.0, 4.0, 1.0, 10.0, 10.0, 6.0, 8.0, 8.0, 8.0, 7.0, 1.0, 6.0, 5.0, 7.0, 3.0, 2.0, 2.0, 6.0, 8.0, 6.0, 8.0, 9.0, 2.0, 3.0, 7.0, 3.0, 3.0, 5.0, 10.0, 2.0, 7.0, 8.0, 10.0, 3.0, 2.0, 10.0, 8.0, 10.0, 10.0, 6.0, 9.0, 10.0, 4.0, 10.0, 7.0, 6.0, 1.0, 3.0, 8.0, 10.0, 9.0, 9.0, 3.0, 7.0, 1.0, 5.0, 9.0, 4.0, 6.0, 6.0, 10.0, 4.0, 4.0, 1.0, 7.0, 4.0, 10.0, 5.0, 1.0, 4.0, 5.0, 10.0, 2.0, 2.0, 7.0, 1.0, 4.0, 9.0, 3.0, 4.0, 9.0, 9.0, 8.0, 7.0, 9.0, 8.0, 1.0, 4.0, 1.0, 9.0, 9.0, 9.0, 10.0, 3.0, 10.0, 3.0, 10.0, 10.0, 9.0, 9.0, 7.0, 5.0, 10.0, 10.0, 7.0, 2.0, 4.0, 6.0, 7.0, 5.0, 2.0, 6.0, 2.0, 4.0, 8.0, 3.0, 8.0, 8.0, 4.0, 4.0, 9.0, 2.0, 7.0, 8.0, 7.0, 7.0, 7.0, 9.0, 6.0, 3.0, 5.0, 8.0, 7.0, 2.0, 9.0, 1.0, 7.0, 3.0, 3.0, 7.0, 7.0, 10.0, 7.0, 8.0, 1.0, 9.0, 2.0, 8.0, 6.0, 7.0, 4.0, 2.0, 10.0, 1.0, 8.0, 2.0, 9.0, 4.0, 5.0, 1.0, 3.0, 6.0, 5.0, 5.0, 5.0, 10.0, 3.0, 10.0, 6.0, 3.0, 7.0, 5.0, 6.0, 10.0, 10.0, 7.0, 10.0, 1.0, 2.0, 3.0, 8.0, 3.0, 4.0, 3.0, 3.0, 8.0, 8.0, 7.0, 10.0, 3.0, 9.0, 4.0, 9.0, 10.0, 7.0, 6.0, 3.0, 8.0, 8.0, 1.0, 3.0, 1.0, 7.0, 6.0, 5.0, 6.0, 10.0, 9.0, 9.0, 2.0, 4.0, 3.0, 5.0, 5.0, 9.0, 1.0, 10.0, 5.0, 10.0, 3.0, 6.0, 1.0, 5.0, 10.0, 5.0, 9.0, 10.0, 2.0, 4.0, 2.0, 6.0, 2.0, 6.0]
global b_y = 10
global p = [0.718, 0.539, 0.292, 0.039, 0.265, 0.084, 0.363, 0.199, 0.139, 0.557, 0.533, 0.854, 0.507, 0.566, 0.21, 0.516, 0.92, 0.399, 0.719, 0.353, 0.939, 0.354, 0.666, 0.175, 0.116, 0.08, 0.846, 0.058, 0.04, 0.438, 0.397, 0.079, 0.976, 0.812, 0.144, 0.069, 0.441, 0.292, 0.395, 0.194, 0.654, 0.583, 0.205, 0.632, 0.341, 0.825, 0.09, 0.593, 0.829, 0.241, 0.259, 0.616, 0.813, 0.011, 0.379, 0.369, 0.591, 0.125, 0.809, 0.833, 0.756, 0.911, 0.533, 0.512, 0.86, 0.76, 0.259, 0.253, 0.226, 0.821, 0.982, 0.421, 0.577, 0.373, 0.313, 0.93, 0.679, 0.6, 0.68, 0.694, 0.183, 0.374, 0.031, 0.262, 0.31, 0.406, 0.897, 0.55, 0.336, 0.111, 0.778, 0.933, 0.548, 0.688, 0.529, 0.158, 0.385, 0.517, 0.125, 0.874, 0.239, 0.293, 0.619, 0.852, 0.704, 0.454, 0.301, 0.001, 0.078, 0.918, 0.783, 0.182, 0.887, 0.658, 0.842, 0.055, 0.74, 0.008, 0.927, 0.654, 0.696, 0.051, 0.717, 0.335, 0.069, 0.147, 0.152, 0.019, 0.077, 0.354, 0.273, 0.938, 0.756, 0.288, 0.11, 0.403, 0.51, 0.715, 0.506, 0.505, 0.568, 0.592, 0.924, 0.957, 0.537, 0.837, 0.042, 0.753, 0.646, 0.252, 0.248, 0.387, 0.204, 0.659, 0.99, 0.738, 0.683, 0.845, 0.773, 0.84, 0.153, 0.604, 0.194, 0.699, 0.618, 0.045, 0.829, 0.954, 0.62, 0.613, 0.466, 0.795, 0.725, 0.278, 0.352, 0.292, 0.432, 0.159, 0.376, 0.038, 0.877, 0.227, 0.099, 0.191, 0.149, 0.705, 0.17, 0.095, 0.834, 0.941, 0.446, 0.094, 0.936, 0.132, 0.804, 0.702, 0.356, 0.71, 0.325, 0.398, 0.201, 0.754, 0.654, 0.302, 0.859, 0.313, 0.704, 0.34, 0.872, 0.908, 0.584, 0.064, 0.611, 0.939, 0.093, 0.164, 0.297, 0.269, 0.069, 0.7, 0.487, 0.956, 0.526, 0.53, 0.941, 0.105, 0.004, 0.008, 0.976, 0.356, 0.955, 0.812, 0.471, 0.92, 0.398, 0.082, 0.722, 0.654, 0.114, 0.426, 0.134, 0.315, 0.245, 0.467, 0.53, 0.989, 0.911, 0.923, 0.624, 0.438, 0.827]
global q = [0.776, 0.671, 0.461, 0.684, 0.274, 0.381, 0.971, 0.659, 0.475, 0.622, 0.68, 0.985, 0.86, 0.933, 0.932, 0.582, 0.94, 0.628, 0.728, 0.9, 0.951, 0.861, 0.91, 0.815, 0.271, 0.77, 0.972, 0.416, 0.076, 0.961, 0.477, 0.512, 0.999, 0.943, 0.644, 0.077, 0.691, 0.709, 0.576, 0.668, 0.73, 0.953, 0.367, 0.919, 0.363, 0.91, 0.325, 0.874, 0.898, 0.746, 0.447, 0.616, 0.968, 0.932, 0.436, 0.706, 0.718, 0.949, 0.948, 0.833, 0.849, 0.96, 0.689, 0.739, 0.966, 0.891, 0.562, 0.7, 0.532, 0.853, 0.993, 0.441, 0.901, 0.418, 0.695, 0.98, 0.705, 0.735, 0.935, 0.834, 0.57, 0.569, 0.253, 0.909, 0.528, 0.633, 0.993, 0.708, 0.931, 0.481, 0.836, 0.943, 0.667, 0.785, 0.563, 0.579, 0.489, 0.876, 0.517, 0.991, 0.706, 0.611, 0.836, 0.916, 0.767, 0.873, 0.701, 0.065, 0.184, 0.924, 0.912, 0.299, 0.968, 0.946, 0.94, 0.39, 0.811, 0.7, 0.942, 0.76, 0.85, 0.053, 0.944, 0.452, 0.914, 0.815, 0.278, 0.156, 0.182, 0.802, 0.546, 0.954, 0.994, 0.288, 0.638, 0.999, 0.747, 0.751, 0.607, 0.881, 0.596, 0.65, 0.988, 0.991, 0.597, 0.929, 0.825, 0.852, 0.665, 0.499, 0.849, 0.813, 0.953, 0.777, 0.996, 0.832, 0.836, 0.944, 0.938, 0.996, 0.295, 0.872, 0.698, 0.88, 0.931, 0.864, 0.942, 0.98, 0.941, 0.696, 0.839, 0.989, 0.804, 0.96, 0.553, 0.971, 0.507, 0.516, 0.474, 0.261, 0.902, 0.891, 0.645, 0.736, 0.958, 0.735, 0.857, 0.214, 0.961, 0.985, 0.88, 0.562, 0.952, 0.953, 0.896, 0.851, 0.379, 0.899, 0.47, 0.71, 0.444, 0.862, 0.857, 0.822, 0.893, 0.877, 0.758, 0.387, 0.985, 0.938, 0.741, 0.97, 0.868, 0.98, 0.199, 0.733, 0.337, 0.733, 0.751, 0.713, 0.855, 0.971, 0.581, 0.703, 0.999, 0.756, 0.033, 0.023, 0.987, 0.788, 0.968, 0.946, 0.56, 0.943, 0.583, 0.424, 0.947, 0.887, 0.721, 0.984, 0.993, 0.704, 0.544, 0.931, 0.643, 0.99, 0.996, 0.98, 0.663, 0.972, 0.93]
global origin = 1
global destination = 50