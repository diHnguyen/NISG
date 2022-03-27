global arcs = [1 15; 1 34; 1 36; 2 26; 2 38; 2 45; 2 48; 2 56; 2 59; 3 2; 3 11; 3 19; 3 29; 3 31; 3 43; 4 9; 4 30; 4 31; 4 32; 4 40; 4 45; 4 48; 5 3; 5 6; 5 9; 5 24; 5 40; 5 51; 6 20; 6 42; 6 53; 7 12; 7 16; 7 50; 8 7; 8 10; 8 16; 8 23; 8 29; 8 35; 8 55; 8 58; 9 4; 9 17; 9 31; 9 51; 10 5; 10 28; 10 29; 10 43; 10 45; 10 51; 11 10; 11 17; 11 24; 11 29; 11 33; 11 44; 11 52; 11 59; 12 20; 12 22; 12 28; 12 30; 12 31; 12 42; 12 43; 12 45; 13 24; 13 48; 14 9; 14 25; 14 31; 15 6; 15 55; 16 4; 16 13; 16 28; 16 32; 16 33; 16 43; 16 56; 17 4; 17 14; 17 26; 17 41; 17 51; 17 53; 17 58; 18 2; 18 8; 18 17; 18 19; 18 24; 18 43; 18 44; 19 8; 19 14; 19 30; 19 37; 19 39; 19 55; 20 2; 20 5; 20 15; 20 21; 20 42; 20 44; 20 58; 21 30; 21 32; 22 9; 22 16; 22 31; 22 35; 22 44; 22 48; 22 51; 23 2; 23 13; 23 16; 23 29; 23 34; 23 57; 24 10; 24 25; 24 26; 24 38; 24 40; 24 49; 25 30; 25 44; 26 5; 26 6; 26 8; 26 12; 26 13; 26 60; 27 8; 27 11; 27 18; 27 36; 27 58; 28 11; 28 17; 28 22; 28 44; 29 12; 29 41; 30 18; 30 32; 30 45; 31 5; 31 12; 31 20; 31 26; 31 28; 31 29; 31 36; 31 40; 31 55; 31 59; 32 9; 32 18; 32 41; 32 46; 32 47; 32 49; 32 52; 33 23; 33 31; 33 51; 33 55; 34 2; 34 9; 34 12; 34 22; 34 36; 34 42; 34 43; 34 57; 34 60; 35 2; 35 6; 35 31; 35 50; 35 52; 36 6; 36 11; 36 22; 36 29; 36 32; 36 35; 36 47; 36 48; 36 52; 36 55; 37 2; 37 13; 37 18; 37 25; 37 47; 38 12; 38 33; 38 37; 38 39; 38 59; 39 9; 39 30; 39 31; 39 34; 39 35; 39 49; 40 7; 40 9; 40 16; 40 35; 40 42; 40 46; 40 50; 40 52; 40 55; 41 11; 41 26; 41 29; 41 32; 41 43; 41 49; 42 6; 42 11; 42 14; 42 21; 42 26; 42 39; 42 40; 42 46; 42 52; 43 5; 43 12; 43 33; 43 37; 43 47; 43 52; 44 3; 44 18; 44 27; 44 31; 44 34; 44 37; 44 39; 44 43; 44 48; 44 49; 45 47; 46 11; 46 14; 46 19; 46 37; 46 41; 46 53; 46 60; 47 6; 47 9; 47 19; 47 26; 48 11; 48 23; 48 33; 48 49; 49 17; 49 19; 49 21; 49 22; 49 23; 49 25; 49 37; 49 47; 49 55; 50 8; 50 17; 50 58; 50 60; 51 47; 52 7; 52 9; 52 20; 52 50; 52 54; 53 3; 53 11; 53 17; 53 47; 53 51; 54 2; 54 27; 54 34; 54 42; 54 56; 55 16; 55 18; 55 50; 55 52; 55 54; 55 57; 55 60; 56 29; 56 35; 56 47; 56 60; 57 2; 57 7; 57 10; 57 23; 57 30; 57 32; 57 50; 58 3; 58 6; 58 14; 58 16; 58 18; 58 20; 58 23; 58 31; 58 42; 59 2; 59 14; 59 23; 59 24; 59 56; 59 58; 59 60]
global d_x = [8.0, 3.0, 8.0, 7.0, 5.0, 5.0, 1.0, 6.0, 5.0, 8.0, 9.0, 4.0, 9.0, 6.0, 10.0, 5.0, 8.0, 2.0, 10.0, 10.0, 4.0, 9.0, 1.0, 4.0, 5.0, 9.0, 5.0, 3.0, 1.0, 1.0, 4.0, 10.0, 4.0, 6.0, 9.0, 8.0, 4.0, 9.0, 1.0, 1.0, 7.0, 8.0, 6.0, 9.0, 7.0, 10.0, 7.0, 9.0, 4.0, 7.0, 1.0, 8.0, 6.0, 10.0, 2.0, 7.0, 5.0, 3.0, 3.0, 4.0, 5.0, 7.0, 5.0, 9.0, 9.0, 3.0, 1.0, 5.0, 1.0, 1.0, 8.0, 7.0, 3.0, 2.0, 1.0, 8.0, 10.0, 2.0, 8.0, 5.0, 3.0, 9.0, 2.0, 9.0, 10.0, 2.0, 2.0, 6.0, 5.0, 4.0, 1.0, 1.0, 4.0, 3.0, 8.0, 6.0, 3.0, 9.0, 4.0, 8.0, 6.0, 3.0, 4.0, 7.0, 1.0, 9.0, 5.0, 5.0, 5.0, 8.0, 5.0, 2.0, 10.0, 4.0, 10.0, 4.0, 2.0, 2.0, 1.0, 10.0, 3.0, 3.0, 6.0, 10.0, 1.0, 1.0, 10.0, 6.0, 3.0, 6.0, 6.0, 9.0, 2.0, 5.0, 2.0, 6.0, 10.0, 5.0, 7.0, 1.0, 4.0, 8.0, 1.0, 4.0, 4.0, 8.0, 6.0, 7.0, 6.0, 4.0, 10.0, 4.0, 7.0, 3.0, 8.0, 5.0, 5.0, 2.0, 10.0, 1.0, 9.0, 8.0, 3.0, 6.0, 9.0, 8.0, 2.0, 2.0, 10.0, 7.0, 9.0, 6.0, 2.0, 5.0, 6.0, 8.0, 7.0, 5.0, 1.0, 1.0, 6.0, 4.0, 2.0, 8.0, 1.0, 1.0, 7.0, 2.0, 10.0, 10.0, 5.0, 1.0, 4.0, 9.0, 4.0, 9.0, 6.0, 7.0, 7.0, 8.0, 10.0, 10.0, 7.0, 9.0, 4.0, 2.0, 9.0, 4.0, 10.0, 8.0, 8.0, 6.0, 10.0, 1.0, 8.0, 3.0, 8.0, 2.0, 5.0, 5.0, 3.0, 10.0, 3.0, 7.0, 4.0, 6.0, 9.0, 4.0, 4.0, 10.0, 4.0, 7.0, 1.0, 7.0, 5.0, 9.0, 8.0, 9.0, 9.0, 7.0, 8.0, 6.0, 7.0, 9.0, 10.0, 1.0, 2.0, 3.0, 4.0, 10.0, 5.0, 6.0, 2.0, 2.0, 4.0, 6.0, 7.0, 1.0, 7.0, 3.0, 7.0, 8.0, 9.0, 7.0, 6.0, 7.0, 4.0, 6.0, 8.0, 3.0, 3.0, 7.0, 6.0, 7.0, 3.0, 7.0, 5.0, 4.0, 2.0, 2.0, 6.0, 10.0, 5.0, 2.0, 3.0, 2.0, 9.0, 9.0, 2.0, 1.0, 6.0, 6.0, 8.0, 5.0, 9.0, 5.0, 2.0, 8.0, 10.0, 9.0, 4.0, 5.0, 8.0, 1.0, 4.0, 4.0, 1.0, 7.0, 6.0, 6.0, 7.0, 8.0, 5.0, 1.0, 4.0, 10.0, 7.0, 2.0, 7.0, 2.0, 8.0, 6.0, 10.0, 3.0, 9.0, 9.0, 8.0, 8.0, 4.0, 8.0, 2.0, 9.0]
global b_x = 5
global d_y = [5.0, 1.0, 8.0, 7.0, 4.0, 4.0, 4.0, 9.0, 6.0, 9.0, 2.0, 7.0, 2.0, 1.0, 7.0, 4.0, 7.0, 5.0, 8.0, 9.0, 10.0, 5.0, 3.0, 3.0, 5.0, 8.0, 4.0, 6.0, 2.0, 4.0, 6.0, 3.0, 7.0, 9.0, 1.0, 2.0, 9.0, 4.0, 8.0, 7.0, 6.0, 7.0, 2.0, 10.0, 10.0, 2.0, 8.0, 10.0, 1.0, 8.0, 8.0, 2.0, 9.0, 4.0, 4.0, 1.0, 10.0, 6.0, 9.0, 2.0, 9.0, 4.0, 8.0, 5.0, 2.0, 3.0, 1.0, 8.0, 10.0, 2.0, 5.0, 3.0, 10.0, 9.0, 8.0, 3.0, 3.0, 5.0, 9.0, 7.0, 4.0, 10.0, 9.0, 6.0, 9.0, 2.0, 2.0, 5.0, 3.0, 3.0, 3.0, 3.0, 7.0, 3.0, 4.0, 9.0, 1.0, 7.0, 8.0, 6.0, 1.0, 9.0, 9.0, 4.0, 2.0, 5.0, 2.0, 2.0, 1.0, 9.0, 10.0, 6.0, 3.0, 1.0, 10.0, 10.0, 7.0, 5.0, 6.0, 1.0, 8.0, 7.0, 10.0, 7.0, 10.0, 6.0, 5.0, 7.0, 3.0, 2.0, 7.0, 2.0, 4.0, 9.0, 10.0, 1.0, 7.0, 2.0, 10.0, 7.0, 4.0, 6.0, 1.0, 8.0, 10.0, 4.0, 4.0, 7.0, 10.0, 9.0, 3.0, 4.0, 10.0, 1.0, 1.0, 9.0, 4.0, 8.0, 9.0, 8.0, 9.0, 7.0, 7.0, 10.0, 1.0, 1.0, 10.0, 3.0, 1.0, 7.0, 10.0, 6.0, 3.0, 5.0, 1.0, 3.0, 4.0, 8.0, 1.0, 8.0, 9.0, 6.0, 8.0, 7.0, 2.0, 10.0, 5.0, 6.0, 1.0, 5.0, 5.0, 7.0, 5.0, 6.0, 10.0, 5.0, 6.0, 7.0, 8.0, 4.0, 3.0, 1.0, 2.0, 7.0, 4.0, 6.0, 7.0, 2.0, 4.0, 4.0, 5.0, 7.0, 5.0, 1.0, 6.0, 2.0, 6.0, 1.0, 6.0, 7.0, 4.0, 2.0, 1.0, 7.0, 6.0, 5.0, 10.0, 2.0, 8.0, 8.0, 9.0, 5.0, 7.0, 4.0, 5.0, 5.0, 5.0, 10.0, 4.0, 8.0, 3.0, 3.0, 9.0, 9.0, 2.0, 8.0, 5.0, 8.0, 1.0, 7.0, 7.0, 5.0, 10.0, 4.0, 1.0, 3.0, 2.0, 1.0, 7.0, 6.0, 2.0, 4.0, 3.0, 5.0, 5.0, 9.0, 3.0, 3.0, 3.0, 2.0, 4.0, 6.0, 1.0, 6.0, 7.0, 2.0, 8.0, 6.0, 5.0, 2.0, 8.0, 6.0, 4.0, 4.0, 4.0, 8.0, 6.0, 1.0, 5.0, 9.0, 5.0, 10.0, 8.0, 2.0, 7.0, 7.0, 8.0, 3.0, 9.0, 6.0, 8.0, 2.0, 9.0, 10.0, 2.0, 3.0, 8.0, 10.0, 2.0, 3.0, 8.0, 3.0, 3.0, 5.0, 2.0, 10.0, 3.0, 7.0, 10.0, 8.0, 9.0, 9.0, 4.0, 5.0, 2.0, 2.0, 9.0, 10.0, 1.0, 8.0, 10.0, 4.0]
global b_y = 10
global p = [0.049, 0.56, 0.055, 0.271, 0.607, 0.046, 0.271, 0.906, 0.362, 0.162, 0.842, 0.283, 0.64, 0.553, 0.114, 0.292, 0.827, 0.677, 0.643, 0.354, 0.779, 0.3, 0.37, 0.34, 0.855, 0.814, 0.562, 0.915, 0.633, 0.126, 0.089, 0.947, 0.407, 0.417, 0.405, 0.494, 0.767, 0.047, 0.433, 0.804, 0.551, 0.413, 0.176, 0.658, 0.694, 0.542, 0.695, 0.387, 0.914, 0.844, 0.535, 0.584, 0.443, 0.744, 0.514, 0.98, 0.95, 0.971, 0.286, 0.38, 0.537, 0.121, 0.639, 0.926, 0.315, 0.869, 0.857, 0.332, 0.179, 0.274, 0.844, 0.413, 0.043, 0.407, 0.484, 0.948, 0.443, 0.187, 0.197, 0.018, 0.933, 0.26, 0.177, 0.983, 0.945, 0.722, 0.717, 0.687, 0.834, 0.726, 0.988, 0.106, 0.233, 0.355, 0.041, 0.53, 0.237, 0.925, 0.304, 0.153, 0.011, 0.207, 0.052, 0.668, 0.649, 0.382, 0.36, 0.46, 0.282, 0.068, 0.328, 0.18, 0.999, 0.639, 0.301, 0.05, 0.598, 0.931, 0.564, 0.718, 0.761, 0.882, 0.099, 0.872, 0.012, 0.265, 0.211, 0.188, 0.455, 0.123, 0.642, 0.956, 0.121, 0.622, 0.826, 0.024, 0.134, 0.436, 0.929, 0.691, 0.484, 0.875, 0.39, 0.761, 0.814, 0.602, 0.556, 0.35, 0.487, 0.592, 0.712, 0.242, 0.365, 0.039, 0.528, 0.566, 0.223, 0.225, 0.508, 0.531, 0.599, 0.746, 0.395, 0.495, 0.533, 0.8, 0.293, 0.229, 0.765, 0.599, 0.469, 0.697, 0.962, 0.183, 0.542, 0.462, 0.692, 0.189, 0.842, 0.983, 0.013, 0.377, 0.957, 0.481, 0.582, 0.431, 0.432, 0.154, 0.7, 0.273, 0.104, 0.714, 0.113, 0.988, 0.122, 0.149, 0.39, 0.878, 0.922, 0.094, 0.574, 0.925, 0.344, 0.676, 0.199, 0.839, 0.004, 0.63, 0.494, 0.65, 0.173, 0.882, 0.657, 0.541, 0.196, 0.836, 0.734, 0.795, 0.84, 0.276, 0.905, 0.573, 0.099, 0.791, 0.089, 0.214, 0.238, 0.453, 0.375, 0.467, 0.077, 0.594, 0.802, 0.776, 0.941, 0.201, 0.988, 0.574, 0.623, 0.322, 0.726, 0.341, 0.149, 0.075, 0.808, 0.052, 0.584, 0.98, 0.663, 0.91, 0.057, 0.548, 0.894, 0.455, 0.131, 0.941, 0.034, 0.552, 0.063, 0.816, 0.435, 0.232, 0.647, 0.856, 0.384, 0.399, 0.546, 0.146, 0.72, 0.52, 0.488, 0.158, 0.752, 0.098, 0.29, 0.03, 0.531, 0.329, 0.666, 0.18, 0.146, 0.623, 0.959, 0.269, 0.669, 0.755, 0.157, 0.43, 0.285, 0.527, 0.331, 0.561, 0.362, 0.272, 0.961, 0.001, 0.277, 0.471, 0.117, 0.524, 0.949, 0.109, 0.648, 0.884, 0.351, 0.759, 0.487, 0.833, 0.992, 0.019, 0.167, 0.323, 0.999, 0.255, 0.824, 0.463, 0.761, 0.32, 0.978, 0.361, 0.116, 0.723, 0.258, 0.892, 0.89, 0.582, 0.656, 0.591, 0.278, 0.752, 0.125, 0.213]
global q = [0.183, 0.624, 0.838, 0.542, 0.757, 0.693, 0.399, 0.995, 0.424, 0.332, 0.911, 0.666, 0.913, 0.688, 0.305, 0.774, 0.955, 0.916, 0.86, 0.447, 0.863, 0.808, 0.794, 0.986, 0.963, 0.904, 0.701, 0.934, 0.944, 0.449, 0.784, 0.987, 0.757, 0.703, 0.976, 0.893, 0.892, 0.329, 0.631, 0.966, 0.634, 0.748, 0.612, 0.927, 0.82, 0.675, 0.938, 0.481, 0.945, 0.982, 0.559, 0.604, 0.888, 0.781, 0.774, 0.988, 0.959, 0.994, 0.363, 0.876, 0.874, 0.724, 0.927, 0.97, 0.728, 0.921, 0.978, 0.4, 0.857, 0.875, 0.994, 0.73, 0.491, 0.778, 0.85, 0.959, 0.471, 0.457, 0.72, 0.532, 0.947, 0.911, 0.801, 0.983, 0.986, 0.727, 0.991, 0.79, 0.936, 0.871, 0.992, 0.475, 0.969, 0.438, 0.799, 0.964, 0.419, 0.942, 0.744, 0.649, 0.818, 0.335, 0.081, 0.728, 0.907, 0.532, 0.493, 0.832, 0.828, 0.875, 0.345, 0.528, 0.999, 0.64, 0.325, 0.986, 0.778, 0.948, 0.566, 0.805, 0.956, 0.945, 0.314, 0.939, 0.754, 0.375, 0.959, 0.781, 0.752, 0.251, 0.711, 0.97, 0.163, 0.643, 0.971, 0.233, 0.463, 0.796, 0.978, 0.809, 0.738, 0.893, 0.964, 0.89, 0.906, 0.638, 0.594, 0.539, 0.944, 0.906, 0.795, 0.479, 0.786, 0.095, 0.633, 0.773, 0.371, 0.673, 0.829, 0.557, 0.691, 0.827, 0.846, 0.496, 0.861, 0.953, 0.312, 0.843, 0.876, 0.945, 0.622, 0.749, 0.963, 0.626, 0.564, 0.995, 0.859, 0.647, 0.923, 0.988, 0.682, 0.835, 0.995, 0.795, 0.603, 0.671, 0.866, 0.321, 0.891, 0.571, 0.474, 0.759, 0.867, 0.988, 0.304, 0.943, 0.772, 0.912, 0.948, 0.251, 0.678, 0.926, 0.986, 0.729, 0.401, 0.98, 0.81, 0.967, 0.82, 0.998, 0.794, 0.966, 0.659, 0.94, 0.585, 0.926, 0.856, 0.835, 0.939, 0.914, 0.925, 0.766, 0.988, 0.852, 0.731, 0.423, 0.957, 0.67, 0.534, 0.651, 0.713, 0.753, 0.893, 0.887, 0.959, 0.563, 0.99, 0.668, 0.862, 0.701, 0.787, 0.771, 0.482, 0.38, 0.843, 0.737, 0.98, 0.984, 0.872, 0.973, 0.502, 0.664, 0.975, 0.855, 0.742, 0.983, 0.769, 0.986, 0.453, 0.875, 0.905, 0.827, 0.735, 0.971, 0.702, 0.936, 0.845, 0.91, 0.721, 0.837, 0.648, 0.519, 0.805, 0.186, 0.68, 0.986, 0.843, 0.418, 0.863, 0.577, 0.933, 0.97, 0.963, 0.9, 0.785, 0.96, 0.667, 0.76, 0.681, 0.787, 0.726, 0.896, 0.447, 0.597, 0.977, 0.76, 0.482, 0.862, 0.177, 0.65, 0.958, 0.206, 0.778, 0.958, 0.781, 0.857, 0.577, 0.952, 0.996, 0.038, 0.664, 0.706, 0.999, 0.689, 0.973, 0.464, 0.89, 0.791, 0.986, 0.598, 0.959, 0.723, 0.874, 0.972, 0.906, 0.853, 0.799, 0.615, 0.726, 0.844, 0.751, 0.424]
global origin = 1
global destination = 60