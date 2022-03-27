global arcs = [1 17; 1 42; 2 10; 2 13; 2 24; 2 28; 2 32; 2 33; 2 41; 2 43; 3 8; 3 19; 3 23; 3 30; 3 32; 4 6; 4 14; 4 18; 4 20; 4 30; 4 33; 4 53; 4 54; 5 17; 5 23; 5 44; 6 9; 6 16; 6 20; 6 25; 6 26; 6 29; 6 41; 6 42; 6 48; 7 6; 7 10; 7 30; 7 51; 8 20; 8 28; 8 33; 8 52; 9 36; 9 38; 9 43; 9 49; 9 52; 9 59; 10 19; 10 32; 10 35; 10 37; 10 48; 10 57; 11 2; 11 9; 11 29; 11 32; 12 15; 12 32; 12 35; 12 37; 12 40; 12 45; 12 60; 13 15; 13 23; 13 26; 13 37; 13 57; 14 9; 14 16; 14 17; 14 36; 14 59; 15 8; 15 18; 15 44; 15 51; 16 5; 16 11; 16 22; 16 25; 16 27; 16 37; 16 39; 16 41; 16 45; 16 46; 16 49; 16 57; 17 12; 17 24; 17 31; 17 33; 17 37; 17 41; 17 44; 17 55; 18 11; 18 14; 18 16; 18 24; 18 25; 18 27; 18 31; 18 32; 18 33; 18 38; 18 42; 18 46; 18 48; 18 49; 19 9; 19 41; 19 46; 19 53; 19 54; 20 12; 20 32; 20 41; 20 44; 20 45; 21 38; 21 43; 21 46; 22 3; 22 23; 22 28; 22 30; 22 36; 22 46; 22 51; 23 5; 23 20; 23 46; 23 50; 23 57; 24 30; 24 33; 24 42; 24 51; 24 56; 25 4; 25 22; 25 32; 25 54; 26 38; 26 40; 27 5; 27 22; 27 44; 28 3; 28 4; 28 24; 28 46; 28 57; 29 15; 29 19; 29 20; 29 22; 29 26; 29 31; 29 47; 29 53; 30 20; 30 33; 30 44; 30 55; 30 56; 31 4; 31 14; 31 17; 31 29; 31 56; 32 12; 32 13; 32 58; 33 2; 33 13; 33 29; 34 14; 34 19; 34 27; 34 53; 34 56; 35 49; 35 56; 36 5; 36 35; 36 38; 36 50; 37 2; 37 55; 37 56; 38 2; 38 7; 38 10; 38 12; 38 13; 38 23; 38 30; 38 32; 38 49; 38 53; 39 2; 39 20; 39 41; 39 46; 39 48; 40 2; 40 3; 40 4; 40 21; 40 28; 40 35; 40 45; 41 2; 41 20; 41 28; 42 59; 43 3; 43 5; 43 13; 43 20; 43 35; 43 37; 43 46; 44 15; 44 18; 44 29; 44 30; 44 39; 44 43; 45 21; 45 22; 45 47; 45 57; 46 18; 46 19; 46 23; 46 28; 46 47; 46 57; 47 2; 47 9; 47 27; 47 28; 47 29; 47 54; 48 2; 48 7; 48 11; 48 12; 48 25; 48 30; 48 32; 48 43; 48 44; 48 53; 49 9; 49 15; 49 31; 49 50; 50 3; 50 13; 50 18; 50 30; 51 20; 51 31; 51 34; 51 42; 51 47; 51 48; 52 6; 52 18; 52 20; 52 21; 52 37; 52 38; 52 43; 53 19; 53 26; 53 28; 53 52; 54 15; 54 18; 54 40; 54 46; 55 2; 55 10; 55 13; 55 14; 55 19; 55 35; 55 42; 55 52; 56 24; 56 28; 56 48; 56 49; 57 12; 57 15; 57 25; 57 34; 57 37; 57 49; 57 51; 57 53; 58 2; 58 13; 58 29; 58 53; 59 3; 59 4; 59 6; 59 25; 59 55]
global d_x = [9.0, 3.0, 1.0, 10.0, 1.0, 9.0, 10.0, 2.0, 10.0, 10.0, 5.0, 3.0, 1.0, 9.0, 10.0, 8.0, 5.0, 5.0, 7.0, 7.0, 5.0, 10.0, 2.0, 9.0, 7.0, 9.0, 6.0, 3.0, 8.0, 4.0, 3.0, 6.0, 10.0, 1.0, 1.0, 8.0, 1.0, 4.0, 4.0, 5.0, 2.0, 9.0, 5.0, 1.0, 10.0, 2.0, 3.0, 2.0, 8.0, 2.0, 5.0, 8.0, 4.0, 7.0, 1.0, 8.0, 10.0, 2.0, 3.0, 1.0, 10.0, 9.0, 7.0, 6.0, 10.0, 8.0, 2.0, 7.0, 4.0, 3.0, 6.0, 4.0, 8.0, 10.0, 6.0, 6.0, 3.0, 7.0, 2.0, 2.0, 1.0, 9.0, 8.0, 3.0, 4.0, 7.0, 6.0, 4.0, 1.0, 10.0, 3.0, 8.0, 6.0, 3.0, 10.0, 6.0, 7.0, 8.0, 7.0, 9.0, 3.0, 7.0, 4.0, 2.0, 3.0, 6.0, 8.0, 8.0, 5.0, 7.0, 3.0, 9.0, 3.0, 9.0, 10.0, 5.0, 4.0, 3.0, 6.0, 7.0, 7.0, 5.0, 1.0, 2.0, 6.0, 10.0, 9.0, 9.0, 7.0, 8.0, 8.0, 9.0, 9.0, 7.0, 8.0, 6.0, 4.0, 1.0, 1.0, 8.0, 9.0, 9.0, 5.0, 8.0, 5.0, 2.0, 10.0, 9.0, 8.0, 7.0, 5.0, 10.0, 8.0, 1.0, 8.0, 2.0, 3.0, 1.0, 9.0, 8.0, 2.0, 4.0, 4.0, 9.0, 1.0, 5.0, 7.0, 5.0, 4.0, 2.0, 1.0, 6.0, 9.0, 2.0, 7.0, 1.0, 5.0, 3.0, 3.0, 10.0, 2.0, 7.0, 8.0, 9.0, 9.0, 4.0, 9.0, 1.0, 9.0, 9.0, 3.0, 5.0, 10.0, 7.0, 3.0, 6.0, 8.0, 2.0, 2.0, 9.0, 4.0, 5.0, 5.0, 8.0, 9.0, 2.0, 3.0, 7.0, 9.0, 4.0, 2.0, 9.0, 1.0, 3.0, 4.0, 7.0, 10.0, 8.0, 3.0, 3.0, 6.0, 10.0, 8.0, 9.0, 10.0, 7.0, 1.0, 3.0, 2.0, 5.0, 5.0, 3.0, 6.0, 8.0, 9.0, 7.0, 7.0, 1.0, 3.0, 3.0, 4.0, 2.0, 4.0, 9.0, 8.0, 3.0, 4.0, 6.0, 1.0, 4.0, 7.0, 6.0, 3.0, 7.0, 2.0, 7.0, 4.0, 5.0, 2.0, 9.0, 6.0, 6.0, 7.0, 9.0, 5.0, 4.0, 6.0, 2.0, 1.0, 8.0, 1.0, 9.0, 3.0, 1.0, 9.0, 5.0, 7.0, 3.0, 1.0, 2.0, 10.0, 7.0, 6.0, 9.0, 5.0, 10.0, 1.0, 1.0, 10.0, 6.0, 7.0, 8.0, 1.0, 4.0, 3.0, 2.0, 4.0, 5.0, 10.0, 1.0, 7.0, 8.0, 10.0, 7.0, 2.0, 1.0, 8.0, 7.0, 10.0, 5.0, 5.0, 7.0, 2.0, 1.0, 5.0, 9.0, 3.0, 6.0, 9.0]
global b_x = 5
global d_y = [10.0, 4.0, 8.0, 7.0, 3.0, 5.0, 5.0, 1.0, 5.0, 6.0, 3.0, 1.0, 3.0, 9.0, 6.0, 6.0, 2.0, 10.0, 2.0, 10.0, 9.0, 1.0, 7.0, 1.0, 3.0, 3.0, 3.0, 3.0, 9.0, 2.0, 2.0, 4.0, 3.0, 8.0, 9.0, 9.0, 9.0, 1.0, 9.0, 9.0, 5.0, 8.0, 6.0, 1.0, 1.0, 4.0, 3.0, 8.0, 9.0, 8.0, 3.0, 4.0, 6.0, 2.0, 6.0, 6.0, 6.0, 2.0, 2.0, 7.0, 2.0, 9.0, 9.0, 7.0, 6.0, 7.0, 7.0, 8.0, 9.0, 7.0, 7.0, 9.0, 8.0, 5.0, 7.0, 3.0, 5.0, 8.0, 8.0, 10.0, 7.0, 3.0, 10.0, 2.0, 6.0, 7.0, 1.0, 7.0, 9.0, 9.0, 7.0, 10.0, 10.0, 6.0, 5.0, 4.0, 3.0, 5.0, 5.0, 2.0, 4.0, 8.0, 3.0, 1.0, 10.0, 3.0, 2.0, 3.0, 8.0, 4.0, 9.0, 2.0, 3.0, 1.0, 10.0, 1.0, 4.0, 4.0, 1.0, 4.0, 2.0, 9.0, 6.0, 4.0, 6.0, 7.0, 9.0, 5.0, 6.0, 7.0, 4.0, 8.0, 6.0, 2.0, 6.0, 3.0, 7.0, 1.0, 6.0, 10.0, 7.0, 3.0, 6.0, 6.0, 2.0, 1.0, 5.0, 2.0, 3.0, 9.0, 10.0, 8.0, 7.0, 6.0, 5.0, 9.0, 7.0, 5.0, 1.0, 3.0, 3.0, 10.0, 1.0, 7.0, 3.0, 5.0, 2.0, 4.0, 3.0, 1.0, 3.0, 5.0, 9.0, 6.0, 6.0, 4.0, 7.0, 10.0, 3.0, 7.0, 3.0, 2.0, 6.0, 9.0, 7.0, 2.0, 4.0, 2.0, 2.0, 9.0, 8.0, 9.0, 5.0, 1.0, 3.0, 2.0, 9.0, 7.0, 7.0, 4.0, 2.0, 2.0, 7.0, 4.0, 1.0, 4.0, 8.0, 6.0, 4.0, 7.0, 2.0, 9.0, 4.0, 4.0, 3.0, 4.0, 1.0, 10.0, 3.0, 1.0, 8.0, 7.0, 5.0, 2.0, 7.0, 9.0, 1.0, 5.0, 2.0, 4.0, 9.0, 2.0, 10.0, 8.0, 3.0, 7.0, 9.0, 6.0, 3.0, 9.0, 8.0, 5.0, 4.0, 1.0, 7.0, 5.0, 1.0, 3.0, 2.0, 1.0, 1.0, 8.0, 4.0, 8.0, 8.0, 5.0, 6.0, 4.0, 2.0, 3.0, 4.0, 8.0, 1.0, 2.0, 4.0, 10.0, 6.0, 6.0, 10.0, 6.0, 10.0, 10.0, 7.0, 8.0, 9.0, 4.0, 7.0, 7.0, 3.0, 10.0, 7.0, 3.0, 2.0, 5.0, 1.0, 10.0, 10.0, 6.0, 8.0, 8.0, 1.0, 10.0, 4.0, 9.0, 10.0, 4.0, 2.0, 6.0, 2.0, 10.0, 1.0, 9.0, 4.0, 2.0, 1.0, 6.0, 5.0, 6.0, 9.0, 8.0, 6.0, 2.0, 4.0, 2.0, 7.0, 6.0, 1.0, 6.0, 2.0]
global b_y = 10
global p = [0.263, 0.062, 0.784, 0.469, 0.294, 0.608, 0.232, 0.093, 0.101, 0.189, 0.135, 0.312, 0.391, 0.9, 0.957, 0.285, 0.697, 0.208, 0.263, 0.8, 0.34, 0.668, 0.852, 0.02, 0.897, 0.615, 0.588, 0.116, 0.841, 0.92, 0.799, 0.974, 0.721, 0.798, 0.042, 0.035, 0.251, 0.641, 0.744, 0.161, 0.034, 0.12, 0.102, 0.331, 0.309, 0.223, 0.063, 0.548, 0.924, 0.183, 0.703, 0.329, 0.899, 0.178, 0.416, 0.085, 0.988, 0.395, 0.05, 0.329, 0.238, 0.694, 0.779, 0.399, 0.457, 0.882, 0.388, 0.43, 0.56, 0.16, 0.124, 0.954, 0.659, 0.642, 0.756, 0.103, 0.253, 0.998, 0.349, 0.771, 0.213, 0.604, 0.97, 0.59, 0.008, 0.336, 0.302, 0.279, 0.016, 0.756, 0.747, 0.497, 0.041, 0.664, 0.335, 0.995, 0.548, 0.928, 0.747, 0.658, 0.611, 0.859, 0.632, 0.662, 0.726, 0.777, 0.682, 0.82, 0.155, 0.496, 0.486, 0.2, 0.362, 0.706, 0.441, 0.288, 0.695, 0.645, 0.883, 0.78, 0.739, 0.611, 0.9, 0.283, 0.477, 0.493, 0.041, 0.895, 0.252, 0.236, 0.132, 0.341, 0.756, 0.111, 0.519, 0.637, 0.973, 0.533, 0.809, 0.395, 0.339, 0.225, 0.125, 0.607, 0.822, 0.596, 0.403, 0.338, 0.996, 0.171, 0.883, 0.543, 0.226, 0.249, 0.151, 0.615, 0.647, 0.052, 0.045, 0.554, 0.441, 0.989, 0.413, 0.811, 0.396, 0.659, 0.854, 0.533, 0.375, 0.084, 0.583, 0.866, 0.406, 0.287, 0.557, 0.751, 0.551, 0.884, 0.762, 0.137, 0.971, 0.706, 0.528, 0.059, 0.067, 0.559, 0.478, 0.885, 0.798, 0.587, 0.122, 0.643, 0.892, 0.096, 0.123, 0.866, 0.381, 0.114, 0.107, 0.888, 0.594, 0.002, 0.472, 0.605, 0.727, 0.533, 0.65, 0.845, 0.804, 0.521, 0.285, 0.871, 0.484, 0.595, 0.186, 0.937, 0.626, 0.611, 0.413, 0.732, 0.054, 0.523, 0.272, 0.338, 0.523, 0.093, 0.114, 0.732, 0.289, 0.165, 0.934, 0.149, 0.733, 0.901, 0.184, 0.788, 0.369, 0.197, 0.3, 0.91, 0.487, 0.757, 0.605, 0.48, 0.344, 0.126, 0.391, 0.048, 0.914, 0.012, 0.918, 0.075, 0.296, 0.992, 0.898, 0.194, 0.699, 0.265, 0.111, 0.134, 0.831, 0.985, 0.824, 0.813, 0.341, 0.568, 0.805, 0.445, 0.944, 0.673, 0.328, 0.738, 0.927, 0.756, 0.672, 0.119, 0.9, 0.785, 0.611, 0.719, 0.737, 0.172, 0.593, 0.836, 0.253, 0.979, 0.624, 0.823, 0.657, 0.595, 0.059, 0.147, 0.319, 0.463, 0.913, 0.093, 0.197, 0.063, 0.575, 0.496, 0.034, 0.247, 0.089, 0.559, 0.095, 0.915, 0.44, 0.227, 0.967, 0.4, 0.676, 0.469, 0.723, 0.865, 0.722, 0.237, 0.967, 0.825, 0.689]
global q = [0.952, 0.824, 0.794, 0.752, 0.616, 0.724, 0.286, 0.941, 0.333, 0.632, 0.928, 0.721, 0.823, 0.946, 0.995, 0.954, 0.701, 0.633, 0.679, 0.941, 0.427, 0.823, 0.912, 0.7, 0.947, 0.617, 0.756, 0.759, 0.974, 0.981, 0.922, 0.975, 0.934, 0.992, 0.729, 0.983, 0.625, 0.82, 0.806, 0.173, 0.455, 0.683, 0.443, 0.86, 0.484, 0.895, 0.318, 0.8, 0.951, 0.821, 0.785, 0.782, 0.954, 0.714, 0.584, 0.842, 0.988, 0.487, 0.817, 0.498, 0.276, 0.857, 0.97, 0.636, 0.482, 0.892, 0.679, 0.535, 0.735, 0.71, 0.169, 0.974, 0.683, 0.886, 0.954, 0.134, 0.857, 0.998, 0.516, 0.92, 0.224, 0.91, 0.991, 0.86, 0.784, 0.723, 0.473, 0.369, 0.064, 0.96, 0.984, 0.816, 0.732, 0.83, 0.634, 0.997, 0.986, 0.977, 0.989, 0.697, 0.805, 0.958, 0.884, 0.891, 0.908, 0.815, 0.799, 0.995, 0.626, 0.596, 0.977, 0.984, 0.94, 0.82, 0.74, 0.304, 0.886, 0.709, 0.965, 0.918, 0.95, 0.82, 0.943, 0.872, 0.869, 0.956, 0.532, 0.915, 0.458, 0.898, 0.365, 0.498, 0.877, 0.806, 0.859, 0.983, 0.995, 0.535, 0.974, 0.888, 0.813, 0.968, 0.816, 0.981, 0.918, 0.747, 0.922, 0.436, 0.999, 0.418, 0.919, 0.572, 0.514, 0.753, 0.755, 0.795, 0.881, 0.807, 0.794, 0.772, 0.807, 0.989, 0.836, 0.94, 0.647, 0.672, 0.887, 0.953, 0.578, 0.593, 0.706, 0.988, 0.698, 0.551, 0.89, 0.981, 0.879, 0.885, 0.784, 0.215, 0.973, 0.756, 0.669, 0.104, 0.925, 0.745, 0.653, 0.973, 0.844, 0.635, 0.145, 0.923, 0.964, 0.776, 0.488, 0.896, 0.962, 0.828, 0.887, 0.971, 0.951, 0.782, 0.473, 0.829, 0.733, 0.925, 0.906, 0.929, 0.874, 0.62, 0.529, 0.942, 0.526, 0.609, 0.748, 0.988, 0.949, 0.889, 0.43, 0.856, 0.485, 0.914, 0.604, 0.823, 0.951, 0.545, 0.538, 0.994, 0.52, 0.924, 0.976, 0.966, 0.951, 0.948, 0.614, 0.931, 0.851, 0.909, 0.43, 0.932, 0.495, 0.985, 0.891, 0.71, 0.581, 0.311, 0.607, 0.77, 0.996, 0.836, 0.981, 0.515, 0.433, 0.996, 0.9, 0.613, 0.83, 0.513, 0.571, 0.865, 0.947, 0.988, 0.874, 0.84, 0.579, 0.625, 0.97, 0.621, 0.976, 0.738, 0.565, 0.894, 0.971, 0.944, 0.994, 0.827, 0.958, 0.96, 0.976, 0.812, 0.791, 0.628, 0.888, 0.84, 0.821, 0.991, 0.822, 0.84, 0.844, 0.648, 0.333, 0.471, 0.946, 0.827, 0.915, 0.188, 0.507, 0.918, 0.908, 0.747, 0.335, 0.606, 0.974, 0.849, 0.732, 0.981, 0.748, 0.911, 0.972, 0.876, 0.846, 0.736, 0.889, 0.925, 0.755, 0.771, 0.969, 0.88, 0.786]
global origin = 1
global destination = 60