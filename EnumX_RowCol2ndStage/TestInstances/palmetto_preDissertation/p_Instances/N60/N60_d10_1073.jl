global arcs = [1 3; 1 23; 2 5; 2 24; 2 26; 2 29; 2 43; 2 50; 3 2; 3 21; 3 29; 3 31; 3 37; 3 56; 4 14; 4 30; 4 48; 4 54; 5 12; 5 15; 5 20; 5 25; 5 36; 5 39; 5 49; 5 55; 6 12; 6 18; 6 35; 6 48; 6 59; 7 19; 7 37; 7 40; 7 49; 7 60; 8 9; 8 26; 8 30; 8 52; 9 47; 10 4; 10 27; 10 31; 10 35; 10 36; 10 57; 11 5; 11 15; 11 20; 11 26; 11 47; 11 57; 11 58; 12 8; 12 19; 12 21; 12 34; 12 45; 12 58; 13 9; 13 10; 13 12; 13 25; 13 29; 13 43; 13 53; 14 4; 14 9; 14 16; 14 32; 14 45; 14 49; 14 57; 15 20; 15 22; 15 31; 15 48; 15 51; 15 58; 16 8; 16 10; 16 21; 16 29; 17 9; 17 10; 17 26; 17 44; 18 3; 18 21; 18 51; 18 54; 19 34; 19 35; 19 36; 19 37; 20 23; 20 27; 20 36; 20 38; 20 48; 20 54; 20 58; 21 8; 21 10; 21 19; 21 24; 21 29; 21 59; 22 4; 22 9; 22 29; 22 39; 22 41; 22 46; 23 5; 23 15; 23 55; 24 17; 24 18; 24 22; 24 25; 24 26; 24 27; 24 28; 24 35; 24 41; 25 14; 25 31; 25 33; 25 43; 25 49; 26 2; 26 4; 26 6; 26 9; 26 12; 26 22; 26 28; 26 38; 26 40; 26 41; 26 47; 26 60; 27 12; 27 22; 27 45; 27 53; 27 56; 28 18; 28 22; 28 23; 28 27; 28 30; 28 35; 28 50; 28 51; 29 16; 29 39; 29 47; 29 50; 29 55; 30 18; 30 28; 30 36; 30 53; 31 3; 31 10; 31 16; 31 30; 31 34; 31 40; 31 53; 32 8; 32 39; 32 42; 32 52; 32 53; 32 60; 33 7; 33 23; 33 52; 34 2; 34 17; 34 29; 34 50; 34 59; 35 5; 35 6; 35 11; 35 36; 35 51; 35 53; 36 3; 36 4; 36 9; 36 18; 36 21; 36 25; 36 26; 36 32; 36 43; 36 55; 37 12; 37 25; 37 32; 37 35; 37 41; 37 50; 37 53; 37 55; 38 14; 38 16; 38 19; 38 28; 38 49; 38 53; 39 7; 39 13; 39 16; 39 28; 39 31; 39 40; 39 46; 39 49; 40 6; 40 22; 40 33; 40 53; 41 19; 41 34; 41 40; 41 43; 41 44; 41 48; 41 52; 41 58; 42 9; 42 13; 42 21; 42 22; 42 29; 42 53; 42 54; 43 2; 43 34; 43 48; 44 5; 44 13; 44 28; 44 39; 44 50; 44 58; 45 15; 45 22; 45 44; 45 47; 45 52; 45 55; 45 58; 46 39; 46 55; 46 59; 47 20; 47 22; 47 26; 47 48; 47 51; 47 56; 48 3; 48 20; 48 26; 48 27; 48 54; 48 60; 49 2; 49 14; 49 39; 50 9; 50 25; 50 33; 50 34; 50 35; 51 18; 51 42; 51 57; 52 2; 52 11; 52 13; 52 17; 52 20; 52 42; 52 55; 52 60; 53 28; 53 34; 53 35; 53 55; 53 56; 54 11; 54 12; 54 14; 54 22; 54 40; 54 48; 55 3; 55 28; 55 39; 55 40; 56 23; 56 46; 56 48; 56 55; 57 3; 57 26; 57 43; 57 44; 57 51; 58 7; 58 25; 58 39; 58 46; 58 57; 58 59; 58 60; 59 5; 59 10; 59 15; 59 19; 59 20; 59 25; 59 35; 59 44; 59 60]
global d_x = [8.0, 7.0, 6.0, 2.0, 7.0, 1.0, 6.0, 10.0, 4.0, 4.0, 6.0, 10.0, 2.0, 4.0, 6.0, 1.0, 9.0, 1.0, 7.0, 3.0, 1.0, 8.0, 1.0, 7.0, 4.0, 4.0, 1.0, 1.0, 10.0, 3.0, 2.0, 3.0, 7.0, 9.0, 1.0, 4.0, 6.0, 4.0, 3.0, 9.0, 7.0, 3.0, 8.0, 8.0, 10.0, 10.0, 3.0, 7.0, 6.0, 5.0, 8.0, 1.0, 8.0, 1.0, 3.0, 10.0, 9.0, 4.0, 9.0, 6.0, 10.0, 10.0, 1.0, 6.0, 10.0, 8.0, 5.0, 5.0, 1.0, 2.0, 5.0, 3.0, 3.0, 6.0, 8.0, 10.0, 3.0, 8.0, 7.0, 1.0, 6.0, 6.0, 6.0, 4.0, 8.0, 1.0, 1.0, 10.0, 10.0, 4.0, 3.0, 1.0, 9.0, 5.0, 9.0, 2.0, 2.0, 10.0, 5.0, 10.0, 8.0, 2.0, 2.0, 10.0, 9.0, 6.0, 10.0, 7.0, 9.0, 2.0, 3.0, 7.0, 5.0, 2.0, 9.0, 9.0, 2.0, 8.0, 10.0, 4.0, 1.0, 8.0, 6.0, 9.0, 5.0, 9.0, 1.0, 5.0, 10.0, 5.0, 1.0, 6.0, 5.0, 9.0, 10.0, 9.0, 9.0, 4.0, 1.0, 5.0, 5.0, 7.0, 3.0, 8.0, 2.0, 1.0, 8.0, 7.0, 2.0, 5.0, 9.0, 8.0, 3.0, 2.0, 2.0, 8.0, 1.0, 5.0, 3.0, 3.0, 5.0, 10.0, 5.0, 7.0, 10.0, 6.0, 7.0, 1.0, 6.0, 7.0, 9.0, 3.0, 8.0, 5.0, 9.0, 7.0, 1.0, 4.0, 2.0, 8.0, 3.0, 10.0, 10.0, 1.0, 5.0, 9.0, 5.0, 7.0, 8.0, 4.0, 6.0, 7.0, 7.0, 3.0, 5.0, 5.0, 9.0, 7.0, 6.0, 7.0, 3.0, 4.0, 1.0, 7.0, 4.0, 1.0, 7.0, 4.0, 8.0, 3.0, 7.0, 5.0, 2.0, 1.0, 8.0, 10.0, 5.0, 1.0, 9.0, 3.0, 9.0, 3.0, 6.0, 9.0, 8.0, 2.0, 6.0, 9.0, 1.0, 10.0, 9.0, 2.0, 8.0, 3.0, 3.0, 6.0, 8.0, 2.0, 6.0, 4.0, 9.0, 2.0, 7.0, 5.0, 10.0, 3.0, 2.0, 6.0, 10.0, 3.0, 6.0, 7.0, 4.0, 4.0, 6.0, 5.0, 3.0, 7.0, 6.0, 6.0, 3.0, 7.0, 9.0, 3.0, 6.0, 8.0, 6.0, 4.0, 4.0, 6.0, 8.0, 8.0, 6.0, 2.0, 6.0, 6.0, 10.0, 9.0, 10.0, 2.0, 1.0, 1.0, 4.0, 2.0, 8.0, 5.0, 7.0, 4.0, 4.0, 7.0, 7.0, 7.0, 8.0, 9.0, 8.0, 2.0, 1.0, 1.0, 4.0, 4.0, 1.0, 2.0, 10.0, 10.0, 6.0, 7.0, 8.0, 4.0, 5.0, 6.0, 8.0, 1.0, 3.0, 1.0, 8.0, 1.0, 1.0, 8.0, 4.0, 2.0, 6.0, 3.0, 6.0, 6.0, 7.0, 1.0, 9.0, 10.0, 10.0, 4.0, 8.0, 1.0, 3.0, 5.0]
global b_x = 5
global d_y = [4.0, 7.0, 6.0, 5.0, 2.0, 2.0, 1.0, 10.0, 9.0, 4.0, 5.0, 2.0, 5.0, 10.0, 4.0, 9.0, 3.0, 10.0, 10.0, 8.0, 4.0, 3.0, 2.0, 7.0, 7.0, 2.0, 3.0, 6.0, 4.0, 8.0, 6.0, 2.0, 9.0, 7.0, 8.0, 5.0, 10.0, 5.0, 5.0, 1.0, 9.0, 5.0, 9.0, 1.0, 10.0, 3.0, 4.0, 5.0, 9.0, 1.0, 9.0, 3.0, 7.0, 6.0, 8.0, 5.0, 9.0, 9.0, 7.0, 2.0, 10.0, 10.0, 1.0, 5.0, 4.0, 6.0, 4.0, 6.0, 2.0, 4.0, 7.0, 10.0, 3.0, 9.0, 5.0, 1.0, 7.0, 8.0, 9.0, 4.0, 2.0, 8.0, 6.0, 7.0, 6.0, 8.0, 2.0, 4.0, 9.0, 5.0, 10.0, 8.0, 8.0, 6.0, 2.0, 10.0, 3.0, 5.0, 1.0, 7.0, 9.0, 9.0, 10.0, 2.0, 7.0, 2.0, 6.0, 3.0, 1.0, 1.0, 6.0, 10.0, 7.0, 2.0, 6.0, 8.0, 4.0, 6.0, 4.0, 6.0, 5.0, 4.0, 5.0, 6.0, 8.0, 9.0, 5.0, 8.0, 4.0, 3.0, 10.0, 6.0, 8.0, 3.0, 10.0, 2.0, 8.0, 1.0, 4.0, 3.0, 3.0, 4.0, 8.0, 7.0, 3.0, 8.0, 3.0, 10.0, 5.0, 1.0, 1.0, 6.0, 10.0, 7.0, 10.0, 6.0, 3.0, 10.0, 8.0, 1.0, 4.0, 10.0, 3.0, 10.0, 2.0, 8.0, 3.0, 5.0, 3.0, 8.0, 3.0, 1.0, 10.0, 3.0, 3.0, 10.0, 8.0, 4.0, 10.0, 6.0, 4.0, 3.0, 5.0, 7.0, 10.0, 3.0, 1.0, 1.0, 7.0, 2.0, 3.0, 3.0, 7.0, 4.0, 4.0, 4.0, 6.0, 4.0, 7.0, 4.0, 9.0, 2.0, 4.0, 10.0, 6.0, 10.0, 1.0, 5.0, 1.0, 9.0, 6.0, 8.0, 7.0, 6.0, 1.0, 7.0, 10.0, 8.0, 6.0, 8.0, 3.0, 5.0, 8.0, 2.0, 9.0, 6.0, 8.0, 10.0, 2.0, 2.0, 1.0, 6.0, 10.0, 3.0, 3.0, 2.0, 7.0, 7.0, 5.0, 2.0, 7.0, 3.0, 4.0, 4.0, 4.0, 7.0, 8.0, 6.0, 7.0, 7.0, 7.0, 8.0, 9.0, 4.0, 2.0, 5.0, 2.0, 1.0, 9.0, 10.0, 7.0, 8.0, 5.0, 7.0, 2.0, 4.0, 6.0, 10.0, 10.0, 7.0, 10.0, 9.0, 1.0, 7.0, 4.0, 7.0, 3.0, 1.0, 6.0, 4.0, 3.0, 4.0, 1.0, 6.0, 2.0, 4.0, 1.0, 4.0, 6.0, 1.0, 4.0, 10.0, 6.0, 9.0, 2.0, 3.0, 5.0, 7.0, 6.0, 6.0, 3.0, 10.0, 9.0, 9.0, 6.0, 5.0, 5.0, 8.0, 2.0, 5.0, 9.0, 7.0, 3.0, 6.0, 8.0, 4.0, 1.0, 7.0, 9.0, 1.0, 3.0, 8.0, 10.0, 3.0, 9.0, 7.0, 3.0, 4.0, 8.0, 4.0, 3.0, 7.0, 10.0, 2.0]
global b_y = 10
global p = [0.788, 0.835, 0.126, 0.606, 0.26, 0.826, 0.486, 0.091, 0.358, 0.087, 0.523, 0.073, 0.753, 0.51, 0.966, 0.955, 0.745, 0.211, 0.479, 0.488, 0.119, 0.426, 0.879, 0.323, 0.282, 0.536, 0.002, 0.776, 0.415, 0.28, 0.561, 0.998, 0.549, 0.887, 0.482, 0.059, 0.246, 0.661, 0.096, 0.347, 0.647, 0.18, 0.597, 0.418, 0.922, 0.358, 0.526, 0.788, 0.033, 0.063, 0.593, 0.873, 0.193, 0.395, 0.301, 0.379, 0.777, 0.746, 0.293, 0.888, 0.736, 0.119, 0.733, 0.165, 0.317, 0.068, 0.049, 0.558, 0.722, 0.843, 0.486, 0.783, 0.353, 0.815, 0.943, 0.973, 0.761, 0.823, 0.41, 0.212, 0.386, 0.599, 0.681, 0.236, 0.923, 0.081, 0.517, 0.394, 0.879, 0.198, 0.508, 0.074, 0.976, 0.986, 0.058, 0.941, 0.113, 0.2, 0.268, 0.617, 0.15, 0.612, 0.66, 0.789, 0.837, 0.463, 0.672, 0.11, 0.308, 0.661, 0.247, 0.116, 0.771, 0.906, 0.825, 0.902, 0.375, 0.759, 0.758, 0.573, 0.915, 0.983, 0.543, 0.257, 0.636, 0.335, 0.357, 0.219, 0.083, 0.812, 0.475, 0.424, 0.384, 0.217, 0.447, 0.141, 0.506, 0.178, 0.939, 0.977, 0.264, 0.624, 0.998, 0.954, 0.991, 0.287, 0.279, 0.498, 0.805, 0.58, 0.472, 0.052, 0.443, 0.163, 0.904, 0.414, 0.815, 0.129, 0.493, 0.848, 0.192, 0.61, 0.387, 0.564, 0.573, 0.947, 0.388, 0.698, 0.148, 0.933, 0.114, 0.27, 0.899, 0.638, 0.341, 0.336, 0.363, 0.585, 0.507, 0.223, 0.211, 0.633, 0.858, 0.139, 0.036, 0.59, 0.793, 0.441, 0.566, 0.459, 0.456, 0.972, 0.3, 0.76, 0.047, 0.503, 0.671, 0.608, 0.212, 0.061, 0.119, 0.463, 0.741, 0.014, 0.925, 0.263, 0.822, 0.127, 0.57, 0.574, 0.133, 0.791, 0.21, 0.455, 0.668, 0.941, 0.959, 0.558, 0.306, 0.615, 0.894, 0.21, 0.489, 0.169, 0.016, 0.549, 0.309, 0.937, 0.158, 0.815, 0.709, 0.572, 0.447, 0.062, 0.095, 0.327, 0.837, 0.575, 0.68, 0.515, 0.088, 0.251, 0.748, 0.815, 0.331, 0.043, 0.865, 0.445, 0.384, 0.101, 0.003, 0.068, 0.799, 0.485, 0.675, 0.871, 0.652, 0.399, 0.264, 0.06, 0.523, 0.944, 0.386, 0.707, 0.919, 0.532, 0.21, 0.273, 0.872, 0.701, 0.77, 0.07, 0.374, 0.983, 0.299, 0.64, 0.606, 0.08, 0.587, 0.572, 0.466, 0.026, 0.198, 0.773, 0.635, 0.243, 0.652, 0.876, 0.013, 0.312, 0.634, 0.355, 0.782, 0.49, 0.596, 0.141, 0.575, 0.063, 0.424, 0.104, 0.887, 0.446, 0.824, 0.939, 0.281, 0.104, 0.497, 0.756, 0.566, 0.686, 0.159, 0.12, 0.318, 0.556, 0.643, 0.171, 0.299, 0.954, 0.449, 0.331, 0.971, 0.844, 0.503, 0.718, 0.663, 0.828, 0.038, 0.972, 0.499, 0.943, 0.325, 0.729, 0.573, 0.994]
global q = [0.893, 0.896, 0.198, 0.698, 0.543, 0.838, 0.883, 0.772, 0.625, 0.801, 0.784, 0.976, 0.82, 0.983, 0.981, 0.99, 0.97, 0.526, 0.912, 0.675, 0.581, 0.834, 0.925, 0.541, 0.705, 0.99, 0.464, 0.935, 0.564, 0.346, 0.673, 0.999, 0.877, 0.922, 0.724, 0.414, 0.28, 0.668, 0.329, 0.657, 0.914, 0.352, 0.677, 0.767, 0.954, 0.862, 0.941, 0.969, 0.815, 0.122, 0.851, 0.905, 0.224, 0.828, 0.418, 0.662, 0.847, 0.972, 0.35, 0.961, 0.847, 0.521, 0.879, 0.943, 0.485, 0.081, 0.502, 0.719, 0.795, 0.87, 0.611, 0.905, 0.543, 0.975, 0.97, 0.993, 0.865, 0.834, 0.509, 0.68, 0.54, 0.943, 0.745, 0.496, 0.953, 0.433, 0.647, 0.974, 0.896, 0.535, 0.614, 0.929, 0.998, 0.986, 0.669, 0.958, 0.282, 0.262, 0.757, 0.675, 0.225, 0.646, 0.971, 0.846, 0.984, 0.565, 0.914, 0.919, 0.613, 0.745, 0.559, 0.691, 0.798, 0.917, 0.915, 0.943, 0.6, 0.765, 0.779, 0.947, 0.984, 0.99, 0.89, 0.64, 0.8, 0.848, 0.5, 0.775, 0.924, 0.998, 0.584, 0.881, 0.542, 0.995, 0.528, 0.726, 0.624, 0.981, 0.963, 0.982, 0.711, 0.777, 0.999, 0.979, 0.997, 0.487, 0.333, 0.711, 0.89, 0.675, 0.622, 0.707, 0.925, 0.498, 0.914, 0.513, 0.953, 0.151, 0.749, 0.944, 0.425, 0.872, 0.429, 0.85, 0.875, 0.974, 0.825, 0.962, 0.425, 0.946, 0.44, 0.56, 0.952, 0.675, 0.884, 0.806, 0.936, 0.714, 0.719, 0.233, 0.752, 0.713, 0.949, 0.488, 0.998, 0.827, 0.975, 0.797, 0.66, 0.573, 0.983, 0.987, 0.717, 0.784, 0.264, 0.611, 0.74, 0.752, 0.37, 0.303, 0.784, 0.913, 0.787, 0.826, 0.997, 0.285, 0.915, 0.822, 0.989, 0.729, 0.159, 0.807, 0.493, 0.949, 0.979, 0.973, 0.977, 0.615, 0.306, 0.988, 0.948, 0.859, 0.954, 0.771, 0.454, 0.942, 0.694, 0.979, 0.61, 0.858, 0.769, 0.696, 0.978, 0.513, 0.835, 0.657, 0.976, 0.674, 0.746, 0.682, 0.983, 0.57, 0.888, 0.921, 0.474, 0.882, 0.935, 0.503, 0.719, 0.498, 0.848, 0.449, 0.917, 0.616, 0.942, 0.952, 0.992, 0.675, 0.617, 0.172, 0.821, 0.961, 0.467, 0.761, 0.998, 0.965, 0.381, 0.817, 0.887, 0.925, 0.871, 0.322, 0.934, 0.983, 0.866, 0.997, 0.782, 0.404, 0.729, 0.821, 0.806, 0.554, 0.626, 0.965, 0.822, 0.739, 0.705, 0.876, 0.101, 0.438, 0.846, 0.609, 0.962, 0.734, 0.671, 0.362, 0.665, 0.687, 0.447, 0.359, 0.906, 0.487, 0.893, 0.978, 0.649, 0.368, 0.809, 0.833, 0.816, 0.999, 0.731, 0.712, 0.723, 0.824, 0.979, 0.494, 0.525, 0.956, 0.467, 0.548, 0.979, 0.919, 0.569, 0.905, 0.706, 0.883, 0.131, 0.999, 0.598, 0.95, 0.986, 0.766, 0.822, 0.996]
global origin = 1
global destination = 60