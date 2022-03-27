global arcs = [1 4; 1 6; 1 14; 1 16; 1 20; 1 22; 1 24; 1 28; 1 49; 2 11; 2 19; 2 40; 2 42; 2 49; 2 57; 3 9; 3 10; 3 13; 3 23; 3 26; 3 39; 3 49; 4 3; 4 13; 4 16; 4 50; 5 11; 5 13; 5 14; 5 24; 5 30; 5 39; 5 41; 5 53; 6 9; 6 11; 6 16; 6 19; 6 20; 6 27; 6 34; 6 36; 6 40; 6 52; 6 57; 7 11; 7 21; 7 25; 7 29; 7 50; 8 28; 8 29; 8 46; 9 5; 9 17; 9 43; 9 60; 10 18; 10 19; 10 57; 11 2; 11 9; 11 29; 11 32; 11 36; 11 43; 11 44; 11 45; 11 54; 11 60; 12 17; 12 21; 12 26; 12 45; 13 19; 13 23; 13 27; 13 54; 13 59; 14 11; 14 15; 14 16; 15 9; 15 14; 15 16; 15 19; 15 29; 15 34; 16 22; 16 34; 17 15; 17 24; 17 30; 17 32; 17 34; 17 47; 17 53; 18 17; 18 21; 18 28; 19 3; 19 4; 19 5; 19 8; 19 10; 19 15; 20 4; 20 17; 20 23; 20 27; 20 30; 20 35; 20 47; 20 48; 21 38; 21 45; 22 9; 22 13; 22 24; 22 38; 22 46; 22 54; 23 5; 23 28; 23 35; 23 37; 23 46; 23 56; 24 15; 24 17; 24 22; 24 27; 24 29; 25 6; 25 26; 25 34; 25 42; 25 46; 25 60; 26 3; 26 17; 26 28; 26 43; 26 57; 26 60; 27 3; 27 9; 27 10; 27 13; 27 23; 27 28; 27 31; 27 43; 27 47; 28 5; 28 12; 28 17; 28 20; 28 36; 28 37; 28 48; 28 58; 29 16; 29 22; 29 25; 29 36; 29 53; 29 55; 29 57; 30 4; 30 5; 30 19; 30 39; 30 43; 30 50; 30 56; 31 11; 31 17; 31 28; 31 34; 31 42; 31 47; 31 54; 31 60; 32 8; 32 23; 32 37; 32 45; 32 52; 33 2; 33 9; 33 35; 33 44; 34 25; 34 37; 34 38; 34 43; 34 52; 34 53; 35 16; 35 27; 35 36; 35 40; 35 41; 36 2; 36 18; 36 20; 36 27; 36 29; 36 40; 36 57; 36 59; 37 53; 37 60; 38 24; 38 32; 38 36; 38 46; 38 48; 38 58; 39 14; 39 20; 39 48; 39 51; 39 53; 40 2; 40 6; 40 12; 40 17; 40 24; 40 29; 40 37; 40 39; 40 57; 41 17; 41 22; 41 24; 41 27; 41 30; 41 35; 41 51; 41 55; 41 57; 42 6; 42 7; 42 14; 42 16; 42 25; 42 51; 43 3; 43 19; 43 23; 43 29; 43 45; 44 33; 44 38; 45 10; 45 14; 45 28; 45 38; 45 43; 45 47; 46 7; 46 23; 46 25; 46 27; 46 48; 46 50; 46 51; 47 10; 47 24; 47 27; 47 35; 47 51; 48 26; 48 27; 48 28; 48 44; 49 41; 49 46; 50 13; 50 29; 50 37; 50 57; 51 26; 51 43; 52 9; 52 12; 52 29; 52 35; 52 45; 52 58; 53 26; 53 29; 53 45; 54 21; 54 26; 54 41; 54 43; 54 48; 54 52; 55 12; 55 35; 55 37; 55 39; 55 41; 55 47; 56 3; 56 4; 56 15; 56 18; 56 39; 56 59; 57 4; 57 12; 57 13; 57 19; 57 43; 57 48; 57 60; 58 5; 58 23; 58 38; 58 44; 58 55; 59 19; 59 37; 59 38; 59 40; 59 49; 59 55; 59 58]
global d_x = [10.0, 1.0, 8.0, 10.0, 6.0, 7.0, 4.0, 5.0, 6.0, 4.0, 3.0, 1.0, 2.0, 8.0, 10.0, 2.0, 10.0, 9.0, 2.0, 6.0, 7.0, 9.0, 5.0, 2.0, 3.0, 9.0, 2.0, 9.0, 3.0, 3.0, 5.0, 4.0, 10.0, 2.0, 2.0, 1.0, 2.0, 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 2.0, 1.0, 5.0, 2.0, 4.0, 3.0, 2.0, 10.0, 9.0, 3.0, 4.0, 1.0, 10.0, 5.0, 6.0, 10.0, 6.0, 9.0, 5.0, 4.0, 9.0, 3.0, 9.0, 4.0, 8.0, 8.0, 9.0, 9.0, 1.0, 6.0, 7.0, 8.0, 2.0, 3.0, 5.0, 9.0, 10.0, 8.0, 6.0, 4.0, 7.0, 5.0, 10.0, 6.0, 1.0, 1.0, 2.0, 10.0, 9.0, 9.0, 5.0, 1.0, 4.0, 5.0, 6.0, 7.0, 2.0, 5.0, 8.0, 5.0, 1.0, 8.0, 10.0, 10.0, 8.0, 7.0, 8.0, 9.0, 7.0, 6.0, 9.0, 7.0, 10.0, 2.0, 2.0, 7.0, 6.0, 3.0, 10.0, 5.0, 6.0, 7.0, 5.0, 8.0, 7.0, 10.0, 8.0, 6.0, 3.0, 1.0, 9.0, 4.0, 2.0, 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 9.0, 9.0, 1.0, 7.0, 6.0, 6.0, 5.0, 4.0, 9.0, 2.0, 3.0, 4.0, 10.0, 5.0, 2.0, 10.0, 3.0, 1.0, 10.0, 1.0, 4.0, 6.0, 4.0, 4.0, 1.0, 3.0, 1.0, 4.0, 4.0, 1.0, 7.0, 10.0, 5.0, 8.0, 8.0, 5.0, 5.0, 2.0, 9.0, 6.0, 2.0, 8.0, 3.0, 5.0, 6.0, 10.0, 10.0, 6.0, 10.0, 9.0, 8.0, 1.0, 5.0, 4.0, 6.0, 7.0, 7.0, 4.0, 6.0, 4.0, 10.0, 3.0, 9.0, 10.0, 3.0, 7.0, 5.0, 3.0, 1.0, 6.0, 9.0, 1.0, 8.0, 9.0, 8.0, 2.0, 8.0, 2.0, 8.0, 10.0, 10.0, 3.0, 6.0, 2.0, 9.0, 5.0, 4.0, 9.0, 6.0, 6.0, 5.0, 10.0, 9.0, 9.0, 3.0, 9.0, 7.0, 1.0, 5.0, 10.0, 2.0, 2.0, 6.0, 7.0, 9.0, 1.0, 6.0, 9.0, 7.0, 10.0, 2.0, 7.0, 5.0, 3.0, 6.0, 1.0, 3.0, 1.0, 10.0, 10.0, 1.0, 8.0, 7.0, 6.0, 10.0, 3.0, 4.0, 8.0, 9.0, 7.0, 2.0, 9.0, 9.0, 4.0, 4.0, 8.0, 2.0, 4.0, 6.0, 4.0, 10.0, 1.0, 3.0, 5.0, 3.0, 8.0, 5.0, 10.0, 8.0, 9.0, 8.0, 8.0, 9.0, 9.0, 10.0, 3.0, 9.0, 10.0, 8.0, 7.0, 6.0, 10.0, 1.0, 6.0, 2.0, 7.0, 2.0, 6.0, 4.0, 5.0, 3.0, 8.0, 9.0, 3.0, 6.0, 10.0, 7.0, 2.0, 1.0, 4.0, 2.0, 9.0, 5.0, 9.0, 9.0, 5.0, 8.0]
global b_x = 5
global d_y = [9.0, 4.0, 8.0, 7.0, 6.0, 10.0, 10.0, 4.0, 4.0, 6.0, 4.0, 2.0, 8.0, 7.0, 1.0, 5.0, 7.0, 2.0, 2.0, 4.0, 4.0, 2.0, 2.0, 1.0, 8.0, 2.0, 7.0, 6.0, 8.0, 5.0, 10.0, 7.0, 6.0, 8.0, 4.0, 4.0, 4.0, 1.0, 6.0, 3.0, 7.0, 8.0, 10.0, 2.0, 7.0, 6.0, 5.0, 7.0, 2.0, 1.0, 10.0, 9.0, 4.0, 2.0, 8.0, 6.0, 7.0, 2.0, 2.0, 5.0, 10.0, 5.0, 9.0, 8.0, 9.0, 10.0, 10.0, 1.0, 5.0, 7.0, 9.0, 9.0, 3.0, 1.0, 3.0, 5.0, 10.0, 3.0, 3.0, 9.0, 9.0, 8.0, 7.0, 2.0, 2.0, 9.0, 1.0, 4.0, 2.0, 1.0, 10.0, 8.0, 3.0, 8.0, 9.0, 7.0, 7.0, 6.0, 2.0, 2.0, 8.0, 10.0, 8.0, 5.0, 6.0, 1.0, 7.0, 8.0, 2.0, 6.0, 1.0, 9.0, 1.0, 6.0, 2.0, 10.0, 4.0, 3.0, 2.0, 8.0, 2.0, 6.0, 6.0, 10.0, 9.0, 9.0, 3.0, 6.0, 1.0, 3.0, 4.0, 8.0, 4.0, 6.0, 2.0, 5.0, 1.0, 6.0, 4.0, 2.0, 5.0, 7.0, 5.0, 1.0, 5.0, 5.0, 4.0, 10.0, 3.0, 2.0, 6.0, 4.0, 2.0, 3.0, 3.0, 3.0, 9.0, 9.0, 7.0, 10.0, 9.0, 6.0, 5.0, 1.0, 8.0, 2.0, 6.0, 1.0, 9.0, 3.0, 2.0, 6.0, 2.0, 4.0, 10.0, 3.0, 3.0, 6.0, 8.0, 2.0, 8.0, 5.0, 1.0, 5.0, 6.0, 7.0, 2.0, 5.0, 2.0, 7.0, 6.0, 8.0, 7.0, 3.0, 10.0, 1.0, 3.0, 8.0, 9.0, 7.0, 6.0, 10.0, 9.0, 1.0, 1.0, 5.0, 7.0, 4.0, 7.0, 4.0, 8.0, 9.0, 3.0, 10.0, 4.0, 2.0, 3.0, 2.0, 1.0, 8.0, 8.0, 5.0, 8.0, 1.0, 6.0, 3.0, 5.0, 7.0, 4.0, 7.0, 1.0, 4.0, 5.0, 2.0, 9.0, 10.0, 3.0, 6.0, 4.0, 8.0, 3.0, 3.0, 5.0, 2.0, 3.0, 5.0, 8.0, 6.0, 4.0, 8.0, 2.0, 4.0, 9.0, 10.0, 9.0, 9.0, 5.0, 6.0, 10.0, 10.0, 9.0, 9.0, 10.0, 5.0, 9.0, 10.0, 1.0, 1.0, 7.0, 5.0, 4.0, 5.0, 2.0, 1.0, 4.0, 3.0, 3.0, 6.0, 5.0, 7.0, 6.0, 2.0, 1.0, 5.0, 2.0, 1.0, 9.0, 1.0, 10.0, 7.0, 3.0, 3.0, 5.0, 5.0, 3.0, 6.0, 1.0, 8.0, 9.0, 2.0, 8.0, 5.0, 9.0, 3.0, 5.0, 3.0, 3.0, 3.0, 5.0, 4.0, 2.0, 1.0, 9.0, 5.0, 4.0, 6.0, 6.0, 1.0, 2.0, 10.0, 9.0, 10.0, 10.0, 1.0, 8.0, 3.0, 10.0, 9.0, 1.0, 8.0, 4.0, 3.0]
global b_y = 10
global p = [0.092, 0.91, 0.225, 0.753, 0.964, 0.638, 0.95, 0.036, 0.681, 0.717, 0.512, 0.744, 0.373, 0.601, 0.352, 0.71, 0.868, 0.933, 0.075, 0.922, 0.572, 0.404, 0.163, 0.701, 0.187, 0.484, 0.072, 0.341, 0.223, 0.575, 0.188, 0.854, 0.787, 0.776, 0.272, 0.835, 0.335, 0.128, 0.139, 0.319, 0.042, 0.952, 0.291, 0.608, 0.77, 0.091, 0.139, 0.235, 0.45, 0.587, 0.175, 0.373, 0.008, 0.49, 0.103, 0.063, 0.68, 0.303, 0.875, 0.524, 0.818, 0.361, 0.45, 0.092, 0.267, 0.265, 0.118, 0.634, 0.613, 0.536, 0.759, 0.126, 0.224, 0.779, 0.32, 0.107, 0.463, 0.301, 0.217, 0.22, 0.681, 0.221, 0.276, 0.171, 0.831, 0.9, 0.71, 0.235, 0.607, 0.521, 0.404, 0.953, 0.062, 0.673, 0.78, 0.375, 0.625, 0.712, 0.447, 0.382, 0.076, 0.864, 0.603, 0.371, 0.112, 0.435, 0.623, 0.057, 0.763, 0.624, 0.472, 0.303, 0.777, 0.959, 0.734, 0.685, 0.137, 0.72, 0.82, 0.571, 0.124, 0.004, 0.339, 0.967, 0.806, 0.643, 0.636, 0.361, 0.318, 0.413, 0.395, 0.977, 0.332, 0.914, 0.195, 0.104, 0.908, 0.85, 0.423, 0.447, 0.477, 0.942, 0.088, 0.879, 0.178, 0.03, 0.907, 0.904, 0.879, 0.464, 0.817, 0.373, 0.013, 0.294, 0.73, 0.132, 0.103, 0.31, 0.98, 0.019, 0.923, 0.826, 0.071, 0.457, 0.362, 0.455, 0.047, 0.986, 0.309, 0.957, 0.886, 0.361, 0.28, 0.404, 0.836, 0.807, 0.84, 0.544, 0.036, 0.448, 0.561, 0.343, 0.796, 0.847, 0.903, 0.408, 0.242, 0.285, 0.969, 0.58, 0.947, 0.086, 0.814, 0.274, 0.412, 0.627, 0.041, 0.644, 0.149, 0.645, 0.589, 0.776, 0.533, 0.377, 0.547, 0.208, 0.962, 0.728, 0.044, 0.648, 0.26, 0.798, 0.961, 0.897, 0.281, 0.779, 0.339, 0.843, 0.998, 0.166, 0.499, 0.212, 0.308, 0.042, 0.655, 0.265, 0.918, 0.566, 0.692, 0.661, 0.376, 0.545, 0.36, 0.439, 0.443, 0.013, 0.781, 0.463, 0.764, 0.782, 0.868, 0.361, 0.838, 0.913, 0.694, 0.48, 0.193, 0.415, 0.115, 0.493, 0.946, 0.44, 0.668, 0.457, 0.214, 0.584, 0.371, 0.107, 0.041, 0.649, 0.251, 0.002, 0.357, 0.618, 0.801, 0.888, 0.511, 0.585, 0.772, 0.838, 0.894, 0.109, 0.439, 0.37, 0.671, 0.169, 0.289, 0.892, 0.638, 0.372, 0.327, 0.861, 0.281, 0.22, 0.655, 0.563, 0.57, 0.016, 0.788, 0.748, 0.036, 0.691, 0.54, 0.404, 0.894, 0.86, 0.626, 0.62, 0.636, 0.291, 0.338, 0.588, 0.601, 0.404, 0.522, 0.089, 0.082, 0.14, 0.473, 0.55, 0.698, 0.517, 0.512, 0.906, 0.61, 0.604, 0.561, 0.045, 0.613, 0.493, 0.003, 0.73, 0.534, 0.888, 0.225, 0.101, 0.922, 0.596, 0.78, 0.257, 0.897, 0.047]
global q = [0.921, 0.974, 0.559, 0.795, 0.985, 0.939, 0.956, 0.172, 0.944, 0.786, 0.906, 0.779, 0.737, 0.925, 0.576, 0.896, 0.878, 0.973, 0.771, 0.968, 0.803, 0.55, 0.513, 0.957, 0.534, 0.51, 0.418, 0.925, 0.591, 0.621, 0.385, 0.894, 0.937, 0.863, 0.656, 0.985, 0.856, 0.733, 0.214, 0.728, 0.949, 0.966, 0.295, 0.842, 0.967, 0.321, 0.53, 0.851, 0.768, 0.733, 0.298, 0.582, 0.39, 0.952, 0.916, 0.588, 0.723, 0.835, 0.996, 0.629, 0.992, 0.922, 0.527, 0.886, 0.542, 0.333, 0.392, 0.688, 0.762, 0.612, 0.791, 0.363, 0.815, 0.819, 0.544, 0.605, 0.564, 0.43, 0.321, 0.391, 0.892, 0.533, 0.67, 0.842, 0.919, 0.946, 0.893, 0.246, 0.734, 0.86, 0.777, 0.985, 0.319, 0.95, 0.892, 0.775, 0.747, 0.729, 0.47, 0.788, 0.559, 0.991, 0.877, 0.518, 0.92, 0.729, 0.88, 0.834, 0.993, 0.849, 0.699, 0.652, 0.873, 0.962, 0.92, 0.91, 0.202, 0.759, 0.993, 0.957, 0.93, 0.811, 0.571, 0.993, 0.949, 0.734, 0.978, 0.832, 0.541, 0.598, 0.994, 0.984, 0.674, 0.926, 0.791, 0.71, 0.961, 0.864, 0.693, 0.61, 0.596, 0.948, 0.152, 0.901, 0.627, 0.787, 0.919, 0.928, 0.966, 0.572, 0.933, 0.56, 0.859, 0.485, 0.988, 0.296, 0.911, 0.899, 0.991, 0.732, 0.936, 0.907, 0.214, 0.995, 0.964, 0.504, 0.398, 0.993, 0.98, 0.993, 0.889, 0.796, 0.675, 0.604, 0.959, 0.967, 0.988, 0.611, 0.447, 0.747, 0.987, 0.899, 0.831, 0.951, 0.931, 0.417, 0.732, 0.333, 0.982, 0.749, 0.998, 0.632, 0.892, 0.574, 0.438, 0.726, 0.553, 0.883, 0.791, 0.841, 0.969, 0.812, 0.74, 0.615, 0.991, 0.212, 0.987, 0.73, 0.928, 0.834, 0.439, 0.94, 0.976, 0.925, 0.514, 0.88, 0.92, 0.93, 0.998, 0.587, 0.763, 0.521, 0.542, 0.062, 0.836, 0.553, 0.949, 0.81, 0.829, 0.997, 0.637, 0.876, 0.521, 0.991, 0.814, 0.348, 0.845, 0.773, 0.997, 0.983, 0.949, 0.624, 0.846, 0.946, 0.907, 0.663, 0.227, 0.713, 0.234, 0.628, 0.967, 0.475, 0.938, 0.567, 0.416, 0.83, 0.523, 0.573, 0.178, 0.852, 0.93, 0.629, 0.706, 0.772, 0.872, 0.969, 0.923, 0.619, 0.948, 0.917, 0.953, 0.307, 0.696, 0.694, 0.913, 0.227, 0.818, 0.894, 0.825, 0.394, 0.771, 0.933, 0.436, 0.783, 0.677, 0.592, 0.603, 0.698, 0.92, 0.927, 0.6, 0.965, 0.927, 0.6, 0.916, 0.919, 0.858, 0.811, 0.816, 0.332, 0.432, 0.88, 0.848, 0.438, 0.848, 0.973, 0.716, 0.246, 0.68, 0.985, 0.74, 0.82, 0.578, 0.917, 0.698, 0.693, 0.808, 0.649, 0.913, 0.519, 0.473, 0.895, 0.838, 0.977, 0.279, 0.313, 0.96, 0.632, 0.964, 0.386, 0.903, 0.985]
global origin = 1
global destination = 60