global arcs = [1 32; 1 48; 2 8; 2 10; 2 14; 2 19; 2 21; 2 40; 2 44; 2 47; 2 49; 3 37; 4 5; 4 6; 4 13; 4 15; 4 17; 4 26; 4 35; 5 22; 5 30; 5 33; 5 52; 5 54; 5 57; 6 3; 6 5; 6 27; 7 12; 7 13; 7 14; 7 19; 7 24; 7 25; 7 37; 7 39; 7 44; 7 46; 7 51; 8 20; 8 48; 9 7; 9 11; 9 16; 9 32; 9 35; 9 36; 9 38; 10 14; 10 28; 10 36; 10 44; 10 55; 11 7; 11 8; 11 20; 11 22; 11 35; 11 59; 12 2; 12 10; 12 20; 12 21; 12 44; 12 52; 13 29; 13 52; 13 55; 14 5; 14 12; 14 15; 14 16; 14 25; 14 29; 14 56; 14 59; 15 2; 15 4; 15 6; 15 10; 16 2; 16 3; 16 19; 16 44; 16 55; 17 3; 17 16; 17 33; 17 36; 17 49; 18 11; 18 13; 18 25; 18 35; 18 36; 18 51; 18 54; 19 12; 19 26; 19 28; 19 30; 19 35; 19 36; 19 47; 20 2; 20 9; 20 34; 20 37; 20 49; 20 59; 21 9; 21 11; 21 48; 22 14; 22 19; 22 25; 22 26; 22 52; 23 6; 23 32; 24 7; 24 10; 24 44; 24 45; 24 46; 25 10; 25 15; 25 20; 25 40; 25 42; 26 11; 26 32; 26 36; 26 43; 26 56; 27 6; 27 11; 27 24; 27 26; 27 31; 27 40; 27 46; 27 50; 27 54; 28 9; 28 10; 28 22; 28 31; 28 32; 28 36; 28 47; 28 48; 28 57; 28 58; 29 32; 29 34; 29 53; 30 4; 30 10; 30 13; 30 14; 30 25; 30 52; 31 25; 31 30; 31 40; 31 41; 31 53; 32 3; 32 6; 32 22; 32 23; 32 49; 33 11; 33 14; 33 17; 33 32; 33 37; 33 38; 33 39; 33 41; 33 42; 34 4; 34 22; 35 11; 35 17; 35 26; 35 30; 35 43; 35 49; 35 59; 36 2; 36 6; 36 11; 36 12; 36 16; 36 17; 36 30; 36 35; 37 5; 37 41; 37 51; 37 55; 37 58; 38 13; 38 37; 38 41; 38 49; 38 50; 39 3; 39 18; 39 22; 39 33; 39 41; 39 43; 39 54; 39 58; 39 59; 40 7; 40 9; 40 12; 40 15; 40 23; 40 33; 40 35; 40 37; 40 57; 41 19; 41 22; 41 52; 41 56; 42 2; 42 13; 42 28; 42 31; 42 36; 42 60; 43 12; 43 18; 43 40; 43 48; 43 52; 44 6; 44 11; 44 40; 44 42; 44 57; 44 58; 44 59; 45 4; 45 12; 45 20; 45 24; 45 40; 45 46; 45 50; 45 55; 46 4; 46 33; 46 50; 46 53; 46 60; 47 13; 47 34; 47 39; 47 51; 47 60; 48 8; 48 16; 48 25; 48 35; 48 50; 48 60; 49 3; 49 8; 49 15; 49 38; 49 42; 49 52; 49 55; 50 10; 50 16; 50 21; 50 32; 50 33; 50 40; 50 43; 51 14; 51 21; 51 26; 51 28; 51 42; 51 43; 51 48; 51 52; 51 59; 52 3; 52 27; 52 39; 52 57; 53 5; 53 9; 53 21; 53 37; 53 46; 53 49; 54 24; 54 25; 54 30; 54 35; 54 46; 54 55; 55 10; 55 25; 55 33; 55 47; 55 50; 56 17; 56 22; 56 30; 56 43; 56 53; 57 22; 57 23; 57 28; 58 3; 58 29; 58 46; 59 6; 59 12; 59 16; 59 30; 59 32; 59 36; 59 39]
global d_x = [2.0, 7.0, 10.0, 10.0, 10.0, 8.0, 3.0, 7.0, 1.0, 7.0, 1.0, 3.0, 3.0, 9.0, 8.0, 8.0, 10.0, 6.0, 3.0, 9.0, 1.0, 3.0, 3.0, 6.0, 6.0, 5.0, 4.0, 5.0, 8.0, 2.0, 5.0, 1.0, 4.0, 1.0, 3.0, 3.0, 4.0, 6.0, 10.0, 4.0, 7.0, 3.0, 8.0, 10.0, 5.0, 5.0, 10.0, 3.0, 10.0, 10.0, 10.0, 5.0, 3.0, 1.0, 2.0, 7.0, 8.0, 3.0, 2.0, 2.0, 8.0, 9.0, 5.0, 8.0, 2.0, 7.0, 4.0, 5.0, 8.0, 9.0, 7.0, 7.0, 5.0, 6.0, 2.0, 7.0, 2.0, 3.0, 4.0, 4.0, 7.0, 1.0, 2.0, 6.0, 9.0, 4.0, 9.0, 10.0, 3.0, 10.0, 7.0, 4.0, 1.0, 7.0, 3.0, 2.0, 3.0, 2.0, 10.0, 2.0, 2.0, 1.0, 6.0, 4.0, 8.0, 9.0, 9.0, 7.0, 10.0, 9.0, 6.0, 6.0, 8.0, 7.0, 4.0, 5.0, 10.0, 6.0, 7.0, 5.0, 5.0, 1.0, 8.0, 10.0, 6.0, 2.0, 2.0, 2.0, 6.0, 7.0, 6.0, 8.0, 9.0, 4.0, 8.0, 6.0, 8.0, 6.0, 9.0, 7.0, 5.0, 2.0, 7.0, 3.0, 6.0, 5.0, 7.0, 3.0, 10.0, 7.0, 10.0, 2.0, 4.0, 8.0, 6.0, 1.0, 8.0, 4.0, 2.0, 10.0, 4.0, 5.0, 1.0, 1.0, 8.0, 8.0, 8.0, 10.0, 5.0, 5.0, 3.0, 7.0, 6.0, 3.0, 7.0, 6.0, 10.0, 9.0, 2.0, 4.0, 2.0, 6.0, 8.0, 9.0, 3.0, 3.0, 9.0, 3.0, 8.0, 3.0, 7.0, 2.0, 2.0, 2.0, 8.0, 4.0, 7.0, 3.0, 8.0, 7.0, 2.0, 8.0, 4.0, 1.0, 2.0, 3.0, 2.0, 5.0, 9.0, 5.0, 2.0, 4.0, 8.0, 5.0, 1.0, 9.0, 5.0, 4.0, 3.0, 2.0, 6.0, 10.0, 10.0, 10.0, 9.0, 4.0, 8.0, 5.0, 1.0, 9.0, 4.0, 1.0, 2.0, 9.0, 6.0, 1.0, 9.0, 2.0, 9.0, 5.0, 9.0, 9.0, 4.0, 2.0, 3.0, 9.0, 2.0, 10.0, 9.0, 1.0, 6.0, 3.0, 8.0, 4.0, 5.0, 10.0, 4.0, 1.0, 2.0, 7.0, 2.0, 2.0, 6.0, 9.0, 2.0, 5.0, 6.0, 10.0, 10.0, 9.0, 6.0, 8.0, 6.0, 10.0, 3.0, 8.0, 7.0, 2.0, 2.0, 2.0, 3.0, 8.0, 8.0, 4.0, 2.0, 3.0, 6.0, 3.0, 1.0, 8.0, 7.0, 9.0, 6.0, 3.0, 10.0, 9.0, 1.0, 6.0, 4.0, 6.0, 2.0, 5.0, 2.0, 9.0, 5.0, 6.0, 4.0, 3.0, 7.0, 3.0, 10.0, 10.0, 8.0, 6.0, 2.0, 10.0, 7.0, 10.0, 9.0, 4.0, 8.0, 8.0, 9.0, 3.0, 8.0, 9.0, 1.0, 9.0, 10.0, 9.0, 7.0, 7.0, 3.0, 1.0, 3.0]
global b_x = 5
global d_y = [6.0, 6.0, 3.0, 9.0, 10.0, 9.0, 1.0, 9.0, 9.0, 2.0, 1.0, 6.0, 10.0, 4.0, 1.0, 10.0, 1.0, 1.0, 9.0, 9.0, 6.0, 2.0, 4.0, 1.0, 3.0, 5.0, 2.0, 10.0, 5.0, 8.0, 4.0, 9.0, 4.0, 5.0, 5.0, 1.0, 6.0, 9.0, 3.0, 5.0, 4.0, 2.0, 8.0, 9.0, 8.0, 1.0, 6.0, 6.0, 2.0, 6.0, 3.0, 1.0, 7.0, 9.0, 9.0, 4.0, 3.0, 9.0, 6.0, 2.0, 6.0, 2.0, 7.0, 6.0, 7.0, 2.0, 9.0, 8.0, 5.0, 3.0, 5.0, 8.0, 1.0, 8.0, 6.0, 9.0, 7.0, 2.0, 2.0, 9.0, 5.0, 7.0, 6.0, 5.0, 7.0, 1.0, 2.0, 10.0, 10.0, 4.0, 1.0, 4.0, 1.0, 4.0, 7.0, 1.0, 4.0, 2.0, 2.0, 5.0, 4.0, 4.0, 6.0, 5.0, 1.0, 4.0, 3.0, 10.0, 4.0, 9.0, 10.0, 8.0, 3.0, 10.0, 5.0, 2.0, 6.0, 10.0, 10.0, 5.0, 2.0, 9.0, 1.0, 7.0, 6.0, 8.0, 8.0, 6.0, 8.0, 9.0, 8.0, 10.0, 6.0, 9.0, 6.0, 1.0, 2.0, 3.0, 9.0, 8.0, 7.0, 4.0, 8.0, 3.0, 1.0, 1.0, 9.0, 6.0, 2.0, 2.0, 2.0, 9.0, 5.0, 5.0, 4.0, 1.0, 2.0, 2.0, 9.0, 4.0, 5.0, 6.0, 3.0, 10.0, 5.0, 4.0, 5.0, 8.0, 2.0, 10.0, 1.0, 5.0, 6.0, 8.0, 1.0, 5.0, 3.0, 3.0, 1.0, 2.0, 5.0, 2.0, 4.0, 4.0, 7.0, 2.0, 6.0, 4.0, 10.0, 3.0, 9.0, 6.0, 4.0, 1.0, 7.0, 6.0, 3.0, 9.0, 7.0, 5.0, 6.0, 3.0, 7.0, 9.0, 10.0, 9.0, 5.0, 1.0, 6.0, 4.0, 6.0, 2.0, 7.0, 2.0, 1.0, 4.0, 9.0, 3.0, 8.0, 5.0, 9.0, 6.0, 4.0, 6.0, 3.0, 3.0, 4.0, 3.0, 10.0, 10.0, 1.0, 6.0, 9.0, 3.0, 5.0, 3.0, 9.0, 6.0, 6.0, 1.0, 10.0, 7.0, 7.0, 1.0, 1.0, 8.0, 1.0, 8.0, 2.0, 4.0, 9.0, 4.0, 5.0, 9.0, 4.0, 5.0, 10.0, 6.0, 3.0, 8.0, 5.0, 9.0, 9.0, 1.0, 7.0, 2.0, 6.0, 9.0, 10.0, 3.0, 2.0, 6.0, 4.0, 3.0, 9.0, 8.0, 9.0, 10.0, 9.0, 1.0, 9.0, 9.0, 4.0, 3.0, 5.0, 7.0, 2.0, 1.0, 6.0, 9.0, 7.0, 9.0, 6.0, 7.0, 6.0, 7.0, 1.0, 9.0, 4.0, 7.0, 8.0, 7.0, 8.0, 10.0, 10.0, 2.0, 6.0, 3.0, 3.0, 8.0, 6.0, 7.0, 5.0, 6.0, 5.0, 8.0, 4.0, 7.0, 3.0, 6.0, 4.0, 4.0, 10.0, 10.0, 5.0, 4.0, 4.0, 2.0, 6.0, 7.0, 3.0, 7.0, 4.0, 7.0, 6.0]
global b_y = 10
global p = [0.618, 0.126, 0.32, 0.424, 0.813, 0.509, 0.92, 0.563, 0.142, 0.163, 0.63, 0.321, 0.042, 0.398, 0.191, 0.74, 0.381, 0.965, 0.117, 0.522, 0.702, 0.135, 0.518, 0.789, 0.495, 0.327, 0.888, 0.19, 0.189, 0.998, 0.284, 0.661, 0.072, 0.364, 0.3, 0.73, 0.458, 0.298, 0.783, 0.485, 0.281, 0.159, 0.154, 0.762, 0.016, 0.768, 0.584, 0.752, 0.712, 0.131, 0.943, 0.733, 0.815, 0.433, 0.907, 0.353, 0.315, 0.915, 0.971, 0.074, 0.482, 0.134, 0.757, 0.767, 0.393, 0.439, 0.525, 0.691, 0.348, 0.922, 0.135, 0.37, 0.113, 0.671, 0.726, 0.704, 0.301, 0.643, 0.581, 0.813, 0.69, 0.518, 0.29, 0.2, 0.799, 0.32, 0.899, 0.789, 0.096, 0.37, 0.525, 0.194, 0.152, 0.455, 0.314, 0.358, 0.495, 0.903, 0.546, 0.597, 0.644, 0.939, 0.852, 0.366, 0.814, 0.973, 0.306, 0.249, 0.092, 0.474, 0.841, 0.683, 0.93, 0.391, 0.708, 0.402, 0.028, 0.229, 0.48, 0.434, 0.164, 0.537, 0.963, 0.268, 0.891, 0.883, 0.295, 0.279, 0.034, 0.105, 0.632, 0.393, 0.531, 0.845, 0.818, 0.021, 0.094, 0.898, 0.75, 0.772, 0.102, 0.67, 0.877, 0.6, 0.19, 0.525, 0.968, 0.781, 0.51, 0.239, 0.005, 0.16, 0.997, 0.194, 0.223, 0.798, 0.898, 0.512, 0.992, 0.905, 0.349, 0.324, 0.197, 0.748, 0.912, 0.548, 0.961, 0.287, 0.637, 0.616, 0.425, 0.553, 0.588, 0.901, 0.95, 0.058, 0.556, 0.887, 0.99, 0.667, 0.181, 0.37, 0.445, 0.649, 0.308, 0.835, 0.686, 0.365, 0.885, 0.938, 0.154, 0.007, 0.566, 0.362, 0.469, 0.128, 0.835, 0.371, 0.819, 0.803, 0.135, 0.225, 0.643, 0.498, 0.155, 0.971, 0.298, 0.147, 0.541, 0.277, 0.014, 0.681, 0.987, 0.213, 0.113, 0.463, 0.123, 0.302, 0.301, 0.96, 0.825, 0.886, 0.876, 0.941, 0.739, 0.553, 0.125, 0.268, 0.402, 0.488, 0.884, 0.111, 0.429, 0.07, 0.416, 0.003, 0.833, 0.75, 0.51, 0.185, 0.825, 0.693, 0.134, 0.969, 0.978, 0.858, 0.018, 0.118, 0.64, 0.802, 0.693, 0.973, 0.436, 0.971, 0.311, 0.045, 0.854, 0.132, 0.397, 0.632, 0.325, 0.084, 0.548, 0.686, 0.077, 0.776, 0.86, 0.451, 0.935, 0.013, 0.269, 0.704, 0.773, 0.827, 0.652, 0.717, 0.62, 0.112, 0.206, 0.866, 0.659, 0.928, 0.3, 0.006, 0.731, 0.787, 0.024, 0.119, 0.729, 0.173, 0.959, 0.474, 0.266, 0.772, 0.631, 0.152, 0.998, 0.143, 0.812, 0.625, 0.464, 0.905, 0.708, 0.068, 0.733, 0.284, 0.017, 0.156, 0.076, 0.301, 0.528, 0.76, 0.711, 0.336, 0.371, 0.541, 0.368, 0.057, 0.964, 0.062, 0.522, 0.75, 0.057, 0.838, 0.73, 0.498, 0.524, 0.327, 0.976, 0.349, 0.236, 0.391, 0.083, 0.737, 0.439]
global q = [0.73, 0.604, 0.683, 0.991, 0.814, 0.551, 0.986, 0.701, 0.547, 0.887, 0.685, 0.366, 0.205, 0.515, 0.956, 0.876, 0.833, 0.98, 0.658, 0.726, 0.805, 0.342, 0.828, 0.814, 0.944, 0.656, 0.912, 0.425, 0.691, 0.999, 0.466, 0.934, 0.389, 0.463, 0.953, 0.916, 0.641, 0.405, 0.934, 0.62, 0.46, 0.29, 0.307, 0.995, 0.137, 0.962, 0.663, 0.977, 0.991, 0.784, 0.979, 0.858, 0.951, 0.764, 0.998, 0.357, 0.831, 0.949, 0.987, 0.864, 0.78, 0.31, 0.761, 0.805, 0.922, 0.49, 0.597, 0.839, 0.993, 0.96, 0.862, 0.741, 0.241, 0.876, 0.986, 0.776, 0.984, 0.885, 0.805, 0.869, 0.98, 0.766, 0.665, 0.876, 0.836, 0.464, 0.97, 0.833, 0.894, 0.522, 0.956, 0.436, 0.615, 0.961, 0.466, 0.724, 0.722, 0.989, 0.637, 0.662, 0.736, 0.991, 0.916, 0.64, 0.906, 0.988, 0.988, 0.827, 0.382, 0.985, 0.856, 0.88, 0.962, 0.826, 0.726, 0.663, 0.817, 0.456, 0.759, 0.461, 0.911, 0.86, 0.999, 0.356, 0.905, 0.98, 0.474, 0.771, 0.592, 0.758, 0.908, 0.551, 0.929, 0.853, 0.846, 0.65, 0.733, 0.948, 0.876, 0.923, 0.483, 0.913, 0.931, 0.92, 0.86, 0.619, 0.989, 0.997, 0.726, 0.374, 0.995, 0.638, 0.997, 0.529, 0.796, 0.918, 0.994, 0.604, 0.997, 0.915, 0.475, 0.519, 0.314, 0.956, 0.934, 0.988, 0.996, 0.916, 0.905, 0.909, 0.788, 0.665, 0.914, 0.938, 0.979, 0.363, 0.888, 0.959, 0.995, 0.691, 0.53, 0.873, 0.996, 0.917, 0.851, 0.911, 0.812, 0.478, 0.964, 0.956, 0.725, 0.303, 0.836, 0.497, 0.606, 0.986, 0.956, 0.505, 0.821, 0.815, 0.717, 0.941, 0.912, 0.934, 0.193, 0.979, 0.654, 0.651, 0.787, 0.484, 0.099, 0.983, 0.991, 0.314, 0.278, 0.908, 0.928, 0.965, 0.314, 0.967, 0.87, 0.963, 0.99, 0.945, 0.835, 0.797, 0.29, 0.863, 0.745, 0.567, 0.934, 0.97, 0.554, 0.66, 0.823, 0.009, 0.964, 0.95, 0.769, 0.342, 0.943, 0.719, 0.862, 0.99, 0.986, 0.934, 0.847, 0.189, 0.74, 0.925, 0.913, 0.976, 0.58, 0.995, 0.447, 0.098, 0.947, 0.913, 0.982, 0.944, 0.913, 0.986, 0.555, 0.956, 0.596, 0.951, 0.939, 0.795, 0.968, 0.536, 0.698, 0.943, 0.798, 0.951, 0.906, 0.937, 0.821, 0.466, 0.891, 0.912, 0.755, 0.973, 0.934, 0.353, 0.771, 0.833, 0.405, 0.25, 0.752, 0.861, 0.978, 0.593, 0.714, 0.916, 0.784, 0.424, 0.998, 0.684, 0.878, 0.697, 0.861, 0.948, 0.827, 0.575, 0.817, 0.452, 0.289, 0.763, 0.723, 0.531, 0.586, 0.924, 0.826, 0.784, 0.559, 0.563, 0.548, 0.161, 0.969, 0.657, 0.948, 0.897, 0.425, 0.839, 0.902, 0.537, 0.586, 0.47, 0.988, 0.74, 0.386, 0.771, 0.093, 0.942, 0.498]
global origin = 1
global destination = 60