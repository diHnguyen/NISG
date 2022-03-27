global arcs = [1 10; 1 31; 1 35; 1 38; 2 18; 2 23; 2 28; 2 30; 2 33; 2 47; 3 14; 3 29; 3 32; 3 43; 4 8; 4 15; 4 33; 4 34; 4 39; 4 42; 4 45; 4 48; 4 49; 5 4; 5 17; 5 29; 5 34; 6 5; 6 18; 6 31; 6 37; 6 50; 7 8; 7 14; 7 25; 7 27; 7 33; 8 6; 8 29; 8 33; 8 39; 8 42; 9 10; 9 16; 9 23; 9 40; 9 50; 10 9; 10 13; 10 25; 10 42; 11 5; 11 13; 11 15; 11 21; 11 38; 11 41; 11 44; 12 3; 12 4; 12 10; 12 11; 13 28; 13 29; 13 35; 13 41; 13 46; 14 32; 14 37; 15 12; 15 13; 15 35; 15 37; 15 46; 15 49; 16 4; 16 6; 16 11; 16 18; 16 42; 16 45; 17 30; 17 46; 18 29; 18 33; 18 44; 18 47; 19 14; 19 15; 19 31; 19 34; 19 39; 19 43; 20 5; 20 12; 20 26; 21 13; 21 34; 21 43; 22 11; 22 17; 23 3; 23 9; 23 36; 23 44; 24 6; 24 8; 24 16; 24 25; 25 4; 25 6; 25 7; 25 8; 25 13; 25 34; 25 36; 25 40; 26 24; 26 47; 27 23; 27 31; 27 45; 28 14; 28 18; 28 20; 28 24; 28 25; 28 35; 28 37; 28 41; 28 47; 29 7; 29 17; 29 40; 29 50; 30 16; 30 19; 30 33; 30 38; 30 44; 30 48; 31 2; 31 24; 31 37; 31 42; 32 10; 32 25; 32 40; 33 14; 33 17; 33 18; 33 28; 33 36; 33 37; 33 43; 33 45; 34 20; 34 25; 34 30; 34 33; 34 45; 34 47; 34 49; 35 6; 35 20; 35 37; 35 42; 35 50; 36 4; 36 22; 36 48; 36 50; 37 7; 37 15; 37 30; 37 33; 37 36; 37 39; 37 50; 38 7; 38 23; 38 43; 38 46; 39 26; 39 33; 39 40; 39 50; 40 5; 40 7; 40 43; 40 48; 41 7; 41 32; 41 48; 42 13; 42 41; 43 11; 43 25; 43 30; 43 39; 43 48; 44 13; 44 16; 44 27; 44 39; 44 40; 44 42; 44 46; 44 47; 44 48; 45 7; 45 10; 45 13; 45 16; 45 18; 46 17; 46 19; 46 25; 46 33; 47 4; 47 11; 47 34; 47 46; 48 2; 48 21; 48 29; 48 37; 48 40; 49 6; 49 20; 49 43; 49 45]
global d_x = [4.0, 1.0, 4.0, 5.0, 4.0, 8.0, 9.0, 1.0, 7.0, 2.0, 3.0, 10.0, 1.0, 3.0, 6.0, 5.0, 8.0, 10.0, 8.0, 9.0, 6.0, 3.0, 6.0, 10.0, 5.0, 2.0, 8.0, 5.0, 7.0, 10.0, 6.0, 5.0, 3.0, 3.0, 9.0, 7.0, 10.0, 7.0, 10.0, 4.0, 7.0, 4.0, 9.0, 8.0, 1.0, 5.0, 6.0, 3.0, 3.0, 2.0, 3.0, 5.0, 7.0, 5.0, 3.0, 5.0, 8.0, 5.0, 7.0, 8.0, 4.0, 9.0, 10.0, 4.0, 6.0, 4.0, 2.0, 3.0, 6.0, 4.0, 1.0, 10.0, 2.0, 1.0, 9.0, 1.0, 9.0, 6.0, 7.0, 3.0, 2.0, 10.0, 4.0, 8.0, 8.0, 2.0, 9.0, 5.0, 2.0, 6.0, 1.0, 6.0, 4.0, 5.0, 1.0, 9.0, 3.0, 6.0, 3.0, 1.0, 6.0, 3.0, 5.0, 2.0, 5.0, 8.0, 4.0, 1.0, 2.0, 2.0, 7.0, 4.0, 5.0, 1.0, 6.0, 6.0, 1.0, 7.0, 2.0, 3.0, 8.0, 5.0, 1.0, 5.0, 1.0, 10.0, 8.0, 9.0, 8.0, 8.0, 4.0, 8.0, 1.0, 1.0, 2.0, 7.0, 10.0, 1.0, 1.0, 9.0, 7.0, 1.0, 4.0, 7.0, 5.0, 8.0, 2.0, 2.0, 8.0, 7.0, 6.0, 5.0, 9.0, 9.0, 4.0, 6.0, 8.0, 9.0, 7.0, 8.0, 6.0, 7.0, 4.0, 4.0, 6.0, 2.0, 1.0, 1.0, 10.0, 4.0, 2.0, 6.0, 1.0, 8.0, 2.0, 5.0, 1.0, 3.0, 8.0, 4.0, 8.0, 3.0, 8.0, 3.0, 10.0, 5.0, 4.0, 5.0, 4.0, 8.0, 10.0, 5.0, 7.0, 2.0, 2.0, 2.0, 1.0, 5.0, 10.0, 3.0, 7.0, 1.0, 8.0, 5.0, 2.0, 7.0, 9.0, 7.0, 3.0, 1.0, 10.0, 8.0, 3.0, 9.0, 6.0, 4.0, 6.0, 9.0, 5.0, 1.0, 5.0, 10.0, 10.0, 9.0, 8.0, 8.0, 8.0, 2.0, 1.0, 9.0, 3.0, 4.0]
global b_x = 5
global d_y = [3.0, 4.0, 2.0, 7.0, 7.0, 9.0, 7.0, 9.0, 3.0, 6.0, 2.0, 2.0, 6.0, 8.0, 3.0, 6.0, 6.0, 10.0, 7.0, 10.0, 1.0, 2.0, 8.0, 10.0, 1.0, 4.0, 3.0, 7.0, 6.0, 3.0, 4.0, 4.0, 1.0, 8.0, 3.0, 7.0, 5.0, 9.0, 2.0, 6.0, 3.0, 8.0, 3.0, 4.0, 6.0, 2.0, 3.0, 4.0, 8.0, 7.0, 5.0, 5.0, 1.0, 8.0, 10.0, 10.0, 6.0, 3.0, 2.0, 9.0, 4.0, 1.0, 1.0, 1.0, 10.0, 7.0, 8.0, 3.0, 2.0, 10.0, 2.0, 3.0, 1.0, 2.0, 5.0, 5.0, 5.0, 6.0, 5.0, 4.0, 8.0, 6.0, 10.0, 4.0, 9.0, 10.0, 8.0, 5.0, 8.0, 1.0, 7.0, 7.0, 8.0, 10.0, 9.0, 2.0, 5.0, 5.0, 3.0, 2.0, 1.0, 10.0, 1.0, 9.0, 4.0, 9.0, 4.0, 3.0, 1.0, 7.0, 5.0, 7.0, 7.0, 7.0, 5.0, 5.0, 2.0, 10.0, 6.0, 6.0, 7.0, 6.0, 8.0, 8.0, 10.0, 8.0, 4.0, 6.0, 4.0, 1.0, 4.0, 6.0, 6.0, 2.0, 9.0, 4.0, 9.0, 6.0, 6.0, 4.0, 5.0, 6.0, 1.0, 3.0, 7.0, 4.0, 1.0, 1.0, 6.0, 4.0, 10.0, 3.0, 7.0, 5.0, 2.0, 1.0, 4.0, 4.0, 5.0, 10.0, 5.0, 3.0, 6.0, 5.0, 7.0, 7.0, 4.0, 5.0, 2.0, 5.0, 9.0, 4.0, 4.0, 5.0, 1.0, 7.0, 1.0, 5.0, 4.0, 10.0, 2.0, 4.0, 8.0, 4.0, 1.0, 4.0, 9.0, 1.0, 1.0, 2.0, 3.0, 7.0, 6.0, 9.0, 6.0, 3.0, 7.0, 6.0, 10.0, 1.0, 9.0, 6.0, 2.0, 4.0, 4.0, 5.0, 9.0, 2.0, 2.0, 2.0, 3.0, 7.0, 7.0, 6.0, 1.0, 6.0, 6.0, 1.0, 10.0, 10.0, 7.0, 9.0, 4.0, 6.0, 5.0, 8.0, 1.0, 1.0, 10.0, 4.0, 1.0, 4.0]
global b_y = 10
global p = [0.682, 0.174, 0.904, 0.744, 0.558, 0.367, 0.329, 0.895, 0.384, 0.259, 0.983, 0.001, 0.344, 0.908, 0.502, 0.32, 0.205, 0.81, 0.179, 0.86, 0.593, 0.924, 0.285, 0.406, 0.7, 0.363, 0.229, 0.466, 0.015, 0.694, 0.299, 0.385, 0.207, 0.681, 0.827, 0.407, 0.443, 0.749, 0.662, 0.569, 0.771, 0.266, 0.574, 0.009, 0.664, 0.814, 0.696, 0.61, 0.162, 0.434, 0.097, 0.242, 0.92, 0.494, 0.259, 0.078, 0.327, 0.404, 0.251, 0.166, 0.81, 0.96, 0.067, 0.795, 0.553, 0.364, 0.139, 0.72, 0.916, 0.498, 0.731, 0.195, 0.639, 0.87, 0.5, 0.133, 0.452, 0.618, 0.423, 0.076, 0.435, 0.744, 0.991, 0.186, 0.295, 0.688, 0.823, 0.38, 0.266, 0.676, 0.882, 0.425, 0.968, 0.156, 0.601, 0.063, 0.368, 0.309, 0.92, 0.97, 0.911, 0.431, 0.124, 0.606, 0.431, 0.65, 0.829, 0.191, 0.306, 0.245, 0.389, 0.428, 0.648, 0.039, 0.493, 0.056, 0.512, 0.758, 0.192, 0.92, 0.715, 0.654, 0.311, 0.488, 0.433, 0.083, 0.88, 0.887, 0.543, 0.87, 0.36, 0.914, 0.146, 0.216, 0.975, 0.979, 0.491, 0.995, 0.366, 0.589, 0.157, 0.989, 0.268, 0.027, 0.566, 0.539, 0.035, 0.66, 0.653, 0.348, 0.019, 0.604, 0.396, 0.421, 0.17, 0.083, 0.908, 0.247, 0.192, 0.3, 0.194, 0.93, 0.288, 0.716, 0.722, 0.123, 0.246, 0.012, 0.514, 0.79, 0.577, 0.102, 0.689, 0.438, 0.39, 0.825, 0.039, 0.481, 0.909, 0.384, 0.163, 0.637, 0.58, 0.768, 0.914, 0.672, 0.837, 0.455, 0.013, 0.893, 0.529, 0.933, 0.457, 0.156, 0.555, 0.91, 0.57, 0.097, 0.867, 0.931, 0.756, 0.864, 0.456, 0.437, 0.616, 0.732, 0.654, 0.849, 0.259, 0.739, 0.46, 0.993, 0.544, 0.632, 0.863, 0.575, 0.986, 0.585, 0.546, 0.694, 0.005, 0.913, 0.645, 0.9, 0.667, 0.595, 0.129, 0.112, 0.785, 0.41, 0.424, 0.566]
global q = [0.748, 0.235, 0.923, 0.874, 0.87, 0.674, 0.81, 0.994, 0.469, 0.312, 0.993, 0.171, 0.951, 0.972, 0.833, 0.338, 0.506, 0.864, 0.984, 0.989, 0.744, 0.939, 0.866, 0.994, 0.926, 0.752, 0.717, 0.899, 0.043, 0.735, 0.672, 0.476, 0.807, 0.728, 0.999, 0.449, 0.705, 0.786, 0.736, 0.594, 0.925, 0.672, 0.909, 0.678, 0.809, 0.986, 0.995, 0.85, 0.748, 0.743, 0.274, 0.678, 0.936, 0.578, 0.293, 0.29, 0.578, 0.807, 0.505, 0.952, 0.961, 0.988, 0.377, 0.933, 0.844, 0.874, 0.476, 0.857, 0.982, 0.499, 0.795, 0.834, 0.962, 0.972, 0.629, 0.881, 0.924, 0.958, 0.741, 0.23, 0.716, 0.945, 0.992, 0.45, 0.562, 0.947, 0.941, 0.55, 0.389, 0.836, 0.969, 0.908, 0.984, 0.382, 0.649, 0.171, 0.468, 0.402, 0.986, 0.997, 0.953, 0.854, 0.974, 0.802, 0.982, 0.76, 0.939, 0.705, 0.487, 0.501, 0.915, 0.623, 0.878, 0.45, 0.743, 0.592, 0.987, 0.981, 0.325, 0.984, 0.856, 0.837, 0.426, 0.545, 0.474, 0.212, 0.971, 0.983, 0.724, 0.943, 0.55, 0.927, 0.673, 0.324, 0.99, 0.991, 0.687, 0.996, 0.566, 0.825, 0.221, 0.997, 0.993, 0.752, 0.857, 0.721, 0.991, 0.886, 0.943, 0.802, 0.102, 0.75, 0.952, 0.924, 0.408, 0.93, 0.931, 0.68, 0.729, 0.633, 0.833, 0.94, 0.319, 0.926, 0.743, 0.444, 0.558, 0.286, 0.8, 0.966, 0.7, 0.626, 0.705, 0.78, 0.742, 0.88, 0.231, 0.898, 0.945, 0.612, 0.829, 0.644, 0.668, 0.779, 0.932, 0.928, 0.901, 0.867, 0.39, 0.933, 0.644, 0.95, 0.845, 0.399, 0.938, 0.994, 0.63, 0.922, 0.899, 0.965, 0.865, 0.981, 0.672, 0.992, 0.737, 0.835, 0.981, 0.89, 0.43, 0.906, 0.731, 0.998, 0.834, 0.834, 0.96, 0.916, 0.986, 0.754, 0.904, 0.967, 0.222, 0.97, 0.975, 0.966, 0.971, 0.66, 0.689, 0.124, 0.808, 0.453, 0.537, 0.763]
global origin = 1
global destination = 50