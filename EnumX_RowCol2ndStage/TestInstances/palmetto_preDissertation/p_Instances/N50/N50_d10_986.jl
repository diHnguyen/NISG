global arcs = [1 7; 1 8; 1 12; 1 13; 1 29; 1 34; 1 44; 1 48; 2 29; 2 38; 2 41; 2 49; 3 9; 3 26; 3 34; 3 37; 4 7; 4 15; 4 17; 4 20; 4 22; 4 28; 4 36; 4 39; 4 44; 5 11; 5 15; 5 41; 6 2; 6 9; 6 18; 6 28; 6 31; 6 39; 6 44; 7 3; 7 19; 7 28; 7 34; 7 37; 7 41; 8 12; 8 19; 8 28; 8 36; 8 37; 8 39; 8 50; 9 12; 9 17; 9 21; 9 33; 9 38; 9 39; 10 42; 11 26; 11 47; 11 48; 12 4; 12 11; 12 48; 13 6; 13 7; 13 12; 13 15; 13 34; 13 40; 14 11; 14 16; 14 21; 14 43; 14 47; 14 48; 15 7; 15 9; 15 18; 15 38; 16 9; 16 17; 16 24; 16 43; 16 49; 17 3; 17 35; 17 44; 18 7; 18 23; 18 26; 18 32; 19 2; 19 24; 20 4; 20 15; 20 17; 20 35; 21 14; 21 16; 21 45; 22 38; 22 39; 22 40; 22 48; 23 6; 23 35; 23 38; 23 45; 23 47; 24 4; 24 6; 24 19; 24 22; 24 43; 24 44; 24 47; 25 12; 25 16; 25 17; 25 19; 25 41; 25 43; 25 46; 26 11; 26 13; 26 24; 26 36; 27 4; 27 10; 27 17; 27 44; 28 4; 28 14; 28 37; 28 39; 28 41; 29 33; 29 46; 30 4; 30 7; 30 11; 30 12; 30 14; 30 19; 30 44; 30 46; 30 47; 31 3; 31 17; 31 24; 31 28; 31 39; 31 47; 32 8; 32 37; 33 2; 33 17; 33 20; 33 42; 34 23; 34 25; 35 18; 35 36; 35 46; 36 32; 36 40; 36 41; 36 46; 36 47; 37 5; 37 21; 37 28; 37 30; 37 45; 37 46; 38 8; 38 22; 38 31; 38 40; 38 49; 39 4; 39 12; 39 16; 39 29; 39 45; 39 48; 40 16; 40 19; 40 37; 40 49; 41 20; 41 24; 42 3; 42 4; 42 5; 42 17; 42 34; 42 45; 42 48; 43 10; 43 11; 43 15; 43 27; 43 37; 43 39; 43 44; 43 50; 44 6; 44 8; 44 9; 44 10; 44 32; 44 39; 45 20; 45 21; 45 31; 46 19; 47 45; 48 8; 48 32; 48 37; 48 38; 48 41; 48 49; 49 17; 49 40; 49 44]
global d_x = [7.0, 4.0, 2.0, 8.0, 1.0, 5.0, 1.0, 5.0, 5.0, 10.0, 5.0, 3.0, 1.0, 4.0, 9.0, 1.0, 10.0, 9.0, 10.0, 10.0, 6.0, 1.0, 5.0, 4.0, 4.0, 8.0, 6.0, 6.0, 6.0, 9.0, 7.0, 7.0, 6.0, 7.0, 2.0, 7.0, 5.0, 8.0, 5.0, 9.0, 4.0, 6.0, 3.0, 2.0, 2.0, 5.0, 8.0, 8.0, 4.0, 5.0, 2.0, 8.0, 1.0, 10.0, 5.0, 5.0, 2.0, 10.0, 3.0, 4.0, 8.0, 8.0, 9.0, 1.0, 3.0, 5.0, 4.0, 3.0, 9.0, 3.0, 4.0, 2.0, 6.0, 3.0, 7.0, 6.0, 6.0, 8.0, 3.0, 8.0, 8.0, 9.0, 5.0, 4.0, 8.0, 8.0, 10.0, 3.0, 6.0, 9.0, 7.0, 1.0, 5.0, 7.0, 10.0, 2.0, 4.0, 10.0, 2.0, 4.0, 2.0, 3.0, 5.0, 10.0, 1.0, 9.0, 7.0, 8.0, 9.0, 8.0, 10.0, 7.0, 9.0, 1.0, 10.0, 9.0, 8.0, 3.0, 7.0, 5.0, 6.0, 6.0, 10.0, 4.0, 5.0, 4.0, 4.0, 8.0, 2.0, 7.0, 9.0, 1.0, 10.0, 6.0, 10.0, 1.0, 6.0, 1.0, 10.0, 8.0, 2.0, 8.0, 6.0, 8.0, 9.0, 2.0, 10.0, 7.0, 7.0, 8.0, 6.0, 2.0, 8.0, 3.0, 10.0, 2.0, 7.0, 8.0, 7.0, 2.0, 10.0, 7.0, 5.0, 1.0, 7.0, 1.0, 10.0, 10.0, 3.0, 8.0, 10.0, 2.0, 4.0, 6.0, 9.0, 8.0, 1.0, 1.0, 2.0, 7.0, 5.0, 5.0, 1.0, 6.0, 5.0, 4.0, 6.0, 5.0, 4.0, 9.0, 6.0, 1.0, 6.0, 9.0, 2.0, 9.0, 8.0, 2.0, 1.0, 9.0, 9.0, 3.0, 6.0, 7.0, 4.0, 2.0, 7.0, 1.0, 7.0, 8.0, 3.0, 9.0, 9.0, 3.0, 9.0, 10.0, 6.0, 9.0, 1.0, 10.0, 2.0, 3.0, 8.0, 10.0, 4.0]
global b_x = 5
global d_y = [2.0, 7.0, 10.0, 2.0, 6.0, 2.0, 6.0, 6.0, 5.0, 4.0, 2.0, 7.0, 9.0, 3.0, 3.0, 7.0, 1.0, 5.0, 1.0, 3.0, 3.0, 1.0, 5.0, 10.0, 2.0, 4.0, 2.0, 7.0, 6.0, 1.0, 8.0, 7.0, 8.0, 1.0, 4.0, 3.0, 3.0, 1.0, 7.0, 10.0, 4.0, 1.0, 6.0, 4.0, 3.0, 8.0, 6.0, 6.0, 4.0, 3.0, 7.0, 8.0, 4.0, 4.0, 5.0, 4.0, 8.0, 9.0, 6.0, 10.0, 3.0, 3.0, 7.0, 8.0, 10.0, 4.0, 1.0, 2.0, 6.0, 4.0, 8.0, 6.0, 7.0, 2.0, 3.0, 5.0, 5.0, 9.0, 3.0, 2.0, 2.0, 4.0, 2.0, 5.0, 3.0, 4.0, 2.0, 7.0, 4.0, 9.0, 3.0, 4.0, 10.0, 7.0, 4.0, 10.0, 10.0, 9.0, 5.0, 5.0, 4.0, 1.0, 6.0, 9.0, 2.0, 5.0, 6.0, 4.0, 4.0, 9.0, 2.0, 3.0, 1.0, 1.0, 2.0, 7.0, 9.0, 5.0, 5.0, 5.0, 3.0, 10.0, 2.0, 8.0, 3.0, 4.0, 10.0, 5.0, 8.0, 10.0, 6.0, 5.0, 2.0, 5.0, 9.0, 1.0, 8.0, 5.0, 7.0, 7.0, 10.0, 7.0, 3.0, 4.0, 7.0, 3.0, 10.0, 7.0, 4.0, 2.0, 2.0, 9.0, 3.0, 10.0, 9.0, 4.0, 5.0, 7.0, 5.0, 2.0, 7.0, 10.0, 9.0, 5.0, 2.0, 8.0, 10.0, 10.0, 4.0, 3.0, 7.0, 6.0, 9.0, 9.0, 1.0, 5.0, 3.0, 8.0, 8.0, 8.0, 5.0, 2.0, 6.0, 2.0, 4.0, 9.0, 2.0, 9.0, 3.0, 4.0, 1.0, 4.0, 4.0, 8.0, 3.0, 8.0, 3.0, 3.0, 6.0, 10.0, 5.0, 3.0, 9.0, 9.0, 5.0, 1.0, 8.0, 9.0, 1.0, 8.0, 9.0, 5.0, 7.0, 4.0, 5.0, 8.0, 5.0, 2.0, 5.0, 2.0, 5.0, 4.0, 1.0, 8.0, 2.0]
global b_y = 10
global p = [0.163, 0.839, 0.533, 0.667, 0.692, 0.456, 0.459, 0.403, 0.622, 0.496, 0.474, 0.828, 0.387, 0.103, 0.274, 0.508, 0.732, 0.946, 0.963, 0.432, 0.163, 0.507, 0.157, 0.906, 0.044, 0.603, 0.749, 0.05, 0.417, 0.661, 0.854, 0.685, 0.846, 0.712, 0.193, 0.612, 0.2, 0.693, 0.429, 0.963, 0.245, 0.945, 0.004, 0.481, 0.992, 0.609, 0.133, 0.231, 0.433, 0.01, 0.255, 0.793, 0.037, 0.294, 0.385, 0.487, 0.347, 0.16, 0.324, 0.078, 0.21, 0.346, 0.142, 0.087, 0.525, 0.09, 0.715, 0.858, 0.114, 0.712, 0.537, 0.893, 0.285, 0.753, 0.426, 0.475, 0.702, 0.244, 0.678, 0.11, 0.777, 0.192, 0.063, 0.329, 0.232, 0.645, 0.635, 0.848, 0.908, 0.538, 0.6, 0.49, 0.857, 0.44, 0.786, 0.07, 0.768, 0.829, 0.483, 0.258, 0.15, 0.432, 0.233, 0.695, 0.556, 0.748, 0.068, 0.115, 0.341, 0.365, 0.116, 0.37, 0.164, 0.975, 0.158, 0.032, 0.324, 0.076, 0.124, 0.89, 0.598, 0.561, 0.717, 0.647, 0.013, 0.711, 0.632, 0.69, 0.317, 0.914, 0.399, 0.366, 0.293, 0.374, 0.454, 0.389, 0.943, 0.749, 0.257, 0.133, 0.147, 0.927, 0.726, 0.546, 0.532, 0.128, 0.232, 0.3, 0.237, 0.611, 0.112, 0.386, 0.567, 0.68, 0.739, 0.194, 0.358, 0.061, 0.047, 0.492, 0.67, 0.582, 0.724, 0.999, 0.993, 0.309, 0.486, 0.26, 0.319, 0.504, 0.868, 0.286, 0.118, 0.692, 0.907, 0.281, 0.111, 0.026, 0.69, 0.471, 0.451, 0.592, 0.088, 0.735, 0.585, 0.108, 0.991, 0.137, 0.655, 0.622, 0.358, 0.158, 0.969, 0.89, 0.762, 0.945, 0.598, 0.286, 0.969, 0.018, 0.826, 0.716, 0.321, 0.233, 0.391, 0.808, 0.888, 0.303, 0.089, 0.405, 0.609, 0.049, 0.45, 0.453, 0.482, 0.438, 0.3, 0.976, 0.858, 0.814, 0.216, 0.693, 0.919, 0.186, 0.055]
global q = [0.934, 0.937, 0.714, 0.699, 0.745, 0.8, 0.523, 0.473, 0.949, 0.824, 0.532, 0.927, 0.994, 0.324, 0.411, 0.788, 0.987, 0.989, 0.99, 0.461, 0.971, 0.583, 0.234, 0.947, 0.242, 0.923, 0.905, 0.904, 0.434, 0.892, 0.901, 0.824, 0.982, 0.762, 0.644, 0.829, 0.271, 0.738, 0.522, 0.973, 0.997, 0.971, 0.835, 0.652, 0.998, 0.855, 0.167, 0.759, 0.818, 0.966, 0.558, 0.833, 0.699, 0.471, 0.86, 0.607, 0.67, 0.702, 0.893, 0.852, 0.423, 0.857, 0.606, 0.386, 0.651, 0.582, 0.992, 0.911, 0.43, 0.989, 0.877, 0.907, 0.894, 0.901, 0.621, 0.937, 0.82, 0.369, 0.739, 0.994, 0.984, 0.942, 0.384, 0.492, 0.324, 0.7, 0.99, 0.887, 0.916, 0.859, 0.965, 0.88, 0.872, 0.866, 0.901, 0.309, 0.959, 0.933, 0.737, 0.543, 0.66, 0.602, 0.891, 0.887, 0.605, 0.872, 0.59, 0.242, 0.425, 0.844, 0.898, 0.596, 0.609, 0.977, 0.875, 0.107, 0.645, 0.384, 0.975, 0.992, 0.746, 0.933, 0.942, 0.772, 0.44, 0.903, 0.994, 0.943, 0.985, 0.928, 0.664, 0.787, 0.671, 0.989, 0.615, 0.952, 0.984, 0.962, 0.499, 0.445, 0.469, 0.993, 0.849, 0.723, 0.904, 0.471, 0.344, 0.791, 0.992, 0.879, 0.527, 0.551, 0.857, 0.903, 0.988, 0.593, 0.959, 0.35, 0.082, 0.595, 0.869, 0.784, 0.817, 0.999, 0.995, 0.802, 0.586, 0.91, 0.407, 0.798, 0.881, 0.588, 0.945, 0.855, 0.966, 0.637, 0.283, 0.937, 0.746, 0.889, 0.58, 0.96, 0.7, 0.989, 0.705, 0.748, 0.991, 0.68, 0.853, 0.847, 0.478, 0.337, 0.997, 0.918, 0.988, 0.975, 0.915, 0.562, 0.969, 0.545, 0.93, 0.727, 0.962, 0.83, 0.512, 0.953, 0.963, 0.31, 0.759, 0.526, 0.892, 0.299, 0.744, 0.464, 0.88, 0.644, 0.827, 0.99, 0.876, 0.998, 0.744, 0.701, 0.932, 0.263, 0.776]
global origin = 1
global destination = 50