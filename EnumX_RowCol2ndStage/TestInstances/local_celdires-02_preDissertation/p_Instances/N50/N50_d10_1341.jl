global arcs = [1 8; 1 19; 1 29; 1 33; 1 41; 1 46; 2 7; 2 9; 2 19; 2 43; 2 46; 3 31; 3 32; 3 42; 4 2; 4 10; 4 27; 4 34; 5 3; 5 10; 5 27; 5 31; 5 32; 5 50; 6 3; 6 17; 6 20; 6 27; 6 28; 6 32; 6 41; 7 3; 7 9; 7 11; 7 34; 7 37; 7 46; 7 47; 7 50; 8 2; 8 5; 8 7; 8 11; 8 14; 8 16; 8 17; 8 26; 8 38; 8 42; 8 43; 9 7; 9 12; 9 30; 9 37; 9 50; 10 17; 10 23; 10 26; 10 30; 11 3; 11 17; 11 29; 11 31; 11 33; 11 35; 12 5; 12 11; 12 21; 12 28; 12 50; 13 11; 13 12; 13 17; 13 31; 13 35; 13 36; 14 5; 14 15; 14 23; 14 38; 14 43; 14 45; 14 48; 14 49; 15 7; 15 16; 15 20; 15 30; 15 34; 15 36; 15 44; 15 46; 16 2; 16 15; 16 27; 16 28; 16 30; 16 37; 16 46; 16 49; 16 50; 17 13; 17 44; 17 50; 18 15; 18 42; 18 43; 18 46; 19 32; 19 47; 20 2; 20 3; 20 18; 20 27; 20 38; 21 4; 21 8; 21 11; 21 15; 21 36; 22 49; 23 11; 23 30; 24 25; 24 48; 24 49; 25 20; 25 26; 25 29; 25 43; 25 44; 26 2; 26 4; 26 12; 26 37; 27 5; 27 15; 27 16; 27 18; 27 30; 27 44; 27 47; 28 7; 28 13; 28 26; 28 29; 28 31; 28 40; 28 47; 28 49; 29 2; 29 15; 29 26; 29 28; 29 34; 29 42; 29 44; 30 5; 30 7; 30 17; 30 23; 30 34; 31 2; 31 10; 31 18; 31 20; 31 23; 31 27; 31 35; 32 2; 32 27; 32 37; 32 41; 33 4; 33 15; 33 19; 33 26; 33 30; 33 39; 34 4; 34 28; 34 41; 34 44; 34 49; 35 10; 35 18; 35 20; 35 24; 35 34; 35 39; 35 44; 35 45; 36 3; 36 5; 36 29; 36 46; 37 4; 37 8; 37 22; 37 41; 37 44; 38 25; 38 32; 38 41; 39 13; 39 17; 39 21; 39 32; 40 20; 40 24; 40 27; 40 34; 40 45; 41 20; 41 46; 42 2; 42 4; 42 16; 42 31; 42 46; 42 47; 42 48; 43 6; 43 9; 43 12; 43 33; 44 2; 44 17; 44 19; 44 21; 44 26; 44 30; 44 34; 45 3; 45 11; 45 21; 45 28; 45 42; 46 4; 46 11; 46 30; 46 36; 46 40; 46 41; 46 48; 47 2; 47 8; 47 15; 47 24; 47 36; 47 46; 48 3; 48 7; 48 25; 48 27; 48 35; 48 45; 49 14; 49 25; 49 39; 49 40]
global d_x = [6.0, 7.0, 3.0, 5.0, 7.0, 7.0, 5.0, 6.0, 9.0, 6.0, 4.0, 4.0, 3.0, 6.0, 6.0, 1.0, 1.0, 4.0, 5.0, 9.0, 4.0, 6.0, 2.0, 3.0, 3.0, 9.0, 2.0, 10.0, 7.0, 5.0, 8.0, 10.0, 8.0, 7.0, 3.0, 4.0, 3.0, 6.0, 3.0, 9.0, 8.0, 9.0, 5.0, 5.0, 1.0, 8.0, 5.0, 7.0, 6.0, 5.0, 8.0, 1.0, 1.0, 10.0, 3.0, 4.0, 5.0, 8.0, 10.0, 2.0, 3.0, 2.0, 5.0, 10.0, 9.0, 4.0, 4.0, 8.0, 4.0, 8.0, 9.0, 10.0, 10.0, 1.0, 9.0, 7.0, 8.0, 6.0, 4.0, 3.0, 10.0, 10.0, 7.0, 2.0, 6.0, 7.0, 3.0, 8.0, 10.0, 6.0, 3.0, 8.0, 7.0, 5.0, 10.0, 9.0, 3.0, 10.0, 6.0, 4.0, 7.0, 1.0, 1.0, 6.0, 7.0, 5.0, 2.0, 9.0, 3.0, 7.0, 10.0, 4.0, 1.0, 4.0, 3.0, 10.0, 9.0, 4.0, 8.0, 2.0, 8.0, 8.0, 4.0, 7.0, 7.0, 1.0, 3.0, 8.0, 1.0, 10.0, 1.0, 2.0, 10.0, 5.0, 7.0, 5.0, 6.0, 5.0, 7.0, 8.0, 1.0, 4.0, 1.0, 6.0, 6.0, 9.0, 4.0, 10.0, 4.0, 9.0, 3.0, 2.0, 4.0, 10.0, 6.0, 5.0, 10.0, 7.0, 3.0, 5.0, 10.0, 9.0, 1.0, 2.0, 7.0, 5.0, 7.0, 3.0, 5.0, 8.0, 10.0, 8.0, 2.0, 10.0, 6.0, 9.0, 9.0, 10.0, 9.0, 8.0, 5.0, 4.0, 4.0, 9.0, 8.0, 4.0, 8.0, 7.0, 5.0, 1.0, 6.0, 10.0, 9.0, 3.0, 10.0, 5.0, 1.0, 10.0, 8.0, 7.0, 2.0, 10.0, 4.0, 6.0, 3.0, 8.0, 5.0, 10.0, 5.0, 2.0, 7.0, 10.0, 6.0, 2.0, 2.0, 9.0, 10.0, 9.0, 6.0, 4.0, 6.0, 6.0, 4.0, 2.0, 2.0, 8.0, 10.0, 7.0, 5.0, 8.0, 3.0, 5.0, 7.0, 2.0, 9.0, 5.0, 9.0, 8.0, 1.0, 1.0, 10.0, 5.0, 9.0, 1.0, 6.0, 4.0, 7.0, 1.0, 9.0, 7.0, 8.0, 8.0, 4.0, 6.0, 9.0, 2.0, 8.0, 3.0, 2.0, 1.0, 10.0]
global b_x = 5
global d_y = [2.0, 7.0, 7.0, 8.0, 2.0, 9.0, 2.0, 7.0, 7.0, 5.0, 10.0, 2.0, 8.0, 2.0, 7.0, 10.0, 8.0, 9.0, 5.0, 1.0, 6.0, 4.0, 4.0, 3.0, 8.0, 10.0, 7.0, 3.0, 6.0, 3.0, 10.0, 6.0, 2.0, 4.0, 2.0, 5.0, 6.0, 5.0, 3.0, 4.0, 9.0, 3.0, 8.0, 3.0, 3.0, 7.0, 8.0, 10.0, 10.0, 9.0, 5.0, 5.0, 8.0, 9.0, 1.0, 6.0, 3.0, 6.0, 7.0, 7.0, 4.0, 1.0, 9.0, 5.0, 4.0, 9.0, 1.0, 10.0, 4.0, 8.0, 1.0, 5.0, 5.0, 5.0, 3.0, 7.0, 3.0, 9.0, 2.0, 2.0, 3.0, 10.0, 10.0, 9.0, 8.0, 5.0, 7.0, 9.0, 7.0, 8.0, 5.0, 3.0, 2.0, 9.0, 2.0, 10.0, 10.0, 4.0, 5.0, 5.0, 6.0, 8.0, 7.0, 3.0, 6.0, 10.0, 4.0, 10.0, 5.0, 4.0, 1.0, 6.0, 8.0, 4.0, 4.0, 8.0, 1.0, 2.0, 8.0, 4.0, 1.0, 5.0, 7.0, 9.0, 3.0, 9.0, 4.0, 7.0, 5.0, 5.0, 2.0, 4.0, 7.0, 1.0, 5.0, 4.0, 1.0, 3.0, 7.0, 6.0, 6.0, 6.0, 2.0, 2.0, 3.0, 6.0, 9.0, 3.0, 4.0, 6.0, 10.0, 4.0, 2.0, 6.0, 7.0, 3.0, 7.0, 8.0, 9.0, 4.0, 9.0, 10.0, 1.0, 5.0, 4.0, 8.0, 3.0, 7.0, 8.0, 4.0, 5.0, 7.0, 7.0, 4.0, 8.0, 4.0, 1.0, 2.0, 3.0, 2.0, 1.0, 6.0, 4.0, 4.0, 4.0, 5.0, 5.0, 4.0, 6.0, 9.0, 5.0, 3.0, 7.0, 2.0, 3.0, 3.0, 6.0, 5.0, 10.0, 1.0, 2.0, 7.0, 4.0, 3.0, 3.0, 8.0, 4.0, 6.0, 2.0, 7.0, 10.0, 10.0, 7.0, 10.0, 8.0, 1.0, 3.0, 1.0, 8.0, 3.0, 2.0, 1.0, 1.0, 5.0, 5.0, 2.0, 10.0, 9.0, 5.0, 4.0, 7.0, 1.0, 10.0, 9.0, 4.0, 8.0, 6.0, 10.0, 2.0, 2.0, 7.0, 6.0, 2.0, 8.0, 7.0, 6.0, 5.0, 6.0, 2.0, 1.0, 6.0, 10.0, 4.0, 2.0, 7.0, 8.0, 7.0, 5.0, 7.0, 9.0, 5.0]
global b_y = 10
global p = [0.2, 0.271, 0.936, 0.788, 0.287, 0.492, 0.077, 0.777, 0.801, 0.07, 0.35, 0.767, 0.042, 0.976, 0.63, 0.741, 0.478, 0.938, 0.337, 0.1, 0.741, 0.213, 0.016, 0.171, 0.847, 0.938, 0.637, 0.082, 0.548, 0.634, 0.488, 0.096, 0.508, 0.493, 0.884, 0.04, 0.166, 0.351, 0.201, 0.31, 0.628, 0.972, 0.102, 0.672, 0.501, 0.788, 0.389, 0.372, 0.003, 0.461, 0.087, 0.619, 0.182, 0.715, 0.921, 0.144, 0.744, 0.282, 0.696, 0.645, 0.76, 0.224, 0.65, 0.114, 0.779, 0.503, 0.136, 0.061, 0.265, 0.483, 0.987, 0.489, 0.397, 0.54, 0.787, 0.457, 0.464, 0.927, 0.865, 0.899, 0.225, 0.352, 0.855, 0.132, 0.785, 0.674, 0.313, 0.832, 0.511, 0.32, 0.023, 0.882, 0.564, 0.562, 0.16, 0.28, 0.594, 0.661, 0.058, 0.355, 0.041, 0.838, 0.545, 0.652, 0.746, 0.152, 0.018, 0.63, 0.861, 0.233, 0.263, 0.917, 0.808, 0.311, 0.351, 0.929, 0.995, 0.912, 0.087, 0.609, 0.08, 0.647, 0.805, 0.471, 0.392, 0.316, 0.954, 0.802, 0.878, 0.69, 0.406, 0.146, 0.181, 0.97, 0.004, 0.506, 0.632, 0.062, 0.723, 0.2, 0.152, 0.749, 0.882, 0.657, 0.907, 0.308, 0.915, 0.488, 0.382, 0.857, 0.884, 0.424, 0.259, 0.323, 0.802, 0.218, 0.295, 0.772, 0.383, 0.904, 0.888, 0.834, 0.854, 0.783, 0.275, 0.579, 0.755, 0.261, 0.429, 0.905, 0.082, 0.408, 0.148, 0.764, 0.705, 0.737, 0.119, 0.08, 0.24, 0.628, 0.362, 0.479, 0.109, 0.188, 0.583, 0.21, 0.676, 0.477, 0.355, 0.449, 0.261, 0.084, 0.14, 0.813, 0.178, 0.224, 0.646, 0.633, 0.263, 0.197, 0.73, 0.213, 0.358, 0.559, 0.756, 0.468, 0.508, 0.841, 0.275, 0.088, 0.914, 0.213, 0.777, 0.09, 0.248, 0.99, 0.309, 0.016, 0.841, 0.147, 0.925, 0.029, 0.809, 0.67, 0.181, 0.809, 0.968, 0.193, 0.815, 0.097, 0.621, 0.245, 0.77, 0.577, 0.538, 0.229, 0.267, 0.672, 0.624, 0.445, 0.676, 0.64, 0.846, 0.167, 0.263, 0.454, 0.031, 0.809, 0.337, 0.242, 0.494, 0.177, 0.83, 0.235, 0.771, 0.749, 0.883, 0.238, 0.779, 0.581, 0.512]
global q = [0.872, 0.782, 0.966, 0.899, 0.643, 0.937, 0.798, 0.887, 0.892, 0.193, 0.876, 0.988, 0.619, 0.993, 0.78, 0.997, 0.713, 0.962, 0.744, 0.69, 0.885, 0.519, 0.236, 0.684, 0.88, 0.968, 0.862, 0.123, 0.769, 0.917, 0.692, 0.681, 0.717, 0.777, 0.895, 0.704, 0.768, 0.457, 0.865, 0.8, 0.709, 0.993, 0.457, 0.803, 0.795, 0.912, 0.917, 0.564, 0.463, 0.774, 0.941, 0.744, 0.76, 0.827, 0.938, 0.985, 0.957, 0.617, 0.973, 0.705, 0.81, 0.53, 0.874, 0.9, 0.806, 0.7, 0.427, 0.908, 0.955, 0.625, 0.989, 0.69, 0.538, 0.665, 0.839, 0.72, 0.541, 0.971, 0.984, 0.927, 0.274, 0.693, 0.902, 0.782, 0.856, 0.995, 0.954, 0.921, 0.557, 0.726, 0.102, 0.998, 0.645, 0.658, 0.595, 0.934, 0.742, 0.896, 0.575, 0.753, 0.12, 0.871, 0.795, 0.764, 0.94, 0.317, 0.745, 0.849, 0.928, 0.503, 0.311, 0.997, 0.868, 0.659, 0.747, 0.965, 0.997, 0.99, 0.239, 0.904, 0.977, 0.69, 0.962, 0.531, 0.557, 0.535, 0.96, 0.984, 0.936, 0.714, 0.922, 0.439, 0.953, 0.995, 0.079, 0.729, 0.71, 0.167, 0.866, 0.649, 0.216, 0.758, 0.971, 0.9, 0.944, 0.802, 0.953, 0.777, 0.916, 0.984, 0.902, 0.803, 0.607, 0.702, 0.83, 0.383, 0.363, 0.948, 0.611, 0.905, 0.893, 0.863, 0.982, 0.968, 0.77, 0.795, 0.758, 0.61, 0.59, 0.908, 0.708, 0.5, 0.728, 0.982, 0.973, 0.889, 0.575, 0.563, 0.52, 0.88, 0.905, 0.782, 0.945, 0.581, 0.694, 0.582, 0.978, 0.7, 0.98, 0.63, 0.971, 0.575, 0.171, 0.857, 0.874, 0.481, 0.821, 0.851, 0.465, 0.759, 0.831, 0.924, 0.478, 0.751, 0.878, 0.78, 0.509, 0.978, 0.53, 0.954, 0.927, 0.564, 0.914, 0.39, 0.878, 0.996, 0.409, 0.926, 0.933, 0.152, 0.932, 0.403, 0.965, 0.813, 0.344, 0.877, 0.99, 0.611, 0.917, 0.954, 0.639, 0.953, 0.941, 0.657, 0.968, 0.238, 0.94, 0.674, 0.85, 0.919, 0.984, 0.841, 0.986, 0.31, 0.707, 0.689, 0.6, 0.915, 0.483, 0.282, 0.757, 0.285, 0.995, 0.278, 0.83, 0.956, 0.977, 0.777, 0.839, 0.959, 0.98]
global origin = 1
global destination = 50