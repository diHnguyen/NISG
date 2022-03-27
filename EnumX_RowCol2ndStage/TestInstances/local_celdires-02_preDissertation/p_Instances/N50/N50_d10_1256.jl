global arcs = [1 2; 1 9; 1 25; 1 30; 1 33; 1 35; 1 44; 1 46; 2 10; 2 12; 2 18; 2 23; 2 39; 2 41; 2 43; 2 46; 3 25; 3 31; 3 32; 3 48; 4 2; 4 12; 4 37; 4 38; 4 47; 4 49; 5 2; 5 8; 5 18; 5 27; 5 35; 5 47; 5 48; 6 13; 6 21; 6 28; 7 6; 7 16; 7 21; 8 3; 8 7; 8 40; 9 4; 9 8; 9 14; 9 18; 9 22; 9 27; 9 29; 9 34; 9 48; 10 21; 10 28; 10 33; 11 6; 11 9; 11 26; 12 4; 12 14; 12 22; 12 42; 12 44; 12 45; 13 11; 13 18; 13 22; 13 25; 14 8; 14 35; 14 49; 15 8; 15 9; 15 13; 15 24; 16 10; 16 22; 16 28; 16 39; 16 46; 16 49; 17 16; 17 20; 17 26; 17 27; 17 38; 17 45; 17 49; 18 7; 18 9; 18 12; 18 31; 18 44; 19 7; 19 37; 19 45; 19 48; 20 5; 20 15; 20 21; 20 23; 20 27; 20 36; 21 15; 21 18; 22 12; 22 39; 22 41; 22 42; 23 7; 23 8; 23 11; 23 35; 24 3; 24 7; 24 11; 24 29; 24 34; 24 38; 24 42; 24 48; 25 7; 25 8; 25 12; 25 15; 25 19; 25 23; 25 28; 25 31; 25 32; 25 46; 26 23; 26 39; 26 42; 26 44; 26 47; 27 14; 27 16; 27 23; 27 26; 27 28; 27 30; 27 34; 27 38; 27 50; 28 4; 28 10; 28 19; 28 25; 28 29; 28 35; 29 48; 30 6; 30 15; 30 20; 30 24; 30 35; 30 40; 30 48; 31 12; 31 16; 31 34; 31 46; 31 50; 32 7; 32 14; 32 19; 32 37; 32 46; 33 11; 33 18; 33 21; 33 34; 33 41; 33 42; 34 9; 34 12; 34 13; 34 30; 34 40; 34 49; 35 25; 35 34; 35 36; 35 44; 35 45; 36 3; 36 7; 36 35; 36 49; 37 21; 37 38; 38 13; 38 44; 38 50; 39 3; 39 13; 39 32; 40 19; 40 27; 40 34; 40 36; 41 22; 41 27; 41 42; 41 45; 42 6; 42 10; 42 21; 42 24; 42 27; 42 28; 42 38; 43 3; 43 20; 43 26; 43 33; 43 41; 44 11; 44 12; 44 45; 45 2; 45 9; 45 20; 45 23; 45 32; 45 43; 45 47; 46 6; 46 9; 46 22; 46 31; 46 36; 46 40; 47 14; 47 17; 47 29; 47 42; 48 5; 48 11; 48 22; 49 2; 49 3; 49 12; 49 28; 49 35]
global d_x = [5.0, 5.0, 7.0, 4.0, 2.0, 5.0, 3.0, 3.0, 9.0, 9.0, 1.0, 10.0, 3.0, 2.0, 3.0, 2.0, 10.0, 10.0, 6.0, 4.0, 1.0, 9.0, 7.0, 3.0, 5.0, 1.0, 6.0, 10.0, 4.0, 2.0, 4.0, 5.0, 1.0, 10.0, 9.0, 7.0, 9.0, 8.0, 6.0, 7.0, 10.0, 4.0, 7.0, 8.0, 1.0, 8.0, 6.0, 8.0, 3.0, 7.0, 8.0, 4.0, 7.0, 4.0, 3.0, 4.0, 5.0, 1.0, 3.0, 7.0, 9.0, 6.0, 6.0, 1.0, 4.0, 3.0, 5.0, 6.0, 1.0, 1.0, 9.0, 8.0, 1.0, 7.0, 4.0, 8.0, 10.0, 6.0, 3.0, 6.0, 2.0, 7.0, 10.0, 3.0, 9.0, 4.0, 3.0, 5.0, 2.0, 2.0, 4.0, 9.0, 2.0, 10.0, 10.0, 9.0, 1.0, 2.0, 1.0, 2.0, 2.0, 5.0, 1.0, 5.0, 1.0, 5.0, 6.0, 5.0, 9.0, 1.0, 6.0, 5.0, 10.0, 7.0, 7.0, 5.0, 3.0, 8.0, 5.0, 7.0, 10.0, 3.0, 10.0, 1.0, 5.0, 6.0, 4.0, 4.0, 10.0, 9.0, 2.0, 4.0, 3.0, 7.0, 2.0, 4.0, 8.0, 2.0, 1.0, 1.0, 2.0, 2.0, 8.0, 6.0, 1.0, 10.0, 2.0, 5.0, 5.0, 10.0, 7.0, 2.0, 7.0, 4.0, 3.0, 9.0, 7.0, 4.0, 9.0, 8.0, 8.0, 3.0, 7.0, 6.0, 3.0, 1.0, 7.0, 8.0, 2.0, 1.0, 6.0, 4.0, 5.0, 4.0, 2.0, 10.0, 6.0, 5.0, 9.0, 9.0, 1.0, 2.0, 7.0, 5.0, 8.0, 1.0, 5.0, 5.0, 9.0, 5.0, 2.0, 2.0, 4.0, 4.0, 7.0, 6.0, 1.0, 7.0, 10.0, 5.0, 6.0, 1.0, 1.0, 2.0, 4.0, 5.0, 6.0, 4.0, 7.0, 10.0, 10.0, 4.0, 7.0, 3.0, 6.0, 1.0, 2.0, 8.0, 1.0, 4.0, 1.0, 6.0, 4.0, 2.0, 3.0, 6.0, 3.0, 2.0, 8.0, 6.0, 10.0, 6.0, 3.0, 1.0, 9.0, 8.0, 5.0, 6.0, 4.0, 1.0, 6.0, 10.0, 10.0, 2.0, 8.0]
global b_x = 5
global d_y = [9.0, 5.0, 4.0, 5.0, 9.0, 10.0, 10.0, 10.0, 10.0, 3.0, 8.0, 5.0, 9.0, 7.0, 10.0, 2.0, 6.0, 3.0, 10.0, 6.0, 10.0, 1.0, 1.0, 10.0, 5.0, 10.0, 7.0, 5.0, 9.0, 4.0, 10.0, 4.0, 2.0, 1.0, 5.0, 9.0, 5.0, 10.0, 10.0, 2.0, 7.0, 6.0, 7.0, 3.0, 10.0, 1.0, 4.0, 4.0, 5.0, 9.0, 3.0, 7.0, 5.0, 7.0, 8.0, 5.0, 1.0, 5.0, 2.0, 6.0, 3.0, 3.0, 10.0, 2.0, 5.0, 6.0, 1.0, 4.0, 8.0, 4.0, 8.0, 3.0, 3.0, 3.0, 6.0, 6.0, 8.0, 3.0, 2.0, 5.0, 2.0, 6.0, 9.0, 1.0, 4.0, 7.0, 1.0, 2.0, 5.0, 3.0, 9.0, 1.0, 1.0, 5.0, 1.0, 7.0, 4.0, 3.0, 9.0, 5.0, 7.0, 7.0, 6.0, 3.0, 6.0, 3.0, 10.0, 1.0, 2.0, 6.0, 4.0, 2.0, 9.0, 5.0, 7.0, 8.0, 7.0, 1.0, 10.0, 6.0, 2.0, 2.0, 2.0, 6.0, 2.0, 5.0, 4.0, 7.0, 3.0, 9.0, 9.0, 7.0, 7.0, 6.0, 10.0, 9.0, 4.0, 4.0, 5.0, 2.0, 9.0, 1.0, 2.0, 3.0, 4.0, 8.0, 4.0, 5.0, 3.0, 9.0, 8.0, 5.0, 1.0, 7.0, 3.0, 3.0, 1.0, 8.0, 8.0, 1.0, 9.0, 9.0, 4.0, 3.0, 5.0, 8.0, 10.0, 8.0, 9.0, 10.0, 7.0, 2.0, 3.0, 9.0, 3.0, 3.0, 6.0, 8.0, 3.0, 5.0, 8.0, 6.0, 5.0, 4.0, 4.0, 4.0, 8.0, 1.0, 2.0, 3.0, 7.0, 2.0, 9.0, 2.0, 10.0, 3.0, 7.0, 4.0, 2.0, 1.0, 1.0, 9.0, 7.0, 9.0, 1.0, 7.0, 9.0, 3.0, 8.0, 2.0, 6.0, 5.0, 6.0, 5.0, 2.0, 10.0, 3.0, 3.0, 4.0, 9.0, 3.0, 3.0, 3.0, 7.0, 9.0, 10.0, 9.0, 10.0, 2.0, 7.0, 6.0, 1.0, 2.0, 2.0, 8.0, 10.0, 1.0, 7.0, 1.0, 7.0, 6.0, 6.0, 1.0, 3.0, 5.0]
global b_y = 10
global p = [0.819, 0.967, 0.911, 0.222, 0.464, 0.422, 0.331, 0.862, 0.893, 0.03, 0.281, 0.076, 0.545, 0.894, 0.93, 0.018, 0.412, 0.578, 0.388, 0.415, 0.505, 0.066, 0.885, 0.554, 0.373, 0.631, 0.649, 0.944, 0.897, 0.856, 0.775, 0.066, 0.3, 0.71, 0.13, 0.099, 0.657, 0.762, 0.149, 0.914, 0.866, 0.99, 0.628, 0.216, 0.398, 0.599, 0.814, 0.253, 0.932, 0.304, 0.252, 0.889, 0.576, 0.01, 0.29, 0.025, 0.216, 0.973, 0.943, 0.641, 0.226, 0.558, 0.175, 0.826, 0.898, 0.765, 0.046, 0.581, 0.686, 0.122, 0.277, 0.945, 0.027, 0.526, 0.479, 0.084, 0.01, 0.229, 0.132, 0.386, 0.557, 0.837, 0.974, 0.117, 0.579, 0.787, 0.294, 0.566, 0.353, 0.17, 0.037, 0.545, 0.463, 0.791, 0.605, 0.841, 0.001, 0.063, 0.664, 0.532, 0.299, 0.569, 0.034, 0.74, 0.512, 0.977, 0.159, 0.384, 0.612, 0.401, 0.216, 0.139, 0.34, 0.107, 0.031, 0.7, 0.418, 0.078, 0.099, 0.431, 0.186, 0.623, 0.632, 0.385, 0.106, 0.641, 0.502, 0.112, 0.008, 0.567, 0.126, 0.017, 0.916, 0.897, 0.712, 0.63, 0.32, 0.471, 0.686, 0.999, 0.189, 0.247, 0.993, 0.224, 0.72, 0.757, 0.687, 0.027, 0.132, 0.377, 0.042, 0.109, 0.458, 0.608, 0.349, 0.511, 0.248, 0.425, 0.461, 0.95, 0.523, 0.619, 0.24, 0.554, 0.613, 0.045, 0.563, 0.009, 0.22, 0.229, 0.493, 0.033, 0.618, 0.198, 0.248, 0.074, 0.886, 0.414, 0.679, 0.779, 0.474, 0.349, 0.204, 0.981, 0.876, 0.845, 0.026, 0.892, 0.715, 0.954, 0.435, 0.371, 0.798, 0.866, 0.354, 0.349, 0.675, 0.819, 0.24, 0.911, 0.32, 0.658, 0.665, 0.304, 0.916, 0.579, 0.624, 0.619, 0.125, 0.199, 0.315, 0.991, 0.487, 0.016, 0.96, 0.037, 0.508, 0.359, 0.381, 0.103, 0.68, 0.208, 0.156, 0.285, 0.722, 0.068, 0.942, 0.098, 0.087, 0.323, 0.003, 0.732, 0.151, 0.979, 0.257, 0.958, 0.435, 0.587, 0.521, 0.315, 0.423, 0.689, 0.187, 0.856, 0.8]
global q = [0.893, 0.969, 0.969, 0.59, 0.724, 0.991, 0.541, 0.939, 0.944, 0.738, 0.955, 0.486, 0.974, 0.983, 0.972, 0.147, 0.893, 0.974, 0.552, 0.516, 0.637, 0.98, 0.983, 0.802, 0.524, 0.721, 0.805, 0.97, 0.916, 0.869, 0.869, 0.492, 0.359, 0.765, 0.982, 0.42, 0.893, 0.962, 0.457, 0.987, 0.94, 0.99, 0.71, 0.57, 0.476, 0.978, 0.938, 0.57, 0.997, 0.65, 0.284, 0.97, 0.68, 0.468, 0.672, 0.269, 0.733, 0.979, 0.999, 0.898, 0.488, 0.637, 0.666, 0.86, 0.912, 0.889, 0.148, 0.74, 0.959, 0.344, 0.703, 0.975, 0.415, 0.665, 0.999, 0.459, 0.676, 0.258, 0.76, 0.679, 0.624, 0.971, 0.982, 0.417, 0.893, 0.96, 0.574, 0.894, 0.976, 0.537, 0.42, 0.991, 0.898, 0.895, 0.812, 0.882, 0.807, 0.211, 0.811, 0.533, 0.42, 0.66, 0.61, 0.944, 0.676, 0.977, 0.504, 0.968, 0.973, 0.586, 0.525, 0.21, 0.849, 0.582, 0.057, 0.951, 0.867, 0.484, 0.494, 0.626, 0.771, 0.907, 0.852, 0.588, 0.236, 0.903, 0.781, 0.641, 0.79, 0.86, 0.682, 0.765, 0.94, 0.92, 0.753, 0.982, 0.978, 0.937, 0.925, 0.999, 0.422, 0.347, 0.997, 0.946, 0.87, 0.861, 0.884, 0.078, 0.514, 0.906, 0.985, 0.256, 0.68, 0.809, 0.822, 0.784, 0.994, 0.924, 0.632, 0.985, 0.584, 0.709, 0.386, 0.886, 0.823, 0.232, 0.772, 0.527, 0.843, 0.454, 0.816, 0.188, 0.655, 0.523, 0.404, 0.721, 0.951, 0.617, 0.93, 0.998, 0.746, 0.354, 0.727, 0.982, 0.933, 0.934, 0.039, 0.938, 0.968, 0.967, 0.728, 0.602, 0.83, 0.885, 0.904, 0.527, 0.748, 0.859, 0.867, 0.921, 0.713, 0.967, 0.863, 0.942, 0.938, 0.776, 0.933, 0.635, 0.984, 0.475, 0.961, 0.998, 0.905, 0.931, 0.978, 0.891, 0.681, 0.923, 0.881, 0.827, 0.89, 0.488, 0.221, 0.348, 0.731, 0.159, 0.98, 0.982, 0.66, 0.693, 0.395, 0.973, 0.38, 0.981, 0.46, 0.968, 0.794, 0.856, 0.875, 0.722, 0.717, 0.984, 0.761, 0.883, 0.823]
global origin = 1
global destination = 50