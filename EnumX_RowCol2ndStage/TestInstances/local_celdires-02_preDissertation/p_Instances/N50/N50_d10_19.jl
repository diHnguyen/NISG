global arcs = [1 5; 1 8; 1 15; 1 16; 1 21; 1 31; 2 8; 2 26; 2 36; 2 42; 2 47; 2 48; 2 49; 3 17; 3 20; 3 25; 3 30; 3 32; 3 34; 4 10; 4 13; 4 25; 4 26; 4 35; 4 48; 5 13; 5 26; 6 4; 6 8; 6 27; 6 32; 6 41; 7 16; 7 20; 7 27; 7 33; 7 36; 7 38; 7 46; 8 11; 8 19; 8 28; 8 43; 8 47; 9 11; 9 16; 9 17; 9 21; 9 22; 9 28; 9 44; 10 7; 10 32; 10 44; 10 48; 11 3; 11 15; 11 20; 11 23; 11 28; 11 31; 11 40; 11 49; 12 2; 12 30; 12 36; 12 39; 13 9; 13 11; 13 29; 13 36; 13 46; 14 2; 14 3; 14 12; 14 19; 14 24; 14 36; 14 37; 14 38; 15 3; 15 7; 15 9; 15 11; 15 20; 15 22; 15 42; 15 44; 15 46; 16 9; 16 26; 16 29; 16 41; 17 9; 17 13; 17 30; 17 36; 17 42; 18 9; 18 33; 18 43; 19 10; 19 14; 19 25; 19 39; 20 13; 20 22; 20 30; 20 34; 20 40; 21 3; 21 11; 21 15; 21 22; 21 25; 21 38; 22 35; 22 45; 23 4; 23 28; 23 35; 23 43; 23 49; 24 17; 24 20; 24 28; 24 30; 25 6; 25 20; 25 23; 25 48; 26 11; 26 19; 26 34; 26 37; 26 41; 26 46; 27 2; 27 18; 27 19; 27 24; 27 30; 27 35; 27 37; 28 4; 28 42; 29 10; 29 16; 29 30; 29 34; 29 41; 29 46; 30 21; 30 27; 31 3; 31 13; 31 16; 31 20; 31 28; 32 6; 32 13; 32 14; 32 19; 32 28; 32 33; 32 39; 32 43; 33 12; 33 18; 33 35; 34 12; 34 17; 34 19; 34 32; 34 35; 35 12; 35 36; 35 48; 36 2; 36 13; 36 30; 36 34; 37 7; 37 12; 37 26; 37 50; 38 13; 38 19; 38 25; 38 30; 38 35; 38 39; 38 41; 38 44; 39 4; 39 7; 39 9; 39 22; 39 23; 39 42; 39 43; 39 45; 39 49; 40 3; 40 4; 40 12; 40 17; 40 20; 40 48; 40 49; 41 7; 41 19; 41 27; 41 34; 41 36; 41 49; 42 3; 42 35; 42 40; 43 18; 43 36; 43 41; 43 45; 43 47; 44 5; 44 32; 44 38; 44 41; 45 6; 45 32; 46 21; 46 35; 46 40; 46 47; 46 48; 47 37; 47 39; 47 46; 47 49; 48 8; 48 9; 48 14; 48 18; 48 23; 48 30; 48 35; 48 47; 49 5; 49 12; 49 14; 49 25; 49 35]
global d_x = [9.0, 9.0, 4.0, 6.0, 7.0, 9.0, 4.0, 5.0, 8.0, 10.0, 7.0, 3.0, 2.0, 8.0, 1.0, 5.0, 3.0, 3.0, 6.0, 10.0, 8.0, 4.0, 2.0, 6.0, 6.0, 2.0, 2.0, 9.0, 9.0, 4.0, 6.0, 5.0, 6.0, 2.0, 1.0, 8.0, 3.0, 4.0, 10.0, 6.0, 1.0, 9.0, 10.0, 5.0, 9.0, 8.0, 10.0, 1.0, 10.0, 5.0, 3.0, 7.0, 3.0, 1.0, 4.0, 8.0, 2.0, 8.0, 3.0, 10.0, 5.0, 9.0, 2.0, 1.0, 4.0, 7.0, 10.0, 9.0, 9.0, 1.0, 10.0, 7.0, 2.0, 1.0, 5.0, 2.0, 5.0, 7.0, 1.0, 5.0, 10.0, 5.0, 6.0, 3.0, 3.0, 9.0, 8.0, 10.0, 2.0, 8.0, 8.0, 2.0, 5.0, 10.0, 1.0, 8.0, 4.0, 3.0, 4.0, 10.0, 1.0, 1.0, 4.0, 7.0, 10.0, 3.0, 5.0, 3.0, 5.0, 5.0, 5.0, 1.0, 1.0, 10.0, 1.0, 6.0, 7.0, 3.0, 1.0, 2.0, 7.0, 2.0, 3.0, 10.0, 3.0, 8.0, 1.0, 1.0, 1.0, 1.0, 7.0, 8.0, 7.0, 5.0, 8.0, 2.0, 6.0, 2.0, 10.0, 3.0, 4.0, 5.0, 9.0, 7.0, 1.0, 1.0, 6.0, 7.0, 4.0, 4.0, 7.0, 3.0, 4.0, 3.0, 8.0, 6.0, 3.0, 4.0, 9.0, 2.0, 3.0, 7.0, 2.0, 1.0, 9.0, 9.0, 8.0, 6.0, 6.0, 1.0, 6.0, 6.0, 8.0, 8.0, 6.0, 5.0, 10.0, 2.0, 7.0, 3.0, 9.0, 6.0, 4.0, 3.0, 7.0, 7.0, 2.0, 10.0, 4.0, 7.0, 2.0, 8.0, 2.0, 9.0, 3.0, 3.0, 3.0, 6.0, 2.0, 1.0, 8.0, 4.0, 2.0, 7.0, 7.0, 9.0, 9.0, 7.0, 2.0, 10.0, 9.0, 8.0, 6.0, 7.0, 7.0, 9.0, 2.0, 2.0, 5.0, 2.0, 1.0, 9.0, 8.0, 7.0, 9.0, 9.0, 6.0, 7.0, 3.0, 9.0, 2.0, 6.0, 3.0, 1.0, 6.0, 3.0, 10.0, 10.0, 8.0, 6.0, 5.0, 1.0, 4.0, 4.0, 9.0, 8.0, 7.0, 1.0, 2.0, 9.0, 8.0, 9.0]
global b_x = 5
global d_y = [9.0, 9.0, 4.0, 6.0, 2.0, 8.0, 9.0, 5.0, 7.0, 7.0, 8.0, 10.0, 4.0, 3.0, 3.0, 4.0, 6.0, 9.0, 10.0, 8.0, 7.0, 8.0, 5.0, 5.0, 9.0, 9.0, 7.0, 9.0, 9.0, 5.0, 2.0, 7.0, 1.0, 5.0, 9.0, 7.0, 5.0, 6.0, 8.0, 2.0, 2.0, 3.0, 1.0, 5.0, 5.0, 10.0, 10.0, 5.0, 4.0, 1.0, 3.0, 2.0, 4.0, 3.0, 4.0, 6.0, 9.0, 5.0, 5.0, 1.0, 3.0, 5.0, 2.0, 9.0, 5.0, 5.0, 8.0, 8.0, 2.0, 9.0, 1.0, 9.0, 5.0, 1.0, 5.0, 4.0, 6.0, 1.0, 4.0, 7.0, 3.0, 8.0, 1.0, 9.0, 5.0, 7.0, 2.0, 9.0, 1.0, 7.0, 4.0, 2.0, 2.0, 6.0, 1.0, 6.0, 8.0, 6.0, 1.0, 9.0, 10.0, 4.0, 6.0, 10.0, 1.0, 10.0, 7.0, 5.0, 7.0, 3.0, 3.0, 5.0, 2.0, 1.0, 10.0, 5.0, 7.0, 3.0, 5.0, 1.0, 9.0, 10.0, 1.0, 7.0, 6.0, 5.0, 3.0, 5.0, 6.0, 10.0, 5.0, 3.0, 8.0, 1.0, 1.0, 1.0, 4.0, 5.0, 5.0, 3.0, 4.0, 9.0, 2.0, 2.0, 9.0, 1.0, 9.0, 6.0, 10.0, 7.0, 3.0, 8.0, 7.0, 7.0, 7.0, 4.0, 10.0, 4.0, 2.0, 9.0, 6.0, 4.0, 7.0, 1.0, 3.0, 4.0, 1.0, 10.0, 6.0, 8.0, 8.0, 1.0, 8.0, 8.0, 2.0, 7.0, 8.0, 9.0, 2.0, 9.0, 10.0, 5.0, 7.0, 9.0, 3.0, 9.0, 1.0, 5.0, 8.0, 8.0, 10.0, 3.0, 10.0, 2.0, 5.0, 5.0, 7.0, 4.0, 1.0, 2.0, 9.0, 7.0, 3.0, 4.0, 1.0, 6.0, 6.0, 5.0, 9.0, 7.0, 10.0, 4.0, 3.0, 9.0, 2.0, 2.0, 1.0, 3.0, 10.0, 10.0, 9.0, 2.0, 7.0, 7.0, 9.0, 5.0, 2.0, 4.0, 6.0, 10.0, 6.0, 7.0, 4.0, 3.0, 6.0, 6.0, 2.0, 10.0, 9.0, 7.0, 2.0, 5.0, 10.0, 8.0, 8.0, 1.0, 8.0, 6.0, 2.0, 10.0, 6.0, 1.0]
global b_y = 10
global p = [0.636, 0.703, 0.794, 0.162, 0.442, 0.948, 0.253, 0.47, 0.18, 0.435, 0.778, 0.996, 0.245, 0.568, 0.306, 0.866, 0.794, 0.784, 0.851, 0.373, 0.209, 0.296, 0.515, 0.113, 0.455, 0.06, 0.412, 0.94, 0.801, 0.01, 0.66, 0.09, 0.44, 0.394, 0.002, 0.232, 0.846, 0.231, 0.79, 0.788, 0.094, 0.813, 0.502, 0.337, 0.113, 0.744, 0.696, 0.067, 0.651, 0.496, 0.283, 0.477, 0.677, 0.109, 0.551, 0.741, 0.942, 0.279, 0.234, 0.742, 0.572, 0.13, 0.842, 0.114, 0.415, 0.674, 0.438, 0.756, 0.536, 0.842, 0.732, 0.079, 0.984, 0.328, 0.949, 0.806, 0.516, 0.633, 0.176, 0.764, 0.383, 0.213, 0.004, 0.656, 0.448, 0.927, 0.168, 0.372, 0.266, 0.997, 0.161, 0.892, 0.491, 0.836, 0.862, 0.818, 0.586, 0.062, 0.489, 0.479, 0.714, 0.403, 0.38, 0.295, 0.407, 0.538, 0.875, 0.814, 0.129, 0.962, 0.839, 0.929, 0.765, 0.126, 0.021, 0.026, 0.213, 0.258, 0.181, 0.948, 0.861, 0.018, 0.604, 0.687, 0.389, 0.727, 0.759, 0.394, 0.548, 0.309, 0.404, 0.403, 0.936, 0.004, 0.173, 0.846, 0.87, 0.086, 0.059, 0.847, 0.57, 0.024, 0.803, 0.127, 0.038, 0.066, 0.759, 0.279, 0.905, 0.347, 0.075, 0.206, 0.908, 0.409, 0.839, 0.67, 0.257, 0.761, 0.003, 0.195, 0.612, 0.391, 0.417, 0.852, 0.09, 0.759, 0.076, 0.228, 0.11, 0.018, 0.707, 0.759, 0.097, 0.63, 0.933, 0.373, 0.064, 0.678, 0.883, 0.642, 0.846, 0.938, 0.655, 0.444, 0.199, 0.406, 0.603, 0.855, 0.747, 0.392, 0.959, 0.599, 0.179, 0.34, 0.359, 0.497, 0.129, 0.147, 0.784, 0.706, 0.366, 0.342, 0.829, 0.887, 0.322, 0.53, 0.918, 0.913, 0.272, 0.397, 0.363, 0.676, 0.65, 0.935, 0.553, 0.668, 0.874, 0.673, 0.627, 0.545, 0.485, 0.846, 0.305, 0.112, 0.287, 0.985, 0.161, 0.956, 0.96, 0.74, 0.074, 0.443, 0.545, 0.356, 0.592, 0.145, 0.024, 0.222, 0.146, 0.563, 0.79, 0.648, 0.808, 0.299, 0.901, 0.671, 0.232, 0.216, 0.988, 0.935, 0.557, 0.164]
global q = [0.782, 0.824, 0.811, 0.597, 0.953, 0.991, 0.391, 0.514, 0.555, 0.85, 0.801, 0.999, 0.535, 0.899, 0.613, 0.947, 0.833, 0.812, 0.892, 0.715, 0.586, 0.438, 0.59, 0.567, 0.956, 0.796, 0.503, 0.976, 0.962, 0.817, 0.817, 0.223, 0.805, 0.84, 0.018, 0.416, 0.981, 0.386, 0.925, 0.832, 0.955, 0.929, 0.522, 0.679, 0.335, 0.847, 0.918, 0.073, 0.971, 0.972, 0.304, 0.721, 0.953, 0.944, 0.98, 0.875, 0.945, 0.842, 0.815, 0.748, 0.717, 0.499, 0.989, 0.33, 0.768, 0.752, 0.817, 0.884, 0.83, 0.904, 0.749, 0.561, 0.989, 0.862, 0.979, 0.944, 0.799, 0.69, 0.741, 0.812, 0.735, 0.404, 0.765, 0.846, 0.78, 0.937, 0.699, 0.752, 0.374, 0.997, 0.97, 0.954, 0.773, 0.958, 0.899, 0.989, 0.732, 0.71, 0.675, 0.974, 0.966, 0.913, 0.999, 0.705, 0.727, 0.587, 0.98, 0.865, 0.71, 0.993, 0.839, 0.975, 0.799, 0.932, 0.807, 0.941, 0.259, 0.37, 0.603, 0.959, 0.943, 0.515, 0.638, 0.77, 0.931, 0.87, 0.958, 0.407, 0.65, 0.949, 0.917, 0.879, 0.988, 0.293, 0.26, 0.956, 0.883, 0.242, 0.069, 0.938, 0.775, 0.574, 0.819, 0.242, 0.77, 0.725, 0.845, 0.41, 0.971, 0.574, 0.34, 0.808, 0.914, 0.553, 0.98, 0.739, 0.627, 0.986, 0.889, 0.857, 0.716, 0.401, 0.528, 0.975, 0.83, 0.776, 0.107, 0.384, 0.974, 0.892, 0.963, 0.954, 0.949, 0.851, 0.971, 0.385, 0.739, 0.947, 0.912, 0.89, 0.981, 0.947, 0.859, 0.844, 0.814, 0.789, 0.693, 0.959, 0.927, 0.407, 0.964, 0.963, 0.685, 0.803, 0.679, 0.789, 0.663, 0.266, 0.997, 0.865, 0.402, 0.887, 0.889, 0.968, 0.497, 0.982, 0.924, 0.95, 0.919, 0.982, 0.966, 0.805, 0.651, 0.941, 0.693, 0.894, 0.893, 0.676, 0.811, 0.634, 0.888, 0.913, 0.6, 0.135, 0.317, 0.998, 0.816, 0.961, 0.983, 0.964, 0.658, 0.95, 0.835, 0.751, 0.824, 0.145, 0.921, 0.782, 0.928, 0.916, 0.964, 0.674, 0.883, 0.647, 0.974, 0.705, 0.555, 0.632, 0.993, 0.954, 0.816, 0.422]
global origin = 1
global destination = 50