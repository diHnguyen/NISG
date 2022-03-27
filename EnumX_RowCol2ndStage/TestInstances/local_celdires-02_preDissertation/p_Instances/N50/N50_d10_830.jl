global arcs = [1 3; 1 5; 1 15; 1 17; 1 23; 1 44; 2 4; 2 9; 2 14; 2 19; 2 22; 2 35; 2 47; 3 21; 3 31; 3 32; 3 35; 3 40; 3 46; 4 8; 4 21; 4 26; 4 27; 4 35; 4 49; 5 13; 5 23; 5 25; 5 30; 6 8; 6 9; 6 21; 6 35; 6 44; 6 48; 7 25; 7 43; 7 46; 8 2; 8 4; 8 5; 8 7; 8 23; 8 29; 8 31; 8 34; 8 39; 9 19; 9 40; 9 44; 10 7; 10 18; 10 21; 10 23; 10 27; 10 30; 11 6; 11 8; 11 9; 11 21; 11 30; 11 46; 12 15; 12 20; 12 47; 12 50; 13 3; 13 19; 13 27; 13 45; 14 34; 14 37; 14 38; 15 7; 15 18; 15 31; 16 14; 16 26; 16 29; 16 46; 17 2; 17 10; 17 22; 17 25; 17 37; 17 39; 18 13; 18 31; 18 36; 18 49; 19 11; 19 29; 19 48; 19 50; 20 11; 20 32; 20 35; 20 36; 21 13; 21 15; 21 34; 21 37; 21 45; 21 47; 22 3; 22 13; 22 20; 22 37; 23 22; 23 29; 23 37; 23 48; 23 50; 24 6; 24 13; 24 14; 24 33; 25 5; 25 14; 25 15; 25 20; 25 30; 25 31; 26 4; 26 16; 26 24; 26 27; 26 31; 26 38; 27 4; 27 11; 27 16; 27 19; 27 37; 27 45; 28 9; 28 25; 28 48; 29 21; 29 24; 29 25; 30 15; 30 38; 30 50; 31 7; 31 25; 31 45; 31 48; 31 49; 32 2; 32 12; 32 17; 32 46; 33 5; 33 16; 33 17; 34 16; 34 30; 34 37; 35 4; 35 20; 35 23; 35 27; 35 28; 35 31; 35 32; 35 50; 36 2; 36 13; 36 34; 37 7; 37 18; 37 25; 38 7; 38 26; 38 30; 38 46; 38 50; 39 4; 39 10; 39 13; 39 17; 39 31; 39 41; 39 50; 40 23; 40 45; 40 50; 41 15; 42 15; 42 24; 42 35; 42 37; 43 12; 43 20; 43 23; 43 44; 44 15; 44 26; 44 28; 44 34; 45 11; 45 16; 45 20; 45 34; 46 11; 46 14; 46 28; 46 39; 46 49; 47 42; 47 46; 48 8; 48 11; 48 14; 48 37; 48 45; 49 27; 49 37; 49 39; 49 47; 49 48]
global d_x = [10.0, 8.0, 2.0, 2.0, 2.0, 6.0, 5.0, 7.0, 6.0, 10.0, 8.0, 1.0, 8.0, 7.0, 8.0, 9.0, 4.0, 7.0, 10.0, 3.0, 6.0, 1.0, 1.0, 6.0, 6.0, 5.0, 6.0, 3.0, 10.0, 4.0, 1.0, 9.0, 6.0, 1.0, 2.0, 3.0, 3.0, 10.0, 2.0, 1.0, 10.0, 5.0, 6.0, 9.0, 9.0, 4.0, 6.0, 10.0, 10.0, 3.0, 5.0, 7.0, 7.0, 10.0, 2.0, 10.0, 10.0, 9.0, 10.0, 7.0, 6.0, 2.0, 7.0, 8.0, 1.0, 2.0, 5.0, 1.0, 4.0, 7.0, 7.0, 10.0, 5.0, 5.0, 8.0, 2.0, 8.0, 1.0, 9.0, 1.0, 3.0, 2.0, 6.0, 4.0, 3.0, 6.0, 1.0, 3.0, 8.0, 10.0, 1.0, 7.0, 5.0, 10.0, 8.0, 5.0, 1.0, 2.0, 4.0, 2.0, 3.0, 5.0, 4.0, 8.0, 5.0, 1.0, 4.0, 1.0, 7.0, 1.0, 9.0, 4.0, 4.0, 6.0, 6.0, 6.0, 8.0, 9.0, 6.0, 2.0, 9.0, 5.0, 7.0, 5.0, 3.0, 3.0, 8.0, 9.0, 2.0, 8.0, 4.0, 4.0, 7.0, 6.0, 5.0, 10.0, 4.0, 6.0, 5.0, 5.0, 9.0, 8.0, 1.0, 1.0, 5.0, 1.0, 9.0, 10.0, 8.0, 2.0, 7.0, 10.0, 10.0, 2.0, 9.0, 4.0, 8.0, 1.0, 8.0, 10.0, 1.0, 3.0, 2.0, 8.0, 6.0, 8.0, 10.0, 5.0, 9.0, 1.0, 4.0, 1.0, 1.0, 7.0, 6.0, 5.0, 2.0, 5.0, 2.0, 4.0, 6.0, 8.0, 7.0, 9.0, 3.0, 7.0, 3.0, 1.0, 3.0, 7.0, 5.0, 3.0, 3.0, 9.0, 4.0, 3.0, 4.0, 5.0, 9.0, 10.0, 1.0, 5.0, 5.0, 3.0, 7.0, 8.0, 7.0, 8.0, 2.0, 5.0, 7.0, 9.0, 8.0, 5.0, 8.0, 7.0, 6.0, 6.0, 5.0, 8.0, 6.0, 1.0]
global b_x = 5
global d_y = [6.0, 7.0, 7.0, 8.0, 4.0, 9.0, 7.0, 10.0, 7.0, 6.0, 3.0, 3.0, 7.0, 4.0, 8.0, 6.0, 7.0, 2.0, 3.0, 9.0, 3.0, 1.0, 4.0, 9.0, 2.0, 5.0, 1.0, 5.0, 4.0, 10.0, 5.0, 10.0, 10.0, 4.0, 6.0, 10.0, 9.0, 3.0, 5.0, 3.0, 1.0, 3.0, 1.0, 2.0, 4.0, 2.0, 4.0, 8.0, 4.0, 9.0, 1.0, 7.0, 6.0, 3.0, 7.0, 2.0, 3.0, 5.0, 2.0, 3.0, 7.0, 4.0, 8.0, 4.0, 5.0, 2.0, 3.0, 4.0, 1.0, 8.0, 9.0, 10.0, 5.0, 9.0, 6.0, 2.0, 9.0, 9.0, 6.0, 8.0, 9.0, 5.0, 10.0, 1.0, 8.0, 4.0, 2.0, 2.0, 7.0, 4.0, 7.0, 2.0, 8.0, 5.0, 7.0, 5.0, 4.0, 1.0, 8.0, 4.0, 7.0, 4.0, 7.0, 10.0, 1.0, 10.0, 5.0, 10.0, 8.0, 4.0, 9.0, 2.0, 3.0, 8.0, 5.0, 2.0, 2.0, 6.0, 2.0, 2.0, 1.0, 9.0, 8.0, 4.0, 7.0, 9.0, 7.0, 7.0, 10.0, 6.0, 5.0, 6.0, 4.0, 3.0, 7.0, 2.0, 2.0, 10.0, 3.0, 8.0, 5.0, 7.0, 10.0, 8.0, 8.0, 8.0, 7.0, 5.0, 2.0, 2.0, 6.0, 7.0, 9.0, 6.0, 8.0, 4.0, 2.0, 1.0, 3.0, 1.0, 3.0, 5.0, 10.0, 3.0, 5.0, 1.0, 7.0, 8.0, 9.0, 4.0, 3.0, 6.0, 2.0, 2.0, 4.0, 8.0, 2.0, 6.0, 1.0, 4.0, 7.0, 7.0, 9.0, 4.0, 9.0, 9.0, 7.0, 10.0, 5.0, 1.0, 5.0, 2.0, 4.0, 7.0, 7.0, 4.0, 9.0, 7.0, 2.0, 6.0, 1.0, 4.0, 10.0, 9.0, 7.0, 7.0, 5.0, 3.0, 6.0, 6.0, 10.0, 6.0, 4.0, 1.0, 10.0, 5.0, 4.0, 3.0, 9.0, 3.0, 5.0, 9.0]
global b_y = 10
global p = [0.143, 0.359, 0.745, 0.895, 0.33, 0.728, 0.237, 0.795, 0.624, 0.291, 0.786, 0.144, 0.922, 0.041, 0.845, 0.502, 0.989, 0.807, 0.756, 0.994, 0.52, 0.816, 0.52, 0.876, 0.676, 0.39, 0.132, 0.263, 0.95, 0.133, 0.63, 0.076, 0.191, 0.562, 0.295, 0.782, 0.939, 0.16, 0.197, 0.717, 0.86, 0.669, 0.01, 0.082, 0.352, 0.676, 0.452, 0.055, 0.548, 0.894, 0.232, 0.219, 0.491, 0.141, 0.081, 0.253, 0.501, 0.125, 0.623, 0.895, 0.045, 0.773, 0.754, 0.841, 0.033, 0.288, 0.395, 0.845, 0.904, 0.791, 0.967, 0.339, 0.635, 0.157, 0.049, 0.341, 0.856, 0.862, 0.626, 0.284, 0.397, 0.125, 0.231, 0.728, 0.412, 0.335, 0.623, 0.127, 0.654, 0.534, 0.486, 0.242, 0.785, 0.175, 0.069, 0.093, 0.692, 0.099, 0.26, 0.753, 0.823, 0.371, 0.566, 0.309, 0.776, 0.146, 0.951, 0.905, 0.773, 0.225, 0.855, 0.473, 0.548, 0.008, 0.657, 0.947, 0.335, 0.489, 0.113, 0.001, 0.616, 0.782, 0.832, 0.739, 0.135, 0.207, 0.384, 0.817, 0.666, 0.078, 0.762, 0.315, 0.997, 0.895, 0.415, 0.59, 0.074, 0.716, 0.253, 0.738, 0.673, 0.96, 0.707, 0.834, 0.592, 0.949, 0.978, 0.294, 0.018, 0.411, 0.161, 0.907, 0.984, 0.865, 0.114, 0.575, 0.627, 0.765, 0.954, 0.864, 0.771, 0.546, 0.033, 0.933, 0.642, 0.968, 0.198, 0.043, 0.584, 0.073, 0.398, 0.874, 0.449, 0.071, 0.203, 0.616, 0.914, 0.923, 0.944, 0.383, 0.814, 0.756, 0.366, 0.003, 0.665, 0.657, 0.507, 0.923, 0.261, 0.594, 0.33, 0.738, 0.629, 0.518, 0.671, 0.346, 0.942, 0.173, 0.3, 0.63, 0.907, 0.782, 0.938, 0.988, 0.945, 0.351, 0.641, 0.6, 0.232, 0.493, 0.717, 0.957, 0.979, 0.46, 0.612, 0.931, 0.981, 0.53, 0.406, 0.245, 0.763, 0.084]
global q = [0.269, 0.817, 0.952, 0.924, 0.628, 0.783, 0.96, 0.863, 0.828, 0.964, 0.82, 0.221, 0.954, 0.645, 0.914, 0.77, 0.997, 0.894, 0.995, 0.995, 0.679, 0.876, 0.74, 0.941, 0.822, 0.931, 0.711, 0.758, 0.968, 0.849, 0.861, 0.54, 0.958, 0.765, 0.669, 0.894, 0.944, 0.232, 0.64, 0.919, 0.998, 0.692, 0.037, 0.696, 0.943, 0.92, 0.545, 0.281, 0.566, 0.974, 0.774, 0.229, 0.764, 0.73, 0.349, 0.412, 0.544, 0.158, 0.89, 0.98, 0.353, 0.92, 0.85, 0.968, 0.05, 0.633, 0.544, 0.974, 0.936, 0.828, 0.998, 0.681, 0.723, 0.934, 0.98, 0.924, 0.964, 0.956, 0.723, 0.888, 0.446, 0.882, 0.856, 0.975, 0.775, 0.39, 0.738, 0.892, 0.689, 0.707, 0.595, 0.764, 0.883, 0.303, 0.472, 0.948, 0.972, 0.204, 0.622, 0.888, 0.988, 0.469, 0.873, 0.321, 0.979, 0.857, 0.958, 0.987, 0.811, 0.601, 0.899, 0.487, 0.746, 0.853, 0.734, 0.986, 0.933, 0.902, 0.502, 0.368, 0.728, 0.811, 0.835, 0.923, 0.418, 0.892, 0.892, 0.891, 0.961, 0.846, 0.851, 0.811, 0.997, 0.934, 0.847, 0.698, 0.369, 0.837, 0.928, 0.808, 0.714, 0.994, 0.932, 0.984, 0.779, 0.972, 0.984, 0.387, 0.673, 0.477, 0.752, 0.979, 0.995, 0.998, 0.217, 0.719, 0.932, 0.803, 0.981, 0.956, 0.91, 0.626, 0.624, 0.981, 0.995, 0.991, 0.498, 0.276, 0.965, 0.311, 0.577, 0.992, 0.478, 0.715, 0.519, 0.851, 0.972, 0.982, 0.95, 0.733, 0.92, 0.877, 0.66, 0.634, 0.982, 0.723, 0.92, 0.979, 0.524, 0.916, 0.544, 0.894, 0.635, 0.621, 0.763, 0.408, 0.967, 0.617, 0.441, 0.644, 0.948, 0.934, 0.958, 0.988, 0.972, 0.356, 0.891, 0.984, 0.237, 0.845, 0.782, 0.973, 0.984, 0.636, 0.621, 0.989, 0.996, 0.776, 0.675, 0.772, 0.943, 0.175]
global origin = 1
global destination = 50