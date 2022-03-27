global arcs = [1 4; 1 15; 1 17; 1 22; 1 24; 1 39; 1 47; 2 8; 2 17; 2 27; 2 44; 2 48; 3 5; 3 19; 3 21; 3 27; 3 33; 3 34; 3 36; 3 42; 3 46; 4 5; 4 7; 4 46; 5 15; 5 44; 5 46; 6 16; 6 28; 6 34; 7 19; 7 28; 7 36; 7 37; 7 44; 7 48; 8 3; 8 20; 9 4; 9 15; 9 33; 9 42; 9 45; 10 37; 11 16; 11 20; 11 21; 11 28; 11 39; 12 13; 12 18; 12 33; 12 46; 13 2; 13 14; 13 17; 13 22; 13 29; 13 42; 14 29; 15 17; 15 21; 15 28; 15 29; 15 35; 16 3; 16 9; 16 41; 17 5; 17 22; 17 37; 17 46; 18 2; 18 10; 19 9; 19 13; 19 36; 19 44; 20 13; 20 49; 21 8; 21 24; 21 27; 21 38; 22 4; 22 26; 22 31; 23 11; 23 15; 23 17; 23 18; 23 32; 24 11; 24 20; 24 25; 24 38; 24 39; 24 47; 24 49; 25 5; 25 14; 25 27; 25 42; 25 43; 25 44; 25 46; 26 16; 26 17; 26 47; 26 50; 27 4; 27 8; 27 10; 27 15; 27 16; 27 40; 27 49; 28 5; 28 8; 28 34; 28 39; 28 43; 28 48; 28 49; 29 2; 29 6; 29 23; 29 24; 29 25; 30 15; 30 33; 30 35; 30 39; 31 19; 31 21; 31 34; 31 48; 32 9; 32 10; 32 18; 32 33; 32 35; 32 42; 32 44; 33 7; 33 19; 33 21; 33 25; 33 41; 33 44; 34 10; 34 44; 35 21; 35 47; 36 23; 36 29; 36 35; 36 43; 37 3; 37 5; 37 15; 37 36; 37 42; 37 50; 38 11; 38 29; 38 41; 38 44; 38 49; 39 9; 39 12; 39 19; 39 26; 39 49; 40 2; 40 13; 40 19; 40 25; 40 39; 41 26; 41 30; 41 50; 42 11; 42 17; 42 28; 42 34; 42 37; 43 21; 43 27; 43 31; 43 46; 44 2; 44 8; 44 29; 44 37; 45 7; 45 20; 45 25; 45 43; 45 44; 46 11; 46 15; 46 30; 46 36; 46 41; 46 43; 46 44; 46 48; 46 50; 47 6; 47 11; 47 21; 47 24; 47 35; 47 42; 48 6; 48 18; 48 25; 48 33; 48 49; 49 16; 49 20; 49 22; 49 28; 49 35; 49 36; 49 40; 49 42]
global d_x = [10.0, 2.0, 1.0, 8.0, 4.0, 8.0, 9.0, 3.0, 3.0, 9.0, 3.0, 2.0, 10.0, 1.0, 4.0, 8.0, 3.0, 5.0, 2.0, 6.0, 2.0, 10.0, 9.0, 10.0, 2.0, 4.0, 5.0, 3.0, 3.0, 1.0, 4.0, 5.0, 9.0, 10.0, 4.0, 4.0, 9.0, 8.0, 3.0, 1.0, 4.0, 9.0, 10.0, 9.0, 7.0, 8.0, 9.0, 10.0, 5.0, 4.0, 3.0, 4.0, 5.0, 7.0, 1.0, 5.0, 4.0, 3.0, 6.0, 9.0, 3.0, 6.0, 9.0, 4.0, 6.0, 4.0, 3.0, 9.0, 2.0, 4.0, 2.0, 7.0, 10.0, 7.0, 2.0, 8.0, 4.0, 6.0, 7.0, 7.0, 5.0, 2.0, 4.0, 3.0, 6.0, 7.0, 3.0, 9.0, 4.0, 5.0, 9.0, 4.0, 3.0, 9.0, 9.0, 9.0, 7.0, 3.0, 5.0, 5.0, 6.0, 9.0, 4.0, 4.0, 6.0, 10.0, 9.0, 10.0, 9.0, 6.0, 9.0, 9.0, 2.0, 8.0, 3.0, 8.0, 2.0, 9.0, 8.0, 8.0, 7.0, 10.0, 3.0, 3.0, 7.0, 2.0, 10.0, 8.0, 10.0, 3.0, 7.0, 8.0, 2.0, 9.0, 8.0, 8.0, 10.0, 5.0, 10.0, 9.0, 9.0, 1.0, 7.0, 7.0, 7.0, 1.0, 5.0, 9.0, 9.0, 1.0, 4.0, 6.0, 4.0, 5.0, 3.0, 7.0, 9.0, 4.0, 7.0, 7.0, 6.0, 7.0, 9.0, 8.0, 10.0, 6.0, 10.0, 1.0, 7.0, 8.0, 8.0, 2.0, 3.0, 5.0, 2.0, 8.0, 5.0, 4.0, 10.0, 4.0, 5.0, 9.0, 4.0, 5.0, 8.0, 9.0, 6.0, 6.0, 6.0, 9.0, 9.0, 3.0, 9.0, 5.0, 1.0, 5.0, 5.0, 6.0, 3.0, 9.0, 7.0, 4.0, 8.0, 3.0, 10.0, 1.0, 4.0, 2.0, 5.0, 9.0, 10.0, 1.0, 9.0, 6.0, 10.0, 1.0, 7.0, 7.0, 4.0, 4.0, 4.0, 6.0, 2.0, 6.0, 4.0, 9.0, 8.0, 1.0]
global b_x = 5
global d_y = [2.0, 6.0, 9.0, 6.0, 1.0, 2.0, 4.0, 1.0, 6.0, 7.0, 6.0, 6.0, 4.0, 8.0, 10.0, 4.0, 1.0, 3.0, 2.0, 9.0, 3.0, 5.0, 7.0, 9.0, 1.0, 4.0, 8.0, 5.0, 8.0, 8.0, 1.0, 9.0, 6.0, 6.0, 4.0, 1.0, 1.0, 1.0, 7.0, 7.0, 2.0, 9.0, 9.0, 7.0, 9.0, 5.0, 1.0, 4.0, 8.0, 8.0, 2.0, 3.0, 1.0, 3.0, 2.0, 4.0, 8.0, 4.0, 9.0, 9.0, 3.0, 9.0, 2.0, 7.0, 8.0, 5.0, 4.0, 8.0, 8.0, 5.0, 2.0, 7.0, 7.0, 6.0, 6.0, 3.0, 5.0, 5.0, 1.0, 7.0, 3.0, 6.0, 9.0, 9.0, 7.0, 6.0, 10.0, 6.0, 2.0, 8.0, 8.0, 7.0, 10.0, 8.0, 6.0, 9.0, 6.0, 8.0, 3.0, 8.0, 2.0, 1.0, 8.0, 5.0, 4.0, 5.0, 5.0, 6.0, 4.0, 4.0, 1.0, 7.0, 4.0, 2.0, 10.0, 3.0, 3.0, 2.0, 3.0, 9.0, 4.0, 1.0, 5.0, 6.0, 6.0, 9.0, 8.0, 4.0, 1.0, 1.0, 3.0, 10.0, 1.0, 3.0, 2.0, 5.0, 6.0, 6.0, 8.0, 1.0, 6.0, 6.0, 8.0, 9.0, 4.0, 10.0, 4.0, 4.0, 6.0, 9.0, 7.0, 5.0, 10.0, 5.0, 7.0, 9.0, 1.0, 10.0, 7.0, 3.0, 8.0, 4.0, 7.0, 7.0, 8.0, 10.0, 1.0, 2.0, 6.0, 9.0, 1.0, 10.0, 6.0, 1.0, 5.0, 7.0, 10.0, 10.0, 2.0, 1.0, 7.0, 3.0, 10.0, 1.0, 5.0, 1.0, 7.0, 10.0, 10.0, 7.0, 7.0, 8.0, 5.0, 3.0, 10.0, 10.0, 8.0, 6.0, 6.0, 10.0, 9.0, 8.0, 5.0, 5.0, 8.0, 8.0, 4.0, 1.0, 7.0, 7.0, 8.0, 2.0, 5.0, 3.0, 10.0, 7.0, 3.0, 8.0, 6.0, 1.0, 3.0, 7.0, 10.0, 3.0, 1.0, 3.0, 3.0, 8.0]
global b_y = 10
global p = [0.98, 0.275, 0.599, 0.527, 0.387, 0.796, 0.348, 0.101, 0.057, 0.717, 0.15, 0.839, 0.613, 0.925, 0.971, 0.217, 0.972, 0.331, 0.905, 0.548, 0.315, 0.741, 0.6, 0.308, 0.178, 0.944, 0.582, 0.426, 0.784, 0.213, 0.922, 0.256, 0.017, 0.72, 0.51, 0.03, 0.874, 0.215, 0.351, 0.577, 0.834, 0.536, 0.645, 0.579, 0.225, 0.442, 0.685, 0.956, 0.205, 0.369, 0.269, 0.803, 0.161, 0.445, 0.48, 0.27, 0.718, 0.258, 0.353, 0.215, 0.075, 0.548, 0.244, 0.084, 0.755, 0.954, 0.476, 0.869, 0.752, 0.079, 0.319, 0.167, 0.115, 0.049, 0.698, 0.426, 0.808, 0.888, 0.104, 0.363, 0.52, 0.007, 0.957, 0.993, 0.898, 0.345, 0.165, 0.275, 0.68, 0.515, 0.81, 0.234, 0.786, 0.576, 0.546, 0.829, 0.177, 0.717, 0.839, 0.15, 0.953, 0.09, 0.142, 0.658, 0.624, 0.946, 0.14, 0.75, 0.052, 0.176, 0.477, 0.953, 0.972, 0.798, 0.335, 0.877, 0.135, 0.452, 0.234, 0.305, 0.471, 0.823, 0.18, 0.12, 0.319, 0.769, 0.404, 0.197, 0.149, 0.686, 0.708, 0.967, 0.68, 0.48, 0.507, 0.028, 0.261, 0.244, 0.814, 0.885, 0.557, 0.188, 0.636, 0.432, 0.587, 0.58, 0.135, 0.426, 0.522, 0.965, 0.019, 0.28, 0.645, 0.03, 0.485, 0.627, 0.462, 0.22, 0.838, 0.298, 0.321, 0.721, 0.792, 0.151, 0.181, 0.652, 0.265, 0.151, 0.964, 0.238, 0.404, 0.982, 0.848, 0.012, 0.325, 0.41, 0.81, 0.208, 0.896, 0.52, 0.222, 0.653, 0.062, 0.849, 0.025, 0.072, 0.288, 0.511, 0.441, 0.995, 0.6, 0.333, 0.615, 0.064, 0.568, 0.024, 0.489, 0.466, 0.971, 0.344, 0.342, 0.099, 0.689, 0.873, 0.477, 0.147, 0.685, 0.765, 0.394, 0.262, 0.214, 0.527, 0.107, 0.583, 0.01, 0.825, 0.519, 0.713, 0.794, 0.811, 0.805, 0.974, 0.516, 0.66, 0.848, 0.181, 0.039, 0.117]
global q = [0.984, 0.924, 0.705, 0.561, 0.45, 0.965, 0.991, 0.84, 0.067, 0.872, 0.861, 0.956, 0.973, 0.946, 0.974, 0.799, 0.979, 0.546, 0.964, 0.657, 0.806, 0.788, 0.859, 0.725, 0.791, 0.995, 0.602, 0.456, 0.833, 0.781, 0.932, 0.433, 0.741, 0.74, 0.849, 0.964, 0.993, 0.469, 0.477, 0.682, 0.921, 0.728, 0.772, 0.872, 0.331, 0.963, 0.981, 0.995, 0.471, 0.717, 0.751, 0.859, 0.271, 0.595, 0.88, 0.415, 0.752, 0.64, 0.79, 0.359, 0.565, 0.602, 0.263, 0.916, 0.809, 0.966, 0.607, 0.945, 0.759, 0.385, 0.357, 0.415, 0.713, 0.586, 0.844, 0.82, 0.96, 0.909, 0.191, 0.629, 0.953, 0.096, 0.985, 0.995, 0.969, 0.936, 0.259, 0.564, 0.929, 0.924, 0.887, 0.846, 0.966, 0.878, 0.792, 0.962, 0.865, 0.864, 0.843, 0.825, 0.961, 0.471, 0.348, 0.794, 0.714, 0.992, 0.651, 0.891, 0.852, 0.886, 0.757, 0.989, 0.988, 0.966, 0.829, 0.929, 0.951, 0.753, 0.506, 0.742, 0.896, 0.875, 0.566, 0.465, 0.497, 0.78, 0.676, 0.545, 0.208, 0.765, 0.799, 0.97, 0.79, 0.607, 0.983, 0.276, 0.761, 0.93, 0.912, 0.922, 0.69, 0.327, 0.809, 0.475, 0.957, 0.663, 0.823, 0.505, 0.851, 0.993, 0.491, 0.989, 0.87, 0.051, 0.588, 0.965, 0.487, 0.531, 0.968, 0.596, 0.763, 0.933, 0.822, 0.795, 0.751, 0.962, 0.406, 0.608, 0.973, 0.668, 0.615, 0.989, 0.981, 0.23, 0.744, 0.75, 0.917, 0.248, 0.928, 0.893, 0.736, 0.975, 0.989, 0.854, 0.459, 0.667, 0.53, 0.635, 0.477, 0.995, 0.887, 0.723, 0.907, 0.812, 0.845, 0.129, 0.555, 0.894, 0.982, 0.73, 0.724, 0.633, 0.725, 0.996, 0.572, 0.172, 0.913, 0.991, 0.458, 0.955, 0.532, 0.818, 0.538, 0.712, 0.258, 0.876, 0.901, 0.982, 0.947, 0.824, 0.857, 0.984, 0.566, 0.887, 0.905, 0.59, 0.479, 0.405]
global origin = 1
global destination = 50