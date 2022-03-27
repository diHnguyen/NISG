global arcs = [1 3; 1 34; 1 38; 1 41; 2 8; 2 17; 2 21; 2 29; 2 48; 3 11; 3 14; 3 32; 3 35; 4 3; 4 8; 4 10; 4 14; 4 25; 4 26; 4 28; 5 17; 5 46; 6 3; 6 4; 6 18; 6 21; 6 23; 6 25; 6 50; 7 4; 7 8; 7 12; 7 40; 7 44; 8 6; 8 10; 8 14; 8 30; 8 35; 8 46; 9 4; 9 28; 10 7; 10 33; 11 7; 11 15; 11 30; 11 33; 12 13; 12 14; 12 20; 12 23; 12 46; 12 50; 13 21; 13 27; 13 32; 13 33; 14 9; 14 13; 14 32; 14 36; 14 49; 15 12; 15 16; 15 17; 15 22; 15 34; 15 40; 15 45; 15 48; 16 33; 16 35; 16 45; 17 2; 17 22; 17 47; 17 50; 18 9; 18 40; 18 48; 19 4; 19 6; 19 11; 19 12; 19 14; 19 16; 19 21; 19 35; 20 4; 20 7; 20 33; 20 47; 20 48; 21 7; 22 6; 22 15; 22 19; 22 20; 23 18; 23 22; 23 27; 23 31; 23 39; 24 3; 24 7; 24 8; 24 11; 24 13; 24 17; 24 26; 24 27; 24 30; 24 48; 24 50; 25 12; 25 16; 25 33; 25 35; 25 45; 25 50; 26 9; 26 18; 26 21; 26 46; 27 2; 27 7; 27 15; 27 32; 27 42; 27 43; 27 46; 28 21; 28 38; 29 42; 29 46; 30 2; 30 7; 30 23; 30 28; 30 35; 31 17; 31 25; 31 26; 31 30; 31 46; 31 50; 32 15; 32 27; 32 37; 32 39; 33 3; 33 19; 33 20; 33 24; 33 46; 34 8; 34 50; 35 10; 35 18; 35 36; 35 43; 36 3; 36 17; 36 20; 36 23; 36 39; 36 44; 36 45; 37 14; 37 16; 37 45; 38 31; 38 42; 39 2; 39 26; 39 35; 39 37; 40 10; 40 24; 40 38; 40 42; 41 3; 41 5; 41 18; 42 13; 42 32; 42 36; 42 43; 43 5; 43 14; 43 16; 43 18; 43 32; 43 44; 44 28; 44 31; 44 36; 44 43; 45 2; 45 21; 45 22; 45 29; 45 40; 46 15; 46 16; 46 23; 46 34; 47 10; 47 19; 47 20; 47 25; 47 31; 47 32; 48 12; 48 16; 48 25; 48 28; 48 40; 49 7; 49 9; 49 10; 49 29; 49 34; 49 40]
global d_x = [3.0, 5.0, 2.0, 8.0, 1.0, 5.0, 4.0, 5.0, 10.0, 6.0, 6.0, 2.0, 2.0, 8.0, 4.0, 1.0, 6.0, 3.0, 4.0, 6.0, 6.0, 5.0, 6.0, 9.0, 10.0, 2.0, 8.0, 2.0, 5.0, 9.0, 2.0, 8.0, 10.0, 7.0, 3.0, 2.0, 5.0, 9.0, 9.0, 4.0, 8.0, 9.0, 7.0, 5.0, 10.0, 7.0, 4.0, 1.0, 2.0, 6.0, 6.0, 5.0, 6.0, 5.0, 2.0, 8.0, 4.0, 4.0, 3.0, 10.0, 7.0, 7.0, 10.0, 6.0, 2.0, 4.0, 4.0, 5.0, 1.0, 6.0, 8.0, 4.0, 7.0, 9.0, 7.0, 3.0, 3.0, 4.0, 4.0, 9.0, 1.0, 2.0, 10.0, 1.0, 6.0, 6.0, 5.0, 8.0, 6.0, 10.0, 4.0, 3.0, 7.0, 10.0, 4.0, 6.0, 3.0, 4.0, 10.0, 5.0, 3.0, 3.0, 4.0, 7.0, 7.0, 5.0, 2.0, 5.0, 1.0, 10.0, 4.0, 8.0, 7.0, 8.0, 4.0, 4.0, 3.0, 1.0, 9.0, 1.0, 1.0, 8.0, 3.0, 10.0, 10.0, 1.0, 9.0, 10.0, 3.0, 8.0, 7.0, 3.0, 3.0, 10.0, 2.0, 5.0, 5.0, 3.0, 2.0, 4.0, 6.0, 2.0, 10.0, 2.0, 3.0, 3.0, 6.0, 5.0, 8.0, 9.0, 3.0, 10.0, 3.0, 6.0, 1.0, 7.0, 1.0, 5.0, 5.0, 9.0, 6.0, 9.0, 7.0, 3.0, 7.0, 10.0, 4.0, 4.0, 5.0, 8.0, 7.0, 5.0, 9.0, 10.0, 5.0, 7.0, 2.0, 7.0, 6.0, 10.0, 5.0, 6.0, 9.0, 3.0, 3.0, 6.0, 1.0, 3.0, 4.0, 2.0, 9.0, 5.0, 1.0, 6.0, 5.0, 9.0, 7.0, 3.0, 10.0, 3.0, 3.0, 9.0, 2.0, 5.0, 4.0, 8.0, 3.0, 7.0, 5.0, 9.0, 3.0, 1.0, 8.0, 7.0, 7.0, 1.0, 1.0, 10.0, 2.0, 4.0, 8.0, 5.0, 9.0, 8.0, 1.0]
global b_x = 5
global d_y = [8.0, 6.0, 1.0, 9.0, 4.0, 10.0, 3.0, 9.0, 4.0, 9.0, 9.0, 9.0, 9.0, 3.0, 2.0, 2.0, 1.0, 3.0, 5.0, 7.0, 10.0, 7.0, 5.0, 5.0, 7.0, 1.0, 7.0, 7.0, 7.0, 6.0, 6.0, 7.0, 3.0, 2.0, 9.0, 9.0, 4.0, 10.0, 3.0, 5.0, 10.0, 10.0, 6.0, 4.0, 7.0, 1.0, 7.0, 2.0, 8.0, 6.0, 3.0, 7.0, 10.0, 10.0, 6.0, 7.0, 9.0, 9.0, 3.0, 7.0, 5.0, 4.0, 6.0, 5.0, 9.0, 7.0, 7.0, 7.0, 4.0, 7.0, 4.0, 6.0, 2.0, 5.0, 5.0, 5.0, 5.0, 7.0, 9.0, 7.0, 6.0, 10.0, 5.0, 10.0, 3.0, 5.0, 6.0, 1.0, 8.0, 4.0, 1.0, 6.0, 1.0, 5.0, 3.0, 3.0, 8.0, 8.0, 10.0, 6.0, 7.0, 6.0, 6.0, 10.0, 4.0, 3.0, 10.0, 10.0, 2.0, 5.0, 7.0, 1.0, 5.0, 1.0, 10.0, 5.0, 1.0, 6.0, 7.0, 7.0, 1.0, 2.0, 3.0, 9.0, 5.0, 4.0, 2.0, 9.0, 3.0, 9.0, 9.0, 7.0, 6.0, 6.0, 3.0, 7.0, 1.0, 2.0, 6.0, 3.0, 8.0, 8.0, 9.0, 7.0, 2.0, 8.0, 5.0, 4.0, 1.0, 9.0, 8.0, 6.0, 7.0, 1.0, 9.0, 5.0, 8.0, 7.0, 4.0, 4.0, 7.0, 6.0, 7.0, 8.0, 3.0, 6.0, 9.0, 5.0, 4.0, 3.0, 9.0, 7.0, 1.0, 5.0, 3.0, 9.0, 4.0, 7.0, 10.0, 7.0, 9.0, 7.0, 1.0, 6.0, 8.0, 2.0, 8.0, 10.0, 4.0, 5.0, 7.0, 6.0, 9.0, 2.0, 2.0, 8.0, 4.0, 3.0, 7.0, 6.0, 8.0, 10.0, 10.0, 5.0, 1.0, 3.0, 9.0, 6.0, 6.0, 1.0, 6.0, 1.0, 5.0, 7.0, 10.0, 9.0, 4.0, 4.0, 10.0, 4.0, 1.0, 7.0, 4.0, 7.0, 8.0]
global b_y = 10
global p = [0.235, 0.351, 0.863, 0.546, 0.885, 0.373, 0.374, 0.942, 0.027, 0.042, 0.018, 0.09, 0.951, 0.006, 0.175, 0.507, 0.478, 0.703, 0.25, 0.625, 0.272, 0.905, 0.792, 0.624, 0.803, 0.958, 0.469, 0.535, 0.159, 0.943, 0.871, 0.006, 0.543, 0.553, 0.269, 0.559, 0.936, 0.566, 0.068, 0.554, 0.338, 0.517, 0.231, 0.96, 0.14, 0.484, 0.403, 0.115, 0.674, 0.096, 0.486, 0.131, 0.626, 0.183, 0.682, 0.732, 0.216, 0.077, 0.483, 0.751, 0.468, 0.733, 0.396, 0.329, 0.561, 0.413, 0.168, 0.035, 0.794, 0.319, 0.187, 0.183, 0.926, 0.616, 0.611, 0.176, 0.33, 0.707, 0.706, 0.108, 0.997, 0.834, 0.698, 0.066, 0.375, 0.236, 0.288, 0.148, 0.732, 0.771, 0.115, 0.019, 0.507, 0.358, 0.905, 0.777, 0.962, 0.281, 0.076, 0.508, 0.439, 0.197, 0.17, 0.83, 0.877, 0.638, 0.175, 0.655, 0.892, 0.843, 0.443, 0.013, 0.506, 0.951, 0.083, 0.274, 0.771, 0.492, 0.372, 0.999, 0.459, 0.566, 0.917, 0.048, 0.469, 0.529, 0.966, 0.403, 0.217, 0.727, 0.813, 0.196, 0.934, 0.799, 0.835, 0.516, 0.274, 0.913, 0.244, 0.669, 0.403, 0.698, 0.177, 0.99, 0.26, 0.846, 0.137, 0.573, 0.356, 0.725, 0.964, 0.896, 0.749, 0.969, 0.162, 0.979, 0.03, 0.662, 0.234, 0.211, 0.695, 0.262, 0.66, 0.325, 0.77, 0.446, 0.89, 0.473, 0.355, 0.711, 0.622, 0.96, 0.092, 0.045, 0.768, 0.571, 0.492, 0.069, 0.727, 0.445, 0.569, 0.214, 0.751, 0.656, 0.323, 0.728, 0.023, 0.694, 0.437, 0.351, 0.486, 0.263, 0.113, 0.77, 0.68, 0.662, 0.703, 0.883, 0.355, 0.312, 0.632, 0.819, 0.312, 0.211, 0.835, 0.839, 0.847, 0.755, 0.49, 0.159, 0.759, 0.759, 0.181, 0.201, 0.504, 0.813, 0.455, 0.434, 0.627, 0.124, 0.569, 0.168, 0.045, 0.149, 0.072]
global q = [0.935, 0.958, 0.98, 0.658, 0.957, 0.568, 0.684, 0.971, 0.385, 0.646, 0.114, 0.706, 0.981, 0.41, 0.888, 0.99, 0.676, 0.933, 0.453, 0.805, 0.926, 0.908, 0.968, 0.875, 0.97, 0.991, 0.675, 0.702, 0.4, 0.957, 0.988, 0.931, 0.696, 0.724, 0.817, 0.851, 0.953, 0.973, 0.154, 0.995, 0.42, 0.533, 0.989, 0.979, 0.693, 0.789, 0.731, 0.402, 0.858, 0.264, 0.637, 0.659, 0.75, 0.245, 0.742, 0.924, 0.867, 0.493, 0.636, 0.931, 0.55, 0.911, 0.979, 0.758, 0.951, 0.62, 0.333, 0.528, 0.982, 0.896, 0.358, 0.872, 0.937, 0.768, 0.913, 0.683, 0.523, 0.765, 0.75, 0.196, 0.997, 0.886, 0.809, 0.609, 0.732, 0.996, 0.427, 0.768, 0.838, 0.968, 0.752, 0.724, 0.605, 0.706, 0.927, 0.888, 0.998, 0.446, 0.433, 0.519, 0.677, 0.709, 0.478, 0.991, 0.998, 0.906, 0.678, 0.847, 0.986, 0.948, 0.984, 0.289, 0.881, 0.975, 0.226, 0.692, 0.914, 0.511, 0.718, 0.999, 0.83, 0.828, 0.922, 0.664, 0.536, 0.689, 0.982, 0.548, 0.298, 0.928, 0.901, 0.632, 0.995, 0.843, 0.9, 0.673, 0.481, 0.917, 0.495, 0.984, 0.885, 0.871, 0.503, 0.996, 0.762, 0.915, 0.398, 0.87, 0.795, 0.773, 0.998, 0.974, 0.95, 0.975, 0.326, 0.992, 0.914, 0.958, 0.584, 0.93, 0.862, 0.5, 0.909, 0.801, 0.999, 0.536, 0.935, 0.676, 0.488, 0.91, 0.869, 0.995, 0.951, 0.984, 0.94, 0.734, 0.778, 0.278, 0.74, 0.986, 0.772, 0.669, 0.887, 0.97, 0.979, 0.927, 0.635, 0.886, 0.477, 0.977, 0.837, 0.639, 0.259, 0.873, 0.947, 0.679, 0.721, 0.912, 0.625, 0.87, 0.949, 0.836, 0.442, 0.386, 0.838, 0.963, 0.924, 0.866, 0.833, 0.814, 0.935, 0.794, 0.654, 0.345, 0.894, 0.883, 0.56, 0.559, 0.925, 0.222, 0.988, 0.279, 0.863, 0.548, 0.785]
global origin = 1
global destination = 50