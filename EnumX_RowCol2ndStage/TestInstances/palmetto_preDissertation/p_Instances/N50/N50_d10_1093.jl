global arcs = [1 12; 1 15; 1 23; 1 31; 1 36; 1 48; 2 14; 2 18; 2 23; 2 24; 2 35; 2 38; 2 46; 3 39; 3 45; 4 2; 4 13; 4 14; 4 15; 4 21; 4 22; 4 25; 4 33; 4 34; 4 42; 4 43; 4 47; 5 15; 5 17; 5 37; 6 14; 6 16; 6 34; 6 37; 6 38; 7 15; 7 18; 7 32; 7 40; 7 47; 8 6; 8 15; 8 19; 8 20; 8 24; 8 37; 8 41; 8 42; 9 8; 9 19; 9 20; 9 21; 9 30; 9 35; 9 43; 9 46; 9 47; 10 17; 10 44; 10 46; 11 50; 12 3; 12 13; 12 17; 12 32; 12 39; 12 47; 12 50; 13 4; 13 24; 13 45; 14 10; 14 11; 14 16; 15 8; 15 11; 15 18; 15 23; 15 26; 15 42; 16 2; 16 9; 16 27; 16 33; 16 41; 16 42; 17 9; 17 19; 17 25; 17 39; 18 6; 18 21; 18 50; 19 3; 19 4; 19 9; 19 14; 19 25; 19 38; 19 41; 20 2; 20 10; 20 13; 20 17; 20 26; 20 38; 20 47; 20 48; 20 49; 21 12; 21 16; 21 17; 21 18; 21 19; 21 35; 21 41; 21 46; 21 50; 22 9; 22 29; 22 32; 22 39; 23 5; 23 13; 24 8; 24 11; 24 33; 24 36; 24 48; 25 17; 25 20; 25 35; 25 37; 25 40; 25 41; 25 43; 26 33; 27 3; 27 5; 27 20; 27 22; 27 29; 27 35; 27 41; 27 47; 28 3; 28 42; 29 7; 29 8; 29 20; 29 25; 29 27; 29 30; 29 35; 30 2; 30 4; 30 8; 30 19; 30 25; 30 34; 30 49; 31 23; 31 27; 31 33; 32 27; 32 40; 32 47; 33 17; 33 18; 33 37; 33 45; 34 9; 34 15; 34 19; 34 20; 34 35; 34 38; 35 3; 35 13; 35 28; 35 32; 35 36; 35 42; 35 43; 35 44; 35 49; 36 15; 36 33; 36 43; 36 46; 37 21; 37 27; 37 34; 37 47; 37 49; 38 2; 38 16; 38 40; 38 48; 39 11; 39 12; 39 28; 39 41; 40 5; 40 8; 40 9; 40 20; 40 31; 40 36; 40 37; 40 42; 40 46; 41 10; 41 12; 41 33; 42 18; 43 3; 43 24; 43 26; 43 29; 43 39; 43 50; 44 12; 44 19; 44 22; 44 30; 44 35; 45 29; 45 44; 46 8; 46 21; 46 45; 47 10; 47 11; 47 22; 47 27; 47 29; 47 32; 47 35; 47 41; 47 42; 47 45; 47 49; 48 27; 48 47; 49 3; 49 10; 49 13; 49 15; 49 18; 49 20; 49 36; 49 37; 49 38; 49 43]
global d_x = [3.0, 4.0, 9.0, 3.0, 6.0, 5.0, 1.0, 5.0, 7.0, 6.0, 7.0, 7.0, 1.0, 9.0, 5.0, 4.0, 8.0, 7.0, 7.0, 1.0, 7.0, 7.0, 3.0, 9.0, 1.0, 5.0, 2.0, 6.0, 2.0, 7.0, 7.0, 8.0, 6.0, 9.0, 10.0, 2.0, 1.0, 1.0, 7.0, 9.0, 5.0, 4.0, 4.0, 6.0, 5.0, 6.0, 4.0, 7.0, 8.0, 2.0, 10.0, 1.0, 3.0, 6.0, 2.0, 10.0, 3.0, 8.0, 5.0, 4.0, 9.0, 3.0, 7.0, 2.0, 8.0, 8.0, 4.0, 6.0, 10.0, 10.0, 8.0, 8.0, 1.0, 8.0, 4.0, 6.0, 7.0, 9.0, 10.0, 4.0, 10.0, 2.0, 5.0, 7.0, 10.0, 1.0, 8.0, 4.0, 5.0, 6.0, 1.0, 1.0, 4.0, 6.0, 1.0, 8.0, 1.0, 1.0, 2.0, 10.0, 5.0, 10.0, 3.0, 7.0, 2.0, 6.0, 4.0, 5.0, 6.0, 8.0, 5.0, 8.0, 2.0, 9.0, 10.0, 6.0, 5.0, 5.0, 8.0, 7.0, 3.0, 9.0, 1.0, 3.0, 3.0, 10.0, 1.0, 7.0, 9.0, 4.0, 2.0, 4.0, 1.0, 9.0, 9.0, 5.0, 5.0, 9.0, 10.0, 6.0, 1.0, 6.0, 9.0, 8.0, 5.0, 8.0, 3.0, 5.0, 7.0, 10.0, 9.0, 5.0, 3.0, 9.0, 2.0, 8.0, 7.0, 1.0, 5.0, 1.0, 6.0, 5.0, 9.0, 10.0, 5.0, 2.0, 5.0, 4.0, 7.0, 3.0, 3.0, 8.0, 4.0, 5.0, 3.0, 9.0, 8.0, 4.0, 8.0, 9.0, 2.0, 7.0, 7.0, 3.0, 2.0, 9.0, 7.0, 7.0, 3.0, 7.0, 9.0, 3.0, 3.0, 1.0, 2.0, 8.0, 3.0, 8.0, 3.0, 4.0, 3.0, 5.0, 2.0, 3.0, 6.0, 5.0, 1.0, 5.0, 4.0, 10.0, 7.0, 7.0, 6.0, 2.0, 1.0, 8.0, 1.0, 5.0, 8.0, 9.0, 8.0, 5.0, 9.0, 6.0, 6.0, 10.0, 6.0, 4.0, 2.0, 1.0, 2.0, 4.0, 9.0, 6.0, 2.0, 2.0, 9.0, 9.0, 6.0, 10.0, 10.0, 2.0, 10.0, 1.0, 2.0, 8.0, 1.0, 1.0, 5.0, 2.0, 7.0, 8.0, 5.0, 2.0, 8.0]
global b_x = 5
global d_y = [1.0, 6.0, 5.0, 10.0, 8.0, 7.0, 5.0, 5.0, 6.0, 6.0, 4.0, 9.0, 5.0, 3.0, 4.0, 4.0, 3.0, 9.0, 3.0, 3.0, 4.0, 10.0, 3.0, 4.0, 10.0, 1.0, 2.0, 1.0, 9.0, 2.0, 9.0, 7.0, 10.0, 7.0, 8.0, 6.0, 3.0, 4.0, 4.0, 9.0, 8.0, 4.0, 8.0, 3.0, 3.0, 2.0, 1.0, 7.0, 9.0, 2.0, 5.0, 1.0, 3.0, 5.0, 10.0, 9.0, 8.0, 7.0, 4.0, 7.0, 10.0, 6.0, 7.0, 10.0, 3.0, 3.0, 10.0, 7.0, 3.0, 6.0, 2.0, 9.0, 3.0, 2.0, 10.0, 10.0, 1.0, 4.0, 4.0, 10.0, 1.0, 5.0, 9.0, 7.0, 10.0, 9.0, 3.0, 6.0, 9.0, 6.0, 8.0, 8.0, 10.0, 9.0, 3.0, 9.0, 5.0, 4.0, 9.0, 9.0, 6.0, 9.0, 4.0, 9.0, 8.0, 3.0, 5.0, 3.0, 6.0, 1.0, 2.0, 8.0, 6.0, 3.0, 9.0, 2.0, 6.0, 1.0, 2.0, 6.0, 9.0, 1.0, 9.0, 9.0, 5.0, 10.0, 9.0, 10.0, 8.0, 1.0, 6.0, 5.0, 3.0, 5.0, 5.0, 7.0, 1.0, 9.0, 5.0, 7.0, 10.0, 6.0, 5.0, 2.0, 4.0, 6.0, 4.0, 2.0, 7.0, 4.0, 8.0, 3.0, 7.0, 8.0, 5.0, 9.0, 5.0, 9.0, 2.0, 9.0, 1.0, 8.0, 7.0, 6.0, 5.0, 1.0, 1.0, 5.0, 10.0, 2.0, 3.0, 9.0, 9.0, 2.0, 3.0, 7.0, 10.0, 10.0, 5.0, 8.0, 3.0, 6.0, 1.0, 6.0, 8.0, 2.0, 8.0, 6.0, 4.0, 9.0, 5.0, 4.0, 8.0, 2.0, 2.0, 2.0, 5.0, 6.0, 8.0, 5.0, 6.0, 6.0, 3.0, 10.0, 2.0, 9.0, 8.0, 7.0, 7.0, 2.0, 5.0, 10.0, 10.0, 6.0, 3.0, 9.0, 1.0, 6.0, 9.0, 5.0, 5.0, 1.0, 7.0, 3.0, 8.0, 4.0, 1.0, 3.0, 10.0, 8.0, 6.0, 8.0, 9.0, 8.0, 9.0, 2.0, 5.0, 1.0, 8.0, 8.0, 7.0, 10.0, 4.0, 10.0, 5.0, 2.0, 6.0, 2.0, 10.0, 1.0, 2.0, 10.0, 9.0, 6.0, 4.0]
global b_y = 10
global p = [0.939, 0.379, 0.368, 0.351, 0.624, 0.866, 0.165, 0.927, 0.13, 0.038, 0.496, 0.632, 0.961, 0.172, 0.243, 0.021, 0.913, 0.44, 0.656, 0.094, 0.553, 0.489, 0.704, 0.234, 0.202, 0.697, 0.742, 0.541, 0.498, 0.991, 0.524, 0.184, 0.878, 0.21, 0.698, 0.93, 0.757, 0.6, 0.353, 0.502, 0.464, 0.97, 0.026, 0.885, 0.564, 0.31, 0.666, 0.095, 0.211, 0.837, 0.359, 0.521, 0.315, 0.642, 0.71, 0.312, 0.313, 0.47, 0.686, 0.442, 0.811, 0.18, 0.881, 0.254, 0.913, 0.654, 0.854, 0.903, 0.723, 0.737, 0.053, 0.823, 0.453, 0.341, 0.216, 0.788, 0.616, 0.538, 0.538, 0.153, 0.952, 0.919, 0.29, 0.498, 0.406, 0.011, 0.21, 0.182, 0.378, 0.592, 0.97, 0.475, 0.533, 0.15, 0.681, 0.139, 0.877, 0.595, 0.693, 0.768, 0.765, 0.942, 0.703, 0.042, 0.973, 0.437, 0.94, 0.586, 0.582, 0.382, 0.438, 0.795, 0.555, 0.788, 0.419, 0.42, 0.237, 0.198, 0.46, 0.946, 0.682, 0.807, 0.727, 0.878, 0.396, 0.196, 0.232, 0.888, 0.174, 0.507, 0.84, 0.637, 0.019, 0.889, 0.336, 0.027, 0.272, 0.111, 0.277, 0.1, 0.482, 0.746, 0.807, 0.856, 0.389, 0.163, 0.929, 0.606, 0.912, 0.838, 0.398, 0.799, 0.711, 0.836, 0.363, 0.364, 0.413, 0.672, 0.644, 0.176, 0.485, 0.412, 0.95, 0.017, 0.168, 0.825, 0.624, 0.82, 0.267, 0.005, 0.406, 0.872, 0.894, 0.192, 0.127, 0.49, 0.749, 0.681, 0.813, 0.539, 0.916, 0.766, 0.411, 0.434, 0.339, 0.938, 0.24, 0.879, 0.618, 0.53, 0.543, 0.252, 0.816, 0.301, 0.353, 0.247, 0.887, 0.278, 0.514, 0.894, 0.038, 0.168, 0.86, 0.094, 0.209, 0.751, 0.788, 0.851, 0.618, 0.615, 0.458, 0.983, 0.787, 0.937, 0.635, 0.884, 0.084, 0.19, 0.943, 0.34, 0.485, 0.967, 0.221, 0.963, 0.233, 0.496, 0.891, 0.176, 0.102, 0.707, 0.342, 0.51, 0.354, 0.328, 0.639, 0.217, 0.994, 0.863, 0.605, 0.32, 0.932, 0.729, 0.397, 0.776, 0.54, 0.339, 0.6, 0.982, 0.022, 0.198, 0.359, 0.665, 0.172, 0.816, 0.995]
global q = [0.954, 0.438, 0.56, 0.432, 0.952, 0.972, 0.618, 0.961, 0.527, 0.544, 0.819, 0.946, 0.991, 0.441, 0.914, 0.512, 0.93, 0.846, 0.791, 0.263, 0.806, 0.626, 0.856, 0.658, 0.425, 0.731, 0.938, 0.639, 0.671, 0.999, 0.744, 0.923, 0.939, 0.374, 0.736, 0.983, 0.804, 0.945, 0.543, 0.509, 0.6, 0.989, 0.244, 0.937, 0.864, 0.941, 0.694, 0.42, 0.966, 0.871, 0.489, 0.788, 0.84, 0.682, 0.869, 0.497, 0.365, 0.613, 0.73, 0.971, 0.989, 0.488, 0.951, 0.601, 0.913, 0.819, 0.894, 0.965, 0.97, 0.817, 0.516, 0.984, 0.741, 0.806, 0.367, 0.874, 0.722, 0.63, 0.707, 0.94, 0.978, 0.948, 0.817, 0.896, 0.61, 0.308, 0.474, 0.942, 0.41, 0.916, 0.98, 0.611, 0.589, 0.357, 0.896, 0.486, 0.944, 0.742, 0.907, 0.785, 0.822, 0.995, 0.918, 0.261, 0.976, 0.738, 0.988, 0.828, 0.974, 0.852, 0.482, 0.825, 0.625, 0.934, 0.442, 0.458, 0.353, 0.717, 0.789, 0.969, 0.854, 0.937, 0.836, 0.97, 0.434, 0.744, 0.395, 0.93, 0.542, 0.759, 0.867, 0.645, 0.025, 0.972, 0.636, 0.783, 0.561, 0.279, 0.705, 0.934, 0.981, 0.994, 0.866, 0.893, 0.859, 0.19, 0.994, 0.85, 0.992, 0.887, 0.817, 0.807, 0.713, 0.984, 0.64, 0.836, 0.658, 0.782, 0.947, 0.188, 0.672, 0.912, 0.99, 0.095, 0.541, 0.867, 0.649, 0.861, 0.328, 0.383, 0.542, 0.915, 0.994, 0.388, 0.995, 0.794, 0.829, 0.888, 0.842, 0.642, 0.935, 0.791, 0.828, 0.465, 0.83, 0.969, 0.554, 0.977, 0.625, 0.713, 0.871, 0.886, 0.873, 0.327, 0.389, 0.514, 0.983, 0.584, 0.965, 0.9, 0.618, 0.986, 0.964, 0.762, 0.462, 0.921, 0.863, 0.905, 0.991, 0.839, 0.823, 0.991, 0.851, 0.951, 0.736, 0.886, 0.507, 0.827, 0.949, 0.495, 0.719, 0.968, 0.794, 0.974, 0.55, 0.91, 0.907, 0.546, 0.181, 0.866, 0.977, 0.569, 0.536, 0.428, 0.977, 0.944, 0.996, 0.891, 0.607, 0.368, 0.984, 0.926, 0.646, 0.963, 0.578, 0.386, 0.689, 0.983, 0.399, 0.794, 0.441, 0.808, 0.505, 0.927, 0.996]
global origin = 1
global destination = 50