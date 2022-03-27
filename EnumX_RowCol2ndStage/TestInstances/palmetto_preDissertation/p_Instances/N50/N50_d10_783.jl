global arcs = [1 7; 1 22; 1 27; 1 32; 2 5; 2 25; 2 27; 2 32; 2 39; 2 40; 2 43; 2 44; 2 50; 3 21; 3 29; 3 36; 4 11; 4 14; 4 25; 4 28; 4 36; 4 38; 5 3; 5 7; 5 13; 6 4; 6 44; 7 5; 7 19; 7 28; 7 31; 8 6; 8 12; 8 24; 8 48; 9 13; 9 24; 9 41; 9 42; 10 2; 10 13; 10 14; 10 15; 10 16; 10 22; 10 40; 11 2; 11 14; 11 15; 11 29; 11 34; 11 41; 12 15; 12 17; 12 19; 12 34; 12 39; 13 5; 13 25; 13 39; 14 3; 14 19; 14 23; 14 26; 14 28; 14 30; 14 35; 14 46; 15 3; 15 17; 15 36; 15 40; 15 44; 15 45; 15 49; 16 8; 16 18; 16 34; 16 36; 16 44; 17 3; 17 7; 17 11; 17 34; 17 36; 17 41; 18 4; 18 7; 18 34; 19 15; 19 18; 19 23; 19 33; 19 45; 20 7; 20 13; 20 37; 20 43; 20 46; 20 49; 21 4; 21 43; 21 48; 22 7; 22 8; 22 17; 22 18; 22 21; 22 48; 23 2; 23 13; 23 15; 23 32; 23 37; 23 47; 24 2; 24 8; 24 9; 24 17; 24 36; 24 37; 24 42; 25 6; 25 7; 25 21; 25 41; 25 49; 26 7; 26 17; 26 20; 26 29; 26 43; 27 8; 27 21; 27 34; 27 42; 28 5; 28 11; 28 25; 28 29; 28 32; 28 45; 28 46; 28 48; 29 4; 29 6; 29 7; 29 8; 29 41; 29 50; 30 14; 31 9; 31 10; 31 41; 31 47; 32 23; 32 36; 33 2; 33 17; 33 18; 33 26; 33 36; 34 3; 34 6; 34 7; 34 11; 34 25; 34 30; 34 31; 34 40; 34 44; 35 4; 35 6; 35 16; 35 31; 35 46; 35 49; 36 9; 36 38; 36 47; 37 19; 37 24; 37 29; 37 31; 37 44; 37 49; 38 20; 38 22; 38 25; 38 40; 39 27; 40 8; 41 16; 41 17; 41 22; 41 46; 41 48; 42 21; 42 22; 42 33; 42 34; 42 38; 42 41; 42 43; 43 3; 43 42; 44 6; 44 11; 44 14; 44 15; 44 19; 44 27; 44 42; 44 50; 45 7; 45 13; 45 21; 45 23; 45 42; 45 50; 46 8; 46 9; 46 10; 46 20; 46 27; 46 28; 46 32; 46 45; 46 50; 47 23; 47 27; 47 30; 47 34; 47 46; 48 6; 48 37; 48 38; 48 40; 49 21; 49 24; 49 31; 49 32]
global d_x = [10.0, 7.0, 4.0, 4.0, 1.0, 4.0, 3.0, 10.0, 9.0, 4.0, 7.0, 5.0, 9.0, 1.0, 5.0, 5.0, 5.0, 5.0, 10.0, 4.0, 1.0, 5.0, 3.0, 8.0, 6.0, 4.0, 10.0, 9.0, 4.0, 8.0, 10.0, 4.0, 5.0, 9.0, 5.0, 5.0, 9.0, 10.0, 9.0, 10.0, 4.0, 4.0, 3.0, 9.0, 8.0, 7.0, 2.0, 10.0, 9.0, 7.0, 7.0, 10.0, 5.0, 10.0, 6.0, 7.0, 2.0, 3.0, 1.0, 4.0, 6.0, 8.0, 9.0, 4.0, 1.0, 7.0, 9.0, 2.0, 5.0, 5.0, 8.0, 3.0, 4.0, 3.0, 2.0, 9.0, 2.0, 4.0, 7.0, 9.0, 5.0, 4.0, 2.0, 1.0, 6.0, 9.0, 7.0, 2.0, 6.0, 8.0, 6.0, 6.0, 2.0, 4.0, 2.0, 4.0, 1.0, 7.0, 4.0, 9.0, 7.0, 9.0, 5.0, 6.0, 6.0, 8.0, 1.0, 2.0, 9.0, 6.0, 1.0, 5.0, 1.0, 1.0, 9.0, 9.0, 8.0, 5.0, 7.0, 9.0, 5.0, 10.0, 10.0, 1.0, 9.0, 3.0, 2.0, 2.0, 10.0, 9.0, 4.0, 1.0, 4.0, 8.0, 8.0, 5.0, 6.0, 6.0, 4.0, 6.0, 7.0, 2.0, 1.0, 6.0, 4.0, 5.0, 3.0, 2.0, 9.0, 10.0, 4.0, 9.0, 4.0, 4.0, 8.0, 2.0, 5.0, 4.0, 3.0, 3.0, 7.0, 8.0, 7.0, 2.0, 3.0, 6.0, 5.0, 5.0, 5.0, 7.0, 5.0, 8.0, 4.0, 2.0, 6.0, 3.0, 1.0, 4.0, 1.0, 6.0, 10.0, 1.0, 10.0, 2.0, 9.0, 2.0, 9.0, 2.0, 10.0, 4.0, 4.0, 1.0, 4.0, 6.0, 1.0, 4.0, 1.0, 4.0, 9.0, 3.0, 8.0, 5.0, 3.0, 9.0, 2.0, 2.0, 3.0, 10.0, 6.0, 8.0, 5.0, 10.0, 2.0, 1.0, 6.0, 4.0, 10.0, 1.0, 5.0, 2.0, 7.0, 4.0, 8.0, 3.0, 4.0, 2.0, 4.0, 4.0, 2.0, 8.0, 5.0, 4.0, 9.0, 7.0, 2.0, 6.0, 4.0, 6.0, 7.0, 5.0, 4.0, 1.0]
global b_x = 5
global d_y = [8.0, 5.0, 6.0, 2.0, 8.0, 10.0, 7.0, 6.0, 8.0, 7.0, 3.0, 1.0, 10.0, 1.0, 1.0, 1.0, 5.0, 9.0, 4.0, 5.0, 1.0, 8.0, 8.0, 3.0, 6.0, 8.0, 2.0, 5.0, 1.0, 8.0, 3.0, 3.0, 7.0, 4.0, 5.0, 2.0, 4.0, 8.0, 7.0, 2.0, 1.0, 7.0, 10.0, 8.0, 4.0, 5.0, 2.0, 2.0, 2.0, 2.0, 8.0, 5.0, 6.0, 4.0, 9.0, 1.0, 1.0, 5.0, 3.0, 3.0, 9.0, 8.0, 7.0, 7.0, 5.0, 6.0, 2.0, 3.0, 5.0, 2.0, 7.0, 8.0, 3.0, 7.0, 2.0, 3.0, 1.0, 9.0, 2.0, 4.0, 10.0, 7.0, 5.0, 3.0, 3.0, 9.0, 9.0, 4.0, 4.0, 4.0, 9.0, 4.0, 8.0, 8.0, 10.0, 5.0, 4.0, 5.0, 5.0, 9.0, 3.0, 5.0, 6.0, 4.0, 8.0, 7.0, 8.0, 10.0, 5.0, 2.0, 4.0, 7.0, 7.0, 6.0, 2.0, 1.0, 2.0, 1.0, 1.0, 3.0, 9.0, 1.0, 8.0, 6.0, 8.0, 8.0, 5.0, 6.0, 1.0, 4.0, 10.0, 10.0, 10.0, 2.0, 2.0, 4.0, 8.0, 10.0, 10.0, 2.0, 5.0, 5.0, 10.0, 2.0, 1.0, 6.0, 2.0, 5.0, 6.0, 1.0, 1.0, 3.0, 2.0, 5.0, 1.0, 3.0, 8.0, 3.0, 6.0, 2.0, 2.0, 4.0, 4.0, 5.0, 4.0, 8.0, 6.0, 9.0, 5.0, 2.0, 3.0, 10.0, 4.0, 7.0, 5.0, 1.0, 1.0, 4.0, 3.0, 6.0, 3.0, 5.0, 6.0, 3.0, 8.0, 5.0, 1.0, 8.0, 5.0, 2.0, 1.0, 1.0, 5.0, 10.0, 8.0, 8.0, 10.0, 9.0, 10.0, 9.0, 6.0, 2.0, 10.0, 5.0, 7.0, 5.0, 7.0, 10.0, 9.0, 5.0, 2.0, 10.0, 1.0, 2.0, 6.0, 1.0, 9.0, 3.0, 2.0, 2.0, 4.0, 5.0, 7.0, 8.0, 7.0, 2.0, 5.0, 6.0, 10.0, 6.0, 8.0, 3.0, 2.0, 2.0, 8.0, 8.0, 8.0, 5.0, 7.0, 9.0, 1.0, 3.0]
global b_y = 10
global p = [0.049, 0.569, 0.23, 0.094, 0.992, 0.758, 0.217, 0.289, 0.858, 0.368, 0.14, 0.697, 0.813, 0.848, 0.579, 0.895, 0.265, 0.435, 0.21, 0.784, 0.604, 0.668, 0.576, 0.518, 0.672, 0.886, 0.678, 0.1, 0.384, 0.621, 0.423, 0.421, 0.553, 0.891, 0.909, 0.247, 0.809, 0.345, 0.464, 0.59, 0.967, 0.504, 0.978, 0.022, 0.658, 0.681, 0.566, 0.775, 0.984, 0.353, 0.971, 0.927, 0.075, 0.134, 0.755, 0.973, 0.83, 0.565, 0.303, 0.167, 0.978, 0.057, 0.01, 0.183, 0.193, 0.903, 0.028, 0.249, 0.05, 0.3, 0.543, 0.147, 0.255, 0.878, 0.086, 0.551, 0.147, 0.549, 0.716, 0.295, 0.361, 0.431, 0.838, 0.594, 0.056, 0.208, 0.718, 0.566, 0.565, 0.877, 0.487, 0.5, 0.08, 0.244, 0.659, 0.889, 0.698, 0.945, 0.22, 0.939, 0.309, 0.583, 0.998, 0.8, 0.547, 0.778, 0.86, 0.373, 0.768, 0.606, 0.025, 0.985, 0.039, 0.375, 0.449, 0.682, 0.649, 0.46, 0.363, 0.539, 0.274, 0.181, 0.167, 0.377, 0.264, 0.595, 0.139, 0.37, 0.003, 0.201, 0.118, 0.607, 0.541, 0.002, 0.338, 0.767, 0.35, 0.281, 0.574, 0.942, 0.077, 0.801, 0.67, 0.695, 0.568, 0.324, 0.05, 0.249, 0.518, 0.749, 0.113, 0.531, 0.914, 0.362, 0.606, 0.807, 0.468, 0.518, 0.697, 0.541, 0.05, 0.75, 0.253, 0.233, 0.264, 0.632, 0.438, 0.555, 0.215, 0.409, 0.846, 0.489, 0.913, 0.746, 0.728, 0.807, 0.411, 0.536, 0.036, 0.799, 0.121, 0.152, 0.493, 0.412, 0.85, 0.89, 0.199, 0.175, 0.139, 0.607, 0.987, 0.118, 0.809, 0.395, 0.356, 0.469, 0.055, 0.774, 0.884, 0.032, 0.842, 0.377, 0.653, 0.364, 0.203, 0.006, 0.934, 0.918, 0.894, 0.672, 0.968, 0.23, 0.019, 0.027, 0.334, 0.541, 0.596, 0.114, 0.83, 0.921, 0.917, 0.895, 0.6, 0.427, 0.306, 0.007, 0.186, 0.158, 0.074, 0.378, 0.927, 0.768, 0.55, 0.071, 0.678, 0.12, 0.065, 0.374, 0.45, 0.359, 0.082, 0.685]
global q = [0.185, 0.73, 0.743, 0.936, 0.992, 0.777, 0.803, 0.958, 0.907, 0.914, 0.966, 0.89, 0.831, 0.931, 0.945, 0.971, 0.553, 0.525, 0.295, 0.88, 0.944, 0.774, 0.695, 0.69, 0.779, 0.937, 0.803, 0.769, 0.692, 0.754, 0.502, 0.449, 0.9, 0.938, 0.91, 0.493, 0.972, 0.912, 0.831, 0.713, 0.973, 0.922, 0.988, 0.296, 0.908, 0.694, 0.57, 0.783, 0.99, 0.406, 0.997, 0.945, 0.72, 0.836, 0.798, 0.992, 0.912, 0.588, 0.348, 0.417, 0.978, 0.339, 0.687, 0.792, 0.282, 0.928, 0.888, 0.537, 0.525, 0.938, 0.799, 0.542, 0.385, 0.961, 0.948, 0.958, 0.174, 0.944, 0.978, 0.595, 0.525, 0.448, 0.95, 0.778, 0.83, 0.276, 0.952, 0.572, 0.584, 0.908, 0.885, 0.93, 0.349, 0.427, 0.849, 0.913, 0.834, 0.949, 0.386, 0.995, 0.355, 0.939, 0.998, 0.873, 0.689, 0.841, 0.886, 0.821, 0.782, 0.658, 0.203, 0.998, 0.28, 0.535, 0.885, 0.795, 0.681, 0.701, 0.605, 0.555, 0.861, 0.575, 0.718, 0.8, 0.43, 0.854, 0.441, 0.494, 0.169, 0.783, 0.831, 0.885, 0.81, 0.852, 0.717, 0.965, 0.443, 0.633, 0.757, 0.996, 0.397, 0.849, 0.941, 0.927, 0.788, 0.512, 0.734, 0.566, 0.551, 0.996, 0.603, 0.58, 0.967, 0.529, 0.896, 0.845, 0.624, 0.651, 0.733, 0.949, 0.113, 0.906, 0.597, 0.942, 0.691, 0.643, 0.711, 0.7, 0.255, 0.57, 0.86, 0.55, 0.974, 0.861, 0.827, 0.842, 0.931, 0.985, 0.916, 0.934, 0.182, 0.458, 0.868, 0.973, 0.986, 0.955, 0.413, 0.803, 0.465, 0.829, 0.993, 0.351, 0.914, 0.516, 0.762, 0.656, 0.509, 0.824, 0.909, 0.865, 0.944, 0.401, 0.789, 0.836, 0.903, 0.812, 0.986, 0.947, 0.998, 0.974, 0.979, 0.973, 0.172, 0.274, 0.736, 0.698, 0.802, 0.733, 0.84, 0.972, 0.967, 0.92, 0.686, 0.488, 0.356, 0.493, 0.263, 0.439, 0.573, 0.884, 0.976, 0.949, 0.888, 0.429, 0.807, 0.406, 0.859, 0.961, 0.891, 0.638, 0.675, 0.735]
global origin = 1
global destination = 50