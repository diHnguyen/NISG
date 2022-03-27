global arcs = [1 4; 1 20; 2 5; 2 10; 2 15; 2 37; 2 49; 3 7; 3 14; 3 16; 3 41; 3 43; 4 8; 4 22; 5 29; 5 37; 6 12; 6 35; 6 40; 6 41; 6 50; 7 22; 7 28; 7 36; 7 39; 8 22; 8 36; 8 41; 8 49; 9 3; 9 6; 9 49; 10 2; 10 13; 10 27; 11 13; 11 18; 11 21; 11 22; 11 32; 11 34; 12 2; 12 7; 12 27; 12 35; 12 41; 13 18; 13 19; 13 25; 13 47; 13 49; 14 3; 14 7; 14 15; 14 28; 14 40; 14 49; 15 3; 15 4; 15 6; 15 16; 15 26; 16 9; 16 25; 16 26; 16 27; 16 30; 16 36; 17 4; 17 27; 17 30; 17 32; 17 43; 17 47; 17 50; 18 17; 18 35; 18 36; 18 40; 18 42; 18 49; 19 4; 19 7; 19 26; 19 41; 19 44; 20 7; 20 12; 20 19; 20 21; 20 24; 20 29; 20 35; 21 2; 21 4; 21 12; 21 23; 21 24; 21 31; 21 33; 21 50; 22 14; 22 25; 22 46; 23 11; 23 39; 24 31; 24 48; 25 6; 25 10; 25 19; 25 29; 25 31; 26 6; 26 7; 26 12; 26 29; 26 38; 26 39; 27 2; 27 6; 27 9; 27 10; 27 23; 27 29; 27 31; 27 35; 27 39; 27 45; 27 47; 28 5; 28 8; 28 30; 28 44; 28 45; 28 47; 29 9; 29 10; 29 20; 29 47; 29 50; 30 8; 30 15; 30 24; 30 32; 30 47; 30 48; 31 6; 31 7; 31 12; 31 18; 31 20; 31 40; 31 44; 32 12; 32 17; 32 22; 32 44; 32 46; 33 7; 33 10; 33 28; 33 32; 33 47; 34 13; 34 16; 34 37; 34 43; 35 9; 35 16; 35 20; 35 26; 35 41; 36 15; 36 16; 36 39; 36 45; 37 5; 37 11; 37 30; 37 32; 37 44; 38 2; 38 11; 38 25; 38 35; 38 44; 38 47; 38 48; 39 7; 39 10; 39 35; 39 37; 40 11; 40 36; 40 37; 40 39; 40 41; 40 44; 41 3; 42 10; 42 16; 42 17; 42 25; 42 30; 42 38; 42 41; 42 45; 43 8; 43 17; 43 31; 44 16; 44 18; 44 20; 44 28; 44 31; 45 12; 45 17; 45 18; 45 30; 45 36; 46 5; 46 31; 46 32; 47 2; 47 8; 47 9; 47 21; 47 30; 47 38; 47 46; 48 7; 48 10; 48 29; 48 31; 48 33; 48 45; 49 43]
global d_x = [6.0, 7.0, 9.0, 8.0, 5.0, 9.0, 6.0, 4.0, 7.0, 9.0, 9.0, 9.0, 4.0, 7.0, 8.0, 8.0, 5.0, 8.0, 5.0, 9.0, 2.0, 7.0, 7.0, 3.0, 8.0, 10.0, 6.0, 1.0, 10.0, 6.0, 10.0, 1.0, 5.0, 8.0, 6.0, 5.0, 4.0, 10.0, 3.0, 5.0, 1.0, 6.0, 5.0, 8.0, 7.0, 2.0, 3.0, 4.0, 2.0, 1.0, 5.0, 5.0, 4.0, 5.0, 2.0, 6.0, 1.0, 2.0, 3.0, 7.0, 10.0, 6.0, 5.0, 4.0, 6.0, 5.0, 7.0, 4.0, 5.0, 1.0, 7.0, 7.0, 8.0, 2.0, 2.0, 3.0, 4.0, 1.0, 10.0, 6.0, 9.0, 7.0, 4.0, 9.0, 6.0, 8.0, 3.0, 2.0, 6.0, 6.0, 6.0, 10.0, 4.0, 6.0, 5.0, 8.0, 7.0, 5.0, 6.0, 1.0, 9.0, 5.0, 7.0, 3.0, 6.0, 10.0, 10.0, 4.0, 7.0, 7.0, 2.0, 7.0, 7.0, 9.0, 10.0, 8.0, 9.0, 5.0, 5.0, 1.0, 6.0, 3.0, 1.0, 10.0, 8.0, 8.0, 2.0, 6.0, 2.0, 2.0, 7.0, 9.0, 2.0, 3.0, 5.0, 6.0, 9.0, 3.0, 10.0, 2.0, 8.0, 7.0, 1.0, 5.0, 3.0, 9.0, 6.0, 3.0, 10.0, 4.0, 2.0, 6.0, 3.0, 3.0, 6.0, 7.0, 4.0, 7.0, 2.0, 10.0, 9.0, 4.0, 3.0, 5.0, 5.0, 7.0, 5.0, 9.0, 4.0, 5.0, 3.0, 10.0, 7.0, 4.0, 4.0, 4.0, 8.0, 9.0, 1.0, 7.0, 6.0, 8.0, 2.0, 6.0, 9.0, 8.0, 9.0, 3.0, 7.0, 2.0, 5.0, 3.0, 1.0, 2.0, 3.0, 10.0, 10.0, 8.0, 5.0, 2.0, 10.0, 2.0, 2.0, 4.0, 2.0, 6.0, 5.0, 6.0, 9.0, 9.0, 2.0, 9.0, 6.0, 2.0, 6.0, 3.0, 8.0, 2.0, 4.0, 9.0, 8.0, 9.0, 3.0, 9.0, 2.0, 5.0, 3.0, 4.0, 10.0, 4.0, 8.0, 2.0, 8.0, 4.0, 4.0, 4.0, 2.0, 7.0]
global b_x = 5
global d_y = [5.0, 10.0, 6.0, 7.0, 8.0, 8.0, 3.0, 7.0, 5.0, 6.0, 7.0, 2.0, 1.0, 2.0, 8.0, 4.0, 10.0, 8.0, 8.0, 3.0, 5.0, 6.0, 8.0, 5.0, 1.0, 2.0, 7.0, 2.0, 3.0, 6.0, 10.0, 5.0, 3.0, 1.0, 3.0, 1.0, 10.0, 4.0, 4.0, 3.0, 1.0, 9.0, 10.0, 7.0, 10.0, 10.0, 8.0, 1.0, 1.0, 1.0, 5.0, 2.0, 7.0, 4.0, 5.0, 5.0, 10.0, 4.0, 7.0, 3.0, 6.0, 4.0, 8.0, 4.0, 3.0, 8.0, 10.0, 9.0, 8.0, 4.0, 2.0, 1.0, 9.0, 4.0, 1.0, 10.0, 8.0, 9.0, 5.0, 8.0, 6.0, 6.0, 4.0, 10.0, 1.0, 10.0, 8.0, 6.0, 9.0, 6.0, 9.0, 3.0, 4.0, 6.0, 10.0, 9.0, 1.0, 6.0, 9.0, 4.0, 9.0, 1.0, 4.0, 2.0, 8.0, 2.0, 9.0, 10.0, 5.0, 9.0, 5.0, 10.0, 4.0, 10.0, 1.0, 10.0, 9.0, 3.0, 8.0, 1.0, 4.0, 8.0, 4.0, 1.0, 2.0, 10.0, 9.0, 7.0, 4.0, 7.0, 6.0, 4.0, 4.0, 1.0, 2.0, 5.0, 7.0, 3.0, 7.0, 2.0, 2.0, 5.0, 6.0, 2.0, 3.0, 6.0, 10.0, 9.0, 3.0, 9.0, 10.0, 5.0, 7.0, 8.0, 3.0, 7.0, 2.0, 4.0, 7.0, 2.0, 3.0, 3.0, 8.0, 6.0, 4.0, 3.0, 8.0, 2.0, 2.0, 7.0, 5.0, 6.0, 6.0, 3.0, 8.0, 7.0, 2.0, 9.0, 2.0, 2.0, 1.0, 9.0, 1.0, 6.0, 9.0, 4.0, 6.0, 9.0, 4.0, 7.0, 10.0, 9.0, 5.0, 10.0, 7.0, 2.0, 2.0, 5.0, 4.0, 3.0, 10.0, 1.0, 7.0, 9.0, 8.0, 4.0, 7.0, 7.0, 7.0, 5.0, 9.0, 10.0, 1.0, 8.0, 7.0, 9.0, 7.0, 4.0, 3.0, 5.0, 9.0, 4.0, 8.0, 4.0, 2.0, 7.0, 2.0, 3.0, 1.0, 7.0, 7.0, 6.0, 7.0, 7.0, 9.0, 6.0, 8.0, 8.0]
global b_y = 10
global p = [0.979, 0.33, 0.549, 0.623, 0.226, 0.061, 0.08, 0.123, 0.844, 0.473, 0.941, 0.692, 0.45, 0.488, 0.052, 0.454, 0.354, 0.297, 0.321, 0.741, 0.697, 0.004, 0.345, 0.841, 0.326, 0.438, 0.632, 0.263, 0.725, 0.197, 0.814, 0.475, 0.05, 0.728, 0.661, 0.979, 0.751, 0.838, 0.607, 0.946, 0.991, 0.981, 0.513, 0.17, 0.47, 0.334, 0.588, 0.3, 0.944, 0.076, 0.786, 0.166, 0.905, 0.521, 0.225, 0.948, 0.439, 0.522, 0.582, 0.36, 0.113, 0.454, 0.797, 0.942, 0.973, 0.08, 0.155, 0.35, 0.134, 0.366, 0.023, 0.697, 0.858, 0.251, 0.27, 0.355, 0.908, 0.175, 0.056, 0.625, 0.742, 0.759, 0.572, 0.293, 0.279, 0.762, 0.684, 0.208, 0.623, 0.885, 0.071, 0.751, 0.157, 0.904, 0.642, 0.982, 0.905, 0.545, 0.801, 0.3, 0.481, 0.324, 0.759, 0.159, 0.686, 0.909, 0.964, 0.735, 0.431, 0.238, 0.667, 0.407, 0.463, 0.228, 0.35, 0.626, 0.092, 0.016, 0.729, 0.138, 0.652, 0.111, 0.38, 0.524, 0.605, 0.535, 0.994, 0.506, 0.369, 0.412, 0.157, 0.829, 0.239, 0.217, 0.177, 0.465, 0.294, 0.289, 0.429, 0.877, 0.993, 0.327, 0.562, 0.583, 0.209, 0.614, 0.744, 0.998, 0.925, 0.178, 0.911, 0.181, 0.842, 0.318, 0.085, 0.278, 0.09, 0.889, 0.428, 0.605, 0.91, 0.333, 0.999, 0.312, 0.989, 0.777, 0.496, 0.079, 0.964, 0.254, 0.947, 0.797, 0.564, 0.033, 0.408, 0.102, 0.285, 0.675, 0.392, 0.988, 0.804, 0.2, 0.585, 0.161, 0.788, 0.137, 0.356, 0.681, 0.337, 0.726, 0.184, 0.341, 0.057, 0.443, 0.22, 0.473, 0.056, 0.207, 0.416, 0.006, 0.521, 0.857, 0.347, 0.482, 0.836, 0.641, 0.869, 0.461, 0.061, 0.672, 0.491, 0.795, 0.76, 0.924, 0.795, 0.853, 0.591, 0.801, 0.28, 0.771, 0.093, 0.021, 0.356, 0.222, 0.015, 0.124, 0.962, 0.172, 0.715, 0.052, 0.359, 0.655, 0.061, 0.032, 0.906, 0.629, 0.087, 0.615]
global q = [0.987, 0.835, 0.653, 0.903, 0.505, 0.501, 0.405, 0.425, 0.859, 0.72, 0.948, 0.878, 0.588, 0.557, 0.719, 0.717, 0.558, 0.936, 0.675, 0.929, 0.782, 0.093, 0.507, 0.971, 0.864, 0.666, 0.94, 0.875, 0.933, 0.34, 0.908, 0.524, 0.624, 0.976, 0.803, 0.992, 0.881, 0.978, 0.936, 0.994, 0.998, 0.982, 0.738, 0.197, 0.565, 0.907, 0.647, 0.604, 0.994, 0.419, 0.941, 0.957, 0.925, 0.883, 0.235, 0.983, 0.866, 0.912, 0.895, 0.895, 0.877, 0.748, 0.952, 0.983, 0.985, 0.819, 0.691, 0.403, 0.143, 0.737, 0.995, 0.942, 0.864, 0.969, 0.786, 0.914, 0.975, 0.35, 0.64, 0.632, 0.747, 0.996, 0.974, 0.701, 0.502, 0.768, 0.714, 0.758, 0.729, 0.907, 0.088, 0.843, 0.603, 0.905, 0.82, 0.983, 0.929, 0.767, 0.919, 0.419, 0.509, 0.867, 0.953, 0.924, 0.79, 0.942, 0.97, 0.932, 0.49, 0.362, 0.765, 0.47, 0.532, 0.928, 0.571, 0.912, 0.993, 0.19, 0.998, 0.463, 0.816, 0.569, 0.762, 0.815, 0.828, 0.86, 0.995, 0.76, 0.651, 0.653, 0.756, 0.97, 0.419, 0.526, 0.839, 0.898, 0.623, 0.93, 0.504, 0.947, 0.993, 0.375, 0.708, 0.981, 0.679, 0.687, 0.866, 0.999, 0.935, 0.297, 0.979, 0.725, 0.976, 0.591, 0.909, 0.713, 0.355, 0.944, 0.447, 0.856, 0.92, 0.416, 0.999, 0.728, 0.991, 0.827, 0.984, 0.111, 0.981, 0.274, 0.978, 0.868, 0.718, 0.833, 0.828, 0.548, 0.469, 0.712, 0.455, 0.991, 0.807, 0.774, 0.954, 0.238, 0.83, 0.533, 0.676, 0.711, 0.376, 0.899, 0.444, 0.415, 0.864, 0.949, 0.763, 0.483, 0.715, 0.922, 0.615, 0.154, 0.891, 0.939, 0.876, 0.872, 0.985, 0.978, 0.96, 0.645, 0.576, 0.815, 0.544, 0.923, 0.945, 0.996, 0.846, 0.862, 0.84, 0.908, 0.725, 0.98, 0.758, 0.82, 0.471, 0.574, 0.726, 0.289, 0.965, 0.725, 0.886, 0.909, 0.481, 0.97, 0.584, 0.945, 0.946, 0.684, 0.412, 0.622]
global origin = 1
global destination = 50