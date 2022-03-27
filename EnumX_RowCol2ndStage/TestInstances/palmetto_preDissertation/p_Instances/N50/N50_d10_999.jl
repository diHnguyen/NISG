global arcs = [1 6; 1 18; 1 40; 1 42; 2 11; 2 33; 2 44; 2 45; 3 4; 3 5; 3 13; 3 22; 3 41; 3 49; 4 31; 4 43; 5 3; 5 26; 6 2; 6 7; 6 10; 6 16; 6 34; 6 37; 6 49; 7 15; 7 21; 7 32; 7 34; 7 35; 7 47; 8 4; 8 10; 8 20; 8 21; 8 22; 8 26; 8 34; 9 21; 9 36; 9 40; 9 42; 9 46; 9 50; 10 2; 10 16; 10 20; 10 40; 10 41; 10 46; 10 49; 11 22; 12 11; 12 22; 12 28; 12 29; 12 47; 13 14; 13 18; 13 20; 13 34; 14 7; 14 9; 14 19; 14 23; 14 26; 14 39; 15 8; 15 17; 15 29; 15 30; 15 43; 16 12; 16 23; 16 29; 16 32; 16 47; 17 22; 17 31; 18 3; 18 5; 18 8; 18 33; 18 34; 18 35; 18 39; 18 42; 18 43; 18 50; 19 2; 19 12; 19 24; 19 47; 20 4; 20 12; 20 15; 20 28; 20 32; 20 46; 21 2; 21 12; 21 24; 22 13; 22 21; 22 26; 22 37; 22 39; 22 44; 23 6; 23 13; 23 17; 23 25; 23 40; 24 6; 24 12; 24 22; 24 34; 24 36; 25 19; 25 29; 25 39; 25 45; 26 14; 26 16; 26 31; 26 41; 26 43; 26 44; 26 50; 27 18; 27 23; 27 32; 27 36; 27 39; 27 41; 27 42; 28 10; 28 12; 28 20; 28 22; 28 33; 28 41; 29 3; 29 5; 29 7; 29 32; 29 39; 29 40; 30 16; 30 17; 30 20; 30 45; 30 46; 31 5; 31 9; 31 11; 31 22; 31 34; 31 39; 32 11; 32 12; 32 15; 32 25; 32 33; 32 35; 32 39; 33 2; 33 7; 33 11; 33 30; 33 36; 33 39; 33 46; 34 10; 34 47; 34 50; 35 8; 35 13; 35 29; 35 43; 35 47; 36 15; 36 16; 36 27; 36 46; 37 7; 37 8; 37 43; 37 47; 38 16; 38 31; 38 40; 39 3; 39 4; 39 13; 39 27; 40 28; 40 48; 40 49; 41 3; 41 10; 41 33; 41 37; 42 8; 42 9; 42 19; 42 39; 43 8; 43 11; 43 13; 43 15; 43 18; 43 20; 43 23; 43 25; 43 30; 43 38; 44 23; 44 25; 44 35; 44 40; 44 50; 45 2; 45 12; 45 14; 45 35; 45 50; 46 20; 46 26; 46 34; 46 47; 47 3; 47 15; 47 21; 47 28; 47 48; 47 49; 48 30; 48 50; 49 4; 49 13; 49 16; 49 24; 49 32; 49 43]
global d_x = [2.0, 5.0, 6.0, 1.0, 5.0, 10.0, 4.0, 7.0, 9.0, 9.0, 10.0, 10.0, 1.0, 4.0, 9.0, 7.0, 3.0, 2.0, 2.0, 6.0, 6.0, 1.0, 10.0, 10.0, 8.0, 2.0, 3.0, 8.0, 4.0, 5.0, 4.0, 9.0, 9.0, 1.0, 5.0, 8.0, 4.0, 10.0, 1.0, 6.0, 1.0, 5.0, 4.0, 7.0, 8.0, 9.0, 7.0, 5.0, 9.0, 1.0, 4.0, 7.0, 1.0, 5.0, 1.0, 6.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 6.0, 2.0, 8.0, 9.0, 5.0, 2.0, 8.0, 7.0, 3.0, 2.0, 6.0, 7.0, 3.0, 3.0, 2.0, 3.0, 2.0, 6.0, 10.0, 3.0, 4.0, 2.0, 8.0, 2.0, 7.0, 4.0, 9.0, 2.0, 3.0, 1.0, 3.0, 9.0, 3.0, 3.0, 8.0, 3.0, 1.0, 8.0, 1.0, 5.0, 9.0, 9.0, 8.0, 9.0, 10.0, 1.0, 9.0, 3.0, 7.0, 10.0, 3.0, 3.0, 6.0, 2.0, 6.0, 4.0, 4.0, 2.0, 9.0, 9.0, 5.0, 7.0, 1.0, 2.0, 1.0, 2.0, 7.0, 10.0, 2.0, 3.0, 6.0, 7.0, 2.0, 3.0, 6.0, 3.0, 7.0, 10.0, 1.0, 7.0, 7.0, 9.0, 8.0, 3.0, 1.0, 1.0, 3.0, 3.0, 1.0, 7.0, 8.0, 9.0, 9.0, 6.0, 1.0, 7.0, 4.0, 5.0, 8.0, 10.0, 5.0, 5.0, 9.0, 10.0, 2.0, 6.0, 5.0, 6.0, 5.0, 7.0, 1.0, 3.0, 10.0, 7.0, 3.0, 8.0, 10.0, 4.0, 1.0, 8.0, 7.0, 10.0, 4.0, 5.0, 5.0, 9.0, 10.0, 10.0, 7.0, 1.0, 1.0, 4.0, 4.0, 8.0, 1.0, 3.0, 7.0, 3.0, 7.0, 1.0, 10.0, 1.0, 8.0, 3.0, 6.0, 10.0, 1.0, 7.0, 10.0, 9.0, 7.0, 5.0, 10.0, 2.0, 1.0, 6.0, 3.0, 1.0, 7.0, 6.0, 3.0, 6.0, 9.0, 4.0, 5.0, 4.0, 8.0, 9.0, 4.0, 3.0, 8.0, 5.0, 3.0, 5.0, 9.0, 5.0, 9.0, 3.0, 8.0, 10.0, 8.0, 7.0, 5.0]
global b_x = 5
global d_y = [1.0, 4.0, 7.0, 5.0, 7.0, 4.0, 6.0, 3.0, 1.0, 2.0, 7.0, 8.0, 1.0, 3.0, 10.0, 8.0, 6.0, 8.0, 6.0, 7.0, 5.0, 6.0, 4.0, 5.0, 2.0, 5.0, 3.0, 5.0, 6.0, 8.0, 10.0, 7.0, 8.0, 8.0, 2.0, 8.0, 5.0, 2.0, 5.0, 6.0, 8.0, 1.0, 10.0, 7.0, 10.0, 6.0, 3.0, 3.0, 7.0, 2.0, 8.0, 3.0, 1.0, 9.0, 9.0, 4.0, 7.0, 3.0, 6.0, 8.0, 7.0, 6.0, 5.0, 8.0, 10.0, 4.0, 2.0, 10.0, 1.0, 4.0, 8.0, 8.0, 8.0, 4.0, 1.0, 4.0, 6.0, 10.0, 1.0, 3.0, 6.0, 3.0, 5.0, 8.0, 1.0, 9.0, 8.0, 6.0, 8.0, 5.0, 6.0, 3.0, 6.0, 10.0, 7.0, 9.0, 2.0, 3.0, 1.0, 3.0, 10.0, 1.0, 6.0, 2.0, 7.0, 4.0, 6.0, 3.0, 9.0, 1.0, 3.0, 3.0, 9.0, 6.0, 8.0, 9.0, 4.0, 8.0, 3.0, 10.0, 3.0, 3.0, 10.0, 5.0, 10.0, 4.0, 5.0, 10.0, 4.0, 4.0, 4.0, 7.0, 8.0, 1.0, 3.0, 2.0, 5.0, 7.0, 1.0, 9.0, 8.0, 4.0, 3.0, 3.0, 10.0, 5.0, 3.0, 6.0, 6.0, 2.0, 3.0, 3.0, 10.0, 8.0, 3.0, 3.0, 2.0, 8.0, 1.0, 7.0, 5.0, 2.0, 4.0, 3.0, 6.0, 2.0, 7.0, 9.0, 5.0, 5.0, 5.0, 6.0, 3.0, 2.0, 2.0, 2.0, 2.0, 4.0, 4.0, 7.0, 7.0, 3.0, 3.0, 7.0, 8.0, 6.0, 5.0, 8.0, 3.0, 3.0, 2.0, 8.0, 3.0, 4.0, 8.0, 3.0, 1.0, 6.0, 6.0, 4.0, 8.0, 4.0, 7.0, 8.0, 7.0, 7.0, 5.0, 4.0, 9.0, 2.0, 10.0, 9.0, 4.0, 1.0, 9.0, 4.0, 6.0, 4.0, 7.0, 3.0, 2.0, 1.0, 6.0, 8.0, 2.0, 5.0, 8.0, 7.0, 4.0, 8.0, 8.0, 8.0, 5.0, 2.0, 2.0, 8.0, 7.0, 8.0, 6.0, 5.0, 7.0, 1.0, 2.0, 1.0, 9.0]
global b_y = 10
global p = [0.05, 0.808, 0.043, 0.602, 0.8, 0.37, 0.812, 0.252, 0.6, 0.28, 0.175, 0.872, 0.412, 0.815, 0.094, 0.453, 0.58, 0.895, 0.662, 0.144, 0.082, 0.882, 0.33, 0.196, 0.932, 0.624, 0.787, 0.246, 0.8, 0.332, 0.601, 0.887, 0.552, 0.297, 0.548, 0.649, 0.108, 0.476, 0.015, 0.228, 0.901, 0.942, 0.446, 0.323, 0.32, 0.243, 0.701, 0.658, 0.131, 0.445, 0.317, 0.471, 0.079, 0.804, 0.259, 0.537, 0.296, 0.073, 0.491, 0.246, 0.782, 0.247, 0.201, 0.506, 0.093, 0.453, 0.883, 0.665, 0.889, 0.445, 0.074, 0.961, 0.942, 0.132, 0.011, 0.518, 0.305, 0.606, 0.216, 0.156, 0.352, 0.369, 0.921, 0.376, 0.339, 0.423, 0.213, 0.152, 0.844, 0.124, 0.955, 0.066, 0.801, 0.48, 0.007, 0.007, 0.717, 0.84, 0.068, 0.21, 0.035, 0.087, 0.561, 0.589, 0.095, 0.703, 0.198, 0.629, 0.312, 0.312, 0.888, 0.174, 0.32, 0.393, 0.755, 0.861, 0.784, 0.319, 0.646, 0.403, 0.784, 0.975, 0.013, 0.807, 0.109, 0.834, 0.048, 0.054, 0.77, 0.758, 0.393, 0.187, 0.18, 0.746, 0.264, 0.177, 0.788, 0.968, 0.705, 0.952, 0.301, 0.629, 0.766, 0.586, 0.074, 0.256, 0.998, 0.065, 0.127, 0.601, 0.988, 0.094, 0.902, 0.345, 0.818, 0.404, 0.991, 0.825, 0.059, 0.648, 0.507, 0.192, 0.39, 0.798, 0.028, 0.073, 0.6, 0.139, 0.211, 0.068, 0.979, 0.876, 0.653, 0.794, 0.755, 0.591, 0.864, 0.663, 0.089, 0.391, 0.487, 0.726, 0.973, 0.626, 0.745, 0.567, 0.762, 0.567, 0.21, 0.348, 0.354, 0.024, 0.377, 0.11, 0.876, 0.469, 0.984, 0.251, 0.811, 0.873, 0.476, 0.912, 0.005, 0.298, 0.295, 0.648, 0.932, 0.619, 0.176, 0.435, 0.895, 0.006, 0.848, 0.333, 0.027, 0.023, 0.269, 0.221, 0.35, 0.44, 0.563, 0.112, 0.369, 0.767, 0.614, 0.313, 0.121, 0.085, 0.819, 0.857, 0.297, 0.023, 0.446, 0.918, 0.793, 0.92, 0.851, 0.206, 0.846, 0.729, 0.134, 0.76, 0.197, 0.363, 0.683]
global q = [0.658, 0.945, 0.589, 0.699, 0.819, 0.656, 0.84, 0.273, 0.858, 0.873, 0.791, 0.92, 0.869, 0.924, 0.314, 0.955, 0.722, 0.973, 0.93, 0.867, 0.12, 0.999, 0.938, 0.709, 0.969, 0.7, 0.974, 0.96, 0.816, 0.976, 0.937, 0.936, 0.787, 0.748, 0.559, 0.872, 0.911, 0.571, 0.328, 0.633, 0.904, 0.949, 0.863, 0.331, 0.541, 0.934, 0.802, 0.815, 0.823, 0.649, 0.796, 0.91, 0.216, 0.812, 0.899, 0.773, 0.529, 0.411, 0.608, 0.753, 0.909, 0.514, 0.872, 0.573, 0.565, 0.903, 0.915, 0.917, 0.938, 0.854, 0.108, 0.967, 0.984, 0.921, 0.081, 0.902, 0.771, 0.824, 0.694, 0.173, 0.505, 0.727, 0.962, 0.823, 0.781, 0.508, 0.969, 0.399, 0.918, 0.574, 0.981, 0.804, 0.931, 0.73, 0.222, 0.235, 0.991, 0.849, 0.502, 0.981, 0.73, 0.988, 0.848, 0.748, 0.716, 0.926, 0.285, 0.672, 0.34, 0.903, 0.99, 0.405, 0.993, 0.591, 0.865, 0.867, 0.862, 0.361, 0.673, 0.991, 0.857, 0.977, 0.657, 0.887, 0.5, 0.92, 0.312, 0.352, 0.9, 0.803, 0.556, 0.511, 0.486, 0.905, 0.445, 0.753, 0.977, 0.983, 0.883, 0.992, 0.696, 0.825, 0.833, 0.788, 0.629, 0.767, 0.998, 0.732, 0.259, 0.927, 0.997, 0.978, 0.956, 0.516, 0.949, 0.603, 0.997, 0.83, 0.437, 0.843, 0.953, 0.193, 0.421, 0.952, 0.966, 0.252, 0.797, 0.77, 0.905, 0.482, 0.99, 0.937, 0.786, 0.995, 0.897, 0.66, 0.921, 0.872, 0.789, 0.946, 0.693, 0.998, 0.979, 0.638, 0.944, 0.653, 0.989, 0.802, 0.317, 0.402, 0.942, 0.109, 0.718, 0.684, 0.985, 0.94, 0.988, 0.586, 0.937, 0.969, 0.975, 0.924, 0.603, 0.569, 0.925, 0.675, 0.959, 0.716, 0.28, 0.987, 0.921, 0.924, 0.984, 0.69, 0.359, 0.74, 0.567, 0.649, 0.945, 0.696, 0.598, 0.517, 0.945, 0.95, 0.956, 0.546, 0.448, 0.503, 0.935, 0.907, 0.426, 0.133, 0.601, 0.933, 0.965, 0.991, 0.866, 0.737, 0.955, 0.937, 0.842, 0.879, 0.268, 0.935, 0.884]
global origin = 1
global destination = 50