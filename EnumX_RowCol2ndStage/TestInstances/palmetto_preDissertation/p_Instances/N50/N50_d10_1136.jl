global arcs = [1 4; 1 25; 1 27; 1 29; 1 39; 1 48; 2 9; 2 22; 2 46; 2 49; 3 13; 3 17; 3 19; 3 36; 4 14; 4 27; 5 11; 5 12; 5 36; 5 46; 6 5; 6 12; 6 14; 6 32; 6 37; 6 50; 7 2; 7 9; 7 37; 7 41; 7 44; 7 45; 8 7; 8 12; 8 15; 8 27; 8 35; 8 38; 8 44; 9 8; 9 13; 9 17; 9 21; 9 22; 9 36; 9 37; 9 38; 9 50; 10 11; 10 13; 10 19; 10 22; 10 28; 10 29; 10 49; 11 3; 11 7; 11 19; 11 38; 12 9; 12 11; 12 13; 12 32; 12 33; 12 34; 12 44; 12 47; 13 5; 13 20; 13 26; 13 29; 13 30; 13 46; 13 48; 14 20; 14 44; 15 17; 15 20; 15 44; 16 6; 16 11; 16 23; 16 26; 16 28; 16 39; 16 48; 17 4; 17 12; 17 20; 17 36; 17 43; 17 44; 17 45; 18 26; 18 31; 18 32; 18 47; 19 3; 19 6; 19 16; 19 31; 19 42; 19 50; 20 31; 21 2; 21 9; 21 31; 21 43; 21 45; 22 3; 22 5; 22 16; 22 29; 22 45; 22 50; 23 7; 23 26; 23 31; 23 50; 24 35; 25 7; 25 33; 26 2; 26 6; 26 10; 26 22; 26 50; 27 9; 27 15; 27 19; 27 21; 27 42; 27 45; 28 7; 28 19; 28 25; 28 33; 28 45; 29 4; 29 19; 29 38; 29 41; 30 17; 30 41; 31 8; 31 10; 31 17; 31 36; 31 38; 31 42; 32 10; 32 18; 32 24; 32 31; 32 33; 32 38; 32 39; 32 42; 32 44; 32 47; 32 48; 33 2; 33 6; 33 13; 33 40; 33 42; 33 45; 34 20; 34 22; 34 24; 34 39; 34 41; 35 7; 35 11; 35 19; 35 40; 35 48; 36 7; 36 14; 36 24; 36 28; 36 38; 36 40; 37 4; 38 6; 38 9; 38 20; 38 32; 38 47; 39 4; 39 5; 39 8; 39 11; 39 17; 40 11; 40 12; 40 17; 40 34; 40 35; 40 41; 40 42; 40 47; 40 48; 40 49; 41 6; 41 10; 41 18; 41 20; 41 31; 41 32; 41 36; 42 37; 43 11; 43 21; 43 46; 44 7; 44 18; 44 21; 44 30; 45 3; 45 6; 45 21; 45 38; 45 40; 45 41; 45 43; 46 4; 46 9; 46 40; 46 43; 47 10; 47 11; 47 13; 47 14; 47 33; 47 35; 48 2; 48 5; 48 27; 48 28; 48 30; 48 40; 49 17; 49 20; 49 29; 49 38]
global d_x = [6.0, 3.0, 3.0, 9.0, 8.0, 3.0, 4.0, 6.0, 10.0, 10.0, 2.0, 8.0, 10.0, 6.0, 10.0, 1.0, 9.0, 4.0, 1.0, 7.0, 6.0, 1.0, 10.0, 5.0, 6.0, 4.0, 3.0, 7.0, 7.0, 1.0, 6.0, 3.0, 1.0, 10.0, 7.0, 3.0, 10.0, 1.0, 9.0, 5.0, 10.0, 7.0, 6.0, 10.0, 6.0, 1.0, 1.0, 4.0, 6.0, 8.0, 4.0, 3.0, 1.0, 7.0, 8.0, 3.0, 1.0, 1.0, 9.0, 9.0, 7.0, 6.0, 3.0, 4.0, 2.0, 9.0, 1.0, 6.0, 10.0, 8.0, 1.0, 6.0, 2.0, 10.0, 3.0, 10.0, 5.0, 8.0, 7.0, 8.0, 10.0, 5.0, 8.0, 7.0, 2.0, 8.0, 9.0, 3.0, 1.0, 7.0, 7.0, 5.0, 5.0, 9.0, 6.0, 9.0, 10.0, 8.0, 1.0, 10.0, 3.0, 4.0, 9.0, 3.0, 3.0, 7.0, 4.0, 7.0, 6.0, 1.0, 9.0, 9.0, 3.0, 10.0, 8.0, 10.0, 5.0, 1.0, 8.0, 10.0, 2.0, 6.0, 4.0, 4.0, 3.0, 6.0, 5.0, 4.0, 4.0, 9.0, 1.0, 9.0, 7.0, 3.0, 9.0, 6.0, 7.0, 9.0, 5.0, 6.0, 7.0, 10.0, 9.0, 8.0, 10.0, 5.0, 2.0, 3.0, 10.0, 9.0, 8.0, 1.0, 6.0, 10.0, 4.0, 5.0, 2.0, 10.0, 1.0, 4.0, 4.0, 1.0, 7.0, 2.0, 10.0, 8.0, 10.0, 9.0, 6.0, 3.0, 5.0, 7.0, 9.0, 2.0, 9.0, 9.0, 3.0, 10.0, 6.0, 7.0, 6.0, 4.0, 7.0, 4.0, 10.0, 7.0, 10.0, 8.0, 9.0, 5.0, 4.0, 9.0, 9.0, 8.0, 4.0, 6.0, 6.0, 1.0, 8.0, 4.0, 10.0, 1.0, 2.0, 6.0, 10.0, 2.0, 1.0, 10.0, 9.0, 1.0, 8.0, 10.0, 1.0, 5.0, 5.0, 1.0, 10.0, 3.0, 7.0, 6.0, 1.0, 2.0, 10.0, 7.0, 2.0, 4.0, 9.0, 8.0, 6.0, 7.0, 8.0, 3.0, 7.0, 3.0, 6.0, 10.0, 2.0, 5.0, 2.0, 2.0, 7.0, 8.0, 6.0, 9.0, 10.0, 8.0]
global b_x = 5
global d_y = [7.0, 9.0, 9.0, 9.0, 2.0, 10.0, 7.0, 9.0, 9.0, 1.0, 3.0, 2.0, 7.0, 2.0, 8.0, 2.0, 5.0, 4.0, 3.0, 5.0, 9.0, 1.0, 2.0, 8.0, 7.0, 4.0, 10.0, 8.0, 9.0, 7.0, 10.0, 3.0, 1.0, 7.0, 2.0, 3.0, 2.0, 6.0, 10.0, 2.0, 9.0, 9.0, 8.0, 7.0, 8.0, 4.0, 7.0, 6.0, 1.0, 8.0, 1.0, 4.0, 6.0, 8.0, 9.0, 1.0, 10.0, 1.0, 6.0, 8.0, 6.0, 2.0, 8.0, 2.0, 3.0, 10.0, 9.0, 2.0, 3.0, 4.0, 9.0, 10.0, 7.0, 4.0, 10.0, 6.0, 7.0, 3.0, 8.0, 3.0, 3.0, 9.0, 10.0, 2.0, 6.0, 3.0, 2.0, 8.0, 4.0, 9.0, 5.0, 7.0, 8.0, 3.0, 3.0, 4.0, 5.0, 2.0, 2.0, 4.0, 10.0, 10.0, 7.0, 7.0, 2.0, 9.0, 10.0, 3.0, 3.0, 6.0, 7.0, 1.0, 6.0, 7.0, 9.0, 1.0, 8.0, 9.0, 9.0, 1.0, 3.0, 9.0, 3.0, 9.0, 9.0, 5.0, 3.0, 5.0, 6.0, 5.0, 8.0, 6.0, 7.0, 3.0, 6.0, 10.0, 8.0, 6.0, 2.0, 3.0, 1.0, 1.0, 8.0, 10.0, 5.0, 8.0, 5.0, 2.0, 10.0, 5.0, 2.0, 9.0, 2.0, 7.0, 9.0, 2.0, 9.0, 8.0, 8.0, 1.0, 4.0, 4.0, 4.0, 8.0, 10.0, 6.0, 6.0, 4.0, 9.0, 4.0, 4.0, 7.0, 10.0, 2.0, 3.0, 10.0, 4.0, 10.0, 7.0, 3.0, 8.0, 10.0, 8.0, 4.0, 5.0, 2.0, 10.0, 7.0, 2.0, 7.0, 3.0, 6.0, 4.0, 3.0, 6.0, 3.0, 4.0, 10.0, 9.0, 5.0, 4.0, 6.0, 7.0, 8.0, 5.0, 5.0, 3.0, 10.0, 3.0, 4.0, 9.0, 9.0, 7.0, 3.0, 6.0, 9.0, 4.0, 7.0, 4.0, 6.0, 4.0, 10.0, 3.0, 4.0, 5.0, 6.0, 3.0, 2.0, 9.0, 2.0, 2.0, 10.0, 9.0, 4.0, 6.0, 4.0, 9.0, 2.0, 4.0, 2.0, 2.0, 3.0, 2.0, 9.0, 3.0, 2.0]
global b_y = 10
global p = [0.894, 0.148, 0.293, 0.945, 0.28, 0.163, 0.93, 0.853, 0.466, 0.521, 0.11, 0.018, 0.086, 0.531, 0.443, 0.816, 0.237, 0.928, 0.492, 0.538, 0.837, 0.692, 0.774, 0.331, 0.292, 0.61, 0.854, 0.319, 0.304, 0.068, 0.98, 0.628, 0.113, 0.545, 0.291, 0.064, 0.785, 0.078, 0.67, 0.934, 0.441, 0.697, 0.717, 0.098, 0.005, 0.258, 0.567, 0.888, 0.075, 0.356, 0.886, 0.933, 0.236, 0.951, 0.045, 0.284, 0.23, 0.416, 0.612, 0.723, 0.272, 0.008, 0.112, 0.963, 0.157, 0.154, 0.886, 0.641, 0.64, 0.052, 0.7, 0.969, 0.599, 0.33, 0.335, 0.054, 0.577, 0.862, 0.23, 0.025, 0.2, 0.516, 0.176, 0.497, 0.084, 0.917, 0.304, 0.525, 0.814, 0.903, 0.668, 0.399, 0.308, 0.276, 0.95, 0.136, 0.087, 0.6, 0.015, 0.982, 0.011, 0.446, 0.05, 0.361, 0.119, 0.245, 0.227, 0.849, 0.803, 0.738, 0.931, 0.996, 0.898, 0.693, 0.59, 0.198, 0.53, 0.263, 0.189, 0.985, 0.563, 0.913, 0.104, 0.38, 0.205, 0.877, 0.867, 0.559, 0.085, 0.7, 0.384, 0.723, 0.408, 0.721, 0.966, 0.335, 0.822, 0.33, 0.679, 0.32, 0.213, 0.691, 0.072, 0.831, 0.282, 0.833, 0.279, 0.039, 0.961, 0.841, 0.129, 0.879, 0.164, 0.565, 0.707, 0.614, 0.494, 0.44, 0.682, 0.135, 0.647, 0.492, 0.42, 0.267, 0.792, 0.571, 0.523, 0.407, 0.101, 0.532, 0.519, 0.61, 0.026, 0.248, 0.481, 0.868, 0.401, 0.35, 0.468, 0.855, 0.577, 0.734, 0.211, 0.812, 0.727, 0.4, 0.41, 0.666, 0.241, 0.724, 0.459, 0.425, 0.877, 0.624, 0.579, 0.382, 0.893, 0.536, 0.565, 0.717, 0.158, 0.859, 0.118, 0.903, 0.701, 0.905, 0.295, 0.224, 0.179, 0.302, 0.546, 0.503, 0.178, 0.057, 0.301, 0.363, 0.42, 0.91, 0.426, 0.372, 0.003, 0.5, 0.854, 0.94, 0.403, 0.763, 0.614, 0.575, 0.783, 0.809, 0.412, 0.738, 0.767, 0.578, 0.962, 0.807, 0.34, 0.244, 0.224, 0.936, 0.336, 0.81, 0.959, 0.891, 0.238, 0.428]
global q = [0.996, 0.68, 0.513, 0.948, 0.475, 0.449, 0.996, 0.859, 0.546, 0.544, 0.807, 0.517, 0.991, 0.915, 0.553, 0.842, 0.345, 0.993, 0.943, 0.972, 0.906, 0.847, 0.833, 0.488, 0.331, 0.724, 0.868, 0.561, 0.826, 0.267, 0.987, 0.804, 0.395, 0.87, 0.745, 0.49, 0.815, 0.218, 0.913, 0.992, 0.959, 0.775, 0.736, 0.834, 0.038, 0.796, 0.816, 0.991, 0.845, 0.655, 0.981, 0.944, 0.527, 0.967, 0.745, 0.536, 0.339, 0.527, 0.698, 0.916, 0.445, 0.251, 0.732, 0.975, 0.862, 0.976, 0.889, 0.739, 0.728, 0.732, 0.998, 0.993, 0.836, 0.493, 0.655, 0.074, 0.737, 0.979, 0.347, 0.378, 0.781, 0.85, 0.301, 0.811, 0.902, 0.992, 0.361, 0.665, 0.966, 0.914, 0.912, 0.644, 0.863, 0.706, 0.98, 0.451, 0.502, 0.615, 0.43, 0.996, 0.185, 0.93, 0.144, 0.605, 0.539, 0.834, 0.577, 0.876, 0.888, 0.928, 0.96, 0.999, 0.932, 0.796, 0.799, 0.205, 0.905, 0.954, 0.632, 0.997, 0.76, 0.993, 0.917, 0.887, 0.545, 0.9, 0.932, 0.589, 0.348, 0.936, 0.742, 0.899, 0.655, 0.73, 0.998, 0.942, 0.975, 0.981, 0.777, 0.944, 0.936, 0.741, 0.074, 0.939, 0.63, 0.862, 0.862, 0.369, 0.99, 0.9, 0.372, 0.999, 0.399, 0.649, 0.783, 0.735, 0.842, 0.902, 0.726, 0.42, 0.965, 0.769, 0.859, 0.293, 0.918, 0.634, 0.523, 0.733, 0.416, 0.713, 0.563, 0.98, 0.26, 0.665, 0.956, 0.894, 0.726, 0.629, 0.911, 0.921, 0.822, 0.971, 0.884, 0.917, 0.804, 0.903, 0.588, 0.757, 0.335, 0.725, 0.497, 0.445, 0.888, 0.848, 0.623, 0.474, 0.926, 0.928, 0.77, 0.929, 0.386, 0.978, 0.516, 0.903, 0.799, 0.918, 0.379, 0.293, 0.33, 0.617, 0.683, 0.899, 0.313, 0.223, 0.401, 0.429, 0.541, 0.968, 0.927, 0.441, 0.042, 0.552, 0.89, 0.986, 0.472, 0.989, 0.764, 0.905, 0.809, 0.944, 0.423, 0.901, 0.828, 0.944, 0.963, 0.869, 0.368, 0.368, 0.469, 0.938, 0.473, 0.837, 0.985, 0.929, 0.62, 0.602]
global origin = 1
global destination = 50