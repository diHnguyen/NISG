global arcs = [1 11; 1 23; 1 41; 1 44; 2 10; 3 2; 3 5; 3 23; 3 27; 3 44; 4 5; 4 28; 4 34; 4 47; 5 12; 5 32; 5 40; 5 47; 6 2; 6 12; 6 16; 6 19; 6 43; 6 47; 7 2; 7 9; 7 10; 7 13; 7 22; 7 29; 7 45; 8 7; 8 12; 8 18; 9 3; 9 12; 9 19; 9 30; 9 33; 9 50; 10 31; 10 38; 11 2; 11 5; 11 18; 11 19; 11 31; 11 36; 11 49; 12 6; 12 11; 12 46; 13 5; 13 18; 13 31; 13 35; 14 10; 14 21; 14 29; 14 48; 14 50; 15 9; 15 19; 15 23; 15 44; 16 21; 16 46; 16 47; 16 50; 17 18; 17 26; 17 46; 18 11; 18 20; 18 21; 18 27; 18 35; 18 43; 19 3; 19 6; 19 15; 19 26; 19 35; 19 42; 19 45; 20 16; 20 19; 21 7; 21 41; 21 42; 21 45; 22 4; 22 13; 22 30; 22 39; 22 41; 22 49; 23 3; 23 11; 23 13; 23 19; 24 9; 24 16; 24 21; 24 27; 24 33; 25 3; 25 22; 25 23; 25 39; 25 46; 26 37; 26 42; 26 49; 27 44; 27 45; 27 48; 28 8; 28 12; 28 13; 28 21; 28 33; 28 34; 28 35; 28 37; 28 47; 28 50; 29 6; 29 45; 29 47; 30 3; 30 24; 30 25; 30 45; 31 2; 31 11; 31 34; 31 45; 32 12; 32 15; 32 17; 32 18; 32 29; 32 33; 32 36; 32 38; 32 39; 32 43; 32 47; 33 13; 33 24; 33 26; 33 36; 34 5; 34 21; 34 25; 34 26; 35 2; 35 3; 35 4; 35 5; 35 24; 35 25; 35 40; 35 42; 36 3; 36 6; 36 29; 36 33; 36 43; 37 7; 37 10; 37 43; 37 45; 38 14; 38 28; 38 41; 38 45; 38 50; 39 6; 39 43; 39 49; 40 9; 40 17; 41 5; 41 25; 41 31; 41 34; 41 37; 42 33; 42 49; 43 13; 43 38; 43 50; 44 2; 44 41; 45 8; 45 14; 45 16; 45 41; 45 50; 46 5; 46 8; 46 23; 46 34; 46 35; 46 36; 46 37; 46 41; 47 4; 47 10; 47 12; 47 17; 47 38; 48 12; 48 14; 48 23; 48 45; 49 8; 49 24; 49 32; 49 43; 49 48]
global d_x = [8.0, 5.0, 1.0, 6.0, 9.0, 1.0, 8.0, 6.0, 3.0, 2.0, 9.0, 7.0, 9.0, 4.0, 9.0, 7.0, 1.0, 1.0, 8.0, 7.0, 5.0, 9.0, 10.0, 8.0, 8.0, 5.0, 6.0, 3.0, 8.0, 8.0, 9.0, 2.0, 2.0, 9.0, 6.0, 3.0, 3.0, 3.0, 9.0, 1.0, 7.0, 1.0, 8.0, 5.0, 7.0, 1.0, 1.0, 4.0, 7.0, 5.0, 5.0, 2.0, 2.0, 9.0, 7.0, 10.0, 2.0, 3.0, 3.0, 9.0, 9.0, 4.0, 3.0, 4.0, 10.0, 2.0, 7.0, 1.0, 9.0, 4.0, 2.0, 7.0, 10.0, 10.0, 2.0, 9.0, 6.0, 2.0, 1.0, 3.0, 7.0, 5.0, 10.0, 2.0, 8.0, 2.0, 4.0, 6.0, 1.0, 10.0, 5.0, 8.0, 4.0, 3.0, 9.0, 2.0, 1.0, 1.0, 7.0, 9.0, 2.0, 5.0, 7.0, 7.0, 3.0, 1.0, 4.0, 9.0, 4.0, 1.0, 6.0, 7.0, 8.0, 8.0, 3.0, 7.0, 3.0, 6.0, 5.0, 7.0, 8.0, 7.0, 1.0, 1.0, 10.0, 3.0, 1.0, 9.0, 10.0, 5.0, 10.0, 10.0, 4.0, 3.0, 8.0, 7.0, 9.0, 5.0, 3.0, 7.0, 1.0, 1.0, 3.0, 6.0, 3.0, 2.0, 3.0, 7.0, 3.0, 3.0, 6.0, 4.0, 8.0, 7.0, 6.0, 10.0, 8.0, 6.0, 4.0, 6.0, 2.0, 6.0, 10.0, 6.0, 7.0, 2.0, 1.0, 4.0, 9.0, 10.0, 9.0, 9.0, 4.0, 5.0, 2.0, 3.0, 1.0, 5.0, 2.0, 8.0, 9.0, 4.0, 3.0, 9.0, 3.0, 1.0, 1.0, 7.0, 6.0, 2.0, 7.0, 3.0, 6.0, 4.0, 3.0, 1.0, 4.0, 6.0, 4.0, 7.0, 3.0, 9.0, 4.0, 10.0, 3.0, 7.0, 2.0, 1.0, 7.0, 2.0, 5.0, 7.0, 6.0, 2.0, 8.0, 1.0, 10.0, 10.0, 1.0, 7.0, 5.0, 3.0, 10.0]
global b_x = 5
global d_y = [9.0, 6.0, 7.0, 5.0, 4.0, 5.0, 6.0, 5.0, 5.0, 7.0, 4.0, 3.0, 1.0, 6.0, 6.0, 8.0, 6.0, 6.0, 6.0, 6.0, 4.0, 1.0, 7.0, 9.0, 10.0, 2.0, 3.0, 7.0, 2.0, 4.0, 8.0, 10.0, 2.0, 8.0, 10.0, 8.0, 10.0, 7.0, 10.0, 9.0, 2.0, 7.0, 8.0, 6.0, 6.0, 4.0, 6.0, 7.0, 4.0, 10.0, 3.0, 9.0, 6.0, 10.0, 10.0, 10.0, 5.0, 3.0, 10.0, 6.0, 1.0, 7.0, 8.0, 9.0, 4.0, 6.0, 3.0, 2.0, 9.0, 7.0, 3.0, 10.0, 2.0, 2.0, 10.0, 9.0, 6.0, 10.0, 1.0, 4.0, 2.0, 4.0, 3.0, 10.0, 8.0, 3.0, 8.0, 4.0, 10.0, 3.0, 7.0, 10.0, 8.0, 4.0, 1.0, 9.0, 2.0, 1.0, 8.0, 2.0, 6.0, 9.0, 8.0, 9.0, 9.0, 8.0, 2.0, 8.0, 7.0, 1.0, 7.0, 1.0, 8.0, 3.0, 10.0, 8.0, 7.0, 5.0, 9.0, 10.0, 10.0, 4.0, 4.0, 1.0, 2.0, 2.0, 9.0, 4.0, 3.0, 8.0, 10.0, 8.0, 6.0, 3.0, 1.0, 6.0, 1.0, 1.0, 1.0, 2.0, 8.0, 4.0, 3.0, 4.0, 5.0, 8.0, 5.0, 8.0, 5.0, 1.0, 5.0, 8.0, 2.0, 7.0, 2.0, 2.0, 4.0, 5.0, 2.0, 5.0, 9.0, 7.0, 9.0, 9.0, 3.0, 4.0, 6.0, 1.0, 3.0, 6.0, 1.0, 3.0, 4.0, 2.0, 6.0, 2.0, 10.0, 10.0, 9.0, 9.0, 10.0, 2.0, 8.0, 5.0, 3.0, 4.0, 8.0, 5.0, 1.0, 9.0, 6.0, 10.0, 3.0, 8.0, 2.0, 10.0, 9.0, 10.0, 5.0, 2.0, 1.0, 3.0, 5.0, 8.0, 4.0, 5.0, 2.0, 4.0, 7.0, 4.0, 9.0, 6.0, 2.0, 2.0, 10.0, 6.0, 5.0, 7.0, 9.0, 4.0, 2.0, 4.0, 4.0]
global b_y = 10
global p = [0.443, 0.002, 0.725, 0.245, 0.962, 0.956, 0.454, 0.425, 0.206, 0.075, 0.957, 0.434, 0.604, 0.105, 0.834, 0.457, 0.561, 0.201, 0.333, 0.144, 0.111, 0.162, 0.335, 0.012, 0.391, 0.082, 0.93, 0.703, 0.335, 0.331, 0.143, 0.947, 0.373, 0.645, 0.389, 0.156, 0.537, 0.228, 0.211, 0.169, 0.211, 0.347, 0.077, 0.776, 0.665, 0.895, 0.044, 0.027, 0.485, 0.541, 0.78, 0.503, 0.149, 0.373, 0.521, 0.69, 0.303, 0.707, 0.106, 0.041, 0.057, 0.949, 0.7, 0.057, 0.377, 0.63, 0.152, 0.936, 0.157, 0.02, 0.562, 0.917, 0.238, 0.224, 0.34, 0.069, 0.551, 0.434, 0.338, 0.701, 0.919, 0.85, 0.116, 0.136, 0.024, 0.502, 0.662, 0.127, 0.559, 0.22, 0.997, 0.494, 0.3, 0.392, 0.595, 0.008, 0.866, 0.154, 0.796, 0.962, 0.269, 0.842, 0.895, 0.078, 0.817, 0.321, 0.079, 0.7, 0.731, 0.449, 0.901, 0.679, 0.754, 0.851, 0.438, 0.607, 0.263, 0.051, 0.235, 0.517, 0.87, 0.54, 0.602, 0.277, 0.755, 0.1, 0.687, 0.783, 0.309, 0.229, 0.314, 0.428, 0.847, 0.24, 0.962, 0.756, 0.246, 0.063, 0.383, 0.779, 0.29, 0.537, 0.606, 0.299, 0.127, 0.294, 0.319, 0.223, 0.235, 0.408, 0.382, 0.362, 0.057, 0.376, 0.981, 0.491, 0.279, 0.531, 0.894, 0.921, 0.843, 0.614, 0.691, 0.413, 0.727, 0.193, 0.709, 0.527, 0.667, 0.097, 0.102, 0.425, 0.668, 0.722, 0.312, 0.529, 0.158, 0.984, 0.748, 0.713, 0.343, 0.835, 0.61, 0.399, 0.9, 0.537, 0.139, 0.259, 0.394, 0.307, 0.073, 0.651, 0.185, 0.717, 0.273, 0.042, 0.407, 0.605, 0.497, 0.665, 0.539, 0.376, 0.892, 0.434, 0.225, 0.651, 0.596, 0.394, 0.731, 0.45, 0.781, 0.703, 0.984, 0.923, 0.12, 0.93, 0.41, 0.7, 0.075, 0.681, 0.082, 0.679, 0.314]
global q = [0.861, 0.38, 0.81, 0.649, 0.962, 0.997, 0.697, 0.514, 0.424, 0.816, 0.992, 0.905, 0.776, 0.124, 0.882, 0.887, 0.617, 0.952, 0.873, 0.351, 0.85, 0.934, 0.796, 0.417, 0.414, 0.417, 0.948, 0.899, 0.498, 0.714, 0.479, 0.952, 0.812, 0.923, 0.999, 0.631, 0.837, 0.566, 0.36, 0.765, 0.92, 0.512, 0.355, 0.842, 0.852, 0.915, 0.812, 0.377, 0.565, 0.616, 0.782, 0.987, 0.192, 0.864, 0.969, 0.873, 0.641, 0.734, 0.208, 0.454, 0.295, 0.98, 0.738, 0.735, 0.66, 0.701, 0.689, 0.966, 0.605, 0.629, 0.687, 0.932, 0.833, 0.736, 0.945, 0.653, 0.87, 0.989, 0.375, 0.744, 0.945, 0.986, 0.363, 0.621, 0.762, 0.826, 0.756, 0.62, 0.729, 0.279, 0.998, 0.545, 0.311, 0.853, 0.714, 0.372, 0.872, 0.962, 0.929, 0.963, 0.328, 0.86, 0.978, 0.263, 0.898, 0.591, 0.543, 0.83, 0.891, 0.669, 0.997, 0.911, 0.856, 0.894, 0.724, 0.869, 0.793, 0.079, 0.674, 0.775, 0.973, 0.677, 0.808, 0.526, 0.86, 0.525, 0.959, 0.824, 0.494, 0.746, 0.464, 0.463, 0.942, 0.519, 0.997, 0.926, 0.411, 0.502, 0.891, 0.929, 0.607, 0.832, 0.862, 0.456, 0.632, 0.545, 0.382, 0.358, 0.266, 0.783, 0.728, 0.94, 0.961, 0.454, 0.988, 0.951, 0.701, 0.767, 0.998, 0.986, 0.987, 0.752, 0.899, 0.573, 0.877, 0.706, 0.922, 0.745, 0.672, 0.373, 0.541, 0.611, 0.809, 0.748, 0.947, 0.619, 0.408, 0.989, 0.937, 0.734, 0.943, 0.883, 0.838, 0.824, 0.968, 0.647, 0.204, 0.66, 0.418, 0.432, 0.244, 0.835, 0.824, 0.747, 0.506, 0.929, 0.672, 0.633, 0.968, 0.748, 0.942, 0.785, 0.965, 0.62, 0.508, 0.841, 0.989, 0.705, 0.967, 0.974, 0.807, 0.882, 0.999, 0.93, 0.435, 0.95, 0.541, 0.87, 0.939, 0.707, 0.687, 0.731, 0.681]
global origin = 1
global destination = 50