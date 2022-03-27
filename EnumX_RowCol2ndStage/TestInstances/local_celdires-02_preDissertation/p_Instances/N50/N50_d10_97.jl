global arcs = [1 7; 1 22; 1 31; 1 35; 1 36; 1 40; 1 41; 2 7; 2 12; 2 15; 2 16; 2 22; 2 36; 2 39; 2 42; 2 43; 2 48; 3 8; 3 10; 3 21; 3 26; 3 31; 3 41; 3 48; 3 50; 4 2; 4 13; 4 17; 4 49; 5 4; 5 11; 5 16; 5 35; 6 10; 6 14; 6 49; 7 15; 7 29; 7 40; 7 41; 7 48; 8 3; 8 4; 8 13; 8 23; 8 31; 8 32; 9 3; 9 18; 9 20; 9 37; 9 45; 9 48; 10 2; 10 13; 10 25; 10 46; 11 7; 11 12; 11 15; 11 18; 11 24; 11 26; 11 28; 11 34; 11 35; 11 41; 11 45; 12 19; 12 27; 12 30; 12 33; 12 47; 13 26; 13 27; 13 35; 13 41; 13 43; 13 46; 14 2; 14 5; 14 7; 14 16; 14 28; 14 35; 15 6; 15 13; 15 17; 15 20; 15 27; 15 34; 15 43; 15 45; 16 10; 16 21; 16 26; 16 28; 16 36; 16 42; 16 43; 16 50; 17 2; 17 19; 17 22; 17 28; 17 36; 17 42; 17 43; 17 44; 18 3; 18 9; 18 40; 19 10; 19 12; 19 14; 19 31; 20 5; 20 7; 20 12; 20 16; 20 34; 21 3; 21 8; 21 13; 21 14; 21 27; 21 29; 21 33; 21 34; 21 50; 22 14; 22 16; 22 28; 22 43; 22 48; 23 2; 23 3; 23 13; 23 19; 23 21; 23 22; 23 28; 23 39; 23 47; 23 48; 23 50; 24 25; 24 35; 25 36; 25 44; 26 10; 26 19; 26 22; 26 25; 26 33; 26 42; 26 50; 27 13; 27 20; 27 30; 27 34; 27 44; 28 2; 28 23; 28 34; 28 44; 29 6; 29 12; 29 14; 29 21; 29 33; 29 35; 30 7; 30 44; 31 13; 32 11; 32 16; 32 38; 32 41; 33 10; 33 31; 33 48; 34 12; 34 13; 34 20; 34 26; 34 27; 34 38; 35 26; 35 30; 35 33; 35 38; 35 40; 36 11; 36 16; 36 27; 36 28; 36 29; 37 2; 37 6; 37 15; 37 38; 37 40; 38 9; 38 14; 38 19; 38 36; 38 37; 39 2; 39 13; 39 21; 39 26; 39 27; 39 46; 40 5; 40 6; 40 36; 41 10; 41 13; 41 20; 41 21; 41 43; 42 3; 42 19; 42 46; 42 47; 43 2; 43 9; 43 20; 43 27; 43 34; 43 50; 44 9; 44 31; 44 42; 44 43; 45 10; 45 15; 45 41; 45 43; 46 28; 46 34; 46 35; 46 36; 46 41; 47 6; 47 13; 47 16; 47 22; 48 9; 48 17; 48 26; 49 2; 49 17; 49 22]
global d_x = [2.0, 7.0, 1.0, 10.0, 8.0, 8.0, 2.0, 6.0, 10.0, 5.0, 4.0, 7.0, 1.0, 5.0, 9.0, 5.0, 9.0, 9.0, 2.0, 1.0, 7.0, 9.0, 3.0, 6.0, 10.0, 10.0, 2.0, 3.0, 1.0, 4.0, 1.0, 2.0, 3.0, 5.0, 2.0, 7.0, 8.0, 5.0, 8.0, 7.0, 3.0, 1.0, 5.0, 7.0, 8.0, 5.0, 7.0, 10.0, 3.0, 6.0, 10.0, 5.0, 9.0, 8.0, 3.0, 1.0, 7.0, 10.0, 7.0, 4.0, 8.0, 8.0, 3.0, 5.0, 6.0, 8.0, 5.0, 7.0, 3.0, 5.0, 1.0, 10.0, 5.0, 3.0, 1.0, 5.0, 6.0, 1.0, 5.0, 3.0, 10.0, 2.0, 10.0, 3.0, 6.0, 2.0, 2.0, 8.0, 8.0, 1.0, 2.0, 10.0, 2.0, 3.0, 6.0, 7.0, 2.0, 4.0, 1.0, 1.0, 4.0, 3.0, 8.0, 8.0, 9.0, 7.0, 9.0, 3.0, 10.0, 2.0, 1.0, 2.0, 8.0, 7.0, 8.0, 4.0, 3.0, 9.0, 8.0, 3.0, 10.0, 7.0, 9.0, 9.0, 7.0, 3.0, 7.0, 3.0, 2.0, 5.0, 4.0, 9.0, 1.0, 1.0, 7.0, 10.0, 7.0, 3.0, 10.0, 5.0, 6.0, 10.0, 1.0, 9.0, 5.0, 1.0, 10.0, 7.0, 5.0, 1.0, 8.0, 6.0, 9.0, 8.0, 10.0, 6.0, 5.0, 5.0, 5.0, 8.0, 1.0, 9.0, 6.0, 8.0, 4.0, 3.0, 9.0, 2.0, 6.0, 10.0, 3.0, 8.0, 8.0, 9.0, 4.0, 9.0, 9.0, 8.0, 8.0, 7.0, 1.0, 9.0, 10.0, 8.0, 7.0, 9.0, 4.0, 10.0, 5.0, 7.0, 1.0, 8.0, 5.0, 9.0, 6.0, 1.0, 1.0, 5.0, 1.0, 8.0, 5.0, 2.0, 4.0, 7.0, 9.0, 4.0, 9.0, 4.0, 6.0, 10.0, 3.0, 2.0, 1.0, 7.0, 1.0, 2.0, 7.0, 1.0, 6.0, 5.0, 10.0, 8.0, 6.0, 1.0, 10.0, 4.0, 9.0, 6.0, 6.0, 5.0, 7.0, 1.0, 2.0, 5.0, 2.0, 8.0, 7.0, 1.0, 6.0, 10.0, 4.0, 5.0, 2.0, 2.0, 3.0, 8.0, 6.0, 3.0, 10.0, 7.0, 4.0, 4.0, 2.0, 8.0, 10.0]
global b_x = 5
global d_y = [10.0, 6.0, 2.0, 3.0, 10.0, 4.0, 8.0, 2.0, 9.0, 10.0, 3.0, 9.0, 10.0, 6.0, 1.0, 10.0, 2.0, 6.0, 5.0, 7.0, 10.0, 4.0, 4.0, 1.0, 3.0, 5.0, 9.0, 6.0, 4.0, 9.0, 9.0, 8.0, 5.0, 9.0, 10.0, 9.0, 1.0, 3.0, 8.0, 9.0, 1.0, 4.0, 3.0, 3.0, 10.0, 4.0, 8.0, 4.0, 4.0, 6.0, 9.0, 9.0, 9.0, 1.0, 4.0, 8.0, 3.0, 5.0, 9.0, 8.0, 5.0, 1.0, 2.0, 1.0, 9.0, 3.0, 7.0, 1.0, 7.0, 10.0, 10.0, 4.0, 1.0, 3.0, 10.0, 6.0, 2.0, 9.0, 1.0, 3.0, 7.0, 4.0, 6.0, 8.0, 6.0, 3.0, 5.0, 7.0, 5.0, 6.0, 10.0, 5.0, 3.0, 7.0, 9.0, 8.0, 4.0, 5.0, 5.0, 5.0, 5.0, 2.0, 10.0, 5.0, 10.0, 8.0, 4.0, 10.0, 9.0, 5.0, 6.0, 3.0, 2.0, 4.0, 3.0, 4.0, 4.0, 1.0, 3.0, 9.0, 9.0, 2.0, 7.0, 3.0, 5.0, 2.0, 4.0, 4.0, 2.0, 5.0, 1.0, 8.0, 9.0, 2.0, 8.0, 5.0, 3.0, 7.0, 9.0, 5.0, 8.0, 10.0, 9.0, 7.0, 7.0, 1.0, 1.0, 9.0, 9.0, 9.0, 4.0, 7.0, 7.0, 4.0, 6.0, 4.0, 9.0, 8.0, 10.0, 3.0, 6.0, 4.0, 5.0, 8.0, 2.0, 7.0, 2.0, 6.0, 7.0, 2.0, 9.0, 10.0, 4.0, 1.0, 4.0, 2.0, 10.0, 1.0, 3.0, 1.0, 5.0, 6.0, 2.0, 3.0, 9.0, 8.0, 8.0, 4.0, 4.0, 6.0, 3.0, 1.0, 4.0, 1.0, 9.0, 3.0, 6.0, 7.0, 2.0, 10.0, 4.0, 4.0, 3.0, 6.0, 4.0, 1.0, 9.0, 4.0, 8.0, 8.0, 1.0, 2.0, 8.0, 8.0, 4.0, 6.0, 6.0, 10.0, 4.0, 9.0, 10.0, 1.0, 7.0, 5.0, 7.0, 3.0, 2.0, 10.0, 6.0, 8.0, 7.0, 5.0, 3.0, 2.0, 10.0, 7.0, 7.0, 4.0, 2.0, 4.0, 6.0, 5.0, 8.0, 9.0, 7.0, 10.0, 8.0, 4.0, 9.0, 5.0, 10.0, 6.0, 1.0, 9.0, 7.0]
global b_y = 10
global p = [0.887, 0.055, 0.523, 0.135, 0.1, 0.016, 0.908, 0.197, 0.035, 0.786, 0.184, 0.703, 0.134, 0.19, 0.171, 0.468, 0.085, 0.803, 0.727, 0.308, 0.112, 0.726, 0.716, 0.695, 0.003, 0.747, 0.389, 0.66, 0.92, 0.735, 0.944, 0.789, 0.891, 0.522, 0.455, 0.09, 0.055, 0.453, 0.509, 0.856, 0.168, 0.354, 0.063, 0.089, 0.937, 0.159, 0.248, 0.456, 0.563, 0.76, 0.065, 0.375, 0.454, 0.996, 0.232, 0.865, 0.327, 0.299, 0.292, 0.796, 0.999, 0.312, 0.808, 0.693, 0.091, 0.298, 0.954, 0.961, 0.124, 0.512, 0.354, 0.428, 0.227, 0.892, 0.236, 0.59, 0.157, 0.257, 0.966, 0.667, 0.081, 0.937, 0.575, 0.496, 0.321, 0.458, 0.21, 0.056, 0.445, 0.852, 0.602, 0.102, 0.999, 0.409, 0.13, 0.741, 0.822, 0.619, 0.099, 0.528, 0.575, 0.466, 0.46, 0.379, 0.87, 0.127, 0.312, 0.737, 0.541, 0.985, 0.09, 0.798, 0.774, 0.251, 0.447, 0.322, 0.398, 0.694, 0.144, 0.681, 0.597, 0.227, 0.254, 0.184, 0.646, 0.333, 0.189, 0.413, 0.428, 0.531, 0.781, 0.718, 0.49, 0.334, 0.12, 0.424, 0.045, 0.978, 0.914, 0.883, 0.266, 0.24, 0.589, 0.529, 0.674, 0.121, 0.503, 0.64, 0.079, 0.242, 0.071, 0.094, 0.938, 0.056, 0.935, 0.558, 0.013, 0.52, 0.887, 0.95, 0.261, 0.469, 0.098, 0.878, 0.629, 0.14, 0.042, 0.895, 0.339, 0.658, 0.69, 0.027, 0.215, 0.026, 0.686, 0.148, 0.47, 0.768, 0.545, 0.166, 0.756, 0.739, 0.994, 0.861, 0.89, 0.583, 0.794, 0.102, 0.042, 0.323, 0.561, 0.805, 0.856, 0.244, 0.203, 0.695, 0.662, 0.039, 0.269, 0.908, 0.2, 0.945, 0.298, 0.991, 0.655, 0.084, 0.563, 0.269, 0.159, 0.138, 0.623, 0.655, 0.702, 0.089, 0.356, 0.916, 0.006, 0.836, 0.204, 0.369, 0.527, 0.762, 0.034, 0.094, 0.329, 0.994, 0.196, 0.899, 0.47, 0.061, 0.865, 0.747, 0.354, 0.774, 0.63, 0.506, 0.501, 0.376, 0.708, 0.265, 0.79, 0.006, 0.904, 0.543, 0.063, 0.801, 0.796, 0.978, 0.499, 0.968, 0.277, 0.6, 0.977, 0.179, 0.239]
global q = [0.905, 0.204, 0.585, 0.679, 0.372, 0.075, 0.908, 0.915, 0.422, 0.839, 0.436, 0.837, 0.369, 0.664, 0.251, 0.914, 0.998, 0.864, 0.767, 0.574, 0.647, 0.897, 0.97, 0.85, 0.276, 0.755, 0.763, 0.972, 0.924, 0.818, 0.952, 0.953, 0.939, 0.544, 0.768, 0.472, 0.346, 0.554, 0.518, 0.879, 0.255, 0.454, 0.342, 0.778, 0.981, 0.816, 0.254, 0.999, 0.607, 0.984, 0.105, 0.81, 0.614, 0.997, 0.283, 0.926, 0.708, 0.334, 0.945, 0.903, 0.999, 0.577, 0.992, 0.755, 0.117, 0.411, 0.975, 0.99, 0.933, 0.628, 0.723, 0.652, 0.806, 0.974, 0.795, 0.976, 0.251, 0.578, 0.984, 0.758, 0.568, 0.991, 0.64, 0.809, 0.551, 0.871, 0.788, 0.311, 0.563, 0.883, 0.875, 0.139, 0.999, 0.773, 0.767, 0.953, 0.952, 0.663, 0.934, 0.727, 0.765, 0.963, 0.654, 0.948, 0.958, 0.794, 0.356, 0.972, 0.705, 0.997, 0.939, 0.918, 0.828, 0.545, 0.448, 0.403, 0.887, 0.818, 0.81, 0.694, 0.826, 0.963, 0.63, 0.761, 0.687, 0.931, 0.752, 0.522, 0.548, 0.867, 0.856, 0.978, 0.897, 0.703, 0.501, 0.47, 0.614, 0.991, 0.939, 0.972, 0.486, 0.721, 0.967, 0.723, 0.968, 0.272, 0.839, 0.746, 0.557, 0.341, 0.161, 0.499, 0.951, 0.525, 0.997, 0.67, 0.489, 0.951, 0.98, 0.977, 0.536, 0.686, 0.388, 0.98, 0.753, 0.332, 0.58, 0.989, 0.858, 0.697, 0.886, 0.828, 0.97, 0.114, 0.722, 0.157, 0.664, 0.984, 0.615, 0.404, 0.76, 0.826, 0.994, 0.911, 0.993, 0.883, 0.806, 0.696, 0.061, 0.513, 0.879, 0.939, 0.91, 0.463, 0.714, 0.93, 0.744, 0.249, 0.539, 0.984, 0.307, 0.967, 0.318, 0.998, 0.849, 0.784, 0.812, 0.288, 0.838, 0.57, 0.883, 0.811, 0.83, 0.238, 0.432, 0.972, 0.49, 0.869, 0.391, 0.967, 0.693, 0.931, 0.685, 0.645, 0.702, 0.994, 0.888, 0.934, 0.789, 0.377, 0.943, 0.893, 0.544, 0.882, 0.861, 0.895, 0.648, 0.79, 0.798, 0.512, 0.987, 0.548, 0.951, 0.872, 0.359, 0.809, 0.837, 0.999, 0.822, 0.987, 0.387, 0.603, 0.99, 0.71, 0.421]
global origin = 1
global destination = 50