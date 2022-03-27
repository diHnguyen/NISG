global arcs = [1 23; 1 37; 1 50; 2 9; 2 21; 2 46; 3 2; 3 4; 3 8; 3 28; 3 29; 4 7; 4 15; 4 17; 4 28; 4 38; 4 40; 5 10; 5 12; 5 17; 5 19; 5 20; 5 23; 5 39; 5 44; 6 10; 6 13; 6 38; 6 41; 6 43; 6 46; 6 49; 7 2; 7 13; 7 23; 7 24; 7 39; 7 46; 8 16; 9 3; 9 26; 9 35; 9 36; 9 39; 9 41; 9 43; 10 20; 10 23; 10 25; 10 27; 10 31; 10 45; 10 46; 11 5; 11 8; 11 9; 11 47; 11 49; 12 14; 12 19; 12 22; 12 27; 12 40; 13 10; 13 38; 13 40; 13 46; 13 47; 14 4; 14 6; 15 20; 15 22; 15 41; 16 11; 16 15; 16 29; 16 34; 16 37; 16 39; 16 40; 16 42; 17 3; 17 6; 17 14; 17 18; 17 21; 17 25; 17 26; 17 45; 17 47; 18 9; 18 10; 18 45; 18 47; 19 4; 19 12; 19 20; 19 44; 19 45; 19 49; 19 50; 20 18; 20 32; 21 8; 21 11; 21 16; 21 29; 22 4; 22 31; 22 36; 22 41; 22 48; 23 15; 23 28; 23 30; 23 35; 23 39; 23 47; 23 48; 24 27; 25 4; 25 23; 25 30; 25 41; 25 44; 25 48; 26 2; 26 3; 26 6; 26 20; 26 22; 26 39; 26 40; 27 3; 27 36; 27 38; 27 49; 28 23; 28 29; 29 25; 29 47; 30 3; 30 23; 30 24; 30 27; 30 31; 30 37; 30 42; 30 44; 31 13; 31 37; 32 6; 32 9; 32 22; 32 26; 33 15; 33 43; 34 7; 34 13; 34 17; 34 18; 34 31; 34 45; 35 5; 35 7; 35 11; 35 39; 35 45; 36 6; 36 22; 36 27; 36 28; 36 48; 37 14; 37 39; 37 47; 37 49; 38 20; 38 22; 38 25; 38 36; 39 2; 39 4; 39 13; 39 33; 39 37; 39 48; 39 49; 40 18; 40 22; 40 32; 40 34; 40 45; 41 5; 41 7; 41 25; 41 37; 41 50; 42 18; 42 29; 42 30; 42 38; 43 10; 43 45; 44 11; 44 26; 44 30; 45 13; 45 15; 45 17; 45 22; 45 29; 45 44; 46 7; 46 13; 46 28; 46 48; 47 8; 47 25; 47 36; 47 50; 48 41; 48 47; 49 9; 49 16; 49 23; 49 37; 49 47; 49 48]
global d_x = [7.0, 7.0, 8.0, 5.0, 6.0, 5.0, 8.0, 3.0, 8.0, 7.0, 10.0, 5.0, 10.0, 2.0, 4.0, 8.0, 1.0, 5.0, 3.0, 5.0, 6.0, 9.0, 8.0, 10.0, 2.0, 3.0, 1.0, 8.0, 9.0, 3.0, 10.0, 5.0, 7.0, 1.0, 3.0, 3.0, 9.0, 4.0, 4.0, 7.0, 1.0, 5.0, 2.0, 4.0, 7.0, 6.0, 2.0, 9.0, 7.0, 4.0, 2.0, 10.0, 9.0, 4.0, 5.0, 1.0, 4.0, 6.0, 1.0, 9.0, 3.0, 7.0, 4.0, 7.0, 3.0, 9.0, 6.0, 3.0, 2.0, 7.0, 3.0, 5.0, 4.0, 7.0, 5.0, 4.0, 1.0, 4.0, 1.0, 1.0, 8.0, 3.0, 1.0, 3.0, 1.0, 1.0, 6.0, 7.0, 1.0, 4.0, 1.0, 5.0, 4.0, 4.0, 8.0, 4.0, 8.0, 4.0, 4.0, 7.0, 2.0, 1.0, 7.0, 10.0, 5.0, 2.0, 5.0, 5.0, 8.0, 4.0, 8.0, 6.0, 5.0, 6.0, 5.0, 1.0, 4.0, 5.0, 10.0, 7.0, 6.0, 10.0, 4.0, 6.0, 6.0, 5.0, 9.0, 4.0, 4.0, 3.0, 8.0, 1.0, 8.0, 7.0, 4.0, 2.0, 7.0, 1.0, 1.0, 1.0, 3.0, 7.0, 3.0, 6.0, 9.0, 7.0, 9.0, 4.0, 4.0, 3.0, 4.0, 7.0, 9.0, 2.0, 2.0, 7.0, 4.0, 7.0, 1.0, 5.0, 4.0, 7.0, 2.0, 1.0, 9.0, 2.0, 2.0, 3.0, 9.0, 1.0, 8.0, 2.0, 5.0, 9.0, 9.0, 5.0, 10.0, 5.0, 5.0, 7.0, 5.0, 4.0, 1.0, 9.0, 6.0, 10.0, 1.0, 8.0, 5.0, 8.0, 3.0, 8.0, 10.0, 9.0, 3.0, 8.0, 7.0, 4.0, 10.0, 3.0, 8.0, 10.0, 7.0, 6.0, 2.0, 7.0, 7.0, 2.0, 10.0, 2.0, 6.0, 2.0, 7.0, 3.0, 4.0, 5.0, 8.0, 5.0, 2.0, 8.0, 4.0, 5.0, 1.0, 4.0, 4.0, 8.0, 2.0, 6.0, 9.0]
global b_x = 5
global d_y = [2.0, 6.0, 1.0, 7.0, 5.0, 2.0, 8.0, 6.0, 7.0, 4.0, 5.0, 6.0, 9.0, 6.0, 6.0, 3.0, 3.0, 4.0, 7.0, 1.0, 8.0, 9.0, 2.0, 5.0, 8.0, 2.0, 2.0, 1.0, 8.0, 9.0, 8.0, 6.0, 5.0, 5.0, 3.0, 5.0, 4.0, 6.0, 3.0, 5.0, 10.0, 7.0, 3.0, 6.0, 1.0, 10.0, 7.0, 4.0, 2.0, 6.0, 7.0, 2.0, 9.0, 3.0, 2.0, 7.0, 6.0, 6.0, 3.0, 1.0, 2.0, 10.0, 7.0, 1.0, 10.0, 6.0, 10.0, 1.0, 7.0, 10.0, 5.0, 10.0, 5.0, 7.0, 10.0, 1.0, 10.0, 6.0, 8.0, 1.0, 6.0, 5.0, 6.0, 7.0, 7.0, 8.0, 1.0, 5.0, 6.0, 6.0, 5.0, 9.0, 2.0, 6.0, 6.0, 3.0, 3.0, 1.0, 10.0, 9.0, 10.0, 10.0, 4.0, 4.0, 4.0, 1.0, 10.0, 7.0, 7.0, 4.0, 9.0, 4.0, 3.0, 3.0, 1.0, 3.0, 1.0, 7.0, 7.0, 9.0, 10.0, 3.0, 2.0, 7.0, 7.0, 3.0, 7.0, 8.0, 10.0, 5.0, 3.0, 3.0, 1.0, 3.0, 1.0, 5.0, 3.0, 8.0, 7.0, 2.0, 6.0, 6.0, 7.0, 10.0, 10.0, 5.0, 2.0, 1.0, 7.0, 7.0, 6.0, 4.0, 8.0, 2.0, 9.0, 7.0, 8.0, 4.0, 2.0, 10.0, 9.0, 3.0, 8.0, 3.0, 2.0, 3.0, 3.0, 9.0, 8.0, 10.0, 5.0, 7.0, 9.0, 8.0, 5.0, 4.0, 6.0, 10.0, 9.0, 3.0, 2.0, 7.0, 3.0, 5.0, 5.0, 2.0, 8.0, 5.0, 7.0, 6.0, 3.0, 3.0, 5.0, 5.0, 10.0, 8.0, 6.0, 5.0, 4.0, 3.0, 10.0, 5.0, 4.0, 5.0, 6.0, 7.0, 9.0, 9.0, 7.0, 8.0, 2.0, 2.0, 6.0, 5.0, 6.0, 9.0, 7.0, 10.0, 6.0, 5.0, 5.0, 4.0, 5.0, 5.0, 10.0, 1.0, 2.0, 8.0, 1.0]
global b_y = 10
global p = [0.329, 0.141, 0.627, 0.442, 0.759, 0.355, 0.303, 0.718, 0.322, 0.136, 0.871, 0.316, 0.837, 0.216, 0.548, 0.186, 0.464, 0.398, 0.576, 0.35, 0.467, 0.543, 0.235, 0.376, 0.826, 0.365, 0.721, 0.382, 0.073, 0.487, 0.251, 0.628, 0.317, 0.147, 0.037, 0.78, 0.363, 0.304, 0.49, 0.902, 0.341, 0.411, 0.629, 0.479, 0.088, 0.538, 0.067, 0.632, 0.035, 0.548, 0.764, 0.638, 0.207, 0.38, 0.083, 0.905, 0.693, 0.77, 0.032, 0.381, 0.624, 0.878, 0.016, 0.061, 0.287, 0.741, 0.615, 0.128, 0.903, 0.829, 0.574, 0.46, 0.594, 0.293, 0.624, 0.916, 0.01, 0.431, 0.111, 0.467, 0.756, 0.167, 0.136, 0.025, 0.715, 0.303, 0.604, 0.742, 0.66, 0.127, 0.575, 0.988, 0.867, 0.672, 0.985, 0.629, 0.044, 0.521, 0.285, 0.439, 0.209, 0.564, 0.776, 0.264, 0.176, 0.944, 0.828, 0.868, 0.233, 0.837, 0.108, 0.872, 0.091, 0.064, 0.461, 0.459, 0.127, 0.05, 0.614, 0.927, 0.652, 0.052, 0.579, 0.983, 0.666, 0.272, 0.212, 0.03, 0.042, 0.726, 0.479, 0.009, 0.973, 0.34, 0.535, 0.989, 0.451, 0.351, 0.778, 0.757, 0.424, 0.447, 0.838, 0.701, 0.456, 0.591, 0.877, 0.1, 0.046, 0.664, 0.469, 0.275, 0.799, 0.804, 0.693, 0.223, 0.111, 0.691, 0.996, 0.211, 0.157, 0.771, 0.855, 0.471, 0.587, 0.209, 0.191, 0.75, 0.114, 0.141, 0.016, 0.181, 0.438, 0.062, 0.3, 0.365, 0.529, 0.025, 0.548, 0.449, 0.345, 0.106, 0.168, 0.057, 0.16, 0.932, 0.777, 0.149, 0.24, 0.707, 0.662, 0.392, 0.066, 0.552, 0.373, 0.788, 0.153, 0.128, 0.885, 0.514, 0.716, 0.098, 0.035, 0.771, 0.55, 0.416, 0.887, 0.119, 0.178, 0.234, 0.825, 0.108, 0.317, 0.62, 0.498, 0.534, 0.127, 0.714, 0.028, 0.008, 0.74, 0.593, 0.671, 0.303, 0.782, 0.858, 0.337, 0.472, 0.886]
global q = [0.964, 0.354, 0.634, 0.484, 0.828, 0.838, 0.831, 0.869, 0.919, 0.662, 0.919, 0.614, 0.926, 0.56, 0.731, 0.984, 0.771, 0.512, 0.809, 0.785, 0.499, 0.873, 0.823, 0.985, 0.862, 0.667, 0.888, 0.96, 0.921, 0.797, 0.827, 0.777, 0.797, 0.863, 0.962, 0.937, 0.76, 0.849, 0.934, 0.921, 0.379, 0.514, 0.917, 0.517, 0.949, 0.817, 0.224, 0.849, 0.185, 0.748, 0.877, 0.651, 0.783, 0.523, 0.231, 0.94, 0.846, 0.875, 0.658, 0.59, 0.971, 0.915, 0.747, 0.471, 0.513, 0.941, 0.803, 0.45, 0.976, 0.907, 0.983, 0.747, 0.911, 0.507, 0.986, 0.93, 0.349, 0.745, 0.48, 0.827, 0.879, 0.789, 0.695, 0.659, 0.77, 0.558, 0.923, 0.952, 0.674, 0.238, 0.593, 0.999, 0.884, 0.913, 0.998, 0.681, 0.732, 0.772, 0.823, 0.646, 0.775, 0.935, 0.814, 0.929, 0.698, 0.965, 0.979, 0.898, 0.956, 0.883, 0.832, 0.927, 0.394, 0.351, 0.788, 0.985, 0.154, 0.728, 0.63, 0.955, 0.969, 0.429, 0.595, 0.993, 0.698, 0.591, 0.96, 0.717, 0.189, 0.8, 0.94, 0.741, 0.996, 0.812, 0.841, 0.994, 0.756, 0.958, 0.924, 0.94, 0.438, 0.977, 0.966, 0.926, 0.934, 0.802, 0.93, 0.301, 0.684, 0.777, 0.515, 0.452, 0.903, 0.886, 0.738, 0.685, 0.467, 0.768, 0.998, 0.448, 0.246, 0.796, 0.922, 0.578, 0.608, 0.96, 0.574, 0.762, 0.833, 0.487, 0.425, 0.869, 0.692, 0.958, 0.568, 0.807, 0.727, 0.492, 0.854, 0.92, 0.921, 0.627, 0.515, 0.841, 0.638, 0.941, 0.864, 0.483, 0.418, 0.735, 0.868, 0.504, 0.644, 0.885, 0.455, 0.904, 0.943, 0.705, 0.961, 0.732, 0.846, 0.806, 0.43, 0.782, 0.597, 0.781, 0.935, 0.894, 0.458, 0.361, 0.881, 0.621, 0.772, 0.628, 0.948, 0.818, 0.27, 0.869, 0.524, 0.975, 0.778, 0.97, 0.891, 0.491, 0.905, 0.977, 0.351, 0.763, 0.93]
global origin = 1
global destination = 50