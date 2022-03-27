global arcs = [1 12; 1 13; 1 20; 1 21; 1 31; 1 47; 2 11; 2 21; 2 25; 2 30; 2 39; 3 2; 3 6; 3 8; 3 16; 3 24; 3 31; 3 42; 4 15; 4 26; 4 40; 4 47; 5 7; 5 39; 5 40; 5 46; 6 13; 6 14; 6 15; 6 19; 6 24; 6 29; 6 45; 7 14; 7 18; 7 25; 7 33; 7 36; 7 45; 8 6; 8 28; 8 30; 8 42; 9 20; 9 22; 9 47; 10 4; 10 22; 11 17; 11 18; 11 46; 11 50; 12 4; 12 14; 12 24; 12 47; 13 17; 13 27; 13 44; 14 20; 14 21; 14 46; 14 50; 15 17; 15 22; 15 29; 15 33; 15 37; 15 38; 15 41; 16 17; 16 29; 16 33; 16 43; 16 44; 16 49; 17 20; 17 26; 17 44; 17 45; 17 46; 18 14; 18 24; 18 41; 19 3; 19 9; 19 24; 19 25; 19 28; 19 38; 19 43; 20 7; 20 13; 20 14; 20 23; 20 50; 21 9; 21 16; 21 19; 21 20; 21 32; 21 35; 21 37; 22 6; 22 14; 22 23; 22 28; 22 30; 22 31; 22 39; 23 5; 23 12; 23 15; 23 22; 23 35; 23 49; 23 50; 24 5; 24 9; 24 16; 24 45; 24 48; 25 8; 25 13; 25 35; 25 38; 25 40; 26 4; 26 8; 26 23; 26 37; 26 47; 27 18; 27 28; 27 36; 27 38; 27 50; 28 2; 28 7; 28 21; 28 29; 29 26; 29 41; 29 44; 30 7; 30 12; 30 28; 30 40; 30 47; 31 5; 31 9; 31 36; 31 46; 32 12; 32 22; 32 25; 32 30; 33 2; 33 5; 33 8; 33 31; 33 45; 34 4; 34 11; 34 15; 34 16; 34 28; 34 44; 34 48; 35 5; 35 16; 35 30; 35 44; 35 49; 36 10; 36 11; 36 22; 36 42; 36 45; 37 28; 37 33; 38 13; 38 14; 38 15; 38 18; 38 27; 38 34; 39 5; 39 9; 39 14; 39 36; 39 40; 39 42; 39 45; 39 48; 40 6; 40 11; 40 33; 41 9; 41 25; 41 26; 41 27; 41 39; 41 43; 41 44; 41 45; 42 37; 42 43; 42 44; 43 5; 43 12; 43 21; 43 27; 43 28; 43 37; 44 5; 44 12; 44 14; 44 27; 44 34; 44 35; 44 36; 44 49; 45 7; 45 9; 45 12; 45 28; 45 30; 45 40; 45 41; 45 42; 46 7; 46 8; 46 11; 46 13; 46 26; 46 27; 46 50; 47 16; 47 21; 47 23; 47 34; 47 37; 47 45; 48 6; 48 9; 48 13; 48 38; 48 45; 48 47; 49 16; 49 17; 49 40; 49 47; 49 48]
global d_x = [4.0, 7.0, 2.0, 7.0, 9.0, 7.0, 4.0, 10.0, 6.0, 9.0, 6.0, 3.0, 6.0, 6.0, 4.0, 10.0, 7.0, 1.0, 10.0, 7.0, 4.0, 10.0, 10.0, 2.0, 4.0, 5.0, 5.0, 3.0, 5.0, 7.0, 3.0, 6.0, 3.0, 1.0, 4.0, 10.0, 3.0, 2.0, 7.0, 4.0, 10.0, 4.0, 2.0, 1.0, 1.0, 10.0, 5.0, 9.0, 1.0, 2.0, 5.0, 10.0, 9.0, 9.0, 6.0, 8.0, 2.0, 3.0, 2.0, 7.0, 1.0, 1.0, 2.0, 8.0, 9.0, 7.0, 9.0, 9.0, 10.0, 9.0, 5.0, 2.0, 1.0, 9.0, 6.0, 7.0, 1.0, 7.0, 9.0, 2.0, 3.0, 2.0, 4.0, 6.0, 7.0, 6.0, 9.0, 4.0, 9.0, 7.0, 1.0, 3.0, 7.0, 2.0, 7.0, 10.0, 3.0, 6.0, 10.0, 1.0, 8.0, 4.0, 7.0, 7.0, 8.0, 1.0, 7.0, 8.0, 4.0, 9.0, 1.0, 2.0, 2.0, 7.0, 9.0, 3.0, 6.0, 1.0, 5.0, 10.0, 6.0, 3.0, 4.0, 6.0, 1.0, 9.0, 2.0, 1.0, 9.0, 4.0, 10.0, 3.0, 9.0, 9.0, 7.0, 10.0, 9.0, 6.0, 6.0, 7.0, 2.0, 7.0, 6.0, 9.0, 6.0, 2.0, 9.0, 8.0, 10.0, 7.0, 1.0, 10.0, 3.0, 2.0, 3.0, 9.0, 3.0, 6.0, 5.0, 5.0, 7.0, 6.0, 2.0, 4.0, 5.0, 3.0, 3.0, 8.0, 3.0, 2.0, 5.0, 2.0, 5.0, 10.0, 10.0, 1.0, 9.0, 5.0, 5.0, 6.0, 4.0, 6.0, 2.0, 3.0, 3.0, 7.0, 6.0, 2.0, 3.0, 5.0, 2.0, 9.0, 2.0, 1.0, 6.0, 5.0, 1.0, 3.0, 1.0, 10.0, 8.0, 5.0, 7.0, 6.0, 3.0, 8.0, 7.0, 7.0, 2.0, 10.0, 3.0, 2.0, 4.0, 8.0, 6.0, 6.0, 5.0, 1.0, 2.0, 6.0, 4.0, 3.0, 2.0, 9.0, 4.0, 2.0, 8.0, 1.0, 2.0, 10.0, 9.0, 8.0, 3.0, 8.0, 10.0, 8.0, 7.0, 2.0, 5.0, 8.0, 2.0, 6.0, 6.0, 1.0, 8.0, 10.0, 3.0, 2.0, 2.0, 8.0, 9.0, 7.0, 6.0, 10.0, 4.0]
global b_x = 5
global d_y = [5.0, 5.0, 2.0, 4.0, 9.0, 2.0, 6.0, 6.0, 5.0, 9.0, 9.0, 7.0, 5.0, 3.0, 2.0, 1.0, 8.0, 6.0, 8.0, 9.0, 4.0, 2.0, 5.0, 5.0, 4.0, 2.0, 3.0, 8.0, 4.0, 10.0, 2.0, 3.0, 8.0, 8.0, 10.0, 6.0, 1.0, 5.0, 8.0, 2.0, 6.0, 5.0, 7.0, 5.0, 1.0, 8.0, 2.0, 8.0, 6.0, 6.0, 9.0, 6.0, 4.0, 2.0, 3.0, 4.0, 10.0, 5.0, 3.0, 6.0, 2.0, 3.0, 1.0, 3.0, 9.0, 8.0, 9.0, 6.0, 2.0, 5.0, 9.0, 5.0, 3.0, 1.0, 2.0, 1.0, 8.0, 10.0, 2.0, 9.0, 4.0, 5.0, 8.0, 2.0, 10.0, 3.0, 9.0, 3.0, 8.0, 8.0, 7.0, 7.0, 5.0, 7.0, 7.0, 1.0, 7.0, 3.0, 4.0, 5.0, 7.0, 2.0, 6.0, 2.0, 2.0, 5.0, 10.0, 1.0, 5.0, 6.0, 8.0, 9.0, 9.0, 9.0, 8.0, 10.0, 2.0, 7.0, 9.0, 8.0, 10.0, 2.0, 7.0, 4.0, 6.0, 8.0, 5.0, 6.0, 1.0, 9.0, 4.0, 7.0, 4.0, 4.0, 10.0, 2.0, 9.0, 9.0, 2.0, 9.0, 8.0, 8.0, 7.0, 3.0, 2.0, 4.0, 3.0, 3.0, 5.0, 9.0, 8.0, 2.0, 4.0, 10.0, 9.0, 2.0, 7.0, 3.0, 4.0, 6.0, 1.0, 1.0, 1.0, 8.0, 5.0, 1.0, 1.0, 2.0, 8.0, 2.0, 3.0, 8.0, 6.0, 4.0, 1.0, 9.0, 8.0, 2.0, 3.0, 3.0, 10.0, 9.0, 5.0, 7.0, 9.0, 7.0, 5.0, 4.0, 8.0, 3.0, 3.0, 4.0, 3.0, 10.0, 4.0, 7.0, 5.0, 9.0, 9.0, 8.0, 7.0, 8.0, 1.0, 7.0, 2.0, 4.0, 8.0, 7.0, 3.0, 2.0, 8.0, 5.0, 3.0, 1.0, 5.0, 10.0, 6.0, 9.0, 8.0, 7.0, 4.0, 7.0, 9.0, 7.0, 10.0, 5.0, 3.0, 10.0, 1.0, 7.0, 7.0, 5.0, 10.0, 1.0, 2.0, 8.0, 8.0, 9.0, 1.0, 1.0, 4.0, 2.0, 1.0, 10.0, 4.0, 2.0, 10.0, 1.0, 7.0, 4.0, 4.0, 4.0, 2.0, 9.0, 8.0]
global b_y = 10
global p = [0.164, 0.923, 0.505, 0.812, 0.803, 0.067, 0.768, 0.746, 0.226, 0.8, 0.151, 0.402, 0.535, 0.021, 0.985, 0.857, 0.046, 0.805, 0.347, 0.693, 0.069, 0.893, 0.461, 0.67, 0.489, 0.927, 0.213, 0.201, 0.528, 0.381, 0.462, 0.763, 0.829, 0.487, 0.283, 0.11, 0.143, 0.108, 0.916, 0.818, 0.628, 0.188, 0.123, 0.488, 0.882, 0.609, 0.932, 0.968, 0.075, 0.066, 0.148, 0.242, 0.434, 0.145, 0.29, 0.495, 0.439, 0.426, 0.213, 0.411, 0.568, 0.752, 0.601, 0.873, 0.917, 0.78, 0.79, 0.546, 0.598, 0.016, 0.463, 0.354, 0.095, 0.138, 0.94, 0.954, 0.264, 0.714, 0.24, 0.716, 0.078, 0.964, 0.925, 0.808, 0.828, 0.216, 0.718, 0.998, 0.58, 0.039, 0.086, 0.621, 0.4, 0.276, 0.104, 0.914, 0.179, 0.462, 0.444, 0.545, 0.167, 0.515, 0.982, 0.878, 0.725, 0.031, 0.813, 0.602, 0.141, 0.423, 0.074, 0.144, 0.774, 0.474, 0.22, 0.729, 0.624, 0.947, 0.105, 0.04, 0.36, 0.112, 0.291, 0.971, 0.033, 0.357, 0.509, 0.188, 0.869, 0.13, 0.689, 0.044, 0.073, 0.627, 0.771, 0.911, 0.957, 0.276, 0.586, 0.595, 0.516, 0.492, 0.392, 0.144, 0.613, 0.799, 0.982, 0.84, 0.952, 0.299, 0.76, 0.396, 0.358, 0.23, 0.433, 0.931, 0.762, 0.714, 0.962, 0.768, 0.817, 0.083, 0.371, 0.594, 0.93, 0.41, 0.618, 0.249, 0.297, 0.124, 0.858, 0.916, 0.325, 0.614, 0.656, 0.267, 0.175, 0.065, 0.422, 0.448, 0.826, 0.161, 0.925, 0.195, 0.687, 0.689, 0.264, 0.267, 0.421, 0.735, 0.941, 0.393, 0.503, 0.895, 0.502, 0.232, 0.575, 0.953, 0.01, 0.6, 0.016, 0.106, 0.971, 0.399, 0.626, 0.553, 0.331, 0.582, 0.546, 0.35, 0.384, 0.054, 0.286, 0.149, 0.91, 0.829, 0.975, 0.454, 0.516, 0.222, 0.774, 0.557, 0.42, 0.745, 0.275, 0.97, 0.805, 0.423, 0.209, 0.998, 0.349, 0.833, 0.568, 0.066, 0.219, 0.973, 0.453, 0.516, 0.785, 0.073, 0.505, 0.552, 0.431, 0.771, 0.681, 0.305, 0.795, 0.762, 0.888, 0.503, 0.637, 0.098, 0.81, 0.239, 0.047]
global q = [0.882, 0.971, 0.577, 0.876, 0.917, 0.3, 0.862, 0.791, 0.759, 0.898, 0.22, 0.916, 0.91, 0.861, 0.994, 0.869, 0.596, 0.97, 0.906, 0.75, 0.471, 0.99, 0.911, 0.749, 0.749, 0.945, 0.286, 0.933, 0.751, 0.505, 0.873, 0.979, 0.982, 0.901, 0.397, 0.276, 0.223, 0.28, 0.92, 0.979, 0.957, 0.413, 0.378, 0.816, 0.939, 0.969, 0.935, 0.998, 0.796, 0.498, 0.509, 0.596, 0.906, 0.37, 0.622, 0.734, 0.613, 0.934, 0.301, 0.639, 0.656, 0.995, 0.951, 0.887, 0.955, 0.805, 0.953, 0.733, 0.728, 0.633, 0.869, 0.501, 0.18, 0.667, 0.962, 0.974, 0.511, 0.763, 0.787, 0.912, 0.182, 0.964, 0.944, 0.816, 0.929, 0.72, 0.989, 0.999, 0.91, 0.607, 0.943, 0.95, 0.588, 0.503, 0.642, 0.997, 0.479, 0.815, 0.719, 0.744, 0.887, 0.987, 0.991, 0.978, 0.867, 0.16, 0.842, 0.757, 0.185, 0.78, 0.274, 0.415, 0.92, 0.849, 0.339, 0.73, 0.914, 0.972, 0.687, 0.331, 0.437, 0.827, 0.762, 0.994, 0.269, 0.66, 0.664, 0.88, 0.996, 0.152, 0.879, 0.689, 0.788, 0.926, 0.903, 0.915, 0.965, 0.982, 0.829, 0.729, 0.698, 0.956, 0.474, 0.332, 0.835, 0.93, 0.999, 0.881, 0.98, 0.749, 0.782, 0.574, 0.926, 0.28, 0.684, 0.971, 0.821, 0.723, 0.989, 0.916, 0.994, 0.672, 0.599, 0.698, 0.985, 0.768, 0.881, 0.641, 0.538, 0.512, 0.979, 0.978, 0.367, 0.921, 0.708, 0.606, 0.971, 0.441, 0.861, 0.92, 0.966, 0.833, 0.984, 0.814, 0.88, 0.74, 0.794, 0.644, 0.499, 0.861, 0.967, 0.889, 0.634, 0.962, 0.767, 0.442, 0.739, 0.956, 0.154, 0.899, 0.416, 0.621, 0.995, 0.728, 0.808, 0.589, 0.41, 0.597, 0.894, 0.481, 0.431, 0.893, 0.822, 0.292, 0.982, 0.905, 0.979, 0.477, 0.945, 0.328, 0.988, 0.762, 0.817, 0.813, 0.452, 0.992, 0.986, 0.901, 0.586, 0.998, 0.627, 0.852, 0.83, 0.554, 0.914, 0.997, 0.85, 0.849, 0.883, 0.332, 0.974, 0.795, 0.789, 0.781, 0.892, 0.975, 0.971, 0.889, 0.9, 0.993, 0.926, 0.745, 0.943, 0.975, 0.197]
global origin = 1
global destination = 50