global arcs = [1 4; 1 22; 1 35; 1 38; 1 43; 1 46; 1 56; 1 60; 2 18; 2 40; 3 2; 3 28; 3 34; 3 35; 4 19; 4 21; 4 29; 4 32; 4 50; 4 52; 5 9; 5 10; 5 23; 5 53; 6 3; 6 13; 6 17; 6 57; 6 58; 7 2; 7 15; 7 19; 7 41; 7 45; 7 47; 7 51; 8 2; 8 9; 8 13; 8 16; 8 35; 8 39; 9 3; 9 12; 9 13; 9 18; 9 19; 9 37; 10 5; 10 7; 10 27; 10 35; 10 46; 10 54; 11 6; 11 7; 11 28; 11 29; 11 50; 11 55; 12 11; 12 26; 12 34; 12 56; 13 7; 13 9; 13 16; 13 19; 13 29; 13 34; 13 38; 13 46; 13 58; 14 19; 15 7; 15 26; 15 36; 15 37; 15 38; 15 41; 15 50; 16 5; 16 7; 16 13; 16 26; 16 48; 16 50; 16 58; 17 3; 17 9; 17 26; 17 34; 17 38; 17 51; 18 6; 18 17; 18 19; 18 24; 18 44; 18 55; 19 3; 19 7; 19 10; 19 25; 19 36; 20 4; 20 16; 20 29; 20 31; 20 43; 20 57; 21 18; 21 41; 21 59; 22 3; 22 40; 22 45; 22 52; 23 22; 23 42; 23 44; 24 3; 24 14; 24 36; 24 41; 24 50; 24 53; 25 4; 25 5; 25 31; 25 45; 26 51; 26 59; 27 4; 27 9; 27 18; 27 26; 27 45; 27 49; 27 60; 28 19; 28 25; 28 30; 28 56; 29 3; 29 7; 29 24; 29 34; 29 51; 30 35; 31 2; 31 6; 31 13; 31 22; 31 33; 31 48; 31 51; 32 2; 32 11; 32 31; 32 37; 32 39; 32 52; 33 8; 33 13; 33 27; 33 43; 34 2; 34 3; 34 21; 34 38; 34 52; 35 45; 35 56; 36 45; 36 56; 37 3; 37 28; 37 39; 38 4; 38 23; 38 56; 39 45; 39 50; 39 53; 39 59; 40 2; 40 32; 40 59; 40 60; 41 34; 41 38; 41 53; 41 55; 41 58; 42 20; 42 23; 42 27; 42 31; 42 37; 42 45; 42 60; 43 4; 43 40; 43 47; 43 48; 44 4; 44 14; 44 19; 44 40; 44 51; 45 4; 45 15; 45 17; 45 34; 45 37; 45 46; 45 55; 45 56; 46 35; 46 38; 46 39; 46 47; 46 52; 47 2; 47 13; 47 16; 47 19; 47 34; 47 37; 47 44; 48 10; 48 28; 48 57; 49 6; 49 24; 49 31; 49 42; 49 48; 50 22; 50 34; 50 35; 50 36; 50 46; 51 12; 51 26; 51 27; 51 42; 51 46; 51 57; 52 4; 52 13; 52 53; 52 56; 53 4; 53 9; 53 12; 53 17; 53 28; 53 47; 53 50; 53 55; 54 6; 54 11; 54 13; 54 31; 54 39; 54 50; 55 2; 55 15; 55 21; 55 30; 55 49; 55 51; 55 53; 55 59; 56 4; 56 7; 56 19; 56 23; 56 34; 56 48; 56 60; 57 2; 57 17; 57 19; 57 27; 57 37; 57 38; 57 47; 57 52; 57 56; 58 5; 58 22; 58 38; 59 7; 59 18; 59 29]
global d_x = [1.0, 5.0, 8.0, 7.0, 4.0, 3.0, 3.0, 8.0, 1.0, 9.0, 9.0, 9.0, 1.0, 8.0, 5.0, 9.0, 1.0, 9.0, 6.0, 3.0, 4.0, 7.0, 7.0, 3.0, 10.0, 6.0, 5.0, 5.0, 3.0, 2.0, 7.0, 1.0, 9.0, 6.0, 5.0, 1.0, 2.0, 3.0, 7.0, 6.0, 5.0, 7.0, 9.0, 1.0, 7.0, 1.0, 4.0, 8.0, 4.0, 6.0, 3.0, 2.0, 3.0, 3.0, 10.0, 10.0, 5.0, 3.0, 7.0, 4.0, 9.0, 8.0, 2.0, 10.0, 3.0, 3.0, 8.0, 9.0, 2.0, 2.0, 4.0, 9.0, 5.0, 7.0, 9.0, 9.0, 5.0, 5.0, 1.0, 9.0, 2.0, 2.0, 9.0, 10.0, 6.0, 3.0, 4.0, 2.0, 2.0, 1.0, 8.0, 9.0, 1.0, 3.0, 2.0, 7.0, 9.0, 9.0, 10.0, 5.0, 10.0, 2.0, 3.0, 6.0, 9.0, 4.0, 1.0, 4.0, 1.0, 2.0, 3.0, 1.0, 1.0, 6.0, 7.0, 4.0, 8.0, 1.0, 3.0, 9.0, 4.0, 4.0, 1.0, 5.0, 6.0, 8.0, 8.0, 5.0, 6.0, 8.0, 5.0, 3.0, 6.0, 9.0, 1.0, 8.0, 6.0, 1.0, 8.0, 8.0, 6.0, 1.0, 5.0, 5.0, 6.0, 1.0, 9.0, 8.0, 2.0, 1.0, 4.0, 8.0, 5.0, 7.0, 5.0, 4.0, 7.0, 8.0, 10.0, 8.0, 6.0, 2.0, 3.0, 7.0, 10.0, 6.0, 1.0, 5.0, 9.0, 7.0, 6.0, 2.0, 8.0, 6.0, 3.0, 7.0, 1.0, 7.0, 8.0, 10.0, 10.0, 4.0, 1.0, 6.0, 4.0, 6.0, 4.0, 10.0, 4.0, 1.0, 10.0, 2.0, 1.0, 6.0, 10.0, 7.0, 4.0, 8.0, 3.0, 2.0, 10.0, 6.0, 2.0, 4.0, 10.0, 2.0, 9.0, 7.0, 9.0, 9.0, 3.0, 10.0, 1.0, 6.0, 2.0, 2.0, 5.0, 3.0, 9.0, 4.0, 10.0, 10.0, 2.0, 1.0, 1.0, 6.0, 2.0, 10.0, 8.0, 3.0, 8.0, 1.0, 2.0, 8.0, 7.0, 2.0, 5.0, 10.0, 5.0, 6.0, 3.0, 1.0, 8.0, 10.0, 6.0, 6.0, 9.0, 7.0, 1.0, 2.0, 5.0, 8.0, 8.0, 10.0, 1.0, 5.0, 8.0, 8.0, 6.0, 9.0, 3.0, 1.0, 5.0, 9.0, 7.0, 6.0, 3.0, 10.0, 4.0, 4.0, 6.0, 9.0, 1.0, 5.0, 1.0, 2.0, 6.0, 6.0, 10.0, 3.0, 1.0, 6.0, 4.0, 6.0, 8.0, 5.0, 2.0, 1.0, 6.0, 7.0, 9.0, 9.0, 3.0, 8.0, 6.0, 3.0, 9.0, 9.0]
global b_x = 5
global d_y = [4.0, 7.0, 8.0, 8.0, 4.0, 1.0, 2.0, 1.0, 4.0, 2.0, 9.0, 8.0, 8.0, 3.0, 4.0, 1.0, 3.0, 1.0, 3.0, 2.0, 7.0, 9.0, 10.0, 8.0, 7.0, 7.0, 9.0, 6.0, 5.0, 1.0, 6.0, 10.0, 8.0, 10.0, 9.0, 2.0, 8.0, 7.0, 2.0, 7.0, 9.0, 4.0, 6.0, 6.0, 4.0, 7.0, 2.0, 9.0, 1.0, 6.0, 6.0, 9.0, 8.0, 6.0, 9.0, 7.0, 1.0, 5.0, 5.0, 6.0, 7.0, 9.0, 8.0, 4.0, 1.0, 6.0, 9.0, 9.0, 5.0, 3.0, 6.0, 10.0, 9.0, 7.0, 10.0, 6.0, 4.0, 6.0, 4.0, 4.0, 10.0, 6.0, 6.0, 6.0, 10.0, 4.0, 7.0, 10.0, 3.0, 8.0, 9.0, 7.0, 2.0, 9.0, 9.0, 7.0, 9.0, 1.0, 1.0, 9.0, 5.0, 1.0, 10.0, 6.0, 2.0, 9.0, 4.0, 3.0, 4.0, 6.0, 8.0, 2.0, 10.0, 5.0, 6.0, 9.0, 9.0, 1.0, 1.0, 4.0, 8.0, 2.0, 4.0, 6.0, 4.0, 9.0, 5.0, 5.0, 1.0, 6.0, 10.0, 5.0, 7.0, 5.0, 1.0, 2.0, 6.0, 5.0, 3.0, 8.0, 9.0, 1.0, 5.0, 9.0, 7.0, 6.0, 6.0, 2.0, 3.0, 6.0, 7.0, 5.0, 9.0, 1.0, 1.0, 9.0, 10.0, 10.0, 2.0, 10.0, 6.0, 2.0, 8.0, 7.0, 5.0, 9.0, 6.0, 8.0, 6.0, 4.0, 6.0, 7.0, 10.0, 10.0, 3.0, 7.0, 5.0, 2.0, 7.0, 3.0, 7.0, 6.0, 2.0, 1.0, 3.0, 4.0, 10.0, 5.0, 8.0, 2.0, 3.0, 10.0, 9.0, 3.0, 10.0, 5.0, 7.0, 2.0, 7.0, 9.0, 7.0, 6.0, 5.0, 5.0, 10.0, 9.0, 6.0, 3.0, 6.0, 7.0, 2.0, 9.0, 4.0, 1.0, 1.0, 3.0, 5.0, 9.0, 1.0, 5.0, 4.0, 6.0, 5.0, 2.0, 5.0, 7.0, 7.0, 2.0, 3.0, 8.0, 3.0, 8.0, 2.0, 6.0, 7.0, 5.0, 3.0, 10.0, 4.0, 9.0, 1.0, 9.0, 2.0, 2.0, 7.0, 3.0, 2.0, 9.0, 8.0, 7.0, 2.0, 1.0, 7.0, 4.0, 10.0, 7.0, 6.0, 9.0, 5.0, 10.0, 6.0, 5.0, 6.0, 9.0, 3.0, 7.0, 9.0, 2.0, 8.0, 5.0, 7.0, 2.0, 10.0, 2.0, 8.0, 9.0, 2.0, 2.0, 4.0, 10.0, 6.0, 5.0, 8.0, 6.0, 2.0, 5.0, 1.0, 10.0, 10.0, 1.0, 6.0, 10.0, 4.0, 4.0, 10.0, 8.0, 1.0, 1.0]
global b_y = 10
global p = [0.916, 0.413, 0.008, 0.977, 0.707, 0.038, 0.091, 0.118, 0.238, 0.015, 0.397, 0.812, 0.714, 0.813, 0.302, 0.49, 0.022, 0.981, 0.031, 0.541, 0.315, 0.745, 0.307, 0.787, 0.478, 0.927, 0.801, 0.668, 0.515, 0.845, 0.575, 0.299, 0.926, 0.902, 0.781, 0.038, 0.614, 0.034, 0.606, 0.75, 0.377, 0.906, 0.343, 0.822, 0.981, 0.824, 0.803, 0.245, 0.844, 0.309, 0.364, 0.187, 0.159, 0.154, 0.84, 0.269, 0.328, 0.353, 0.051, 0.558, 0.645, 0.944, 0.168, 0.383, 0.961, 0.911, 0.452, 0.304, 0.476, 0.558, 0.267, 0.181, 0.613, 0.153, 0.295, 0.748, 0.141, 0.805, 0.81, 0.099, 0.31, 0.552, 0.496, 0.374, 0.995, 0.414, 0.842, 0.772, 0.073, 0.763, 0.761, 0.122, 0.556, 0.586, 0.304, 0.895, 0.418, 0.519, 0.013, 0.975, 0.371, 0.71, 0.992, 0.081, 0.707, 0.078, 0.349, 0.214, 0.04, 0.075, 0.671, 0.081, 0.866, 0.626, 0.957, 0.932, 0.556, 0.909, 0.129, 0.635, 0.196, 0.105, 0.239, 0.585, 0.435, 0.232, 0.924, 0.508, 0.599, 0.545, 0.32, 0.158, 0.228, 0.498, 0.002, 0.666, 0.295, 0.016, 0.251, 0.957, 0.294, 0.598, 0.675, 0.583, 0.219, 0.899, 0.164, 0.788, 0.392, 0.122, 0.036, 0.194, 0.227, 0.751, 0.32, 0.314, 0.801, 0.816, 0.436, 0.83, 0.621, 0.412, 0.101, 0.295, 0.71, 0.695, 0.034, 0.373, 0.687, 0.108, 0.066, 0.097, 0.425, 0.584, 0.838, 0.666, 0.426, 0.813, 0.009, 0.867, 0.209, 0.082, 0.741, 0.74, 0.45, 0.897, 0.361, 0.684, 0.518, 0.36, 0.318, 0.92, 0.123, 0.429, 0.184, 0.874, 0.752, 0.295, 0.088, 0.771, 0.851, 0.879, 0.453, 0.264, 0.382, 0.142, 0.034, 0.607, 0.739, 0.953, 0.326, 0.995, 0.767, 0.676, 0.184, 0.145, 0.521, 0.157, 0.602, 0.838, 0.239, 0.054, 0.071, 0.051, 0.531, 0.146, 0.026, 0.738, 0.448, 0.746, 0.108, 0.79, 0.165, 0.545, 0.46, 0.627, 0.87, 0.265, 0.523, 0.339, 0.268, 0.747, 0.374, 0.922, 0.031, 0.646, 0.378, 0.803, 0.836, 0.102, 0.684, 0.285, 0.623, 0.16, 0.797, 0.332, 0.967, 0.489, 0.397, 0.056, 0.755, 0.49, 0.429, 0.994, 0.323, 0.475, 0.522, 0.852, 0.268, 0.978, 0.936, 0.483, 0.426, 0.75, 0.119, 0.034, 0.492, 0.681, 0.349, 0.209, 0.962, 0.657, 0.671, 0.109, 0.722, 0.348, 0.826, 0.288, 0.093, 0.062, 0.267, 0.221, 0.277, 0.232, 0.317, 0.24, 0.586, 0.052]
global q = [0.928, 0.568, 0.426, 0.986, 0.974, 0.634, 0.244, 0.759, 0.765, 0.14, 0.635, 0.841, 0.865, 0.983, 0.604, 0.978, 0.38, 0.991, 0.497, 0.625, 0.86, 0.868, 0.703, 0.946, 0.608, 0.975, 0.863, 0.823, 0.986, 0.929, 0.985, 0.821, 0.955, 0.993, 0.79, 0.226, 0.699, 0.216, 0.906, 0.972, 0.466, 0.937, 0.767, 0.84, 0.981, 0.916, 0.833, 0.912, 0.907, 0.799, 0.49, 0.751, 0.777, 0.802, 0.949, 0.715, 0.982, 0.971, 0.907, 0.786, 0.879, 0.998, 0.349, 0.618, 0.992, 0.952, 0.859, 0.371, 0.625, 0.935, 0.663, 0.968, 0.957, 0.523, 0.695, 0.935, 0.47, 0.957, 0.847, 0.639, 0.861, 0.62, 0.687, 0.541, 0.996, 0.936, 0.872, 0.854, 0.884, 0.97, 0.968, 0.934, 0.882, 0.793, 0.359, 0.928, 0.661, 0.523, 0.906, 0.976, 0.849, 0.769, 0.997, 0.384, 0.91, 0.092, 0.615, 0.501, 0.981, 0.678, 0.75, 0.554, 0.902, 0.994, 0.957, 0.944, 0.808, 0.919, 0.876, 0.698, 0.924, 0.249, 0.829, 0.815, 0.686, 0.946, 0.965, 0.827, 0.757, 0.939, 0.691, 0.889, 0.982, 0.839, 0.567, 0.983, 0.935, 0.021, 0.928, 0.983, 0.856, 0.762, 0.798, 0.897, 0.883, 0.987, 0.653, 0.945, 0.992, 0.787, 0.329, 0.429, 0.575, 0.885, 0.752, 0.408, 0.833, 0.896, 0.76, 0.871, 0.691, 0.943, 0.923, 0.596, 0.925, 0.87, 0.195, 0.447, 0.927, 0.808, 0.771, 0.724, 0.461, 0.983, 0.847, 0.725, 0.781, 0.94, 0.92, 0.936, 0.885, 0.537, 0.889, 0.829, 0.873, 0.925, 0.439, 0.765, 0.596, 0.421, 0.318, 0.951, 0.461, 0.974, 0.586, 0.887, 0.958, 0.637, 0.689, 0.871, 0.889, 0.985, 0.962, 0.774, 0.676, 0.537, 0.407, 0.905, 0.75, 0.994, 0.411, 0.996, 0.954, 0.723, 0.886, 0.735, 0.882, 0.44, 0.683, 0.902, 0.886, 0.473, 0.457, 0.472, 0.909, 0.658, 0.694, 0.881, 0.816, 0.751, 0.656, 0.993, 0.324, 0.712, 0.758, 0.656, 0.922, 0.9, 0.59, 0.61, 0.907, 0.818, 0.422, 0.97, 0.441, 0.734, 0.393, 0.878, 0.873, 0.872, 0.949, 0.705, 0.682, 0.305, 0.823, 0.822, 0.968, 0.616, 0.619, 0.539, 0.906, 0.681, 0.475, 0.996, 0.962, 0.684, 0.681, 0.959, 0.864, 0.982, 0.995, 0.813, 0.519, 0.895, 0.81, 0.235, 0.766, 0.952, 0.723, 0.995, 0.997, 0.776, 0.801, 0.452, 0.982, 0.545, 0.908, 0.373, 0.434, 0.712, 0.568, 0.294, 0.634, 0.506, 0.848, 0.262, 0.606, 0.599]
global origin = 1
global destination = 60