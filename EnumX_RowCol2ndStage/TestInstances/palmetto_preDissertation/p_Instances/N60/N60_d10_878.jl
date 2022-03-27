global arcs = [1 15; 1 30; 1 36; 1 54; 1 56; 2 4; 2 28; 2 38; 2 43; 2 57; 3 2; 3 7; 3 11; 3 19; 3 38; 3 40; 3 44; 3 52; 3 53; 3 57; 4 7; 4 17; 4 30; 4 39; 5 29; 5 34; 5 37; 5 59; 5 60; 6 11; 6 29; 6 49; 6 56; 7 5; 7 14; 7 17; 7 47; 7 52; 8 7; 8 21; 8 47; 8 52; 8 55; 8 59; 9 13; 9 20; 9 23; 9 41; 9 43; 9 53; 10 15; 10 22; 10 36; 10 45; 10 60; 11 2; 11 18; 11 28; 11 32; 11 52; 11 55; 12 2; 12 4; 12 21; 12 22; 12 30; 12 38; 13 23; 13 30; 13 33; 14 3; 14 21; 14 24; 14 25; 14 29; 14 37; 14 39; 14 47; 14 60; 15 2; 15 3; 15 5; 15 17; 15 47; 15 50; 15 52; 15 58; 16 4; 16 23; 16 29; 16 40; 16 41; 16 44; 17 14; 17 16; 17 29; 17 37; 17 40; 17 56; 17 59; 18 2; 18 7; 18 10; 18 28; 18 33; 19 18; 19 33; 19 37; 20 7; 20 12; 20 25; 20 32; 20 33; 20 45; 21 3; 21 14; 21 20; 21 25; 21 56; 21 58; 22 9; 22 11; 22 51; 22 54; 23 5; 23 6; 23 11; 23 21; 23 24; 23 25; 23 42; 23 44; 23 57; 23 60; 24 19; 24 25; 24 38; 24 49; 25 31; 25 60; 26 2; 26 4; 26 5; 26 6; 26 14; 26 16; 26 18; 26 19; 26 21; 26 44; 27 3; 27 6; 27 10; 27 14; 27 16; 27 17; 27 19; 27 22; 27 32; 27 33; 27 38; 27 53; 27 54; 28 30; 29 13; 29 35; 29 42; 29 44; 30 13; 30 23; 30 43; 30 53; 31 10; 31 30; 31 36; 32 12; 32 21; 32 22; 32 29; 32 40; 32 41; 32 47; 33 2; 33 6; 33 9; 33 22; 33 45; 33 52; 34 2; 34 8; 34 28; 34 29; 34 32; 34 36; 34 48; 34 58; 35 2; 35 23; 35 24; 35 25; 35 30; 35 42; 35 52; 35 59; 36 12; 36 13; 36 18; 36 27; 36 31; 36 56; 37 22; 37 41; 37 52; 38 9; 38 16; 38 37; 38 41; 38 52; 39 23; 39 45; 40 3; 40 18; 40 26; 40 30; 40 32; 40 36; 40 53; 40 54; 41 25; 41 34; 41 50; 41 53; 41 55; 41 60; 42 5; 42 12; 42 13; 42 20; 42 60; 43 6; 43 12; 43 24; 43 50; 43 52; 43 57; 44 19; 44 28; 44 30; 44 34; 44 41; 44 45; 44 55; 45 6; 45 12; 45 28; 45 31; 45 46; 46 6; 46 11; 46 27; 46 47; 47 8; 47 13; 47 48; 48 35; 48 37; 48 42; 48 47; 49 5; 49 18; 49 29; 50 6; 50 16; 50 22; 50 29; 50 49; 51 2; 51 3; 51 14; 51 31; 51 33; 51 44; 51 59; 52 4; 52 16; 52 23; 52 26; 52 35; 52 45; 53 21; 53 32; 53 35; 54 9; 54 15; 54 39; 54 43; 54 46; 55 13; 55 36; 55 40; 55 43; 55 45; 55 47; 55 54; 56 9; 56 16; 56 19; 56 26; 56 27; 56 52; 57 5; 57 39; 57 52; 58 13; 58 20; 58 28; 58 29; 58 32; 58 45; 58 49; 58 52; 58 57; 59 16; 59 20; 59 22; 59 32; 59 34]
global d_x = [4.0, 4.0, 10.0, 2.0, 10.0, 6.0, 6.0, 5.0, 4.0, 2.0, 2.0, 10.0, 10.0, 9.0, 4.0, 7.0, 6.0, 4.0, 7.0, 3.0, 1.0, 10.0, 1.0, 7.0, 3.0, 2.0, 4.0, 10.0, 5.0, 4.0, 1.0, 8.0, 2.0, 3.0, 4.0, 5.0, 4.0, 2.0, 1.0, 4.0, 4.0, 10.0, 4.0, 7.0, 1.0, 5.0, 5.0, 7.0, 3.0, 5.0, 2.0, 5.0, 6.0, 9.0, 1.0, 8.0, 7.0, 6.0, 9.0, 6.0, 8.0, 10.0, 7.0, 4.0, 3.0, 2.0, 8.0, 10.0, 6.0, 4.0, 2.0, 9.0, 6.0, 6.0, 10.0, 2.0, 6.0, 6.0, 8.0, 1.0, 7.0, 4.0, 10.0, 2.0, 5.0, 7.0, 7.0, 10.0, 7.0, 9.0, 1.0, 4.0, 1.0, 8.0, 7.0, 5.0, 8.0, 4.0, 1.0, 10.0, 10.0, 2.0, 1.0, 4.0, 3.0, 10.0, 5.0, 4.0, 2.0, 9.0, 2.0, 1.0, 7.0, 6.0, 10.0, 7.0, 6.0, 8.0, 8.0, 10.0, 2.0, 2.0, 1.0, 6.0, 6.0, 5.0, 6.0, 9.0, 2.0, 9.0, 6.0, 5.0, 3.0, 1.0, 5.0, 1.0, 6.0, 8.0, 3.0, 7.0, 7.0, 9.0, 6.0, 4.0, 10.0, 10.0, 5.0, 10.0, 7.0, 8.0, 6.0, 1.0, 10.0, 4.0, 8.0, 9.0, 2.0, 10.0, 6.0, 4.0, 1.0, 5.0, 6.0, 8.0, 4.0, 1.0, 2.0, 1.0, 5.0, 10.0, 10.0, 3.0, 2.0, 3.0, 2.0, 6.0, 6.0, 9.0, 3.0, 1.0, 6.0, 2.0, 8.0, 5.0, 7.0, 9.0, 4.0, 9.0, 9.0, 10.0, 7.0, 1.0, 5.0, 9.0, 2.0, 9.0, 5.0, 5.0, 3.0, 5.0, 7.0, 9.0, 9.0, 7.0, 5.0, 7.0, 10.0, 10.0, 2.0, 9.0, 4.0, 7.0, 2.0, 8.0, 1.0, 5.0, 3.0, 8.0, 6.0, 1.0, 1.0, 9.0, 1.0, 7.0, 6.0, 8.0, 1.0, 2.0, 9.0, 10.0, 3.0, 5.0, 9.0, 9.0, 7.0, 8.0, 7.0, 9.0, 3.0, 5.0, 2.0, 8.0, 7.0, 4.0, 7.0, 2.0, 7.0, 8.0, 3.0, 9.0, 5.0, 3.0, 7.0, 4.0, 9.0, 5.0, 10.0, 1.0, 5.0, 9.0, 1.0, 5.0, 3.0, 7.0, 2.0, 3.0, 8.0, 3.0, 5.0, 5.0, 6.0, 2.0, 3.0, 4.0, 1.0, 3.0, 7.0, 5.0, 2.0, 9.0, 4.0, 10.0, 4.0, 9.0, 5.0, 6.0, 4.0, 9.0, 4.0, 8.0, 4.0, 4.0, 8.0, 4.0, 8.0, 6.0, 4.0, 8.0, 3.0, 7.0, 7.0, 7.0, 9.0, 6.0, 6.0, 2.0, 8.0, 2.0, 8.0, 3.0, 10.0, 9.0, 10.0, 4.0, 1.0, 5.0, 1.0, 9.0, 1.0, 1.0, 2.0, 10.0, 1.0, 8.0, 3.0, 5.0, 10.0]
global b_x = 5
global d_y = [2.0, 8.0, 5.0, 6.0, 4.0, 5.0, 6.0, 6.0, 7.0, 9.0, 8.0, 7.0, 9.0, 1.0, 10.0, 9.0, 3.0, 9.0, 9.0, 2.0, 6.0, 5.0, 8.0, 2.0, 3.0, 2.0, 2.0, 5.0, 5.0, 9.0, 2.0, 3.0, 10.0, 8.0, 7.0, 3.0, 8.0, 1.0, 6.0, 1.0, 2.0, 2.0, 3.0, 1.0, 3.0, 7.0, 5.0, 1.0, 6.0, 1.0, 6.0, 1.0, 4.0, 4.0, 9.0, 4.0, 10.0, 9.0, 6.0, 3.0, 7.0, 5.0, 6.0, 4.0, 10.0, 5.0, 3.0, 1.0, 7.0, 6.0, 9.0, 2.0, 2.0, 10.0, 9.0, 4.0, 1.0, 3.0, 7.0, 2.0, 7.0, 5.0, 4.0, 5.0, 3.0, 10.0, 4.0, 3.0, 4.0, 8.0, 5.0, 5.0, 4.0, 5.0, 2.0, 4.0, 2.0, 4.0, 9.0, 7.0, 9.0, 5.0, 8.0, 3.0, 5.0, 2.0, 7.0, 4.0, 3.0, 10.0, 10.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 7.0, 2.0, 6.0, 10.0, 3.0, 9.0, 6.0, 10.0, 2.0, 1.0, 2.0, 8.0, 10.0, 10.0, 9.0, 4.0, 4.0, 10.0, 2.0, 6.0, 1.0, 8.0, 10.0, 10.0, 4.0, 9.0, 3.0, 1.0, 3.0, 2.0, 4.0, 10.0, 7.0, 3.0, 4.0, 5.0, 4.0, 4.0, 1.0, 5.0, 8.0, 8.0, 6.0, 2.0, 1.0, 8.0, 7.0, 1.0, 9.0, 2.0, 3.0, 10.0, 6.0, 7.0, 10.0, 4.0, 7.0, 4.0, 9.0, 6.0, 2.0, 3.0, 6.0, 7.0, 9.0, 2.0, 5.0, 7.0, 1.0, 4.0, 10.0, 4.0, 5.0, 4.0, 8.0, 1.0, 5.0, 4.0, 10.0, 3.0, 7.0, 5.0, 10.0, 9.0, 5.0, 9.0, 10.0, 6.0, 5.0, 2.0, 5.0, 3.0, 4.0, 4.0, 9.0, 7.0, 3.0, 3.0, 5.0, 10.0, 6.0, 1.0, 6.0, 7.0, 5.0, 8.0, 3.0, 4.0, 4.0, 9.0, 6.0, 2.0, 6.0, 2.0, 4.0, 8.0, 7.0, 3.0, 8.0, 10.0, 5.0, 2.0, 4.0, 9.0, 7.0, 10.0, 1.0, 3.0, 9.0, 10.0, 1.0, 7.0, 3.0, 6.0, 1.0, 1.0, 8.0, 5.0, 3.0, 8.0, 4.0, 8.0, 5.0, 6.0, 6.0, 10.0, 9.0, 3.0, 10.0, 10.0, 6.0, 8.0, 5.0, 1.0, 2.0, 7.0, 2.0, 8.0, 7.0, 4.0, 6.0, 4.0, 10.0, 1.0, 4.0, 2.0, 3.0, 9.0, 6.0, 10.0, 3.0, 2.0, 7.0, 6.0, 6.0, 1.0, 7.0, 10.0, 7.0, 10.0, 3.0, 2.0, 3.0, 10.0, 6.0, 4.0, 3.0, 8.0, 9.0, 1.0, 4.0, 8.0, 6.0, 8.0, 4.0, 5.0, 6.0, 7.0, 6.0, 4.0, 5.0, 2.0, 8.0, 9.0, 1.0, 2.0, 8.0, 5.0, 2.0, 2.0]
global b_y = 10
global p = [0.263, 0.205, 0.346, 0.939, 0.053, 0.042, 0.154, 0.931, 0.878, 0.018, 0.306, 0.224, 0.103, 0.403, 0.698, 0.141, 0.278, 0.201, 0.515, 0.763, 0.756, 0.453, 0.33, 0.274, 0.272, 0.644, 0.655, 0.273, 0.47, 0.693, 0.996, 0.067, 0.524, 0.72, 0.786, 0.615, 0.25, 0.823, 0.953, 0.733, 0.994, 0.324, 0.38, 0.323, 0.839, 0.635, 0.15, 0.445, 0.286, 0.058, 0.419, 0.67, 0.291, 0.084, 0.05, 0.703, 0.628, 0.422, 0.307, 0.107, 0.357, 0.734, 0.7, 0.657, 0.43, 0.364, 0.125, 0.24, 0.988, 0.9, 0.81, 0.902, 0.392, 0.042, 0.024, 0.18, 0.469, 0.644, 0.36, 0.559, 0.144, 0.735, 0.778, 0.332, 0.27, 0.338, 0.835, 0.788, 0.779, 0.573, 0.25, 0.352, 0.181, 0.584, 0.247, 0.798, 0.268, 0.473, 0.304, 0.51, 0.258, 0.58, 0.797, 0.441, 0.345, 0.824, 0.144, 0.292, 0.219, 0.855, 0.383, 0.726, 0.572, 0.1, 0.149, 0.918, 0.441, 0.146, 0.825, 0.373, 0.202, 0.595, 0.173, 0.918, 0.914, 0.965, 0.638, 0.757, 0.864, 0.911, 0.952, 0.045, 0.834, 0.119, 0.635, 0.6, 0.748, 0.891, 0.75, 0.739, 0.971, 0.504, 0.816, 0.19, 0.641, 0.16, 0.528, 0.984, 0.197, 0.196, 0.966, 0.779, 0.355, 0.4, 0.619, 0.56, 0.625, 0.894, 0.464, 0.902, 0.497, 0.159, 0.491, 0.145, 0.988, 0.365, 0.775, 0.908, 0.904, 0.829, 0.552, 0.33, 0.473, 0.661, 0.611, 0.409, 0.625, 0.314, 0.333, 0.927, 0.835, 0.575, 0.479, 0.07, 0.06, 0.072, 0.636, 0.679, 0.663, 0.094, 0.746, 0.19, 0.584, 0.202, 0.017, 0.803, 0.471, 0.994, 0.175, 0.566, 0.423, 0.074, 0.227, 0.686, 0.426, 0.093, 0.865, 0.808, 0.461, 0.21, 0.113, 0.911, 0.286, 0.475, 0.773, 0.172, 0.558, 0.377, 0.523, 0.332, 0.069, 0.128, 0.477, 0.855, 0.506, 0.889, 0.401, 0.917, 0.466, 0.101, 0.607, 0.454, 0.367, 0.178, 0.165, 0.353, 0.754, 0.183, 0.848, 0.857, 0.068, 0.517, 0.821, 0.722, 0.816, 0.118, 0.989, 0.694, 0.162, 0.415, 0.802, 0.477, 0.399, 0.11, 0.404, 0.199, 0.418, 0.972, 0.861, 0.846, 0.561, 0.731, 0.928, 0.908, 0.392, 0.709, 0.212, 0.766, 0.866, 0.169, 0.096, 0.863, 0.561, 0.509, 0.067, 0.026, 0.162, 0.898, 0.058, 0.765, 0.945, 0.141, 0.478, 0.337, 0.677, 0.133, 0.023, 0.629, 0.485, 0.103, 0.759, 0.704, 0.539, 0.239, 0.05, 0.161, 0.751, 0.74, 0.354, 0.083, 0.683, 0.189, 0.483, 0.092, 0.46, 0.815, 0.906, 0.598, 0.684, 0.398, 0.465, 0.978, 0.505, 0.218, 0.63, 0.918, 0.701, 0.52, 0.767, 0.017, 0.589, 0.122, 0.838, 0.019, 0.648, 0.633, 0.756]
global q = [0.694, 0.948, 0.911, 0.967, 0.208, 0.52, 0.658, 0.985, 0.959, 0.77, 0.513, 0.536, 0.768, 0.859, 0.772, 0.855, 0.789, 0.563, 0.848, 0.841, 0.91, 0.708, 0.553, 0.905, 0.906, 0.726, 0.937, 0.391, 0.837, 0.914, 0.998, 0.871, 0.747, 0.959, 0.848, 0.727, 0.253, 0.924, 0.992, 0.908, 0.997, 0.659, 0.682, 0.479, 0.931, 0.769, 0.899, 0.806, 0.817, 0.559, 0.477, 0.949, 0.593, 0.436, 0.917, 0.853, 0.794, 0.859, 0.799, 0.32, 0.491, 0.956, 0.864, 0.737, 0.484, 0.793, 0.557, 0.513, 0.992, 0.98, 0.964, 0.932, 0.399, 0.123, 0.866, 0.704, 0.625, 0.716, 0.497, 0.571, 0.211, 0.805, 0.85, 0.963, 0.597, 0.478, 0.906, 0.801, 0.932, 0.712, 0.838, 0.382, 0.296, 0.873, 0.956, 0.837, 0.292, 0.747, 0.905, 0.884, 0.322, 0.782, 0.928, 0.689, 0.747, 0.826, 0.474, 0.331, 0.797, 0.975, 0.452, 0.749, 0.778, 0.925, 0.504, 0.932, 0.58, 0.418, 0.851, 0.931, 0.791, 0.623, 0.448, 0.995, 0.972, 0.967, 0.783, 0.991, 0.999, 0.931, 0.97, 0.091, 0.896, 0.883, 0.897, 0.866, 0.917, 0.982, 0.982, 0.761, 0.991, 0.772, 0.94, 0.664, 0.721, 0.459, 0.575, 0.984, 0.252, 0.768, 0.968, 0.968, 0.357, 0.682, 0.849, 0.755, 0.76, 0.902, 0.996, 0.922, 0.789, 0.802, 0.652, 0.321, 0.988, 0.758, 0.889, 0.942, 0.938, 0.845, 0.9, 0.888, 0.768, 0.891, 0.637, 0.804, 0.696, 0.782, 0.983, 0.992, 0.868, 0.936, 0.854, 0.492, 0.155, 0.224, 0.834, 0.727, 0.683, 0.867, 0.86, 0.829, 0.913, 0.79, 0.519, 0.854, 0.541, 0.997, 0.345, 0.861, 0.964, 0.312, 0.771, 0.802, 0.708, 0.228, 0.89, 0.922, 0.817, 0.94, 0.557, 0.954, 0.541, 0.785, 0.91, 0.758, 0.765, 0.976, 0.677, 0.989, 0.782, 0.723, 0.657, 0.988, 0.533, 0.941, 0.828, 0.998, 0.745, 0.754, 0.761, 0.853, 0.927, 0.584, 0.437, 0.427, 0.833, 0.347, 0.993, 0.947, 0.875, 0.584, 0.88, 0.867, 0.827, 0.709, 0.994, 0.984, 0.696, 0.605, 0.833, 0.656, 0.597, 0.96, 0.478, 0.35, 0.665, 0.978, 0.971, 0.908, 0.618, 0.872, 0.967, 0.987, 0.71, 0.781, 0.369, 0.903, 0.929, 0.97, 0.573, 0.878, 0.813, 0.516, 0.124, 0.875, 0.478, 0.949, 0.258, 0.833, 0.993, 0.376, 0.486, 0.453, 0.894, 0.721, 0.027, 0.798, 0.543, 0.74, 0.968, 0.872, 0.562, 0.999, 0.12, 0.464, 0.796, 0.796, 0.748, 0.898, 0.837, 0.195, 0.825, 0.177, 0.726, 0.876, 0.968, 0.612, 0.759, 0.613, 0.751, 0.989, 0.552, 0.687, 0.899, 0.987, 0.97, 0.818, 0.802, 0.18, 0.667, 0.485, 0.914, 0.161, 0.909, 0.917, 0.866]
global origin = 1
global destination = 60