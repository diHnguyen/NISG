global arcs = [1 12; 1 13; 1 18; 1 40; 1 44; 2 7; 2 14; 2 25; 2 35; 2 44; 2 56; 3 6; 3 21; 3 24; 3 28; 3 47; 3 49; 4 3; 4 56; 5 22; 5 32; 5 33; 5 56; 6 3; 6 13; 6 15; 6 32; 6 33; 6 59; 7 5; 7 8; 7 12; 7 14; 7 25; 7 41; 7 56; 8 25; 8 26; 8 30; 8 34; 8 57; 9 3; 9 13; 9 40; 9 57; 9 58; 10 32; 10 33; 10 49; 11 3; 11 5; 11 20; 11 33; 11 57; 11 58; 11 59; 12 14; 12 17; 12 40; 12 41; 12 54; 13 2; 13 12; 13 33; 13 39; 13 57; 14 6; 14 7; 14 18; 14 19; 14 20; 14 40; 14 51; 14 55; 14 58; 15 5; 15 9; 15 55; 15 56; 15 58; 15 59; 16 18; 16 21; 16 26; 16 48; 16 59; 17 6; 17 8; 17 20; 17 44; 17 45; 17 49; 17 54; 18 10; 18 19; 19 9; 19 35; 19 53; 20 31; 20 47; 21 12; 21 17; 21 44; 21 47; 21 53; 21 60; 22 7; 22 8; 22 50; 22 54; 23 7; 23 14; 23 15; 23 20; 23 26; 23 35; 23 41; 23 43; 24 16; 24 17; 24 54; 24 55; 24 56; 25 13; 25 15; 25 16; 25 19; 25 27; 25 29; 26 3; 26 4; 26 8; 26 12; 26 30; 26 33; 26 45; 26 58; 27 14; 27 15; 27 39; 27 43; 27 56; 28 22; 28 29; 28 56; 29 8; 29 10; 29 34; 29 55; 30 11; 30 16; 30 24; 30 32; 30 33; 30 35; 30 38; 30 41; 30 42; 30 54; 30 58; 31 3; 31 15; 31 18; 31 34; 31 35; 31 53; 32 2; 32 16; 32 35; 32 36; 32 47; 32 52; 32 59; 33 9; 33 11; 33 13; 33 14; 33 46; 34 25; 34 32; 34 38; 34 40; 34 58; 35 3; 35 11; 35 18; 35 31; 35 32; 35 42; 35 43; 36 3; 36 6; 36 14; 36 26; 36 46; 36 48; 36 57; 37 24; 37 35; 37 36; 37 39; 37 42; 37 46; 37 60; 38 5; 38 6; 39 6; 39 10; 39 14; 39 15; 39 23; 39 34; 39 44; 39 58; 40 7; 40 13; 40 15; 40 20; 40 57; 41 5; 41 6; 41 18; 42 9; 42 10; 42 44; 42 46; 43 2; 43 12; 43 19; 43 39; 43 55; 43 59; 44 14; 44 18; 44 22; 44 28; 44 30; 44 51; 45 3; 45 6; 45 10; 45 19; 46 36; 46 55; 47 2; 47 14; 47 28; 47 42; 47 48; 47 52; 47 54; 47 55; 48 2; 48 12; 48 28; 48 37; 48 40; 48 58; 49 4; 49 21; 49 27; 49 32; 49 39; 49 41; 49 46; 49 56; 50 9; 50 12; 50 15; 50 17; 50 31; 50 34; 50 42; 50 44; 51 20; 51 30; 51 31; 51 34; 51 36; 52 22; 52 33; 52 40; 52 59; 53 3; 53 6; 53 15; 53 20; 53 24; 53 49; 53 50; 54 5; 54 25; 54 30; 54 36; 54 60; 55 4; 55 11; 55 12; 55 23; 55 29; 55 42; 55 52; 55 54; 56 5; 56 13; 56 20; 56 46; 56 49; 56 50; 57 16; 57 31; 57 49; 57 50; 58 3; 58 6; 58 11; 58 12; 58 26; 58 27; 58 39; 58 46; 58 53; 58 55; 59 4; 59 7; 59 41; 59 53]
global d_x = [10.0, 2.0, 2.0, 3.0, 6.0, 7.0, 2.0, 9.0, 5.0, 1.0, 2.0, 1.0, 7.0, 10.0, 1.0, 8.0, 7.0, 4.0, 1.0, 3.0, 4.0, 9.0, 5.0, 2.0, 5.0, 7.0, 10.0, 9.0, 6.0, 10.0, 9.0, 4.0, 4.0, 2.0, 4.0, 10.0, 4.0, 8.0, 4.0, 5.0, 7.0, 7.0, 3.0, 2.0, 1.0, 3.0, 2.0, 6.0, 1.0, 4.0, 4.0, 6.0, 9.0, 6.0, 5.0, 6.0, 4.0, 1.0, 6.0, 8.0, 1.0, 8.0, 1.0, 9.0, 2.0, 8.0, 8.0, 6.0, 6.0, 5.0, 1.0, 7.0, 1.0, 4.0, 4.0, 7.0, 6.0, 8.0, 9.0, 4.0, 1.0, 5.0, 9.0, 2.0, 5.0, 2.0, 4.0, 6.0, 8.0, 4.0, 10.0, 4.0, 3.0, 1.0, 6.0, 7.0, 9.0, 10.0, 1.0, 4.0, 7.0, 3.0, 1.0, 3.0, 1.0, 3.0, 8.0, 5.0, 5.0, 10.0, 6.0, 1.0, 9.0, 1.0, 4.0, 10.0, 1.0, 6.0, 4.0, 9.0, 7.0, 6.0, 5.0, 10.0, 9.0, 10.0, 2.0, 7.0, 7.0, 3.0, 3.0, 9.0, 3.0, 6.0, 9.0, 1.0, 4.0, 4.0, 7.0, 1.0, 4.0, 2.0, 6.0, 1.0, 3.0, 6.0, 9.0, 9.0, 5.0, 10.0, 10.0, 4.0, 3.0, 10.0, 5.0, 9.0, 5.0, 6.0, 7.0, 3.0, 2.0, 7.0, 4.0, 9.0, 5.0, 4.0, 2.0, 2.0, 4.0, 6.0, 4.0, 8.0, 6.0, 2.0, 3.0, 5.0, 5.0, 10.0, 8.0, 8.0, 9.0, 10.0, 4.0, 8.0, 7.0, 6.0, 9.0, 5.0, 3.0, 5.0, 10.0, 3.0, 9.0, 7.0, 1.0, 10.0, 5.0, 8.0, 10.0, 3.0, 2.0, 6.0, 8.0, 10.0, 3.0, 10.0, 7.0, 10.0, 2.0, 5.0, 6.0, 3.0, 5.0, 8.0, 5.0, 2.0, 6.0, 3.0, 2.0, 4.0, 9.0, 6.0, 2.0, 4.0, 1.0, 7.0, 5.0, 9.0, 9.0, 9.0, 4.0, 8.0, 1.0, 8.0, 8.0, 1.0, 3.0, 7.0, 5.0, 8.0, 7.0, 4.0, 7.0, 5.0, 4.0, 1.0, 3.0, 2.0, 6.0, 10.0, 10.0, 4.0, 1.0, 6.0, 2.0, 2.0, 7.0, 7.0, 10.0, 8.0, 5.0, 1.0, 8.0, 2.0, 5.0, 3.0, 2.0, 8.0, 4.0, 3.0, 7.0, 5.0, 2.0, 8.0, 10.0, 3.0, 8.0, 4.0, 9.0, 9.0, 7.0, 10.0, 6.0, 4.0, 9.0, 9.0, 7.0, 2.0, 9.0, 3.0, 2.0, 5.0, 9.0, 7.0, 10.0, 1.0, 9.0, 8.0, 7.0, 4.0, 3.0, 10.0, 6.0, 2.0, 2.0, 8.0, 5.0, 10.0, 4.0, 4.0, 7.0, 7.0, 4.0, 10.0, 2.0, 8.0, 9.0, 6.0, 5.0, 7.0, 10.0, 10.0, 10.0, 2.0, 1.0, 5.0, 7.0]
global b_x = 5
global d_y = [9.0, 6.0, 3.0, 5.0, 4.0, 8.0, 9.0, 8.0, 4.0, 1.0, 10.0, 5.0, 8.0, 9.0, 2.0, 3.0, 5.0, 6.0, 7.0, 8.0, 3.0, 10.0, 7.0, 1.0, 6.0, 8.0, 8.0, 1.0, 4.0, 2.0, 5.0, 2.0, 1.0, 1.0, 3.0, 5.0, 2.0, 8.0, 10.0, 3.0, 9.0, 1.0, 7.0, 2.0, 4.0, 10.0, 9.0, 4.0, 8.0, 2.0, 8.0, 6.0, 5.0, 4.0, 4.0, 9.0, 9.0, 4.0, 2.0, 9.0, 3.0, 10.0, 3.0, 10.0, 2.0, 6.0, 2.0, 9.0, 1.0, 2.0, 5.0, 8.0, 8.0, 4.0, 3.0, 1.0, 9.0, 10.0, 10.0, 1.0, 8.0, 3.0, 5.0, 7.0, 10.0, 9.0, 1.0, 1.0, 1.0, 5.0, 1.0, 8.0, 9.0, 3.0, 5.0, 8.0, 4.0, 10.0, 4.0, 4.0, 8.0, 1.0, 3.0, 4.0, 1.0, 8.0, 9.0, 5.0, 7.0, 3.0, 1.0, 2.0, 7.0, 4.0, 7.0, 2.0, 1.0, 7.0, 6.0, 8.0, 6.0, 10.0, 6.0, 5.0, 9.0, 6.0, 8.0, 4.0, 9.0, 2.0, 7.0, 4.0, 5.0, 6.0, 5.0, 5.0, 1.0, 3.0, 5.0, 3.0, 7.0, 3.0, 3.0, 2.0, 8.0, 8.0, 6.0, 7.0, 10.0, 9.0, 8.0, 8.0, 3.0, 9.0, 7.0, 9.0, 9.0, 8.0, 1.0, 7.0, 1.0, 9.0, 5.0, 4.0, 9.0, 1.0, 1.0, 4.0, 3.0, 8.0, 8.0, 5.0, 8.0, 5.0, 5.0, 1.0, 1.0, 2.0, 1.0, 7.0, 6.0, 7.0, 4.0, 4.0, 6.0, 9.0, 3.0, 5.0, 9.0, 10.0, 6.0, 6.0, 6.0, 9.0, 10.0, 5.0, 7.0, 5.0, 6.0, 3.0, 10.0, 3.0, 2.0, 7.0, 5.0, 7.0, 3.0, 4.0, 6.0, 4.0, 5.0, 2.0, 6.0, 1.0, 10.0, 7.0, 10.0, 6.0, 10.0, 10.0, 10.0, 3.0, 5.0, 5.0, 2.0, 5.0, 5.0, 1.0, 4.0, 6.0, 8.0, 9.0, 1.0, 3.0, 4.0, 7.0, 3.0, 5.0, 7.0, 10.0, 7.0, 8.0, 1.0, 2.0, 9.0, 9.0, 5.0, 9.0, 5.0, 5.0, 6.0, 2.0, 10.0, 7.0, 4.0, 7.0, 2.0, 5.0, 4.0, 3.0, 6.0, 2.0, 1.0, 1.0, 7.0, 10.0, 10.0, 10.0, 4.0, 7.0, 7.0, 8.0, 5.0, 9.0, 10.0, 3.0, 9.0, 8.0, 1.0, 6.0, 1.0, 6.0, 8.0, 3.0, 10.0, 4.0, 6.0, 6.0, 6.0, 7.0, 9.0, 9.0, 5.0, 6.0, 1.0, 2.0, 10.0, 7.0, 5.0, 3.0, 3.0, 3.0, 1.0, 6.0, 10.0, 5.0, 1.0, 9.0, 8.0, 4.0, 2.0, 5.0, 6.0, 7.0, 1.0, 9.0, 1.0, 2.0, 4.0, 2.0, 2.0, 7.0, 1.0, 3.0, 9.0, 5.0, 4.0]
global b_y = 10
global p = [0.223, 0.029, 0.658, 0.307, 0.433, 0.469, 0.657, 0.206, 0.056, 0.71, 0.092, 0.344, 0.912, 0.866, 0.125, 0.503, 0.646, 0.885, 0.781, 0.627, 0.617, 0.23, 0.579, 0.215, 0.34, 0.953, 0.376, 0.372, 0.62, 0.771, 0.989, 0.113, 0.174, 0.603, 0.24, 0.314, 0.068, 0.4, 0.516, 0.629, 0.902, 0.728, 0.734, 0.27, 0.189, 0.217, 0.202, 0.502, 0.256, 0.13, 0.543, 0.661, 0.48, 0.525, 0.629, 0.892, 0.976, 0.958, 0.27, 0.871, 0.709, 0.566, 0.091, 0.216, 0.007, 0.531, 0.66, 0.88, 0.94, 0.341, 0.586, 0.347, 0.024, 0.802, 0.983, 0.384, 0.224, 0.929, 0.314, 0.998, 0.417, 0.884, 0.759, 0.31, 0.341, 0.45, 0.75, 0.321, 0.257, 0.736, 0.443, 0.877, 0.685, 0.636, 0.382, 0.456, 0.108, 0.03, 0.188, 0.635, 0.496, 0.214, 0.243, 0.286, 0.5, 0.533, 0.8, 0.584, 0.829, 0.623, 0.627, 0.152, 0.854, 0.876, 0.156, 0.802, 0.046, 0.338, 0.451, 0.481, 0.032, 0.166, 0.739, 0.838, 0.021, 0.682, 0.873, 0.159, 0.338, 0.507, 0.241, 0.79, 0.142, 0.96, 0.71, 0.343, 0.954, 0.028, 0.192, 0.146, 0.768, 0.399, 0.048, 0.748, 0.998, 0.145, 0.303, 0.241, 0.186, 0.619, 0.511, 0.096, 0.845, 0.703, 0.461, 0.028, 0.226, 0.447, 0.239, 0.584, 0.244, 0.196, 0.586, 0.185, 0.896, 0.575, 0.105, 0.119, 0.643, 0.556, 0.98, 0.58, 0.167, 0.005, 0.975, 0.254, 0.791, 0.584, 0.055, 0.68, 0.478, 0.486, 0.421, 0.284, 0.127, 0.163, 0.715, 0.438, 0.09, 0.769, 0.597, 0.142, 0.962, 0.337, 0.648, 0.897, 0.044, 0.609, 0.526, 0.367, 0.513, 0.772, 0.171, 0.924, 0.171, 0.508, 0.337, 0.628, 0.974, 0.659, 0.85, 0.481, 0.754, 0.122, 0.137, 0.387, 0.576, 0.996, 0.396, 0.851, 0.145, 0.11, 0.847, 0.857, 0.337, 0.963, 0.233, 0.98, 0.416, 0.708, 0.918, 0.29, 0.684, 0.359, 0.715, 0.273, 0.848, 0.017, 0.82, 0.088, 0.569, 0.079, 0.74, 0.867, 0.386, 0.086, 0.777, 0.213, 0.095, 0.118, 0.295, 0.523, 0.052, 0.058, 0.88, 0.111, 0.323, 0.98, 0.167, 0.77, 0.443, 0.337, 0.266, 0.68, 0.507, 0.838, 0.835, 0.172, 0.001, 0.21, 0.27, 0.344, 0.823, 0.467, 0.834, 0.818, 0.186, 0.561, 0.373, 0.641, 0.264, 0.608, 0.941, 0.345, 0.27, 0.634, 0.948, 0.424, 0.672, 0.804, 0.696, 0.885, 0.321, 0.225, 0.485, 0.586, 0.64, 0.811, 0.467, 0.564, 0.911, 0.608, 0.366, 0.05, 0.724, 0.036, 0.83, 0.491, 0.962, 0.549, 0.687, 0.976, 0.191, 0.179, 0.919, 0.091, 0.911, 0.331, 0.419, 0.496, 0.73, 0.39, 0.306, 0.38, 0.957, 0.42, 0.098]
global q = [0.906, 0.82, 0.976, 0.662, 0.942, 0.484, 0.718, 0.964, 0.11, 0.783, 0.346, 0.421, 0.961, 0.911, 0.93, 0.614, 0.723, 0.904, 0.884, 0.84, 0.93, 0.402, 0.798, 0.434, 0.398, 0.982, 0.805, 0.751, 0.635, 0.781, 0.989, 0.471, 0.398, 0.706, 0.616, 0.734, 0.302, 0.455, 0.595, 0.723, 0.971, 0.938, 0.933, 0.952, 0.715, 0.526, 0.913, 0.914, 0.477, 0.261, 0.819, 0.951, 0.688, 0.767, 0.676, 0.98, 0.991, 0.965, 0.463, 0.967, 0.809, 0.654, 0.285, 0.548, 0.599, 0.751, 0.667, 0.973, 0.973, 0.827, 0.89, 0.503, 0.983, 0.88, 0.988, 0.853, 0.796, 0.995, 0.67, 0.999, 0.63, 0.971, 0.825, 0.541, 0.424, 0.63, 0.883, 0.594, 0.402, 0.848, 0.852, 0.986, 0.878, 0.936, 0.973, 0.572, 0.489, 0.066, 0.888, 0.855, 0.732, 0.805, 0.265, 0.903, 0.633, 0.649, 0.841, 0.678, 0.958, 0.625, 0.979, 0.987, 0.869, 0.896, 0.57, 0.921, 0.665, 0.71, 0.909, 0.599, 0.528, 0.87, 0.958, 0.934, 0.473, 0.775, 0.999, 0.603, 0.842, 0.843, 0.504, 0.999, 0.878, 0.987, 0.891, 0.977, 0.994, 0.303, 0.974, 0.747, 0.872, 0.684, 0.629, 0.749, 0.998, 0.941, 0.62, 0.723, 0.533, 0.726, 0.736, 0.855, 0.941, 0.827, 0.662, 0.536, 0.536, 0.988, 0.534, 0.795, 0.975, 0.755, 0.914, 0.55, 0.95, 0.579, 0.928, 0.928, 0.799, 0.851, 0.982, 0.834, 0.741, 0.171, 0.992, 0.839, 0.871, 0.833, 0.593, 0.819, 0.902, 0.605, 0.982, 0.35, 0.55, 0.303, 0.897, 0.45, 0.985, 0.802, 0.827, 0.302, 0.964, 0.969, 0.731, 0.975, 0.429, 0.91, 0.881, 0.507, 0.712, 0.792, 0.515, 0.987, 0.939, 0.996, 0.554, 0.994, 0.978, 0.762, 0.87, 0.779, 0.892, 0.913, 0.38, 0.701, 0.687, 0.997, 0.573, 0.89, 0.229, 0.938, 0.931, 0.867, 0.785, 0.972, 0.427, 0.996, 0.849, 0.962, 0.97, 0.81, 0.684, 0.508, 0.857, 0.708, 0.912, 0.928, 0.84, 0.74, 0.585, 0.418, 0.777, 0.931, 0.542, 0.271, 0.816, 0.75, 0.612, 0.671, 0.795, 0.905, 0.741, 0.843, 0.94, 0.774, 0.688, 0.98, 0.849, 0.813, 0.795, 0.506, 0.459, 0.748, 0.621, 0.937, 0.835, 0.305, 0.372, 0.884, 0.459, 0.647, 0.913, 0.67, 0.974, 0.847, 0.309, 0.843, 0.547, 0.789, 0.561, 0.814, 0.971, 0.968, 0.832, 0.828, 0.999, 0.486, 0.672, 0.812, 0.733, 0.926, 0.87, 0.652, 0.841, 0.652, 0.927, 0.898, 0.643, 0.962, 0.947, 0.755, 0.571, 0.327, 0.828, 0.19, 0.936, 0.618, 0.986, 0.803, 0.696, 0.982, 0.951, 0.473, 0.954, 0.313, 0.998, 0.665, 0.427, 0.937, 0.884, 0.72, 0.506, 0.469, 0.971, 0.889, 0.612]
global origin = 1
global destination = 60