global arcs = [1 11; 1 17; 1 20; 1 26; 1 32; 1 38; 1 40; 1 41; 1 44; 2 8; 2 17; 2 39; 2 46; 2 52; 3 12; 3 18; 3 23; 3 39; 3 48; 3 53; 3 54; 3 57; 4 5; 4 11; 4 16; 4 34; 4 37; 4 46; 5 20; 5 33; 5 35; 5 42; 5 44; 5 56; 5 58; 6 11; 6 13; 6 16; 6 22; 6 45; 6 52; 7 3; 7 22; 7 42; 7 46; 7 58; 7 59; 8 17; 8 18; 8 25; 8 28; 8 58; 8 60; 9 17; 9 19; 9 25; 9 34; 9 46; 10 7; 10 20; 10 22; 10 26; 10 44; 10 45; 10 56; 11 7; 11 9; 11 26; 11 28; 11 47; 11 49; 11 60; 12 10; 12 35; 12 49; 12 54; 13 17; 13 18; 13 24; 13 31; 13 48; 14 6; 14 13; 14 16; 14 18; 14 26; 14 55; 15 21; 15 52; 15 57; 15 60; 16 2; 16 5; 16 39; 17 4; 17 12; 17 22; 17 27; 17 32; 17 41; 18 2; 18 19; 18 47; 19 7; 19 28; 19 57; 19 59; 20 13; 20 32; 20 42; 20 43; 20 47; 20 59; 21 17; 21 18; 21 23; 21 36; 21 38; 21 42; 21 43; 22 23; 22 41; 22 42; 22 44; 23 30; 23 35; 23 43; 23 57; 24 12; 24 26; 24 48; 24 49; 24 58; 25 10; 25 12; 25 38; 25 49; 25 55; 25 60; 26 6; 26 22; 26 28; 26 37; 26 47; 26 55; 26 57; 27 5; 27 26; 27 48; 28 9; 28 36; 28 43; 28 47; 28 51; 29 15; 29 23; 29 38; 30 8; 30 19; 30 27; 30 28; 30 45; 30 56; 30 58; 30 60; 31 11; 31 21; 31 27; 31 28; 31 38; 31 60; 32 2; 32 22; 32 31; 32 45; 32 58; 33 24; 33 26; 33 28; 33 30; 33 32; 33 34; 33 39; 33 52; 33 55; 34 5; 34 8; 34 25; 34 30; 34 47; 34 49; 34 53; 35 7; 35 24; 35 27; 35 30; 35 33; 35 37; 35 53; 35 57; 35 59; 36 15; 36 41; 37 17; 37 20; 37 23; 37 27; 37 32; 38 8; 38 11; 38 17; 38 19; 38 21; 38 31; 39 43; 39 44; 39 55; 40 26; 40 52; 41 14; 41 16; 41 17; 41 20; 41 30; 41 46; 42 3; 42 5; 42 12; 42 17; 42 25; 42 30; 42 53; 43 21; 43 33; 43 40; 44 5; 44 18; 44 25; 44 28; 45 7; 45 41; 45 46; 45 51; 46 7; 46 14; 46 20; 46 23; 46 24; 46 37; 46 39; 46 40; 46 44; 46 45; 46 47; 47 19; 47 23; 47 24; 47 29; 47 33; 47 53; 48 6; 48 10; 48 30; 48 36; 48 43; 48 45; 48 47; 49 5; 49 22; 49 29; 49 33; 49 39; 49 59; 50 4; 50 27; 50 30; 51 5; 51 8; 51 11; 51 20; 51 22; 51 27; 51 35; 51 47; 51 54; 52 2; 52 30; 52 36; 52 37; 53 13; 53 21; 53 33; 53 40; 54 5; 54 17; 54 34; 54 49; 55 2; 55 20; 55 30; 55 50; 56 35; 56 40; 56 44; 56 55; 57 38; 57 49; 58 5; 58 6; 58 11; 58 21; 58 31; 59 9; 59 20; 59 31; 59 39]
global d_x = [3.0, 3.0, 4.0, 9.0, 3.0, 3.0, 1.0, 1.0, 6.0, 8.0, 9.0, 10.0, 9.0, 1.0, 4.0, 5.0, 7.0, 9.0, 9.0, 6.0, 7.0, 9.0, 8.0, 7.0, 6.0, 6.0, 6.0, 9.0, 7.0, 6.0, 10.0, 10.0, 10.0, 2.0, 6.0, 9.0, 4.0, 1.0, 4.0, 5.0, 4.0, 9.0, 6.0, 7.0, 2.0, 7.0, 8.0, 6.0, 2.0, 10.0, 6.0, 2.0, 2.0, 3.0, 9.0, 6.0, 1.0, 5.0, 8.0, 7.0, 2.0, 2.0, 7.0, 3.0, 10.0, 8.0, 10.0, 4.0, 8.0, 7.0, 6.0, 6.0, 7.0, 10.0, 6.0, 5.0, 8.0, 10.0, 1.0, 3.0, 9.0, 2.0, 8.0, 8.0, 5.0, 6.0, 9.0, 4.0, 10.0, 9.0, 3.0, 3.0, 7.0, 1.0, 2.0, 6.0, 8.0, 6.0, 10.0, 7.0, 7.0, 10.0, 5.0, 5.0, 6.0, 3.0, 6.0, 1.0, 3.0, 5.0, 5.0, 7.0, 8.0, 5.0, 9.0, 6.0, 4.0, 4.0, 3.0, 2.0, 10.0, 9.0, 4.0, 7.0, 3.0, 10.0, 7.0, 2.0, 1.0, 2.0, 7.0, 1.0, 3.0, 5.0, 7.0, 10.0, 2.0, 3.0, 7.0, 5.0, 1.0, 3.0, 4.0, 8.0, 1.0, 9.0, 9.0, 7.0, 2.0, 2.0, 8.0, 10.0, 2.0, 7.0, 8.0, 2.0, 7.0, 6.0, 6.0, 4.0, 3.0, 3.0, 1.0, 5.0, 2.0, 5.0, 4.0, 7.0, 1.0, 9.0, 4.0, 6.0, 5.0, 8.0, 5.0, 2.0, 10.0, 6.0, 8.0, 1.0, 1.0, 6.0, 10.0, 1.0, 10.0, 2.0, 3.0, 4.0, 8.0, 10.0, 7.0, 3.0, 9.0, 6.0, 4.0, 4.0, 2.0, 1.0, 7.0, 9.0, 5.0, 1.0, 7.0, 10.0, 8.0, 4.0, 6.0, 9.0, 7.0, 4.0, 2.0, 7.0, 4.0, 4.0, 1.0, 9.0, 7.0, 1.0, 9.0, 10.0, 1.0, 2.0, 6.0, 3.0, 4.0, 3.0, 7.0, 4.0, 10.0, 10.0, 6.0, 2.0, 10.0, 4.0, 8.0, 7.0, 8.0, 8.0, 2.0, 6.0, 3.0, 7.0, 3.0, 6.0, 4.0, 7.0, 4.0, 9.0, 9.0, 8.0, 9.0, 9.0, 8.0, 4.0, 10.0, 8.0, 7.0, 4.0, 8.0, 2.0, 4.0, 5.0, 5.0, 5.0, 7.0, 5.0, 8.0, 10.0, 9.0, 3.0, 4.0, 7.0, 9.0, 5.0, 7.0, 9.0, 4.0, 5.0, 10.0, 10.0, 4.0, 10.0, 7.0, 9.0, 2.0, 8.0, 3.0, 4.0, 4.0, 1.0, 6.0, 7.0, 6.0, 4.0, 4.0, 6.0, 5.0, 1.0, 8.0, 6.0, 1.0, 1.0, 3.0, 4.0, 7.0, 4.0, 1.0, 10.0, 2.0, 9.0, 7.0, 1.0, 1.0, 8.0, 9.0, 5.0]
global b_x = 5
global d_y = [4.0, 3.0, 5.0, 7.0, 9.0, 6.0, 5.0, 9.0, 7.0, 6.0, 2.0, 8.0, 7.0, 1.0, 7.0, 8.0, 7.0, 4.0, 9.0, 1.0, 4.0, 9.0, 7.0, 8.0, 8.0, 8.0, 10.0, 7.0, 9.0, 1.0, 8.0, 6.0, 9.0, 6.0, 7.0, 2.0, 3.0, 10.0, 9.0, 8.0, 5.0, 3.0, 3.0, 8.0, 1.0, 4.0, 1.0, 10.0, 5.0, 6.0, 6.0, 5.0, 10.0, 9.0, 6.0, 9.0, 6.0, 2.0, 3.0, 8.0, 8.0, 8.0, 10.0, 9.0, 7.0, 2.0, 10.0, 8.0, 1.0, 3.0, 5.0, 7.0, 3.0, 2.0, 2.0, 7.0, 7.0, 6.0, 5.0, 6.0, 3.0, 6.0, 10.0, 7.0, 9.0, 2.0, 7.0, 9.0, 5.0, 9.0, 6.0, 7.0, 2.0, 6.0, 6.0, 8.0, 6.0, 1.0, 1.0, 4.0, 9.0, 2.0, 7.0, 5.0, 8.0, 2.0, 7.0, 2.0, 4.0, 2.0, 4.0, 2.0, 3.0, 5.0, 7.0, 3.0, 9.0, 2.0, 7.0, 9.0, 7.0, 6.0, 2.0, 7.0, 8.0, 3.0, 6.0, 7.0, 2.0, 9.0, 6.0, 6.0, 6.0, 4.0, 9.0, 7.0, 4.0, 2.0, 2.0, 8.0, 7.0, 10.0, 5.0, 4.0, 7.0, 10.0, 3.0, 4.0, 3.0, 4.0, 2.0, 9.0, 6.0, 9.0, 1.0, 2.0, 8.0, 9.0, 4.0, 6.0, 7.0, 9.0, 8.0, 2.0, 5.0, 6.0, 4.0, 3.0, 1.0, 6.0, 10.0, 8.0, 2.0, 10.0, 5.0, 1.0, 1.0, 4.0, 10.0, 6.0, 2.0, 2.0, 7.0, 2.0, 3.0, 3.0, 7.0, 9.0, 9.0, 2.0, 10.0, 8.0, 7.0, 9.0, 5.0, 3.0, 8.0, 1.0, 1.0, 9.0, 3.0, 10.0, 5.0, 3.0, 4.0, 3.0, 9.0, 6.0, 3.0, 5.0, 3.0, 1.0, 10.0, 5.0, 2.0, 1.0, 4.0, 9.0, 7.0, 10.0, 5.0, 5.0, 5.0, 8.0, 3.0, 7.0, 7.0, 10.0, 1.0, 2.0, 5.0, 10.0, 8.0, 8.0, 4.0, 8.0, 5.0, 10.0, 9.0, 1.0, 5.0, 1.0, 1.0, 1.0, 4.0, 7.0, 6.0, 5.0, 1.0, 10.0, 6.0, 7.0, 7.0, 5.0, 9.0, 9.0, 8.0, 3.0, 10.0, 10.0, 4.0, 1.0, 4.0, 6.0, 2.0, 3.0, 4.0, 4.0, 1.0, 5.0, 5.0, 7.0, 10.0, 7.0, 4.0, 1.0, 10.0, 6.0, 9.0, 6.0, 10.0, 3.0, 3.0, 1.0, 6.0, 4.0, 1.0, 2.0, 5.0, 10.0, 10.0, 4.0, 8.0, 2.0, 10.0, 7.0, 6.0, 10.0, 6.0, 9.0, 2.0, 1.0, 3.0, 5.0, 2.0, 2.0, 2.0, 8.0, 1.0, 9.0, 2.0, 10.0, 3.0, 8.0, 2.0, 3.0]
global b_y = 10
global p = [0.186, 0.122, 0.967, 0.923, 0.374, 0.676, 0.211, 0.317, 0.096, 0.668, 0.3, 0.391, 0.357, 0.227, 0.945, 0.135, 0.237, 0.726, 0.615, 0.681, 0.701, 0.493, 0.639, 0.403, 0.988, 0.349, 0.386, 0.577, 0.193, 0.07, 0.495, 0.135, 0.785, 0.425, 0.124, 0.686, 0.369, 0.912, 0.705, 0.04, 0.523, 0.056, 0.445, 0.019, 0.014, 0.452, 0.675, 0.794, 0.222, 0.449, 0.125, 0.537, 0.952, 0.332, 0.64, 0.991, 0.527, 0.149, 0.035, 0.823, 0.155, 0.517, 0.986, 0.165, 0.341, 0.917, 0.349, 0.007, 0.49, 0.082, 0.037, 0.383, 0.622, 0.422, 0.602, 0.107, 0.491, 0.08, 0.198, 0.76, 0.77, 0.395, 0.289, 0.233, 0.515, 0.059, 0.715, 0.243, 0.025, 0.405, 0.322, 0.64, 0.477, 0.388, 0.208, 0.96, 0.377, 0.019, 0.476, 0.461, 0.976, 0.352, 0.693, 0.034, 0.583, 0.279, 0.941, 0.981, 0.283, 0.745, 0.943, 0.367, 0.698, 0.659, 0.947, 0.066, 0.34, 0.02, 0.85, 0.164, 0.047, 0.917, 0.897, 0.843, 0.385, 0.574, 0.215, 0.847, 0.077, 0.974, 0.818, 0.185, 0.912, 0.751, 0.61, 0.488, 0.98, 0.118, 0.641, 0.425, 0.548, 0.52, 0.308, 0.488, 0.979, 0.62, 0.963, 0.866, 0.712, 0.016, 0.349, 0.333, 0.683, 0.179, 0.546, 0.924, 0.381, 0.811, 0.495, 0.877, 0.428, 0.056, 0.247, 0.268, 0.982, 0.53, 0.437, 0.62, 0.168, 0.997, 0.116, 0.768, 0.064, 0.458, 0.9, 0.693, 0.816, 0.367, 0.683, 0.445, 0.093, 0.45, 0.854, 0.019, 0.818, 0.06, 0.916, 0.644, 0.585, 0.999, 0.491, 0.023, 0.188, 0.153, 0.546, 0.418, 0.057, 0.45, 0.099, 0.123, 0.18, 0.202, 0.26, 0.846, 0.292, 0.86, 0.917, 0.041, 0.085, 0.101, 0.963, 0.125, 0.475, 0.683, 0.596, 0.566, 0.147, 0.327, 0.905, 0.628, 0.342, 0.315, 0.414, 0.441, 0.918, 0.003, 0.66, 0.587, 0.887, 0.103, 0.242, 0.279, 0.768, 0.969, 0.426, 0.768, 0.608, 0.454, 0.753, 0.75, 0.883, 0.07, 0.43, 0.178, 0.546, 0.284, 0.544, 0.432, 0.228, 0.504, 0.127, 0.554, 0.076, 0.205, 0.638, 0.437, 0.249, 0.39, 0.293, 0.998, 0.232, 0.372, 0.532, 0.678, 0.558, 0.971, 0.297, 0.041, 0.782, 0.821, 0.982, 0.422, 0.74, 0.915, 0.865, 0.801, 0.597, 0.551, 0.939, 0.849, 0.51, 0.003, 0.78, 0.723, 0.451, 0.819, 0.206, 0.585, 0.094, 0.981, 0.084, 0.456, 0.894, 0.622, 0.095, 0.529, 0.074, 0.477, 0.899, 0.076, 0.259, 0.839, 0.718, 0.704, 0.196, 0.925, 0.137, 0.897, 0.432, 0.514, 0.63, 0.062, 0.341, 0.362, 0.789, 0.185]
global q = [0.868, 0.601, 0.995, 0.94, 0.711, 0.994, 0.567, 0.822, 0.341, 0.778, 0.723, 0.43, 0.685, 0.956, 0.956, 0.873, 0.422, 0.792, 0.867, 0.955, 0.758, 0.965, 0.944, 0.742, 0.991, 0.66, 0.968, 0.625, 0.745, 0.865, 0.685, 0.898, 0.855, 0.651, 0.643, 0.734, 0.44, 0.987, 0.858, 0.068, 0.925, 0.421, 0.627, 0.667, 0.73, 0.592, 0.855, 0.868, 0.548, 0.92, 0.977, 0.788, 0.952, 0.677, 0.894, 0.998, 0.603, 0.447, 0.649, 0.866, 0.244, 0.635, 0.992, 0.427, 0.769, 0.917, 0.584, 0.025, 0.581, 0.275, 0.622, 0.622, 0.999, 0.519, 0.803, 0.644, 0.565, 0.87, 0.713, 0.928, 0.846, 0.487, 0.49, 0.594, 0.84, 0.256, 0.995, 0.687, 0.323, 0.986, 0.933, 0.678, 0.822, 0.576, 0.79, 0.983, 0.609, 0.554, 0.528, 0.724, 0.976, 0.98, 0.891, 0.221, 0.851, 0.372, 0.986, 0.998, 0.693, 0.968, 0.967, 0.42, 0.757, 0.775, 0.979, 0.535, 0.61, 0.229, 0.943, 0.719, 0.205, 0.939, 0.938, 0.924, 0.532, 0.894, 0.802, 0.956, 0.833, 0.992, 0.91, 0.792, 0.969, 0.928, 0.879, 0.693, 0.992, 0.303, 0.903, 0.513, 0.864, 0.944, 0.639, 0.55, 0.999, 0.869, 0.996, 0.938, 0.832, 0.855, 0.924, 0.88, 0.951, 0.809, 0.976, 0.994, 0.711, 0.837, 0.56, 0.942, 0.782, 0.782, 0.649, 0.344, 0.986, 0.746, 0.481, 0.861, 0.33, 0.999, 0.36, 0.782, 0.883, 0.589, 0.962, 0.81, 0.906, 0.671, 0.857, 0.575, 0.155, 0.717, 0.954, 0.894, 0.937, 0.846, 0.918, 0.687, 0.744, 0.999, 0.564, 0.648, 0.95, 0.434, 0.555, 0.876, 0.823, 0.727, 0.929, 0.476, 0.771, 0.306, 0.283, 0.915, 0.991, 0.951, 0.928, 0.505, 0.551, 0.534, 0.995, 0.288, 0.602, 0.957, 0.878, 0.917, 0.812, 0.526, 0.93, 0.741, 0.344, 0.427, 0.936, 0.807, 0.971, 0.409, 0.924, 0.865, 0.98, 0.378, 0.762, 0.918, 0.921, 0.982, 0.507, 0.931, 0.759, 0.983, 0.88, 0.925, 0.988, 0.816, 0.672, 0.808, 0.829, 0.865, 0.615, 0.483, 0.74, 0.505, 0.828, 0.674, 0.531, 0.595, 0.746, 0.674, 0.647, 0.717, 0.966, 0.999, 0.896, 0.395, 0.932, 0.743, 0.606, 0.994, 0.309, 0.333, 0.928, 0.886, 0.983, 0.52, 0.742, 0.992, 0.99, 0.9, 0.863, 0.683, 0.942, 0.857, 0.736, 0.501, 0.877, 0.941, 0.932, 0.923, 0.419, 0.659, 0.316, 0.985, 0.255, 0.596, 0.969, 0.667, 0.205, 0.951, 0.99, 0.895, 0.998, 0.398, 0.835, 0.877, 0.871, 0.784, 0.239, 0.982, 0.467, 0.907, 0.651, 0.779, 0.744, 0.875, 0.867, 0.526, 0.868, 0.308]
global origin = 1
global destination = 60