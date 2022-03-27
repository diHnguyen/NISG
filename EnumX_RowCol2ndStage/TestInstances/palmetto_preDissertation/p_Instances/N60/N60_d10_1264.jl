global arcs = [1 2; 1 4; 1 32; 1 49; 2 3; 2 14; 2 20; 2 37; 2 42; 2 45; 2 60; 3 2; 3 20; 3 24; 3 42; 3 48; 3 50; 3 55; 4 16; 4 40; 4 45; 4 50; 5 8; 5 30; 5 40; 5 59; 6 10; 6 14; 6 28; 6 34; 6 47; 6 57; 7 11; 7 12; 7 13; 7 21; 7 23; 7 46; 7 56; 7 57; 8 3; 8 10; 8 23; 8 24; 8 30; 8 53; 9 11; 9 53; 9 54; 9 55; 10 2; 10 18; 10 28; 10 33; 10 35; 10 40; 10 44; 10 60; 11 6; 11 14; 11 24; 11 35; 11 42; 11 56; 11 58; 12 22; 12 37; 12 48; 12 55; 12 58; 13 43; 13 46; 13 56; 14 5; 14 9; 14 20; 14 35; 14 36; 14 60; 15 4; 15 16; 15 28; 15 41; 15 51; 15 59; 16 21; 16 25; 16 30; 16 49; 16 54; 16 58; 17 9; 17 16; 17 44; 17 57; 17 59; 18 19; 18 42; 19 15; 19 18; 19 32; 19 47; 20 5; 20 11; 20 13; 20 21; 20 33; 20 45; 21 11; 21 23; 21 37; 21 41; 21 60; 22 15; 22 19; 22 27; 22 37; 22 47; 22 54; 23 18; 23 33; 23 42; 23 52; 24 3; 24 12; 24 22; 24 30; 24 44; 24 59; 25 12; 25 35; 25 37; 25 46; 25 60; 26 15; 26 39; 26 49; 27 16; 27 19; 27 21; 27 25; 27 28; 27 33; 27 46; 27 48; 27 52; 28 39; 28 42; 29 9; 29 46; 29 47; 30 13; 30 15; 30 25; 30 46; 30 47; 30 48; 30 50; 30 58; 31 2; 31 7; 31 14; 31 27; 31 47; 31 52; 32 22; 32 33; 32 35; 32 58; 33 3; 33 19; 33 21; 33 36; 34 11; 34 15; 34 22; 34 30; 34 45; 34 52; 34 53; 35 13; 35 18; 35 59; 36 7; 36 22; 36 26; 36 44; 37 4; 37 6; 37 20; 37 40; 37 43; 37 53; 37 60; 38 6; 38 14; 38 21; 38 29; 38 36; 38 49; 39 5; 39 16; 39 56; 40 2; 40 34; 40 53; 40 54; 41 2; 41 7; 41 17; 41 22; 41 28; 41 39; 41 40; 41 43; 41 52; 41 54; 42 37; 42 41; 43 9; 43 19; 43 23; 43 29; 43 34; 44 8; 44 9; 44 22; 44 27; 44 29; 44 30; 44 37; 44 45; 44 51; 45 2; 45 19; 45 24; 45 27; 45 57; 46 3; 46 14; 46 17; 46 21; 46 34; 46 47; 47 4; 47 9; 47 21; 47 34; 47 38; 47 41; 47 44; 47 53; 48 5; 48 17; 48 19; 48 24; 48 29; 48 44; 48 49; 48 53; 49 6; 49 11; 49 21; 49 26; 49 30; 50 27; 50 30; 50 38; 50 53; 50 57; 51 15; 51 16; 51 40; 52 23; 52 34; 52 35; 52 45; 52 56; 52 57; 53 3; 53 17; 53 30; 53 35; 53 48; 54 9; 54 11; 54 13; 54 17; 54 19; 54 21; 54 22; 54 25; 54 31; 54 36; 54 46; 54 47; 55 15; 55 20; 55 23; 55 32; 55 46; 55 48; 56 59; 57 3; 57 58; 58 8; 58 10; 58 14; 58 17; 58 21; 58 36; 58 55; 58 57; 59 13; 59 34; 59 44]
global d_x = [3.0, 4.0, 10.0, 4.0, 8.0, 8.0, 2.0, 10.0, 6.0, 7.0, 5.0, 10.0, 8.0, 10.0, 10.0, 8.0, 1.0, 4.0, 7.0, 2.0, 8.0, 4.0, 5.0, 5.0, 9.0, 1.0, 8.0, 1.0, 7.0, 4.0, 9.0, 9.0, 2.0, 3.0, 8.0, 6.0, 10.0, 7.0, 1.0, 5.0, 2.0, 2.0, 3.0, 5.0, 9.0, 1.0, 6.0, 4.0, 9.0, 9.0, 3.0, 7.0, 2.0, 3.0, 5.0, 10.0, 8.0, 4.0, 8.0, 1.0, 9.0, 5.0, 5.0, 5.0, 2.0, 2.0, 9.0, 10.0, 3.0, 10.0, 8.0, 7.0, 6.0, 4.0, 1.0, 10.0, 3.0, 6.0, 10.0, 10.0, 4.0, 4.0, 1.0, 1.0, 2.0, 8.0, 4.0, 6.0, 7.0, 1.0, 1.0, 10.0, 5.0, 10.0, 7.0, 7.0, 9.0, 8.0, 1.0, 4.0, 10.0, 3.0, 3.0, 7.0, 3.0, 2.0, 3.0, 2.0, 7.0, 6.0, 3.0, 5.0, 2.0, 4.0, 8.0, 1.0, 8.0, 1.0, 6.0, 6.0, 6.0, 8.0, 10.0, 6.0, 8.0, 4.0, 5.0, 10.0, 2.0, 10.0, 10.0, 7.0, 2.0, 6.0, 10.0, 3.0, 6.0, 8.0, 3.0, 6.0, 5.0, 6.0, 3.0, 5.0, 8.0, 4.0, 7.0, 5.0, 7.0, 9.0, 6.0, 2.0, 1.0, 7.0, 2.0, 4.0, 5.0, 1.0, 3.0, 7.0, 7.0, 9.0, 9.0, 9.0, 1.0, 6.0, 10.0, 5.0, 3.0, 10.0, 10.0, 2.0, 8.0, 9.0, 10.0, 2.0, 4.0, 1.0, 3.0, 10.0, 4.0, 8.0, 4.0, 6.0, 3.0, 5.0, 4.0, 3.0, 2.0, 2.0, 4.0, 8.0, 10.0, 10.0, 4.0, 8.0, 8.0, 4.0, 3.0, 3.0, 7.0, 3.0, 2.0, 6.0, 1.0, 2.0, 7.0, 2.0, 7.0, 8.0, 3.0, 3.0, 7.0, 6.0, 1.0, 3.0, 10.0, 2.0, 2.0, 8.0, 5.0, 9.0, 6.0, 3.0, 7.0, 8.0, 10.0, 4.0, 4.0, 5.0, 1.0, 5.0, 9.0, 2.0, 7.0, 1.0, 8.0, 2.0, 9.0, 1.0, 1.0, 9.0, 5.0, 8.0, 5.0, 6.0, 2.0, 7.0, 9.0, 9.0, 6.0, 3.0, 5.0, 10.0, 6.0, 10.0, 3.0, 10.0, 8.0, 1.0, 10.0, 9.0, 8.0, 8.0, 10.0, 4.0, 6.0, 10.0, 1.0, 7.0, 8.0, 3.0, 10.0, 1.0, 10.0, 5.0, 9.0, 3.0, 2.0, 2.0, 6.0, 3.0, 1.0, 10.0, 4.0, 1.0, 6.0, 10.0, 1.0, 10.0, 9.0, 4.0, 4.0, 5.0, 9.0, 4.0, 4.0, 9.0, 9.0, 3.0, 8.0, 1.0, 8.0, 2.0, 8.0, 4.0, 1.0, 6.0, 5.0, 9.0, 5.0, 6.0, 2.0, 8.0, 10.0, 6.0]
global b_x = 5
global d_y = [5.0, 7.0, 8.0, 7.0, 5.0, 9.0, 1.0, 9.0, 5.0, 7.0, 1.0, 9.0, 10.0, 5.0, 4.0, 4.0, 3.0, 7.0, 4.0, 1.0, 5.0, 7.0, 5.0, 9.0, 6.0, 4.0, 10.0, 4.0, 3.0, 7.0, 4.0, 6.0, 7.0, 8.0, 6.0, 9.0, 3.0, 7.0, 2.0, 4.0, 10.0, 8.0, 10.0, 10.0, 4.0, 9.0, 1.0, 7.0, 3.0, 8.0, 3.0, 9.0, 1.0, 10.0, 6.0, 7.0, 9.0, 9.0, 4.0, 10.0, 7.0, 3.0, 8.0, 2.0, 5.0, 1.0, 3.0, 4.0, 6.0, 7.0, 1.0, 9.0, 4.0, 7.0, 7.0, 7.0, 2.0, 1.0, 9.0, 7.0, 2.0, 10.0, 2.0, 3.0, 8.0, 7.0, 7.0, 3.0, 10.0, 9.0, 8.0, 2.0, 5.0, 2.0, 9.0, 10.0, 9.0, 1.0, 7.0, 4.0, 10.0, 9.0, 4.0, 6.0, 7.0, 9.0, 3.0, 2.0, 7.0, 7.0, 7.0, 9.0, 7.0, 4.0, 2.0, 9.0, 4.0, 7.0, 5.0, 3.0, 5.0, 8.0, 3.0, 9.0, 3.0, 5.0, 10.0, 8.0, 4.0, 6.0, 7.0, 10.0, 3.0, 2.0, 1.0, 5.0, 6.0, 7.0, 2.0, 8.0, 7.0, 1.0, 2.0, 2.0, 8.0, 10.0, 6.0, 3.0, 6.0, 7.0, 9.0, 6.0, 4.0, 2.0, 1.0, 10.0, 9.0, 4.0, 7.0, 4.0, 10.0, 9.0, 9.0, 7.0, 9.0, 7.0, 1.0, 2.0, 5.0, 7.0, 4.0, 7.0, 7.0, 4.0, 2.0, 6.0, 5.0, 3.0, 1.0, 8.0, 4.0, 6.0, 10.0, 4.0, 5.0, 10.0, 9.0, 2.0, 5.0, 9.0, 9.0, 5.0, 4.0, 2.0, 10.0, 10.0, 2.0, 7.0, 7.0, 8.0, 2.0, 8.0, 3.0, 3.0, 2.0, 4.0, 5.0, 8.0, 2.0, 9.0, 1.0, 10.0, 5.0, 5.0, 9.0, 6.0, 8.0, 6.0, 10.0, 7.0, 10.0, 3.0, 7.0, 2.0, 1.0, 7.0, 4.0, 8.0, 6.0, 9.0, 9.0, 10.0, 8.0, 8.0, 5.0, 10.0, 9.0, 9.0, 3.0, 6.0, 2.0, 9.0, 1.0, 5.0, 1.0, 4.0, 4.0, 6.0, 3.0, 4.0, 9.0, 6.0, 10.0, 7.0, 7.0, 6.0, 8.0, 7.0, 6.0, 10.0, 2.0, 4.0, 4.0, 10.0, 5.0, 6.0, 9.0, 1.0, 10.0, 7.0, 3.0, 10.0, 9.0, 4.0, 1.0, 5.0, 2.0, 2.0, 4.0, 4.0, 10.0, 5.0, 7.0, 6.0, 10.0, 6.0, 4.0, 4.0, 4.0, 5.0, 5.0, 4.0, 4.0, 2.0, 1.0, 8.0, 9.0, 6.0, 8.0, 4.0, 8.0, 2.0, 7.0, 3.0, 5.0, 3.0, 1.0, 5.0, 8.0, 2.0, 1.0, 8.0, 2.0, 8.0, 1.0, 2.0]
global b_y = 10
global p = [0.491, 0.475, 0.997, 0.976, 0.035, 0.534, 0.323, 0.994, 0.391, 0.166, 0.442, 0.05, 0.904, 0.474, 0.129, 0.336, 0.985, 0.258, 0.285, 0.442, 0.357, 0.601, 0.686, 0.222, 0.857, 0.888, 0.854, 0.945, 0.624, 0.735, 0.908, 0.759, 0.171, 0.095, 0.625, 0.37, 0.388, 0.624, 0.266, 0.309, 0.259, 0.672, 0.468, 0.514, 0.341, 0.135, 0.329, 0.157, 0.355, 0.69, 0.262, 0.915, 0.667, 0.848, 0.239, 0.33, 0.843, 0.172, 0.383, 0.151, 0.305, 0.859, 0.88, 0.372, 0.02, 0.451, 0.781, 0.04, 0.957, 0.323, 0.821, 0.619, 0.067, 0.537, 0.791, 0.71, 0.106, 0.367, 0.658, 0.666, 0.594, 0.607, 0.234, 0.581, 0.266, 0.267, 0.022, 0.762, 0.309, 0.774, 0.356, 0.225, 0.914, 0.329, 0.812, 0.755, 0.553, 0.723, 0.012, 0.422, 0.11, 0.368, 0.87, 0.342, 0.221, 0.191, 0.221, 0.824, 0.143, 0.548, 0.756, 0.889, 0.216, 0.78, 0.465, 0.936, 0.645, 0.198, 0.536, 0.062, 0.906, 0.292, 0.454, 0.78, 0.843, 0.299, 0.292, 0.568, 0.231, 0.922, 0.111, 0.77, 0.015, 0.793, 0.605, 0.572, 0.501, 0.826, 0.611, 0.13, 0.329, 0.979, 0.708, 0.842, 0.048, 0.024, 0.647, 0.517, 0.634, 0.615, 0.561, 0.522, 0.085, 0.669, 0.291, 0.455, 0.451, 0.501, 0.717, 0.955, 0.975, 0.343, 0.782, 0.897, 0.425, 0.61, 0.898, 0.529, 0.314, 0.268, 0.29, 0.244, 0.285, 0.999, 0.135, 0.441, 0.02, 0.744, 0.931, 0.843, 0.355, 0.71, 0.076, 0.648, 0.78, 0.228, 0.549, 0.657, 0.378, 0.99, 0.245, 0.537, 0.571, 0.148, 0.569, 0.4, 0.07, 0.129, 0.813, 0.808, 0.31, 0.274, 0.782, 0.738, 0.86, 0.957, 0.297, 0.374, 0.957, 0.034, 0.996, 0.378, 0.849, 0.451, 0.348, 0.635, 0.116, 0.038, 0.936, 0.822, 0.854, 0.776, 0.786, 0.459, 0.698, 0.591, 0.027, 0.43, 0.48, 0.954, 0.361, 0.303, 0.758, 0.073, 0.912, 0.428, 0.971, 0.124, 0.509, 0.446, 0.193, 0.611, 0.811, 0.467, 0.681, 0.783, 0.056, 0.314, 0.741, 0.73, 0.926, 0.597, 0.186, 0.576, 0.75, 0.351, 0.165, 0.979, 0.397, 0.812, 0.409, 0.86, 0.429, 0.18, 0.859, 0.311, 0.384, 0.21, 0.904, 0.392, 0.967, 0.29, 0.386, 0.74, 0.137, 0.557, 0.613, 0.365, 0.54, 0.614, 0.463, 0.636, 0.889, 0.712, 0.161, 0.389, 0.281, 0.147, 0.746, 0.737, 0.544, 0.953, 0.548, 0.549, 0.179, 0.607, 0.348, 0.227, 0.873, 0.994, 0.916, 0.238, 0.521, 0.666, 0.245, 0.951, 0.139, 0.848, 0.358, 0.581, 0.912, 0.304, 0.465, 0.866, 0.112, 0.985]
global q = [0.697, 0.706, 0.997, 0.979, 0.188, 0.645, 0.548, 0.997, 0.877, 0.562, 0.601, 0.955, 0.931, 0.921, 0.698, 0.974, 0.995, 0.796, 0.594, 0.752, 0.996, 0.967, 0.844, 0.41, 0.99, 0.928, 0.878, 0.998, 0.652, 0.739, 0.984, 0.781, 0.394, 0.48, 0.81, 0.853, 0.944, 0.946, 0.606, 0.997, 0.807, 0.883, 0.611, 0.817, 0.4, 0.648, 0.484, 0.698, 0.412, 0.96, 0.591, 0.989, 0.974, 0.919, 0.722, 0.959, 0.904, 0.532, 0.387, 0.964, 0.957, 0.886, 0.951, 0.832, 0.162, 0.561, 0.881, 0.935, 0.978, 0.768, 0.986, 0.887, 0.499, 0.731, 0.991, 0.937, 0.634, 0.575, 0.713, 0.897, 0.977, 0.778, 0.893, 0.745, 0.404, 0.888, 0.243, 0.973, 0.812, 0.919, 0.673, 0.857, 0.988, 0.915, 0.826, 0.87, 0.592, 0.82, 0.679, 0.446, 0.727, 0.459, 0.876, 0.955, 0.289, 0.869, 0.916, 0.915, 0.623, 0.869, 0.93, 0.907, 0.798, 0.821, 0.768, 0.976, 0.933, 0.731, 0.644, 0.557, 0.962, 0.89, 0.527, 0.808, 0.961, 0.971, 0.443, 0.629, 0.582, 0.951, 0.144, 0.793, 0.422, 0.973, 0.932, 0.626, 0.879, 0.917, 0.946, 0.906, 0.821, 0.996, 0.95, 0.947, 0.566, 0.784, 0.874, 0.613, 0.72, 0.67, 0.602, 0.802, 0.732, 0.692, 0.503, 0.935, 0.803, 0.706, 0.824, 0.977, 0.993, 0.45, 0.817, 0.976, 0.684, 0.81, 0.956, 0.994, 0.917, 0.736, 0.983, 0.928, 0.995, 0.999, 0.653, 0.609, 0.537, 0.936, 0.974, 0.88, 0.511, 0.76, 0.465, 0.865, 0.836, 0.308, 0.635, 0.727, 0.932, 0.991, 0.91, 0.842, 0.804, 0.197, 0.817, 0.713, 0.324, 0.797, 0.827, 0.847, 0.978, 0.882, 0.862, 0.761, 0.989, 0.982, 0.938, 0.961, 0.99, 0.106, 0.997, 0.613, 0.905, 0.671, 0.574, 0.866, 0.121, 0.979, 0.991, 0.884, 0.948, 0.95, 0.869, 0.685, 0.799, 0.71, 0.697, 0.637, 0.914, 0.954, 0.426, 0.441, 0.964, 0.623, 0.94, 0.428, 0.988, 0.743, 0.772, 0.512, 0.541, 0.871, 0.965, 0.612, 0.948, 0.883, 0.592, 0.806, 0.786, 0.845, 0.936, 0.762, 0.504, 0.807, 0.783, 0.826, 0.722, 0.994, 0.972, 0.888, 0.609, 0.877, 0.782, 0.826, 0.893, 0.839, 0.878, 0.468, 0.953, 0.44, 0.972, 0.559, 0.955, 0.752, 0.795, 0.959, 0.907, 0.492, 0.575, 0.772, 0.52, 0.886, 0.913, 0.974, 0.516, 0.71, 0.454, 0.796, 0.853, 0.887, 0.831, 0.956, 0.887, 0.91, 0.867, 0.672, 0.921, 0.532, 0.927, 0.999, 0.935, 0.992, 0.896, 0.829, 0.353, 0.98, 0.751, 0.924, 0.992, 0.934, 0.956, 0.416, 0.941, 0.932, 0.296, 0.999]
global origin = 1
global destination = 60