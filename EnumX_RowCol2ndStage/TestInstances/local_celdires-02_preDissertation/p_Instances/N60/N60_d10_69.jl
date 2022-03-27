global arcs = [1 6; 1 18; 1 42; 1 52; 2 33; 2 35; 2 37; 3 13; 3 24; 3 54; 4 20; 4 34; 4 46; 4 47; 5 3; 5 26; 5 34; 6 20; 6 32; 6 44; 6 47; 6 53; 6 54; 7 5; 7 6; 7 23; 7 48; 8 16; 8 19; 8 20; 8 22; 8 34; 8 40; 8 45; 8 46; 8 55; 9 7; 9 10; 9 21; 9 33; 9 48; 9 54; 10 12; 10 23; 10 24; 10 37; 10 52; 11 7; 11 8; 11 9; 11 10; 11 18; 11 23; 11 31; 11 37; 11 60; 12 9; 12 27; 12 29; 12 40; 12 42; 12 45; 13 2; 13 25; 13 49; 13 58; 14 13; 14 20; 14 29; 14 37; 14 45; 14 48; 14 49; 15 20; 15 33; 15 37; 15 47; 15 49; 16 8; 16 46; 16 55; 17 26; 17 30; 17 39; 17 50; 18 6; 18 14; 18 19; 18 51; 18 55; 19 13; 19 30; 19 31; 19 47; 19 55; 20 2; 20 10; 20 16; 20 34; 20 35; 20 36; 20 37; 21 14; 21 17; 21 19; 21 32; 21 55; 21 59; 22 4; 22 25; 23 11; 23 36; 23 56; 23 58; 24 16; 24 23; 25 19; 25 44; 25 48; 25 56; 26 7; 26 13; 26 17; 26 36; 26 39; 26 40; 26 44; 26 49; 27 9; 27 26; 28 13; 28 18; 28 19; 28 24; 28 39; 28 60; 29 33; 29 50; 29 57; 30 3; 30 15; 30 35; 30 46; 30 54; 30 60; 31 6; 31 16; 31 17; 31 26; 31 55; 32 9; 32 15; 32 27; 32 36; 32 52; 32 53; 32 55; 33 14; 33 25; 33 41; 33 44; 33 45; 33 48; 33 49; 33 56; 34 2; 34 15; 34 16; 34 22; 34 28; 34 33; 34 49; 34 52; 35 19; 35 24; 35 56; 36 19; 36 34; 37 5; 37 12; 37 25; 37 29; 37 30; 37 33; 37 43; 37 46; 38 3; 38 8; 38 9; 38 16; 38 24; 38 25; 38 27; 38 33; 38 41; 38 48; 38 50; 38 52; 38 57; 39 7; 39 26; 39 29; 39 59; 40 10; 40 14; 40 24; 40 29; 40 34; 40 38; 40 42; 40 51; 41 7; 41 20; 41 25; 41 28; 41 31; 41 50; 42 5; 42 30; 42 34; 42 38; 42 46; 42 52; 42 58; 43 7; 43 15; 43 25; 43 26; 43 35; 43 37; 43 54; 44 7; 44 9; 44 10; 44 23; 44 32; 44 47; 44 49; 45 9; 45 12; 45 23; 45 24; 45 43; 46 8; 46 30; 46 35; 46 36; 46 56; 46 57; 47 11; 47 26; 47 31; 47 32; 47 35; 47 56; 47 57; 47 60; 48 10; 48 21; 48 25; 48 28; 48 38; 49 8; 49 12; 49 19; 49 35; 49 40; 49 44; 49 60; 50 37; 50 38; 50 45; 50 47; 51 5; 51 8; 51 10; 51 12; 51 18; 51 39; 52 4; 52 15; 52 26; 52 27; 52 38; 52 40; 53 38; 53 40; 54 33; 54 37; 54 53; 54 56; 55 10; 55 22; 55 40; 55 44; 56 9; 56 20; 56 26; 56 27; 56 35; 56 38; 56 42; 56 50; 56 54; 56 57; 56 59; 57 12; 57 13; 57 23; 57 51; 58 12; 58 15; 58 25; 58 27; 58 30; 58 37; 58 47; 59 23; 59 38; 59 39; 59 50; 59 54]
global d_x = [9.0, 4.0, 9.0, 7.0, 9.0, 2.0, 10.0, 4.0, 2.0, 4.0, 5.0, 8.0, 10.0, 10.0, 5.0, 9.0, 6.0, 3.0, 1.0, 7.0, 7.0, 1.0, 9.0, 1.0, 8.0, 7.0, 8.0, 4.0, 2.0, 6.0, 3.0, 7.0, 3.0, 9.0, 4.0, 7.0, 6.0, 2.0, 5.0, 1.0, 3.0, 6.0, 2.0, 10.0, 10.0, 5.0, 6.0, 8.0, 8.0, 10.0, 6.0, 5.0, 7.0, 1.0, 5.0, 9.0, 8.0, 7.0, 6.0, 5.0, 6.0, 6.0, 7.0, 5.0, 7.0, 6.0, 3.0, 7.0, 8.0, 7.0, 6.0, 2.0, 9.0, 2.0, 4.0, 9.0, 3.0, 4.0, 1.0, 4.0, 7.0, 8.0, 2.0, 4.0, 2.0, 9.0, 10.0, 2.0, 2.0, 7.0, 9.0, 10.0, 5.0, 5.0, 9.0, 6.0, 9.0, 8.0, 2.0, 10.0, 3.0, 8.0, 5.0, 4.0, 2.0, 7.0, 10.0, 2.0, 10.0, 4.0, 8.0, 8.0, 7.0, 10.0, 9.0, 4.0, 3.0, 5.0, 8.0, 4.0, 7.0, 10.0, 8.0, 5.0, 6.0, 8.0, 6.0, 2.0, 7.0, 2.0, 10.0, 6.0, 6.0, 10.0, 2.0, 1.0, 10.0, 9.0, 4.0, 4.0, 2.0, 10.0, 3.0, 5.0, 1.0, 3.0, 1.0, 2.0, 9.0, 2.0, 2.0, 6.0, 2.0, 5.0, 10.0, 10.0, 4.0, 2.0, 2.0, 6.0, 8.0, 8.0, 5.0, 10.0, 9.0, 10.0, 8.0, 1.0, 10.0, 4.0, 10.0, 10.0, 8.0, 10.0, 10.0, 2.0, 9.0, 6.0, 9.0, 5.0, 9.0, 7.0, 3.0, 10.0, 8.0, 2.0, 7.0, 5.0, 7.0, 6.0, 10.0, 8.0, 10.0, 10.0, 9.0, 3.0, 2.0, 1.0, 1.0, 6.0, 2.0, 2.0, 8.0, 6.0, 7.0, 2.0, 8.0, 2.0, 8.0, 8.0, 6.0, 5.0, 9.0, 1.0, 10.0, 4.0, 2.0, 3.0, 4.0, 9.0, 6.0, 4.0, 8.0, 5.0, 2.0, 3.0, 8.0, 4.0, 10.0, 10.0, 7.0, 3.0, 2.0, 2.0, 8.0, 2.0, 1.0, 7.0, 10.0, 2.0, 4.0, 6.0, 8.0, 8.0, 7.0, 9.0, 8.0, 9.0, 8.0, 4.0, 6.0, 6.0, 5.0, 9.0, 1.0, 8.0, 5.0, 1.0, 7.0, 3.0, 5.0, 4.0, 10.0, 5.0, 5.0, 5.0, 7.0, 7.0, 6.0, 8.0, 1.0, 3.0, 7.0, 2.0, 8.0, 4.0, 2.0, 2.0, 5.0, 8.0, 1.0, 5.0, 9.0, 10.0, 5.0, 1.0, 7.0, 5.0, 3.0, 4.0, 9.0, 5.0, 8.0, 2.0, 4.0, 10.0, 6.0, 1.0, 3.0, 2.0, 6.0, 5.0, 10.0, 4.0, 5.0, 3.0, 8.0, 8.0, 3.0, 7.0, 6.0, 9.0, 8.0, 10.0, 5.0, 4.0, 7.0, 8.0, 9.0, 10.0, 7.0, 3.0]
global b_x = 5
global d_y = [2.0, 3.0, 6.0, 7.0, 2.0, 7.0, 1.0, 6.0, 8.0, 4.0, 8.0, 1.0, 7.0, 2.0, 7.0, 5.0, 6.0, 4.0, 6.0, 5.0, 6.0, 2.0, 8.0, 5.0, 5.0, 9.0, 7.0, 10.0, 4.0, 7.0, 10.0, 7.0, 8.0, 10.0, 3.0, 7.0, 3.0, 3.0, 9.0, 10.0, 4.0, 1.0, 9.0, 1.0, 6.0, 10.0, 6.0, 9.0, 4.0, 3.0, 7.0, 10.0, 4.0, 5.0, 4.0, 6.0, 9.0, 10.0, 4.0, 4.0, 9.0, 2.0, 7.0, 10.0, 1.0, 7.0, 4.0, 5.0, 8.0, 9.0, 1.0, 9.0, 6.0, 10.0, 9.0, 1.0, 9.0, 7.0, 1.0, 5.0, 9.0, 4.0, 2.0, 3.0, 3.0, 6.0, 2.0, 2.0, 7.0, 4.0, 4.0, 9.0, 10.0, 1.0, 1.0, 2.0, 8.0, 10.0, 7.0, 5.0, 4.0, 4.0, 4.0, 7.0, 9.0, 2.0, 5.0, 7.0, 8.0, 7.0, 2.0, 2.0, 8.0, 5.0, 9.0, 6.0, 3.0, 6.0, 9.0, 9.0, 3.0, 1.0, 4.0, 1.0, 6.0, 1.0, 7.0, 2.0, 8.0, 8.0, 10.0, 10.0, 4.0, 7.0, 4.0, 8.0, 8.0, 9.0, 10.0, 8.0, 3.0, 1.0, 4.0, 10.0, 1.0, 3.0, 5.0, 9.0, 2.0, 5.0, 6.0, 4.0, 7.0, 4.0, 5.0, 1.0, 3.0, 7.0, 4.0, 1.0, 4.0, 9.0, 6.0, 5.0, 8.0, 5.0, 3.0, 8.0, 4.0, 1.0, 7.0, 5.0, 8.0, 10.0, 8.0, 9.0, 4.0, 8.0, 2.0, 4.0, 7.0, 4.0, 5.0, 7.0, 10.0, 9.0, 6.0, 1.0, 8.0, 2.0, 10.0, 8.0, 4.0, 8.0, 5.0, 2.0, 10.0, 8.0, 9.0, 3.0, 8.0, 3.0, 3.0, 5.0, 1.0, 6.0, 6.0, 10.0, 8.0, 10.0, 8.0, 3.0, 10.0, 3.0, 2.0, 1.0, 6.0, 8.0, 7.0, 8.0, 2.0, 10.0, 3.0, 3.0, 3.0, 6.0, 1.0, 1.0, 7.0, 4.0, 6.0, 3.0, 1.0, 7.0, 8.0, 3.0, 8.0, 6.0, 5.0, 1.0, 10.0, 3.0, 7.0, 5.0, 5.0, 8.0, 5.0, 2.0, 8.0, 7.0, 10.0, 1.0, 5.0, 5.0, 1.0, 2.0, 1.0, 1.0, 5.0, 10.0, 7.0, 3.0, 9.0, 9.0, 10.0, 8.0, 2.0, 8.0, 6.0, 4.0, 8.0, 8.0, 2.0, 9.0, 4.0, 2.0, 9.0, 6.0, 4.0, 8.0, 10.0, 1.0, 9.0, 2.0, 2.0, 3.0, 5.0, 9.0, 2.0, 10.0, 6.0, 6.0, 6.0, 8.0, 5.0, 1.0, 9.0, 5.0, 5.0, 1.0, 5.0, 6.0, 4.0, 2.0, 6.0, 4.0, 4.0, 10.0, 8.0, 9.0, 1.0, 4.0, 2.0, 6.0, 9.0, 4.0, 8.0, 4.0, 5.0, 9.0, 6.0, 8.0]
global b_y = 10
global p = [0.805, 0.744, 0.676, 0.197, 0.061, 0.21, 0.32, 0.064, 0.837, 0.979, 0.307, 0.821, 0.636, 0.368, 0.325, 0.546, 0.432, 0.502, 0.238, 0.83, 0.887, 0.885, 0.899, 0.369, 0.312, 0.826, 0.681, 0.88, 0.714, 0.187, 0.796, 0.295, 0.946, 0.439, 0.416, 0.888, 0.86, 0.187, 0.645, 0.553, 0.757, 0.197, 0.937, 0.704, 0.315, 0.705, 0.723, 0.905, 0.118, 0.948, 0.455, 0.841, 0.548, 0.939, 0.28, 0.109, 0.869, 0.776, 0.226, 0.064, 0.11, 0.138, 0.911, 0.639, 0.582, 0.377, 0.918, 0.923, 0.578, 0.977, 0.662, 0.654, 0.235, 0.649, 0.28, 0.029, 0.225, 0.437, 0.352, 0.84, 0.869, 0.837, 0.938, 0.287, 0.203, 0.227, 0.936, 0.302, 0.915, 0.377, 0.282, 0.935, 0.399, 0.866, 0.332, 0.817, 0.702, 0.913, 0.465, 0.91, 0.107, 0.147, 0.235, 0.425, 0.377, 0.169, 0.447, 0.54, 0.681, 0.98, 0.19, 0.299, 0.612, 0.723, 0.33, 0.9, 0.572, 0.053, 0.029, 0.482, 0.041, 0.939, 0.027, 0.705, 0.154, 0.508, 0.752, 0.611, 0.564, 0.378, 0.85, 0.697, 0.712, 0.068, 0.948, 0.278, 0.365, 0.735, 0.54, 0.445, 0.49, 0.117, 0.737, 0.19, 0.602, 0.468, 0.197, 0.97, 0.1, 0.609, 0.742, 0.333, 0.284, 0.205, 0.612, 0.778, 0.755, 0.074, 0.833, 0.666, 0.244, 0.732, 0.123, 0.825, 0.792, 0.549, 0.129, 0.735, 0.623, 0.011, 0.482, 0.816, 0.51, 0.523, 0.22, 0.806, 0.203, 0.468, 0.277, 0.093, 0.516, 0.054, 0.828, 0.932, 0.878, 0.991, 0.788, 0.428, 0.135, 0.346, 0.107, 0.034, 0.987, 0.594, 0.935, 0.904, 0.428, 0.267, 0.952, 0.7, 0.693, 0.052, 0.858, 0.692, 0.534, 0.446, 0.143, 0.556, 0.393, 0.563, 0.769, 0.196, 0.21, 0.578, 0.286, 0.131, 0.278, 0.031, 0.48, 0.45, 0.833, 0.007, 0.205, 0.074, 0.539, 0.405, 0.917, 0.676, 0.269, 0.645, 0.873, 0.234, 0.056, 0.01, 0.569, 0.382, 0.365, 0.613, 0.583, 0.283, 0.285, 0.756, 0.965, 0.746, 0.03, 0.634, 0.632, 0.118, 0.099, 0.377, 0.427, 0.521, 0.694, 0.131, 0.797, 0.241, 0.351, 0.957, 0.351, 0.071, 0.96, 0.847, 0.293, 0.059, 0.702, 0.957, 0.638, 0.134, 0.152, 0.656, 0.301, 0.653, 0.331, 0.437, 0.991, 0.657, 0.177, 0.911, 0.457, 0.358, 0.939, 0.067, 0.13, 0.881, 0.503, 0.901, 0.796, 0.9, 0.275, 0.937, 0.453, 0.604, 0.157, 0.476, 0.069, 0.557, 0.597, 0.369, 0.212, 0.34, 0.95, 0.16, 0.217, 0.228, 0.706, 0.03, 0.534, 0.474, 0.39, 0.069, 0.606, 0.907, 0.631, 0.129, 0.798, 0.956, 0.666, 0.442, 0.186, 0.813, 0.567, 0.154]
global q = [0.965, 0.763, 0.965, 0.639, 0.677, 0.338, 0.798, 0.315, 0.897, 0.985, 0.899, 0.977, 0.726, 0.576, 0.677, 0.705, 0.995, 0.607, 0.684, 0.85, 0.941, 0.989, 0.96, 0.437, 0.571, 0.862, 0.85, 0.943, 0.949, 0.582, 0.99, 0.462, 0.978, 0.853, 0.53, 0.924, 0.949, 0.644, 0.919, 0.816, 0.962, 0.412, 0.977, 0.971, 0.997, 0.866, 0.965, 0.974, 0.578, 0.963, 0.689, 0.951, 0.75, 0.981, 0.34, 0.585, 0.9, 0.86, 0.691, 0.555, 0.999, 0.47, 0.94, 0.776, 0.617, 0.836, 0.922, 0.965, 0.968, 0.991, 0.823, 0.994, 0.914, 0.703, 0.81, 0.902, 0.993, 0.769, 0.919, 0.998, 0.959, 0.886, 0.943, 0.626, 0.674, 0.593, 0.957, 0.653, 0.915, 0.648, 0.696, 0.957, 0.978, 0.878, 0.79, 0.97, 0.715, 0.999, 0.471, 0.993, 0.382, 0.679, 0.826, 0.937, 0.77, 0.5, 0.591, 0.994, 0.898, 0.983, 0.566, 0.994, 0.806, 0.957, 0.776, 0.958, 0.585, 0.251, 0.344, 0.596, 0.919, 0.973, 0.968, 0.739, 0.963, 0.763, 0.973, 0.923, 0.908, 0.835, 0.974, 0.773, 0.78, 0.622, 0.966, 0.508, 0.947, 0.901, 0.913, 0.961, 0.515, 0.333, 0.737, 0.385, 0.914, 0.499, 0.865, 0.985, 0.555, 0.976, 0.808, 0.424, 0.474, 0.252, 0.887, 0.901, 0.985, 0.383, 0.857, 0.897, 0.535, 0.765, 0.869, 0.895, 0.865, 0.558, 0.781, 0.796, 0.686, 0.455, 0.858, 0.894, 0.7, 0.742, 0.295, 0.928, 0.995, 0.922, 0.896, 0.985, 0.948, 0.453, 0.954, 0.996, 0.986, 0.996, 0.948, 0.625, 0.656, 0.998, 0.174, 0.485, 0.996, 0.689, 0.962, 0.931, 0.859, 0.875, 0.999, 0.88, 0.869, 0.899, 0.975, 0.782, 0.577, 0.78, 0.68, 0.742, 0.924, 0.595, 0.975, 0.774, 0.677, 0.902, 0.422, 0.807, 0.631, 0.641, 0.951, 0.624, 0.937, 0.286, 0.304, 0.219, 0.638, 0.709, 0.988, 0.961, 0.854, 0.802, 0.886, 0.877, 0.346, 0.207, 0.915, 0.473, 0.552, 0.93, 0.751, 0.716, 0.5, 0.952, 0.968, 0.87, 0.719, 0.876, 0.731, 0.985, 0.115, 0.471, 0.519, 0.889, 0.76, 0.874, 0.916, 0.777, 0.431, 0.982, 0.589, 0.286, 0.96, 0.902, 0.593, 0.828, 0.86, 0.992, 0.715, 0.272, 0.681, 0.876, 0.803, 0.808, 0.671, 0.632, 0.995, 0.791, 0.88, 0.922, 0.659, 0.66, 0.945, 0.539, 0.45, 0.883, 0.542, 0.953, 0.811, 0.996, 0.491, 0.974, 0.854, 0.847, 0.448, 0.594, 0.444, 0.833, 0.774, 0.729, 0.735, 0.737, 0.963, 0.489, 0.57, 0.261, 0.929, 0.599, 0.863, 0.481, 0.435, 0.468, 0.782, 0.927, 0.707, 0.175, 0.822, 0.962, 0.858, 0.981, 0.912, 0.95, 0.674, 0.249]
global origin = 1
global destination = 60