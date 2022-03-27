global arcs = [1 18; 1 37; 1 40; 1 52; 2 7; 2 12; 2 23; 2 25; 2 41; 2 49; 3 9; 3 13; 3 17; 3 24; 3 28; 3 40; 3 41; 3 52; 4 18; 4 23; 4 24; 4 38; 4 49; 4 53; 5 17; 5 23; 5 36; 5 51; 6 24; 6 31; 7 10; 7 21; 7 38; 7 42; 8 23; 8 37; 8 60; 9 13; 9 16; 9 18; 9 23; 9 27; 9 33; 10 7; 10 8; 10 22; 10 24; 10 40; 10 55; 11 12; 11 29; 11 30; 11 34; 11 35; 12 3; 12 6; 12 22; 12 25; 12 43; 13 6; 13 8; 13 18; 13 26; 13 32; 13 33; 13 45; 13 46; 13 48; 14 11; 14 16; 14 23; 14 29; 14 39; 14 50; 15 8; 15 14; 15 34; 15 36; 16 11; 16 14; 16 25; 16 37; 16 53; 16 57; 17 5; 17 13; 17 25; 17 30; 17 37; 17 43; 18 21; 18 52; 19 44; 19 53; 19 59; 20 39; 20 51; 20 52; 21 3; 21 11; 21 20; 22 36; 22 42; 22 43; 22 49; 22 50; 22 55; 22 59; 22 60; 23 2; 23 5; 23 10; 23 41; 23 54; 24 21; 24 38; 24 41; 25 2; 25 19; 25 21; 25 35; 25 36; 25 39; 25 59; 26 2; 26 9; 26 11; 26 12; 26 16; 26 21; 26 32; 27 5; 27 16; 27 20; 27 26; 27 38; 27 42; 27 58; 28 3; 28 6; 28 27; 28 56; 29 7; 29 23; 29 25; 29 28; 29 57; 30 32; 30 41; 30 52; 31 28; 31 32; 31 58; 32 9; 32 16; 32 19; 32 22; 32 33; 32 36; 32 51; 32 52; 33 9; 33 15; 33 22; 33 35; 33 47; 33 58; 34 8; 34 15; 34 22; 34 23; 34 33; 35 6; 35 7; 35 12; 35 24; 35 37; 35 45; 35 56; 35 59; 35 60; 36 9; 36 16; 36 29; 36 34; 36 48; 36 52; 37 4; 37 11; 37 12; 37 13; 37 20; 37 23; 37 24; 37 33; 37 34; 37 44; 37 48; 37 49; 38 14; 38 22; 38 39; 38 58; 39 6; 39 8; 39 31; 39 41; 39 45; 39 48; 40 5; 40 18; 40 54; 40 57; 40 60; 41 5; 41 15; 41 36; 41 40; 41 46; 41 49; 42 25; 42 27; 42 28; 42 38; 42 53; 42 56; 43 3; 43 4; 43 8; 43 16; 43 27; 43 30; 43 46; 43 47; 43 53; 43 55; 44 20; 44 33; 44 42; 45 3; 45 31; 45 43; 46 4; 46 7; 46 24; 46 28; 46 30; 46 32; 46 36; 46 58; 47 5; 47 12; 47 35; 47 40; 47 48; 47 49; 47 53; 47 58; 48 11; 48 28; 48 46; 48 59; 49 3; 49 21; 49 35; 49 40; 49 58; 50 25; 50 28; 50 29; 50 31; 50 45; 50 55; 51 9; 51 11; 51 22; 51 30; 51 35; 51 39; 51 44; 51 50; 51 52; 51 56; 52 6; 52 27; 52 39; 52 58; 53 3; 53 42; 53 43; 53 50; 53 52; 53 57; 53 60; 54 6; 54 31; 54 32; 54 35; 54 37; 54 39; 54 43; 54 51; 55 7; 55 12; 55 15; 55 18; 55 29; 55 59; 56 2; 56 37; 56 47; 57 14; 57 19; 57 44; 57 45; 57 51; 57 56; 58 16; 58 21; 58 30; 58 37; 58 43; 58 53; 59 3; 59 5; 59 6; 59 19; 59 29; 59 40; 59 51; 59 56]
global d_x = [8.0, 5.0, 4.0, 1.0, 8.0, 5.0, 4.0, 3.0, 4.0, 4.0, 10.0, 1.0, 7.0, 1.0, 6.0, 4.0, 10.0, 8.0, 9.0, 2.0, 6.0, 1.0, 3.0, 3.0, 7.0, 9.0, 9.0, 1.0, 10.0, 9.0, 5.0, 5.0, 1.0, 6.0, 6.0, 7.0, 5.0, 10.0, 1.0, 5.0, 10.0, 10.0, 5.0, 7.0, 2.0, 1.0, 3.0, 8.0, 5.0, 1.0, 2.0, 1.0, 9.0, 8.0, 5.0, 5.0, 1.0, 6.0, 10.0, 7.0, 3.0, 1.0, 3.0, 7.0, 6.0, 1.0, 10.0, 6.0, 3.0, 9.0, 1.0, 3.0, 9.0, 5.0, 2.0, 2.0, 7.0, 3.0, 10.0, 7.0, 8.0, 5.0, 8.0, 4.0, 1.0, 7.0, 4.0, 2.0, 1.0, 9.0, 3.0, 9.0, 6.0, 4.0, 8.0, 10.0, 1.0, 7.0, 3.0, 5.0, 9.0, 4.0, 7.0, 4.0, 1.0, 4.0, 8.0, 7.0, 10.0, 10.0, 9.0, 9.0, 7.0, 10.0, 1.0, 10.0, 3.0, 6.0, 2.0, 3.0, 4.0, 3.0, 3.0, 4.0, 4.0, 9.0, 4.0, 10.0, 9.0, 5.0, 2.0, 6.0, 4.0, 8.0, 1.0, 1.0, 4.0, 3.0, 7.0, 3.0, 7.0, 3.0, 6.0, 3.0, 2.0, 6.0, 4.0, 4.0, 6.0, 6.0, 10.0, 6.0, 3.0, 3.0, 1.0, 9.0, 7.0, 3.0, 9.0, 2.0, 6.0, 5.0, 5.0, 3.0, 4.0, 10.0, 3.0, 4.0, 8.0, 8.0, 9.0, 1.0, 1.0, 4.0, 6.0, 9.0, 1.0, 3.0, 5.0, 4.0, 5.0, 10.0, 10.0, 10.0, 5.0, 6.0, 2.0, 1.0, 9.0, 8.0, 5.0, 6.0, 10.0, 9.0, 4.0, 2.0, 5.0, 5.0, 2.0, 6.0, 2.0, 2.0, 9.0, 8.0, 9.0, 4.0, 9.0, 3.0, 1.0, 7.0, 5.0, 2.0, 4.0, 8.0, 7.0, 9.0, 9.0, 3.0, 2.0, 9.0, 7.0, 9.0, 4.0, 6.0, 6.0, 7.0, 6.0, 5.0, 1.0, 7.0, 7.0, 2.0, 7.0, 5.0, 4.0, 7.0, 8.0, 1.0, 10.0, 4.0, 9.0, 5.0, 1.0, 2.0, 4.0, 8.0, 5.0, 8.0, 7.0, 10.0, 3.0, 1.0, 10.0, 2.0, 5.0, 3.0, 5.0, 1.0, 8.0, 3.0, 9.0, 7.0, 4.0, 10.0, 2.0, 9.0, 6.0, 7.0, 5.0, 4.0, 7.0, 7.0, 6.0, 7.0, 3.0, 5.0, 10.0, 3.0, 8.0, 10.0, 2.0, 5.0, 9.0, 6.0, 7.0, 4.0, 1.0, 1.0, 3.0, 5.0, 3.0, 9.0, 2.0, 1.0, 8.0, 5.0, 3.0, 8.0, 6.0, 8.0, 1.0, 8.0, 8.0, 3.0, 4.0, 7.0, 6.0, 4.0, 6.0, 4.0, 1.0, 2.0, 6.0, 5.0, 6.0, 8.0, 5.0, 3.0, 3.0, 8.0, 7.0, 9.0, 6.0, 10.0, 5.0, 3.0, 8.0, 2.0, 1.0, 2.0, 3.0]
global b_x = 5
global d_y = [7.0, 2.0, 4.0, 5.0, 9.0, 10.0, 4.0, 5.0, 7.0, 10.0, 3.0, 7.0, 9.0, 6.0, 7.0, 1.0, 5.0, 10.0, 8.0, 4.0, 4.0, 5.0, 6.0, 8.0, 6.0, 6.0, 3.0, 10.0, 8.0, 2.0, 9.0, 3.0, 9.0, 8.0, 7.0, 6.0, 8.0, 1.0, 7.0, 4.0, 1.0, 10.0, 9.0, 1.0, 6.0, 4.0, 6.0, 6.0, 10.0, 2.0, 4.0, 1.0, 3.0, 5.0, 5.0, 2.0, 1.0, 4.0, 4.0, 9.0, 7.0, 5.0, 1.0, 8.0, 5.0, 1.0, 6.0, 8.0, 10.0, 9.0, 9.0, 5.0, 6.0, 4.0, 6.0, 4.0, 1.0, 9.0, 6.0, 1.0, 7.0, 8.0, 2.0, 6.0, 5.0, 5.0, 9.0, 9.0, 5.0, 9.0, 7.0, 2.0, 5.0, 10.0, 3.0, 5.0, 2.0, 4.0, 3.0, 5.0, 9.0, 5.0, 2.0, 10.0, 5.0, 1.0, 8.0, 6.0, 2.0, 4.0, 10.0, 2.0, 4.0, 5.0, 8.0, 4.0, 1.0, 9.0, 3.0, 7.0, 10.0, 6.0, 10.0, 7.0, 2.0, 6.0, 1.0, 7.0, 7.0, 6.0, 10.0, 5.0, 9.0, 1.0, 6.0, 10.0, 7.0, 3.0, 9.0, 9.0, 10.0, 9.0, 10.0, 1.0, 1.0, 4.0, 1.0, 4.0, 1.0, 3.0, 8.0, 1.0, 10.0, 4.0, 2.0, 7.0, 3.0, 5.0, 3.0, 10.0, 1.0, 4.0, 4.0, 7.0, 6.0, 5.0, 2.0, 8.0, 1.0, 3.0, 4.0, 1.0, 8.0, 3.0, 1.0, 2.0, 8.0, 3.0, 2.0, 10.0, 9.0, 4.0, 4.0, 10.0, 2.0, 4.0, 1.0, 5.0, 6.0, 9.0, 7.0, 3.0, 5.0, 5.0, 9.0, 6.0, 5.0, 2.0, 9.0, 10.0, 5.0, 1.0, 7.0, 7.0, 6.0, 4.0, 1.0, 9.0, 10.0, 4.0, 2.0, 10.0, 1.0, 8.0, 4.0, 8.0, 4.0, 9.0, 4.0, 1.0, 2.0, 2.0, 9.0, 9.0, 6.0, 6.0, 2.0, 1.0, 5.0, 7.0, 8.0, 2.0, 5.0, 10.0, 7.0, 2.0, 10.0, 9.0, 1.0, 3.0, 4.0, 5.0, 10.0, 10.0, 3.0, 10.0, 10.0, 2.0, 3.0, 3.0, 3.0, 6.0, 5.0, 1.0, 4.0, 9.0, 10.0, 8.0, 8.0, 10.0, 8.0, 1.0, 3.0, 7.0, 4.0, 2.0, 7.0, 10.0, 8.0, 7.0, 1.0, 3.0, 4.0, 7.0, 2.0, 4.0, 4.0, 3.0, 2.0, 2.0, 8.0, 8.0, 10.0, 1.0, 10.0, 8.0, 10.0, 6.0, 1.0, 10.0, 1.0, 3.0, 10.0, 2.0, 7.0, 5.0, 1.0, 1.0, 8.0, 5.0, 7.0, 6.0, 7.0, 7.0, 4.0, 1.0, 6.0, 5.0, 10.0, 5.0, 8.0, 2.0, 6.0, 9.0, 8.0, 6.0, 4.0, 2.0, 1.0, 8.0, 7.0, 10.0, 9.0, 8.0, 7.0, 4.0, 9.0, 1.0, 1.0, 8.0, 6.0]
global b_y = 10
global p = [0.474, 0.576, 0.741, 0.74, 0.813, 0.093, 0.858, 0.433, 0.883, 0.127, 0.238, 0.084, 0.54, 0.367, 0.715, 0.103, 0.41, 0.804, 0.227, 0.433, 0.465, 0.656, 0.203, 0.908, 0.069, 0.076, 0.967, 0.341, 0.113, 0.607, 0.723, 0.388, 0.04, 0.508, 0.58, 0.163, 0.442, 0.439, 0.474, 0.272, 0.074, 0.86, 0.334, 0.608, 0.647, 0.781, 0.624, 0.313, 0.072, 0.399, 0.11, 0.204, 0.458, 0.428, 0.198, 0.569, 0.838, 0.555, 0.758, 0.82, 0.573, 0.978, 0.126, 0.979, 0.113, 0.503, 0.932, 0.546, 0.157, 0.189, 0.325, 0.551, 0.902, 0.995, 0.453, 0.56, 0.334, 0.977, 0.318, 0.606, 0.788, 0.007, 0.689, 0.184, 0.075, 0.025, 0.224, 0.72, 0.886, 0.594, 0.242, 0.553, 0.706, 0.26, 0.197, 0.791, 0.488, 0.674, 0.051, 0.106, 0.878, 0.866, 0.783, 0.452, 0.079, 0.724, 0.327, 0.863, 0.825, 0.864, 0.998, 0.505, 0.329, 0.565, 0.694, 0.661, 0.973, 0.426, 0.739, 0.563, 0.509, 0.844, 0.288, 0.224, 0.135, 0.053, 0.565, 0.812, 0.904, 0.34, 0.303, 0.217, 0.342, 0.984, 0.57, 0.731, 0.973, 0.186, 0.218, 0.397, 0.956, 0.12, 0.031, 0.297, 0.845, 0.944, 0.605, 0.723, 0.519, 0.057, 0.098, 0.138, 0.847, 0.004, 0.348, 0.056, 0.029, 0.524, 0.026, 0.677, 0.098, 0.697, 0.345, 0.427, 0.743, 0.456, 0.65, 0.535, 0.558, 0.402, 0.269, 0.036, 0.778, 0.037, 0.733, 0.764, 0.128, 0.138, 0.833, 0.04, 0.037, 0.759, 0.101, 0.319, 0.401, 0.062, 0.761, 0.659, 0.991, 0.26, 0.081, 0.298, 0.706, 0.162, 0.914, 0.066, 0.007, 0.593, 0.722, 0.831, 0.402, 0.917, 0.128, 0.497, 0.698, 0.885, 0.676, 0.674, 0.335, 0.799, 0.829, 0.368, 0.45, 0.194, 0.075, 0.409, 0.104, 0.646, 0.754, 0.125, 0.258, 0.515, 0.821, 0.861, 0.263, 0.973, 0.787, 0.517, 0.842, 0.86, 0.582, 0.958, 0.758, 0.854, 0.704, 0.142, 0.616, 0.264, 0.964, 0.756, 0.146, 0.112, 0.994, 0.681, 0.74, 0.724, 0.047, 0.917, 0.647, 0.946, 0.761, 0.553, 0.641, 0.01, 0.749, 0.434, 0.252, 0.657, 0.637, 0.498, 0.533, 0.77, 0.98, 0.471, 0.34, 0.928, 0.133, 0.022, 0.143, 0.415, 0.824, 0.441, 0.213, 0.326, 0.17, 0.588, 0.223, 0.056, 0.64, 0.648, 0.683, 0.394, 0.835, 0.961, 0.45, 0.176, 0.945, 0.871, 0.663, 0.842, 0.637, 0.04, 0.725, 0.981, 0.309, 0.026, 0.696, 0.73, 0.006, 0.097, 0.355, 0.91, 0.195, 0.352, 0.622, 0.529, 0.55, 0.717, 0.607, 0.61, 0.469, 0.753, 0.706, 0.127, 0.315, 0.614, 0.403, 0.262, 0.491, 0.157, 0.612, 0.505, 0.743, 0.387, 0.554, 0.752, 0.326, 0.847, 0.95, 0.175, 0.806]
global q = [0.772, 0.707, 0.743, 0.944, 0.835, 0.99, 0.983, 0.604, 0.967, 0.915, 0.259, 0.834, 0.797, 0.727, 0.883, 0.884, 0.586, 0.938, 0.292, 0.57, 0.671, 0.746, 0.394, 0.911, 0.878, 0.579, 0.989, 0.729, 0.396, 0.608, 0.86, 0.546, 0.221, 0.61, 0.584, 0.509, 0.873, 0.843, 0.733, 0.922, 0.558, 0.903, 0.967, 0.982, 0.739, 0.83, 0.911, 0.729, 0.389, 0.803, 0.91, 0.54, 0.475, 0.906, 0.47, 0.639, 0.894, 0.573, 0.784, 0.842, 0.845, 0.986, 0.868, 0.993, 0.824, 0.567, 0.932, 0.626, 0.44, 0.667, 0.816, 0.741, 0.968, 0.996, 0.996, 0.951, 0.428, 0.981, 0.769, 0.626, 0.885, 0.244, 0.891, 0.677, 0.988, 0.792, 0.474, 0.981, 0.939, 0.831, 0.596, 0.622, 0.75, 0.743, 0.682, 0.868, 0.901, 0.799, 0.945, 0.509, 0.973, 0.972, 0.895, 0.689, 0.117, 0.791, 0.458, 0.966, 0.855, 0.927, 0.998, 0.58, 0.786, 0.785, 0.878, 0.905, 0.977, 0.579, 0.95, 0.664, 0.868, 0.958, 0.601, 0.334, 0.37, 0.198, 0.772, 0.953, 0.927, 0.537, 0.821, 0.344, 0.494, 0.994, 0.605, 0.738, 0.984, 0.607, 0.284, 0.777, 0.983, 0.602, 0.597, 0.769, 0.985, 0.999, 0.855, 0.739, 0.526, 0.942, 0.399, 0.258, 0.861, 0.634, 0.592, 0.568, 0.934, 0.766, 0.993, 0.88, 0.452, 0.792, 0.677, 0.792, 0.784, 0.559, 0.731, 0.777, 0.927, 0.979, 0.811, 0.583, 0.84, 0.065, 0.796, 0.938, 0.137, 0.45, 0.98, 0.217, 0.947, 0.824, 0.143, 0.959, 0.527, 0.655, 0.929, 0.981, 0.999, 0.869, 0.174, 0.909, 0.824, 0.835, 0.98, 0.173, 0.565, 0.898, 0.752, 0.951, 0.557, 0.981, 0.481, 0.723, 0.876, 0.931, 0.682, 0.74, 0.787, 0.875, 0.857, 0.8, 0.83, 0.496, 0.298, 0.741, 0.84, 0.853, 0.782, 0.281, 0.525, 0.724, 0.972, 0.887, 0.808, 0.98, 0.84, 0.559, 0.912, 0.929, 0.607, 0.994, 0.789, 0.995, 0.744, 0.67, 0.632, 0.697, 0.985, 0.944, 0.335, 0.235, 0.997, 0.884, 0.973, 0.868, 0.559, 0.942, 0.667, 0.973, 0.789, 0.867, 0.941, 0.067, 0.895, 0.449, 0.328, 0.865, 0.792, 0.812, 0.692, 0.845, 0.993, 0.551, 0.828, 0.987, 0.531, 0.041, 0.688, 0.851, 0.942, 0.817, 0.238, 0.444, 0.362, 0.788, 0.918, 0.711, 0.776, 0.994, 0.943, 0.762, 0.972, 0.993, 0.88, 0.342, 0.97, 0.998, 0.751, 0.952, 0.875, 0.81, 0.939, 0.997, 0.973, 0.992, 0.711, 0.935, 0.123, 0.892, 0.472, 0.938, 0.46, 0.387, 0.666, 0.699, 0.629, 0.88, 0.657, 0.772, 0.97, 0.835, 0.979, 0.861, 0.967, 0.686, 0.557, 0.816, 0.718, 0.707, 0.954, 0.948, 0.809, 0.887, 0.909, 0.841, 0.44, 0.892, 0.981, 0.556, 0.891]
global origin = 1
global destination = 60