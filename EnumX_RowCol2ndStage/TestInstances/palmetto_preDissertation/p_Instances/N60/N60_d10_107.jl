global arcs = [1 2; 1 5; 1 12; 1 28; 1 48; 2 4; 2 13; 2 14; 2 18; 2 32; 2 38; 2 46; 2 47; 3 4; 3 16; 3 38; 3 45; 3 52; 3 59; 4 5; 4 17; 4 26; 4 32; 4 35; 4 47; 4 56; 5 10; 5 26; 5 42; 5 47; 5 51; 5 60; 6 2; 6 3; 6 41; 6 45; 7 13; 7 33; 7 42; 7 43; 7 44; 8 15; 8 31; 8 34; 8 44; 8 54; 9 2; 9 23; 9 27; 9 36; 9 38; 9 41; 10 3; 10 6; 10 13; 10 25; 10 29; 10 31; 10 32; 10 41; 10 48; 10 51; 11 2; 11 14; 11 22; 11 26; 11 37; 11 46; 12 38; 12 44; 12 59; 13 5; 13 6; 13 17; 13 27; 13 33; 13 38; 13 43; 13 57; 14 21; 14 24; 14 44; 14 45; 14 46; 14 53; 15 3; 15 5; 15 8; 15 9; 15 19; 15 49; 15 53; 15 55; 16 6; 16 19; 16 21; 16 59; 17 2; 17 9; 17 32; 17 53; 18 12; 18 26; 18 33; 18 45; 18 48; 18 50; 18 52; 18 56; 18 58; 19 8; 19 20; 19 27; 19 28; 19 42; 19 47; 19 48; 19 58; 20 3; 20 18; 20 54; 20 58; 21 9; 21 15; 22 33; 22 34; 22 45; 22 55; 23 24; 23 57; 24 17; 24 22; 24 25; 24 40; 24 60; 25 2; 25 16; 25 20; 25 27; 25 30; 25 38; 25 40; 25 41; 25 43; 25 52; 25 53; 26 2; 26 11; 26 42; 26 43; 27 2; 27 9; 27 15; 27 54; 27 57; 28 6; 28 32; 28 38; 28 40; 28 59; 29 6; 29 22; 29 23; 29 32; 29 39; 30 6; 30 34; 30 45; 30 56; 31 8; 31 19; 31 23; 31 40; 31 54; 32 11; 32 17; 32 20; 32 28; 32 45; 33 26; 33 30; 33 51; 34 12; 34 25; 34 26; 34 39; 34 40; 34 59; 35 9; 35 10; 35 12; 35 25; 35 27; 35 46; 35 51; 36 12; 36 27; 36 33; 36 38; 36 51; 37 11; 37 25; 37 28; 37 32; 37 35; 37 47; 38 11; 38 13; 38 17; 38 19; 38 22; 38 23; 38 27; 38 37; 38 55; 38 56; 39 26; 39 34; 39 46; 39 47; 40 10; 40 15; 40 34; 40 39; 40 57; 41 19; 41 22; 41 25; 41 29; 41 47; 41 50; 41 60; 42 14; 42 30; 42 32; 42 37; 42 47; 43 27; 43 35; 43 47; 43 52; 44 9; 44 30; 44 46; 44 47; 45 2; 45 6; 45 10; 45 13; 45 15; 45 38; 45 44; 46 17; 46 30; 46 31; 46 42; 46 55; 47 38; 47 39; 47 44; 47 48; 47 51; 47 53; 47 57; 48 3; 48 30; 48 31; 49 13; 49 35; 49 38; 49 57; 50 21; 50 24; 51 7; 51 19; 51 25; 51 31; 51 48; 52 4; 52 5; 52 21; 52 24; 52 36; 52 40; 52 53; 52 54; 53 13; 53 30; 53 34; 53 39; 53 49; 53 50; 53 54; 53 60; 54 22; 54 23; 54 33; 54 38; 54 43; 54 56; 54 57; 54 58; 55 24; 55 26; 55 30; 55 33; 55 52; 56 3; 56 5; 56 15; 56 22; 56 39; 56 41; 57 9; 57 23; 57 28; 57 31; 57 49; 57 50; 57 58; 58 5; 58 10; 58 25; 58 34; 58 37; 59 39; 59 45; 59 56]
global d_x = [5.0, 4.0, 1.0, 2.0, 6.0, 3.0, 2.0, 10.0, 8.0, 8.0, 3.0, 8.0, 4.0, 6.0, 6.0, 5.0, 1.0, 8.0, 6.0, 1.0, 1.0, 6.0, 7.0, 10.0, 9.0, 3.0, 2.0, 4.0, 4.0, 3.0, 10.0, 3.0, 7.0, 10.0, 8.0, 1.0, 7.0, 10.0, 6.0, 3.0, 6.0, 9.0, 3.0, 2.0, 3.0, 2.0, 7.0, 6.0, 6.0, 5.0, 1.0, 9.0, 5.0, 2.0, 9.0, 7.0, 2.0, 1.0, 6.0, 2.0, 4.0, 8.0, 1.0, 6.0, 5.0, 2.0, 6.0, 8.0, 5.0, 3.0, 3.0, 1.0, 5.0, 5.0, 2.0, 1.0, 4.0, 8.0, 4.0, 1.0, 6.0, 2.0, 2.0, 9.0, 2.0, 8.0, 1.0, 9.0, 8.0, 7.0, 1.0, 10.0, 9.0, 3.0, 8.0, 2.0, 5.0, 6.0, 6.0, 6.0, 4.0, 10.0, 8.0, 7.0, 5.0, 1.0, 2.0, 7.0, 9.0, 3.0, 7.0, 7.0, 4.0, 7.0, 2.0, 6.0, 4.0, 7.0, 3.0, 8.0, 5.0, 2.0, 8.0, 2.0, 2.0, 8.0, 6.0, 10.0, 9.0, 7.0, 6.0, 6.0, 6.0, 2.0, 5.0, 2.0, 3.0, 5.0, 8.0, 1.0, 8.0, 6.0, 9.0, 7.0, 6.0, 9.0, 5.0, 1.0, 5.0, 7.0, 4.0, 9.0, 6.0, 1.0, 1.0, 4.0, 9.0, 1.0, 10.0, 7.0, 10.0, 3.0, 6.0, 10.0, 4.0, 7.0, 2.0, 3.0, 6.0, 9.0, 2.0, 1.0, 9.0, 5.0, 7.0, 3.0, 6.0, 6.0, 3.0, 5.0, 7.0, 8.0, 2.0, 6.0, 10.0, 3.0, 1.0, 8.0, 10.0, 5.0, 10.0, 9.0, 1.0, 6.0, 10.0, 10.0, 9.0, 6.0, 10.0, 5.0, 9.0, 10.0, 2.0, 8.0, 3.0, 9.0, 8.0, 1.0, 1.0, 3.0, 8.0, 1.0, 10.0, 1.0, 7.0, 1.0, 1.0, 3.0, 3.0, 6.0, 6.0, 3.0, 6.0, 6.0, 9.0, 4.0, 5.0, 6.0, 8.0, 3.0, 10.0, 3.0, 6.0, 2.0, 6.0, 8.0, 6.0, 10.0, 8.0, 5.0, 1.0, 9.0, 6.0, 5.0, 2.0, 6.0, 2.0, 1.0, 7.0, 7.0, 3.0, 4.0, 8.0, 7.0, 4.0, 8.0, 5.0, 1.0, 8.0, 3.0, 3.0, 9.0, 3.0, 4.0, 6.0, 2.0, 3.0, 5.0, 1.0, 3.0, 9.0, 7.0, 5.0, 4.0, 9.0, 1.0, 2.0, 2.0, 7.0, 5.0, 10.0, 4.0, 6.0, 1.0, 3.0, 9.0, 7.0, 2.0, 3.0, 3.0, 1.0, 8.0, 1.0, 6.0, 6.0, 5.0, 2.0, 8.0, 7.0, 9.0, 1.0, 1.0, 6.0, 1.0, 3.0, 6.0, 7.0, 4.0, 4.0, 7.0, 8.0, 1.0, 3.0, 5.0, 10.0, 6.0, 10.0, 1.0, 3.0, 7.0, 9.0, 4.0, 9.0, 8.0, 5.0, 4.0, 10.0, 4.0]
global b_x = 5
global d_y = [4.0, 2.0, 2.0, 7.0, 4.0, 6.0, 5.0, 8.0, 10.0, 10.0, 10.0, 8.0, 1.0, 7.0, 1.0, 6.0, 3.0, 7.0, 8.0, 5.0, 5.0, 3.0, 3.0, 1.0, 8.0, 2.0, 10.0, 9.0, 3.0, 6.0, 8.0, 10.0, 4.0, 8.0, 1.0, 2.0, 5.0, 7.0, 7.0, 7.0, 8.0, 10.0, 4.0, 8.0, 7.0, 5.0, 10.0, 10.0, 7.0, 10.0, 2.0, 7.0, 3.0, 4.0, 3.0, 10.0, 3.0, 9.0, 3.0, 8.0, 8.0, 10.0, 7.0, 9.0, 7.0, 1.0, 3.0, 4.0, 9.0, 7.0, 5.0, 1.0, 3.0, 2.0, 4.0, 10.0, 9.0, 9.0, 5.0, 10.0, 5.0, 4.0, 5.0, 4.0, 7.0, 1.0, 5.0, 1.0, 10.0, 6.0, 3.0, 10.0, 1.0, 2.0, 9.0, 8.0, 7.0, 9.0, 10.0, 4.0, 3.0, 7.0, 6.0, 9.0, 2.0, 2.0, 8.0, 5.0, 7.0, 1.0, 6.0, 7.0, 10.0, 4.0, 1.0, 5.0, 7.0, 8.0, 2.0, 6.0, 1.0, 8.0, 1.0, 9.0, 10.0, 5.0, 6.0, 4.0, 8.0, 3.0, 7.0, 3.0, 3.0, 6.0, 1.0, 8.0, 7.0, 8.0, 9.0, 2.0, 1.0, 10.0, 7.0, 5.0, 10.0, 4.0, 2.0, 5.0, 8.0, 1.0, 10.0, 9.0, 6.0, 8.0, 9.0, 5.0, 2.0, 3.0, 2.0, 7.0, 6.0, 5.0, 2.0, 5.0, 10.0, 2.0, 9.0, 5.0, 6.0, 9.0, 7.0, 7.0, 2.0, 2.0, 2.0, 10.0, 6.0, 2.0, 4.0, 4.0, 9.0, 1.0, 3.0, 7.0, 9.0, 6.0, 7.0, 10.0, 5.0, 10.0, 1.0, 5.0, 9.0, 4.0, 5.0, 6.0, 2.0, 5.0, 3.0, 4.0, 8.0, 9.0, 9.0, 4.0, 1.0, 3.0, 5.0, 7.0, 6.0, 2.0, 10.0, 9.0, 8.0, 7.0, 7.0, 7.0, 4.0, 7.0, 7.0, 7.0, 8.0, 9.0, 2.0, 4.0, 4.0, 2.0, 1.0, 2.0, 4.0, 5.0, 5.0, 4.0, 7.0, 8.0, 6.0, 5.0, 6.0, 8.0, 4.0, 2.0, 2.0, 9.0, 6.0, 8.0, 5.0, 6.0, 2.0, 5.0, 5.0, 8.0, 3.0, 5.0, 1.0, 8.0, 6.0, 10.0, 6.0, 2.0, 6.0, 2.0, 9.0, 5.0, 8.0, 9.0, 4.0, 3.0, 10.0, 3.0, 8.0, 2.0, 3.0, 4.0, 6.0, 6.0, 1.0, 8.0, 2.0, 2.0, 1.0, 1.0, 2.0, 8.0, 9.0, 10.0, 8.0, 1.0, 3.0, 2.0, 8.0, 2.0, 3.0, 9.0, 7.0, 7.0, 8.0, 7.0, 8.0, 7.0, 4.0, 5.0, 6.0, 9.0, 8.0, 7.0, 2.0, 6.0, 1.0, 7.0, 4.0, 9.0, 3.0, 10.0, 10.0, 4.0, 1.0, 5.0, 4.0, 6.0, 3.0, 10.0, 1.0, 2.0, 6.0, 3.0, 8.0, 4.0, 3.0, 1.0]
global b_y = 10
global p = [0.022, 0.329, 0.889, 0.982, 0.063, 0.762, 0.749, 0.538, 0.999, 0.813, 0.737, 0.018, 0.101, 0.53, 0.697, 0.547, 0.726, 0.829, 0.093, 0.613, 0.13, 0.507, 0.875, 0.808, 0.307, 0.959, 0.375, 0.755, 0.804, 0.28, 0.804, 0.859, 0.189, 0.025, 0.485, 0.794, 0.72, 0.941, 0.044, 0.762, 0.759, 0.715, 0.643, 0.214, 0.09, 0.336, 0.366, 0.377, 0.252, 0.723, 0.849, 0.888, 0.002, 0.818, 0.957, 0.612, 0.279, 0.732, 0.998, 0.151, 0.942, 0.639, 0.061, 0.781, 0.93, 0.874, 0.783, 0.213, 0.17, 0.983, 0.446, 0.629, 0.571, 0.807, 0.009, 0.347, 0.702, 0.873, 0.029, 0.339, 0.522, 0.804, 0.601, 0.242, 0.539, 0.075, 0.253, 0.306, 0.262, 0.164, 0.375, 0.101, 0.102, 0.327, 0.56, 0.337, 0.484, 0.011, 0.228, 0.598, 0.453, 0.191, 0.907, 0.945, 0.914, 0.674, 0.212, 0.264, 0.878, 0.087, 0.8, 0.507, 0.732, 0.978, 0.097, 0.105, 0.869, 0.54, 0.078, 0.599, 0.571, 0.103, 0.535, 0.827, 0.265, 0.408, 0.706, 0.907, 0.048, 0.164, 0.33, 0.268, 0.809, 0.491, 0.026, 0.225, 0.346, 0.354, 0.049, 0.743, 0.388, 0.296, 0.337, 0.388, 0.861, 0.317, 0.718, 0.087, 0.416, 0.734, 0.244, 0.502, 0.657, 0.525, 0.789, 0.979, 0.676, 0.208, 0.54, 0.083, 0.647, 0.921, 0.692, 0.426, 0.348, 0.597, 0.134, 0.148, 0.085, 0.607, 0.257, 0.127, 0.38, 0.35, 0.079, 0.888, 0.594, 0.916, 0.179, 0.601, 0.669, 0.638, 0.172, 0.38, 0.036, 0.355, 0.967, 0.872, 0.6, 0.359, 0.933, 0.537, 0.118, 0.205, 0.194, 0.643, 0.821, 0.983, 0.04, 0.045, 0.952, 0.634, 0.5, 0.788, 0.054, 0.857, 0.482, 0.487, 0.952, 0.812, 0.259, 0.858, 0.443, 0.433, 0.095, 0.316, 0.164, 0.169, 0.758, 0.444, 0.732, 0.22, 0.978, 0.444, 0.341, 0.672, 0.071, 0.977, 0.869, 0.792, 0.087, 0.002, 0.155, 0.626, 0.165, 0.482, 0.099, 0.872, 0.062, 0.438, 0.429, 0.13, 0.334, 0.367, 0.885, 0.987, 0.303, 0.861, 0.8, 0.538, 0.37, 0.497, 0.157, 0.738, 0.85, 0.97, 0.453, 0.418, 0.192, 0.498, 0.372, 0.332, 0.408, 0.07, 0.594, 0.656, 0.761, 0.835, 0.573, 0.292, 0.233, 0.858, 0.911, 0.606, 0.15, 0.655, 0.127, 0.37, 0.861, 0.205, 0.829, 0.505, 0.22, 0.54, 0.829, 0.306, 0.422, 0.082, 0.809, 0.178, 0.712, 0.801, 0.986, 0.276, 0.618, 0.744, 0.188, 0.724, 0.971, 0.61, 0.452, 0.164, 0.03, 0.33, 0.346, 0.995, 0.813, 0.793, 0.474, 0.598, 0.881, 0.871, 0.315, 0.503, 0.275, 0.027, 0.115, 0.896, 0.47, 0.784, 0.369, 0.124, 0.76, 0.482, 0.752, 0.233, 0.976, 0.824]
global q = [0.645, 0.784, 0.976, 0.999, 0.212, 0.886, 0.805, 0.655, 0.999, 0.973, 0.95, 0.201, 0.977, 0.623, 0.775, 0.687, 0.775, 0.87, 0.205, 0.619, 0.348, 0.8, 0.946, 0.875, 0.363, 0.974, 0.943, 0.892, 0.812, 0.358, 0.942, 0.958, 0.938, 0.628, 0.525, 0.986, 0.979, 0.981, 0.21, 0.918, 0.973, 0.882, 0.918, 0.747, 0.695, 0.685, 0.87, 0.562, 0.958, 0.908, 0.872, 0.973, 0.35, 0.992, 0.968, 0.862, 0.822, 0.877, 0.999, 0.465, 0.972, 0.734, 0.214, 0.831, 0.957, 0.91, 0.822, 0.4, 0.781, 0.992, 0.725, 0.69, 0.889, 0.81, 0.038, 0.832, 0.914, 0.923, 0.451, 0.418, 0.934, 0.828, 0.759, 0.696, 0.629, 0.639, 0.555, 0.867, 0.787, 0.732, 0.647, 0.436, 0.704, 0.579, 0.739, 0.477, 0.807, 0.461, 0.563, 0.929, 0.525, 0.898, 0.995, 0.954, 0.997, 0.845, 0.413, 0.879, 0.898, 0.957, 0.918, 0.535, 0.786, 0.991, 0.505, 0.28, 0.878, 0.744, 0.975, 0.603, 0.759, 0.145, 0.71, 0.952, 0.823, 0.992, 0.796, 0.916, 0.706, 0.402, 0.976, 0.892, 0.942, 0.877, 0.327, 0.801, 0.734, 0.984, 0.16, 0.753, 0.933, 0.672, 0.665, 0.686, 0.999, 0.726, 0.851, 0.348, 0.684, 0.987, 0.576, 0.802, 0.879, 0.595, 0.899, 0.996, 0.898, 0.629, 0.823, 0.927, 0.76, 0.986, 0.813, 0.717, 0.781, 0.932, 0.794, 0.299, 0.899, 0.746, 0.865, 0.603, 0.895, 0.823, 0.593, 0.957, 0.921, 0.939, 0.744, 0.928, 0.93, 0.699, 0.529, 0.974, 0.046, 0.97, 0.978, 0.911, 0.637, 0.97, 0.99, 0.661, 0.359, 0.38, 0.563, 0.719, 0.873, 0.986, 0.571, 0.255, 0.964, 0.916, 0.792, 0.815, 0.129, 0.895, 0.514, 0.963, 0.958, 0.99, 0.841, 0.876, 0.515, 0.731, 0.676, 0.992, 0.813, 0.428, 0.85, 0.695, 0.974, 0.944, 0.999, 0.998, 0.437, 0.728, 0.514, 0.993, 0.937, 0.886, 0.473, 0.978, 0.523, 0.633, 0.51, 0.823, 0.941, 0.897, 0.356, 0.542, 0.471, 0.775, 0.583, 0.765, 0.999, 0.994, 0.433, 0.924, 0.934, 0.87, 0.617, 0.653, 0.812, 0.794, 0.989, 0.975, 0.55, 0.699, 0.802, 0.638, 0.895, 0.522, 0.833, 0.271, 0.935, 0.824, 0.806, 0.999, 0.957, 0.724, 0.629, 0.908, 0.969, 0.975, 0.767, 0.687, 0.337, 0.683, 0.979, 0.375, 0.88, 0.98, 0.555, 0.693, 0.889, 0.48, 0.823, 0.274, 0.894, 0.979, 0.922, 0.973, 0.99, 0.433, 0.624, 0.869, 0.342, 0.902, 0.989, 0.889, 0.595, 0.839, 0.639, 0.539, 0.541, 0.999, 0.887, 0.827, 0.869, 0.625, 0.928, 0.877, 0.333, 0.542, 0.53, 0.356, 0.812, 0.966, 0.922, 0.866, 0.848, 0.473, 0.862, 0.88, 0.911, 0.569, 0.993, 0.868]
global origin = 1
global destination = 60