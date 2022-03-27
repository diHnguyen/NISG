global arcs = [1 8; 1 14; 1 15; 1 39; 1 48; 1 49; 1 52; 1 55; 2 8; 2 15; 2 18; 2 20; 2 23; 2 27; 2 42; 2 53; 2 58; 3 12; 3 21; 3 31; 3 40; 3 52; 3 60; 4 14; 4 26; 4 41; 4 54; 5 2; 5 4; 5 7; 5 44; 5 49; 5 56; 6 23; 6 38; 6 55; 7 9; 7 21; 7 52; 8 3; 8 9; 8 25; 8 49; 8 51; 8 54; 8 60; 9 10; 9 12; 9 26; 9 52; 9 54; 9 55; 9 56; 10 22; 10 23; 10 31; 10 45; 10 58; 11 2; 11 4; 11 27; 11 41; 11 44; 11 45; 12 47; 13 10; 13 15; 13 28; 13 35; 13 36; 13 48; 13 49; 13 53; 14 7; 14 20; 15 21; 15 32; 15 33; 15 36; 16 13; 16 18; 16 27; 16 31; 16 51; 17 4; 17 25; 17 33; 17 43; 17 52; 17 53; 18 5; 18 8; 18 19; 18 27; 18 37; 18 54; 19 5; 19 17; 19 24; 19 35; 19 45; 19 59; 20 10; 20 11; 20 23; 20 30; 20 39; 20 42; 20 44; 20 49; 21 17; 21 47; 22 10; 22 13; 22 19; 22 36; 22 38; 22 43; 23 27; 23 29; 23 46; 23 53; 24 6; 24 14; 24 26; 24 37; 24 40; 24 43; 24 46; 24 49; 24 50; 24 53; 24 58; 25 11; 25 13; 25 14; 25 20; 25 37; 25 54; 25 56; 26 8; 26 24; 26 36; 26 42; 26 50; 26 55; 26 56; 26 57; 27 16; 27 17; 27 34; 27 36; 27 37; 27 40; 27 42; 27 47; 27 56; 28 8; 28 9; 28 25; 29 20; 29 27; 29 38; 29 41; 30 10; 30 16; 30 17; 30 39; 30 47; 31 10; 31 16; 31 21; 31 24; 31 35; 31 43; 31 45; 31 48; 31 55; 31 59; 32 4; 32 29; 32 30; 33 21; 33 22; 33 30; 33 36; 33 46; 33 54; 33 56; 34 13; 34 40; 34 42; 35 4; 35 21; 35 22; 35 39; 35 44; 35 48; 35 51; 35 56; 35 58; 36 2; 36 13; 36 20; 36 48; 36 57; 37 17; 37 26; 37 31; 37 32; 37 33; 37 34; 37 43; 38 27; 38 51; 39 3; 39 6; 39 8; 39 21; 39 47; 39 53; 39 57; 40 3; 40 5; 40 7; 40 36; 40 50; 41 4; 41 7; 41 10; 41 14; 41 25; 41 28; 41 29; 41 31; 41 33; 41 34; 41 39; 41 40; 42 55; 43 2; 43 29; 44 8; 44 30; 44 35; 44 52; 45 33; 45 43; 45 52; 46 9; 46 19; 46 20; 46 23; 46 31; 46 38; 46 48; 46 57; 46 58; 47 3; 47 4; 47 22; 47 28; 47 31; 47 36; 48 13; 48 30; 48 31; 48 47; 48 49; 48 55; 48 60; 49 23; 49 25; 49 51; 49 52; 50 2; 50 33; 50 39; 50 45; 50 48; 50 49; 50 53; 50 57; 51 17; 51 26; 51 33; 51 40; 51 44; 51 48; 51 49; 52 6; 52 8; 52 19; 52 22; 52 23; 52 33; 52 42; 52 56; 53 8; 53 19; 53 32; 53 35; 54 5; 54 23; 54 26; 54 29; 55 16; 55 41; 55 53; 56 13; 56 22; 56 35; 56 60; 57 8; 57 14; 57 31; 57 53; 58 4; 58 8; 58 12; 58 13; 58 24; 58 28; 58 29; 58 30; 58 38; 59 29; 59 56; 59 60]
global d_x = [2.0, 7.0, 3.0, 5.0, 8.0, 5.0, 1.0, 7.0, 8.0, 8.0, 5.0, 4.0, 8.0, 3.0, 2.0, 5.0, 7.0, 6.0, 9.0, 5.0, 2.0, 9.0, 6.0, 6.0, 5.0, 5.0, 5.0, 1.0, 9.0, 9.0, 6.0, 5.0, 10.0, 2.0, 8.0, 4.0, 9.0, 10.0, 2.0, 10.0, 10.0, 10.0, 9.0, 4.0, 8.0, 2.0, 3.0, 5.0, 4.0, 6.0, 2.0, 5.0, 1.0, 10.0, 1.0, 3.0, 7.0, 4.0, 8.0, 4.0, 2.0, 3.0, 6.0, 7.0, 3.0, 10.0, 10.0, 8.0, 1.0, 8.0, 9.0, 3.0, 8.0, 10.0, 2.0, 10.0, 7.0, 2.0, 9.0, 9.0, 10.0, 1.0, 5.0, 4.0, 2.0, 6.0, 8.0, 4.0, 5.0, 8.0, 3.0, 10.0, 8.0, 2.0, 6.0, 10.0, 2.0, 6.0, 2.0, 6.0, 9.0, 10.0, 6.0, 1.0, 8.0, 3.0, 6.0, 9.0, 7.0, 8.0, 2.0, 6.0, 3.0, 9.0, 6.0, 6.0, 4.0, 4.0, 1.0, 9.0, 7.0, 8.0, 10.0, 3.0, 5.0, 6.0, 1.0, 10.0, 6.0, 9.0, 6.0, 1.0, 2.0, 5.0, 1.0, 4.0, 3.0, 1.0, 10.0, 9.0, 6.0, 2.0, 5.0, 5.0, 10.0, 3.0, 2.0, 7.0, 1.0, 9.0, 4.0, 8.0, 10.0, 4.0, 3.0, 4.0, 3.0, 8.0, 3.0, 7.0, 8.0, 3.0, 7.0, 5.0, 4.0, 1.0, 7.0, 10.0, 4.0, 8.0, 9.0, 6.0, 10.0, 4.0, 5.0, 10.0, 3.0, 1.0, 10.0, 4.0, 9.0, 7.0, 1.0, 4.0, 4.0, 2.0, 9.0, 10.0, 7.0, 6.0, 8.0, 1.0, 4.0, 6.0, 2.0, 7.0, 5.0, 4.0, 1.0, 9.0, 2.0, 8.0, 3.0, 4.0, 4.0, 5.0, 10.0, 8.0, 3.0, 5.0, 4.0, 5.0, 2.0, 2.0, 6.0, 9.0, 6.0, 9.0, 5.0, 4.0, 10.0, 9.0, 2.0, 8.0, 4.0, 10.0, 2.0, 3.0, 1.0, 5.0, 4.0, 4.0, 1.0, 2.0, 2.0, 2.0, 6.0, 5.0, 3.0, 3.0, 3.0, 9.0, 4.0, 5.0, 8.0, 1.0, 2.0, 8.0, 10.0, 9.0, 6.0, 6.0, 3.0, 3.0, 5.0, 2.0, 7.0, 8.0, 4.0, 5.0, 10.0, 10.0, 2.0, 3.0, 8.0, 5.0, 2.0, 9.0, 7.0, 10.0, 1.0, 4.0, 7.0, 5.0, 9.0, 5.0, 10.0, 7.0, 8.0, 8.0, 7.0, 7.0, 10.0, 5.0, 3.0, 3.0, 10.0, 7.0, 7.0, 6.0, 3.0, 9.0, 10.0, 1.0, 9.0, 7.0, 2.0, 10.0, 7.0, 3.0, 8.0, 10.0, 10.0, 7.0, 1.0, 10.0, 5.0, 9.0, 7.0, 2.0, 7.0, 3.0, 8.0, 9.0, 10.0, 7.0, 1.0, 3.0, 5.0, 7.0, 9.0, 5.0, 6.0, 7.0, 9.0, 1.0, 9.0, 6.0, 2.0]
global b_x = 5
global d_y = [2.0, 5.0, 5.0, 5.0, 6.0, 2.0, 5.0, 9.0, 7.0, 7.0, 6.0, 6.0, 4.0, 3.0, 3.0, 4.0, 7.0, 2.0, 8.0, 4.0, 10.0, 8.0, 4.0, 5.0, 1.0, 7.0, 10.0, 2.0, 5.0, 3.0, 8.0, 7.0, 4.0, 7.0, 6.0, 9.0, 3.0, 4.0, 3.0, 10.0, 2.0, 7.0, 9.0, 7.0, 9.0, 8.0, 4.0, 8.0, 10.0, 5.0, 6.0, 6.0, 3.0, 3.0, 9.0, 4.0, 7.0, 2.0, 2.0, 9.0, 3.0, 2.0, 2.0, 8.0, 7.0, 7.0, 2.0, 8.0, 6.0, 7.0, 4.0, 9.0, 1.0, 9.0, 8.0, 7.0, 10.0, 2.0, 6.0, 7.0, 2.0, 10.0, 9.0, 8.0, 8.0, 8.0, 5.0, 1.0, 4.0, 4.0, 7.0, 9.0, 4.0, 4.0, 6.0, 8.0, 7.0, 9.0, 5.0, 8.0, 6.0, 5.0, 6.0, 4.0, 4.0, 2.0, 4.0, 10.0, 7.0, 2.0, 1.0, 2.0, 7.0, 2.0, 3.0, 4.0, 4.0, 8.0, 9.0, 3.0, 4.0, 9.0, 1.0, 2.0, 6.0, 3.0, 1.0, 5.0, 9.0, 8.0, 5.0, 7.0, 3.0, 2.0, 4.0, 7.0, 2.0, 1.0, 8.0, 9.0, 10.0, 2.0, 4.0, 8.0, 8.0, 5.0, 9.0, 1.0, 9.0, 8.0, 10.0, 10.0, 5.0, 10.0, 1.0, 7.0, 3.0, 7.0, 5.0, 1.0, 2.0, 9.0, 10.0, 1.0, 4.0, 7.0, 6.0, 4.0, 3.0, 3.0, 5.0, 9.0, 1.0, 4.0, 2.0, 9.0, 6.0, 1.0, 9.0, 6.0, 1.0, 8.0, 6.0, 4.0, 7.0, 3.0, 5.0, 4.0, 3.0, 10.0, 1.0, 6.0, 7.0, 8.0, 5.0, 8.0, 9.0, 5.0, 6.0, 4.0, 4.0, 2.0, 5.0, 8.0, 10.0, 7.0, 5.0, 9.0, 9.0, 3.0, 8.0, 2.0, 1.0, 3.0, 9.0, 9.0, 8.0, 4.0, 9.0, 6.0, 3.0, 6.0, 10.0, 1.0, 9.0, 5.0, 9.0, 1.0, 5.0, 3.0, 10.0, 4.0, 2.0, 6.0, 3.0, 7.0, 9.0, 8.0, 2.0, 1.0, 5.0, 4.0, 6.0, 5.0, 3.0, 4.0, 8.0, 4.0, 1.0, 9.0, 2.0, 5.0, 5.0, 3.0, 5.0, 4.0, 2.0, 3.0, 10.0, 2.0, 1.0, 5.0, 5.0, 4.0, 3.0, 2.0, 7.0, 2.0, 7.0, 1.0, 7.0, 5.0, 10.0, 7.0, 5.0, 9.0, 7.0, 1.0, 4.0, 9.0, 6.0, 3.0, 8.0, 8.0, 10.0, 2.0, 10.0, 1.0, 7.0, 4.0, 4.0, 8.0, 7.0, 10.0, 10.0, 9.0, 8.0, 9.0, 9.0, 8.0, 3.0, 9.0, 9.0, 7.0, 9.0, 3.0, 5.0, 4.0, 2.0, 2.0, 9.0, 1.0, 4.0, 4.0, 4.0, 5.0, 2.0, 4.0, 8.0, 9.0, 3.0, 9.0, 10.0, 9.0, 7.0, 2.0, 8.0, 10.0, 1.0]
global b_y = 10
global p = [0.735, 0.462, 0.899, 0.541, 0.004, 0.648, 0.888, 0.674, 0.939, 0.524, 0.252, 0.371, 0.007, 0.79, 0.703, 0.75, 0.595, 0.429, 0.889, 0.217, 0.6, 0.545, 0.09, 0.116, 0.973, 0.812, 0.756, 0.123, 0.693, 0.918, 0.652, 0.407, 0.31, 0.867, 0.735, 0.843, 0.694, 0.628, 0.606, 0.747, 0.562, 0.272, 0.12, 0.981, 0.946, 0.456, 0.32, 0.458, 0.495, 0.536, 0.836, 0.15, 0.949, 0.609, 0.319, 0.127, 0.683, 0.983, 0.754, 0.517, 0.634, 0.379, 0.981, 0.791, 0.94, 0.631, 0.627, 0.315, 0.601, 0.184, 0.788, 0.621, 0.172, 0.283, 0.027, 0.43, 0.111, 0.436, 0.18, 0.312, 0.075, 0.381, 0.098, 0.149, 0.61, 0.561, 0.705, 0.368, 0.674, 0.87, 0.9, 0.365, 0.738, 0.109, 0.234, 0.863, 0.497, 0.03, 0.36, 0.206, 0.544, 0.174, 0.945, 0.648, 0.342, 0.577, 0.476, 0.565, 0.78, 0.852, 0.859, 0.822, 0.052, 0.351, 0.96, 0.323, 0.454, 0.578, 0.27, 0.969, 0.536, 0.091, 0.484, 0.908, 0.171, 0.636, 0.033, 0.579, 0.746, 0.095, 0.924, 0.006, 0.389, 0.53, 0.11, 0.038, 0.692, 0.677, 0.62, 0.519, 0.266, 0.409, 0.296, 0.936, 0.239, 0.497, 0.436, 0.765, 0.214, 0.798, 0.537, 0.171, 0.943, 0.193, 0.371, 0.504, 0.608, 0.485, 0.949, 0.11, 0.03, 0.623, 0.706, 0.628, 0.126, 0.744, 0.728, 0.287, 0.59, 0.799, 0.426, 0.692, 0.977, 0.788, 0.951, 0.322, 0.639, 0.462, 0.71, 0.065, 0.408, 0.825, 0.118, 0.731, 0.961, 0.691, 0.84, 0.019, 0.664, 0.015, 0.88, 0.208, 0.149, 0.344, 0.575, 0.001, 0.461, 0.126, 0.597, 0.893, 0.589, 0.771, 0.783, 0.457, 0.098, 0.709, 0.685, 0.01, 0.02, 0.677, 0.245, 0.479, 0.313, 0.932, 0.018, 0.08, 0.271, 0.393, 0.556, 0.96, 0.57, 0.817, 0.376, 0.081, 0.803, 0.414, 0.426, 0.087, 0.11, 0.755, 0.571, 0.024, 0.72, 0.472, 0.412, 0.294, 0.806, 0.624, 0.255, 0.618, 0.995, 0.527, 0.815, 0.337, 0.144, 0.449, 0.8, 0.794, 0.305, 0.883, 0.705, 0.046, 0.187, 0.001, 0.295, 0.683, 0.35, 0.755, 0.331, 0.354, 0.361, 0.521, 0.598, 0.234, 0.727, 0.019, 0.963, 0.161, 0.879, 0.614, 0.998, 0.994, 0.849, 0.396, 0.387, 0.401, 0.493, 0.999, 0.574, 0.396, 0.592, 0.142, 0.503, 0.924, 0.44, 0.282, 0.072, 0.179, 0.476, 0.826, 0.888, 0.165, 0.177, 0.457, 0.511, 0.238, 0.592, 0.447, 0.985, 0.582, 0.589, 0.683, 0.66, 0.027, 0.666, 0.493, 0.312, 0.248, 0.805, 0.323, 0.738, 0.133, 0.5, 0.994, 0.119, 0.938, 0.524, 0.413, 0.662, 0.544, 0.505, 0.896, 0.276, 0.693, 0.699, 0.641, 0.982, 0.317, 0.396]
global q = [0.797, 0.92, 0.961, 0.878, 0.32, 0.713, 0.973, 0.745, 0.946, 0.584, 0.688, 0.439, 0.67, 0.848, 0.771, 0.861, 0.823, 0.576, 0.917, 0.476, 0.818, 0.68, 0.244, 0.961, 0.978, 0.864, 0.963, 0.447, 0.844, 0.925, 0.97, 0.599, 0.735, 0.929, 0.872, 0.984, 0.781, 0.653, 0.958, 0.903, 0.647, 0.982, 0.698, 0.984, 0.998, 0.477, 0.516, 0.572, 0.717, 0.907, 0.907, 0.81, 0.999, 0.609, 0.48, 0.836, 0.882, 0.991, 0.768, 0.542, 0.946, 0.734, 0.995, 0.91, 0.992, 0.827, 0.725, 0.948, 0.912, 0.36, 0.809, 0.874, 0.35, 0.471, 0.285, 0.531, 0.573, 0.771, 0.769, 0.886, 0.935, 0.506, 0.375, 0.439, 0.74, 0.672, 0.816, 0.689, 0.804, 0.924, 0.994, 0.459, 0.854, 0.165, 0.668, 0.888, 0.888, 0.586, 0.757, 0.819, 0.815, 0.465, 0.955, 0.89, 0.603, 0.666, 0.867, 0.663, 0.919, 0.868, 0.972, 0.944, 0.773, 0.657, 0.991, 0.637, 0.63, 0.674, 0.851, 0.972, 0.909, 0.783, 0.85, 0.956, 0.399, 0.836, 0.955, 0.652, 0.777, 0.441, 0.997, 0.736, 0.533, 0.821, 0.766, 0.83, 0.751, 0.712, 0.713, 0.904, 0.4, 0.83, 0.639, 0.96, 0.462, 0.584, 0.991, 0.914, 0.508, 0.878, 0.63, 0.64, 0.975, 0.363, 0.985, 0.552, 0.962, 0.637, 0.956, 0.378, 0.636, 0.742, 0.981, 0.967, 0.972, 0.864, 0.909, 0.566, 0.998, 0.946, 0.864, 0.715, 0.987, 0.808, 0.964, 0.476, 0.796, 0.929, 0.846, 0.869, 0.453, 0.943, 0.846, 0.984, 0.999, 0.932, 0.964, 0.303, 0.776, 0.658, 0.887, 0.72, 0.718, 0.95, 0.601, 0.45, 0.96, 0.61, 0.716, 0.991, 0.669, 0.95, 0.906, 0.863, 0.681, 0.92, 0.971, 0.897, 0.35, 0.776, 0.966, 0.485, 0.984, 0.962, 0.425, 0.5, 0.615, 0.793, 0.852, 0.996, 0.732, 0.917, 0.685, 0.351, 0.847, 0.59, 0.994, 0.962, 0.388, 0.972, 0.953, 0.75, 0.997, 0.634, 0.684, 0.557, 0.81, 0.962, 0.603, 0.912, 0.999, 0.7, 0.894, 0.585, 0.715, 0.662, 0.886, 0.881, 0.916, 0.936, 0.888, 0.216, 0.588, 0.525, 0.775, 0.87, 0.492, 0.928, 0.575, 0.879, 0.756, 0.935, 0.827, 0.884, 0.946, 0.553, 0.968, 0.603, 0.925, 0.978, 0.998, 0.998, 0.959, 0.425, 0.398, 0.973, 0.706, 0.999, 0.904, 0.398, 0.782, 0.596, 0.806, 0.967, 0.555, 0.635, 0.477, 0.248, 0.818, 0.935, 0.993, 0.36, 0.483, 0.952, 0.653, 0.3, 0.666, 0.95, 0.996, 0.993, 0.596, 0.759, 0.821, 0.534, 0.667, 0.744, 0.495, 0.336, 0.983, 0.752, 0.745, 0.522, 0.888, 0.999, 0.884, 0.952, 0.676, 0.507, 0.735, 0.889, 0.718, 0.997, 0.373, 0.974, 0.915, 0.672, 0.988, 0.7, 0.407]
global origin = 1
global destination = 60