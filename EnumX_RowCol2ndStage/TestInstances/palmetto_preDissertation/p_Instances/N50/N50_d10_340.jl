global arcs = [1 5; 1 13; 1 32; 1 35; 1 42; 1 48; 1 49; 1 50; 2 6; 2 8; 2 36; 2 39; 3 21; 3 28; 3 35; 4 3; 4 13; 4 24; 4 27; 4 31; 4 47; 4 49; 5 13; 5 20; 5 29; 6 7; 6 8; 6 19; 6 24; 6 29; 6 30; 6 50; 7 2; 7 18; 7 21; 7 24; 7 35; 8 6; 8 10; 8 11; 8 14; 8 31; 8 45; 9 17; 9 18; 9 19; 9 32; 10 3; 10 19; 10 30; 10 33; 10 36; 10 41; 10 43; 10 48; 11 6; 11 20; 11 39; 11 46; 11 47; 12 2; 12 9; 12 18; 12 20; 12 35; 12 36; 12 39; 12 49; 13 8; 13 9; 13 41; 14 12; 14 23; 14 32; 14 33; 14 40; 14 41; 14 45; 15 16; 15 17; 15 20; 15 24; 15 27; 15 29; 15 37; 15 43; 16 15; 16 50; 17 26; 17 43; 18 2; 18 10; 18 12; 18 13; 18 40; 18 43; 18 44; 19 2; 19 5; 19 29; 19 36; 19 39; 19 41; 19 43; 20 14; 20 16; 20 19; 20 25; 20 26; 21 9; 21 11; 21 15; 21 18; 21 27; 22 10; 22 14; 22 17; 22 19; 22 27; 22 36; 22 45; 23 7; 23 10; 23 11; 23 13; 23 17; 23 38; 23 43; 24 2; 24 9; 24 23; 24 31; 24 37; 24 41; 25 3; 25 14; 25 38; 26 4; 26 5; 26 8; 26 19; 26 20; 26 22; 27 28; 27 31; 27 33; 27 42; 28 2; 28 3; 28 5; 28 29; 28 44; 28 45; 29 7; 29 8; 29 14; 29 16; 29 26; 29 27; 29 31; 30 21; 30 32; 30 33; 30 48; 31 9; 31 16; 31 19; 32 4; 32 5; 32 11; 32 12; 32 21; 32 29; 32 39; 32 44; 33 24; 33 36; 33 38; 34 12; 34 16; 34 18; 34 31; 34 35; 35 3; 35 34; 35 37; 35 50; 36 6; 36 12; 36 24; 36 40; 37 3; 37 45; 37 47; 38 10; 38 18; 38 32; 39 6; 39 7; 39 17; 39 22; 39 36; 39 50; 40 13; 40 15; 40 30; 40 38; 41 13; 41 26; 41 39; 42 2; 42 7; 42 9; 42 16; 42 21; 42 36; 42 38; 43 8; 43 11; 43 18; 43 31; 43 36; 43 38; 43 44; 44 23; 44 47; 45 4; 45 36; 45 47; 46 15; 46 22; 46 31; 46 35; 46 39; 46 48; 47 16; 47 25; 48 7; 48 11; 48 16; 48 22; 48 25; 48 33; 48 43; 48 45; 49 44]
global d_x = [6.0, 4.0, 3.0, 10.0, 4.0, 6.0, 2.0, 7.0, 7.0, 8.0, 8.0, 8.0, 7.0, 3.0, 8.0, 8.0, 2.0, 6.0, 10.0, 3.0, 5.0, 1.0, 9.0, 2.0, 5.0, 1.0, 7.0, 5.0, 8.0, 2.0, 10.0, 8.0, 10.0, 8.0, 1.0, 3.0, 3.0, 10.0, 1.0, 3.0, 8.0, 7.0, 10.0, 6.0, 4.0, 8.0, 3.0, 10.0, 1.0, 1.0, 4.0, 5.0, 8.0, 2.0, 2.0, 3.0, 10.0, 2.0, 3.0, 1.0, 2.0, 4.0, 7.0, 8.0, 6.0, 3.0, 3.0, 8.0, 2.0, 9.0, 6.0, 3.0, 10.0, 6.0, 2.0, 6.0, 6.0, 2.0, 3.0, 5.0, 4.0, 9.0, 5.0, 3.0, 9.0, 5.0, 5.0, 5.0, 10.0, 5.0, 4.0, 5.0, 8.0, 7.0, 7.0, 4.0, 8.0, 7.0, 6.0, 4.0, 4.0, 3.0, 5.0, 9.0, 7.0, 6.0, 9.0, 9.0, 1.0, 3.0, 2.0, 4.0, 5.0, 1.0, 6.0, 2.0, 10.0, 10.0, 7.0, 2.0, 5.0, 2.0, 3.0, 6.0, 10.0, 10.0, 4.0, 4.0, 10.0, 8.0, 3.0, 9.0, 2.0, 4.0, 1.0, 6.0, 1.0, 1.0, 9.0, 8.0, 4.0, 8.0, 6.0, 9.0, 8.0, 2.0, 4.0, 7.0, 6.0, 5.0, 1.0, 1.0, 3.0, 2.0, 4.0, 2.0, 7.0, 10.0, 9.0, 10.0, 7.0, 6.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 5.0, 10.0, 7.0, 7.0, 9.0, 3.0, 7.0, 6.0, 5.0, 5.0, 9.0, 8.0, 2.0, 3.0, 9.0, 5.0, 10.0, 2.0, 4.0, 4.0, 1.0, 7.0, 2.0, 5.0, 8.0, 6.0, 7.0, 1.0, 2.0, 7.0, 7.0, 2.0, 2.0, 8.0, 3.0, 10.0, 6.0, 2.0, 10.0, 9.0, 6.0, 7.0, 7.0, 4.0, 3.0, 8.0, 5.0, 8.0, 7.0, 6.0, 1.0, 2.0, 6.0, 10.0, 4.0, 8.0, 5.0, 9.0, 3.0, 6.0, 7.0, 5.0, 3.0, 8.0, 3.0, 1.0, 9.0, 9.0, 8.0, 1.0, 5.0, 8.0, 3.0, 2.0, 9.0, 10.0, 5.0]
global b_x = 5
global d_y = [2.0, 9.0, 10.0, 2.0, 10.0, 10.0, 10.0, 6.0, 3.0, 5.0, 1.0, 6.0, 7.0, 10.0, 6.0, 10.0, 9.0, 9.0, 9.0, 8.0, 5.0, 7.0, 6.0, 10.0, 3.0, 3.0, 7.0, 2.0, 2.0, 6.0, 5.0, 2.0, 9.0, 4.0, 3.0, 10.0, 3.0, 10.0, 3.0, 6.0, 3.0, 3.0, 5.0, 2.0, 2.0, 1.0, 7.0, 3.0, 5.0, 2.0, 4.0, 2.0, 3.0, 2.0, 5.0, 5.0, 6.0, 7.0, 7.0, 1.0, 8.0, 6.0, 5.0, 8.0, 9.0, 5.0, 8.0, 10.0, 5.0, 8.0, 6.0, 9.0, 1.0, 9.0, 4.0, 7.0, 7.0, 5.0, 9.0, 2.0, 3.0, 6.0, 10.0, 8.0, 10.0, 5.0, 7.0, 1.0, 9.0, 2.0, 3.0, 5.0, 6.0, 10.0, 5.0, 10.0, 5.0, 5.0, 6.0, 7.0, 7.0, 10.0, 4.0, 7.0, 7.0, 2.0, 1.0, 2.0, 7.0, 4.0, 10.0, 2.0, 7.0, 4.0, 2.0, 3.0, 10.0, 10.0, 2.0, 2.0, 1.0, 9.0, 4.0, 3.0, 5.0, 6.0, 4.0, 10.0, 7.0, 8.0, 5.0, 2.0, 6.0, 2.0, 9.0, 4.0, 2.0, 7.0, 2.0, 2.0, 2.0, 5.0, 2.0, 10.0, 7.0, 4.0, 10.0, 1.0, 6.0, 4.0, 5.0, 2.0, 2.0, 2.0, 9.0, 10.0, 10.0, 9.0, 3.0, 6.0, 8.0, 8.0, 8.0, 5.0, 7.0, 5.0, 10.0, 4.0, 3.0, 4.0, 7.0, 8.0, 8.0, 1.0, 8.0, 2.0, 1.0, 8.0, 10.0, 8.0, 5.0, 4.0, 10.0, 4.0, 10.0, 10.0, 9.0, 9.0, 6.0, 4.0, 1.0, 8.0, 9.0, 1.0, 8.0, 10.0, 5.0, 4.0, 9.0, 10.0, 8.0, 4.0, 10.0, 6.0, 4.0, 2.0, 2.0, 1.0, 5.0, 6.0, 3.0, 5.0, 2.0, 2.0, 1.0, 8.0, 3.0, 6.0, 10.0, 4.0, 5.0, 10.0, 6.0, 5.0, 8.0, 6.0, 8.0, 7.0, 7.0, 8.0, 9.0, 8.0, 3.0, 4.0, 2.0, 2.0, 5.0, 1.0, 4.0, 5.0, 5.0, 8.0, 5.0, 6.0, 7.0, 10.0]
global b_y = 10
global p = [0.863, 0.838, 0.523, 0.046, 0.987, 0.873, 0.514, 0.316, 0.184, 0.129, 0.791, 0.85, 0.399, 0.039, 0.85, 0.016, 0.778, 0.097, 0.809, 0.4, 0.756, 0.346, 0.766, 0.663, 0.569, 0.05, 0.397, 0.757, 0.389, 0.38, 0.703, 0.421, 0.772, 0.77, 0.332, 0.515, 0.409, 0.774, 0.932, 0.347, 0.832, 0.552, 0.562, 0.194, 0.296, 0.118, 0.35, 0.278, 0.819, 0.559, 0.866, 0.153, 0.2, 0.285, 0.163, 0.316, 0.028, 0.683, 0.54, 0.137, 0.544, 0.989, 0.564, 0.925, 0.296, 0.791, 0.575, 0.439, 0.231, 0.178, 0.824, 0.556, 0.615, 0.542, 0.392, 0.48, 0.673, 0.617, 0.527, 0.707, 0.871, 0.051, 0.42, 0.935, 0.48, 0.536, 0.782, 0.825, 0.906, 0.389, 0.934, 0.506, 0.146, 0.147, 0.09, 0.514, 0.836, 0.263, 0.242, 0.215, 0.332, 0.27, 0.739, 0.837, 0.174, 0.951, 0.746, 0.066, 0.792, 0.765, 0.067, 0.865, 0.91, 0.465, 0.61, 0.023, 0.807, 0.785, 0.676, 0.432, 0.716, 0.318, 0.122, 0.969, 0.289, 0.87, 0.759, 0.323, 0.238, 0.948, 0.606, 0.103, 0.188, 0.07, 0.385, 0.652, 0.715, 0.426, 0.135, 0.274, 0.897, 0.294, 0.498, 0.92, 0.801, 0.384, 0.354, 0.162, 0.069, 0.18, 0.847, 0.463, 0.569, 0.473, 0.085, 0.506, 0.054, 0.895, 0.245, 0.284, 0.203, 0.459, 0.453, 0.701, 0.158, 0.559, 0.769, 0.628, 0.107, 0.448, 0.321, 0.304, 0.095, 0.083, 0.989, 0.297, 0.751, 0.274, 0.696, 0.466, 0.004, 0.999, 0.218, 0.032, 0.207, 0.08, 0.536, 0.83, 0.484, 0.028, 0.764, 0.845, 0.982, 0.911, 0.687, 0.821, 0.664, 0.795, 0.177, 0.131, 0.313, 0.812, 0.34, 0.15, 0.752, 0.58, 0.395, 0.73, 0.445, 0.342, 0.562, 0.872, 0.5, 0.272, 0.362, 0.846, 0.577, 0.758, 0.458, 0.467, 0.783, 0.266, 0.516, 0.391, 0.826, 0.591, 0.652, 0.519, 0.745, 0.29, 0.006, 0.797, 0.4, 0.941, 0.3, 0.771, 0.692, 0.553, 0.914, 0.135, 0.264, 0.27, 0.797, 0.409, 0.426, 0.491]
global q = [0.985, 0.938, 0.789, 0.914, 0.987, 0.967, 0.608, 0.855, 0.221, 0.517, 0.907, 0.955, 0.597, 0.802, 0.879, 0.563, 0.853, 0.508, 0.813, 0.48, 0.913, 0.954, 0.869, 0.86, 0.749, 0.509, 0.61, 0.766, 0.984, 0.413, 0.85, 0.72, 0.821, 0.788, 0.943, 0.966, 0.409, 0.936, 0.954, 0.645, 0.882, 0.606, 0.713, 0.852, 0.509, 0.167, 0.883, 0.899, 0.85, 0.717, 0.993, 0.859, 0.504, 0.893, 0.201, 0.808, 0.139, 0.947, 0.681, 0.465, 0.931, 0.993, 0.987, 0.929, 0.962, 0.858, 0.787, 0.617, 0.348, 0.559, 0.913, 0.948, 0.755, 0.548, 0.653, 0.838, 0.988, 0.666, 0.822, 0.88, 0.882, 0.705, 0.563, 0.955, 0.913, 0.723, 0.846, 0.899, 0.913, 0.852, 0.989, 0.95, 0.929, 0.272, 0.194, 0.641, 0.875, 0.422, 0.265, 0.409, 0.47, 0.874, 0.961, 0.904, 0.392, 0.954, 0.806, 0.699, 0.819, 0.861, 0.608, 0.989, 0.917, 0.917, 0.888, 0.703, 0.822, 0.904, 0.715, 0.475, 0.903, 0.821, 0.729, 0.993, 0.695, 0.982, 0.876, 0.469, 0.863, 0.978, 0.862, 0.493, 0.755, 0.815, 0.624, 0.978, 0.998, 0.859, 0.478, 0.532, 0.941, 0.717, 0.525, 0.936, 0.973, 0.459, 0.425, 0.315, 0.8, 0.834, 0.953, 0.956, 0.715, 0.972, 0.68, 0.743, 0.633, 0.955, 0.829, 0.437, 0.308, 0.518, 0.515, 0.844, 0.524, 0.979, 0.931, 0.79, 0.888, 0.949, 0.743, 0.381, 0.613, 0.289, 0.994, 0.99, 0.784, 0.817, 0.737, 0.929, 0.825, 0.999, 0.877, 0.207, 0.816, 0.62, 0.984, 0.943, 0.853, 0.598, 0.859, 0.879, 0.987, 0.953, 0.922, 0.992, 0.959, 0.852, 0.229, 0.167, 0.464, 0.975, 0.414, 0.194, 0.783, 0.804, 0.443, 0.902, 0.556, 0.599, 0.952, 0.957, 0.886, 0.845, 0.785, 0.859, 0.833, 0.989, 0.48, 0.711, 0.925, 0.667, 0.66, 0.712, 0.84, 0.814, 0.773, 0.839, 0.77, 0.58, 0.705, 0.998, 0.843, 0.942, 0.631, 0.965, 0.761, 0.954, 0.916, 0.793, 0.281, 0.599, 0.837, 0.467, 0.893, 0.923]
global origin = 1
global destination = 50