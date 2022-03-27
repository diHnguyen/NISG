global arcs = [1 22; 1 38; 1 40; 1 47; 1 59; 2 7; 2 8; 2 21; 2 31; 2 37; 2 43; 2 49; 3 8; 3 16; 3 18; 3 24; 3 30; 3 39; 3 40; 4 37; 5 13; 5 16; 5 24; 5 46; 5 49; 5 52; 5 54; 6 4; 6 15; 6 28; 6 30; 6 42; 6 48; 7 34; 7 58; 8 26; 9 12; 9 14; 9 18; 9 32; 9 44; 10 27; 10 41; 10 49; 11 38; 11 44; 12 3; 12 48; 12 50; 13 8; 13 39; 13 48; 13 49; 13 54; 13 56; 14 22; 14 36; 14 44; 14 46; 14 53; 15 17; 15 36; 15 57; 16 38; 16 48; 16 55; 16 57; 17 23; 17 58; 18 2; 18 17; 18 22; 18 39; 18 42; 18 45; 18 60; 19 10; 19 23; 19 31; 19 58; 20 19; 20 24; 20 42; 20 47; 21 9; 21 17; 21 19; 21 47; 22 15; 22 41; 22 60; 23 11; 23 40; 23 44; 24 8; 24 18; 24 28; 24 37; 24 57; 24 59; 25 6; 25 8; 25 16; 25 36; 26 4; 26 28; 26 45; 26 53; 27 8; 27 22; 27 43; 27 53; 28 2; 28 17; 28 25; 28 26; 28 29; 28 32; 28 43; 29 11; 29 19; 29 57; 30 11; 30 12; 30 13; 30 18; 30 21; 30 32; 30 35; 30 57; 31 4; 31 10; 31 36; 31 37; 31 52; 32 2; 32 4; 32 9; 32 45; 32 56; 33 8; 33 23; 33 29; 33 31; 33 49; 34 2; 34 3; 34 11; 34 15; 34 21; 34 40; 34 45; 34 60; 35 7; 35 40; 36 2; 36 7; 36 23; 36 33; 36 35; 36 43; 36 46; 36 48; 37 14; 37 21; 37 35; 37 52; 38 18; 38 28; 38 36; 38 39; 38 48; 39 20; 39 22; 39 27; 39 35; 39 47; 39 52; 39 58; 40 6; 40 16; 40 18; 40 21; 40 26; 40 33; 40 34; 40 38; 40 43; 40 56; 41 13; 41 60; 42 5; 42 30; 42 31; 42 36; 42 49; 42 54; 42 57; 43 7; 43 14; 43 15; 43 30; 43 51; 43 54; 43 56; 44 10; 44 11; 44 12; 44 27; 44 36; 44 46; 44 55; 44 58; 45 4; 45 8; 45 12; 45 16; 45 26; 45 42; 45 46; 46 5; 46 8; 46 9; 46 14; 46 15; 46 17; 46 35; 46 47; 46 53; 46 56; 47 2; 47 5; 47 12; 47 15; 47 37; 47 46; 47 48; 47 51; 48 4; 48 10; 48 13; 48 27; 48 56; 48 57; 49 20; 49 21; 49 26; 49 40; 49 55; 50 3; 50 30; 50 51; 50 60; 51 21; 51 48; 52 13; 52 28; 52 36; 52 46; 52 48; 52 57; 53 3; 53 8; 53 14; 53 18; 53 29; 53 35; 53 38; 54 6; 54 8; 54 10; 54 25; 54 42; 54 43; 54 47; 55 24; 55 32; 55 35; 55 36; 55 37; 55 44; 55 47; 55 54; 55 56; 55 58; 55 59; 56 9; 56 16; 56 20; 56 42; 57 3; 57 19; 57 38; 58 9; 58 33; 58 38; 58 39; 58 42; 58 43; 59 15; 59 17; 59 18; 59 37; 59 45; 59 46; 59 52; 59 57]
global d_x = [6.0, 1.0, 10.0, 7.0, 8.0, 3.0, 1.0, 7.0, 7.0, 1.0, 9.0, 6.0, 6.0, 9.0, 8.0, 8.0, 4.0, 9.0, 5.0, 7.0, 3.0, 6.0, 8.0, 6.0, 5.0, 3.0, 8.0, 3.0, 5.0, 8.0, 3.0, 6.0, 2.0, 7.0, 9.0, 4.0, 4.0, 2.0, 7.0, 2.0, 8.0, 10.0, 4.0, 4.0, 5.0, 2.0, 1.0, 5.0, 8.0, 4.0, 1.0, 8.0, 2.0, 5.0, 2.0, 6.0, 6.0, 8.0, 7.0, 8.0, 5.0, 9.0, 4.0, 8.0, 2.0, 2.0, 10.0, 4.0, 4.0, 2.0, 1.0, 7.0, 9.0, 10.0, 3.0, 7.0, 4.0, 2.0, 2.0, 7.0, 3.0, 8.0, 7.0, 4.0, 7.0, 9.0, 8.0, 4.0, 10.0, 10.0, 7.0, 8.0, 6.0, 8.0, 3.0, 7.0, 10.0, 8.0, 4.0, 3.0, 7.0, 2.0, 6.0, 5.0, 2.0, 1.0, 9.0, 6.0, 3.0, 9.0, 5.0, 6.0, 5.0, 7.0, 3.0, 10.0, 1.0, 8.0, 8.0, 5.0, 4.0, 1.0, 1.0, 3.0, 10.0, 10.0, 2.0, 5.0, 6.0, 10.0, 4.0, 8.0, 5.0, 6.0, 4.0, 2.0, 4.0, 8.0, 10.0, 9.0, 7.0, 5.0, 3.0, 6.0, 5.0, 7.0, 5.0, 8.0, 2.0, 10.0, 6.0, 7.0, 7.0, 8.0, 1.0, 3.0, 8.0, 4.0, 4.0, 6.0, 3.0, 6.0, 2.0, 10.0, 8.0, 5.0, 4.0, 4.0, 10.0, 9.0, 2.0, 3.0, 9.0, 1.0, 9.0, 6.0, 8.0, 7.0, 6.0, 1.0, 1.0, 8.0, 1.0, 3.0, 7.0, 4.0, 2.0, 9.0, 4.0, 2.0, 9.0, 8.0, 6.0, 1.0, 3.0, 8.0, 4.0, 7.0, 10.0, 2.0, 4.0, 3.0, 9.0, 4.0, 7.0, 5.0, 2.0, 8.0, 6.0, 7.0, 7.0, 7.0, 1.0, 2.0, 4.0, 7.0, 9.0, 5.0, 2.0, 8.0, 2.0, 1.0, 5.0, 1.0, 1.0, 8.0, 1.0, 2.0, 5.0, 9.0, 5.0, 1.0, 5.0, 3.0, 3.0, 4.0, 4.0, 1.0, 3.0, 8.0, 3.0, 1.0, 8.0, 9.0, 1.0, 8.0, 4.0, 2.0, 9.0, 9.0, 6.0, 8.0, 1.0, 2.0, 10.0, 5.0, 7.0, 3.0, 7.0, 10.0, 8.0, 1.0, 8.0, 8.0, 4.0, 3.0, 1.0, 9.0, 2.0, 1.0, 9.0, 5.0, 3.0, 9.0, 2.0, 5.0, 4.0, 3.0, 2.0, 4.0, 6.0, 10.0, 8.0, 5.0, 3.0, 5.0, 9.0, 8.0, 8.0, 1.0, 6.0, 9.0, 6.0, 7.0, 6.0, 5.0, 2.0, 10.0, 7.0, 2.0, 5.0, 5.0, 10.0, 9.0, 6.0, 5.0, 3.0]
global b_x = 5
global d_y = [9.0, 8.0, 10.0, 4.0, 5.0, 2.0, 2.0, 4.0, 1.0, 6.0, 8.0, 7.0, 6.0, 10.0, 2.0, 1.0, 10.0, 4.0, 7.0, 9.0, 10.0, 1.0, 7.0, 7.0, 6.0, 8.0, 8.0, 3.0, 1.0, 5.0, 8.0, 8.0, 4.0, 6.0, 5.0, 9.0, 2.0, 2.0, 10.0, 10.0, 10.0, 4.0, 4.0, 6.0, 7.0, 2.0, 8.0, 6.0, 3.0, 10.0, 8.0, 9.0, 4.0, 3.0, 4.0, 8.0, 3.0, 7.0, 4.0, 10.0, 6.0, 7.0, 10.0, 5.0, 6.0, 5.0, 8.0, 8.0, 7.0, 5.0, 10.0, 6.0, 10.0, 6.0, 1.0, 2.0, 3.0, 6.0, 10.0, 2.0, 2.0, 8.0, 6.0, 8.0, 9.0, 9.0, 3.0, 10.0, 5.0, 1.0, 5.0, 4.0, 1.0, 5.0, 1.0, 1.0, 3.0, 6.0, 5.0, 2.0, 6.0, 9.0, 10.0, 10.0, 1.0, 2.0, 10.0, 1.0, 3.0, 3.0, 8.0, 10.0, 8.0, 8.0, 5.0, 1.0, 9.0, 1.0, 3.0, 10.0, 4.0, 3.0, 9.0, 8.0, 4.0, 8.0, 2.0, 9.0, 8.0, 7.0, 10.0, 4.0, 10.0, 9.0, 10.0, 4.0, 10.0, 5.0, 9.0, 5.0, 2.0, 3.0, 2.0, 6.0, 1.0, 4.0, 3.0, 9.0, 1.0, 5.0, 1.0, 10.0, 3.0, 2.0, 3.0, 5.0, 7.0, 2.0, 9.0, 7.0, 1.0, 6.0, 3.0, 6.0, 6.0, 4.0, 5.0, 8.0, 6.0, 9.0, 9.0, 3.0, 3.0, 4.0, 10.0, 5.0, 9.0, 4.0, 7.0, 7.0, 1.0, 9.0, 2.0, 1.0, 1.0, 1.0, 10.0, 1.0, 4.0, 1.0, 2.0, 7.0, 2.0, 3.0, 8.0, 10.0, 10.0, 3.0, 7.0, 5.0, 9.0, 6.0, 5.0, 5.0, 5.0, 2.0, 1.0, 7.0, 4.0, 7.0, 4.0, 9.0, 8.0, 8.0, 4.0, 2.0, 6.0, 3.0, 3.0, 4.0, 8.0, 8.0, 3.0, 7.0, 6.0, 1.0, 1.0, 6.0, 9.0, 1.0, 9.0, 9.0, 3.0, 2.0, 9.0, 3.0, 10.0, 7.0, 5.0, 5.0, 8.0, 1.0, 1.0, 10.0, 9.0, 8.0, 4.0, 8.0, 2.0, 5.0, 5.0, 10.0, 4.0, 1.0, 7.0, 9.0, 10.0, 1.0, 2.0, 3.0, 3.0, 1.0, 3.0, 7.0, 10.0, 8.0, 4.0, 6.0, 1.0, 8.0, 7.0, 8.0, 6.0, 7.0, 9.0, 1.0, 1.0, 10.0, 4.0, 5.0, 1.0, 10.0, 6.0, 3.0, 3.0, 8.0, 7.0, 1.0, 3.0, 1.0, 4.0, 3.0, 10.0, 5.0, 5.0, 9.0, 5.0, 2.0, 3.0, 9.0, 5.0, 2.0, 10.0, 6.0, 6.0, 3.0, 1.0]
global b_y = 10
global p = [0.603, 0.405, 0.464, 0.697, 0.766, 0.107, 0.614, 0.305, 0.933, 0.5, 0.836, 0.096, 0.779, 0.153, 0.955, 0.829, 0.213, 0.955, 0.661, 0.497, 0.954, 0.641, 0.189, 0.192, 0.302, 0.393, 0.726, 0.292, 0.405, 0.605, 0.84, 0.058, 0.467, 0.992, 0.815, 0.716, 0.126, 0.247, 0.084, 0.026, 0.35, 0.479, 0.827, 0.732, 0.222, 0.032, 0.362, 0.572, 0.92, 0.583, 0.989, 0.556, 0.14, 0.791, 0.309, 0.301, 0.619, 0.006, 0.523, 0.211, 0.667, 0.236, 0.212, 0.764, 0.097, 0.596, 0.155, 0.222, 0.584, 0.658, 0.722, 0.569, 0.055, 0.084, 0.873, 0.516, 0.048, 0.959, 0.634, 0.372, 0.102, 0.539, 0.243, 0.488, 0.403, 0.534, 0.629, 0.547, 0.247, 0.088, 0.262, 0.789, 0.022, 0.795, 0.272, 0.888, 0.637, 0.633, 0.889, 0.226, 0.29, 0.666, 0.758, 0.741, 0.535, 0.309, 0.715, 0.24, 0.3, 0.175, 0.395, 0.295, 0.834, 0.522, 0.755, 0.555, 0.624, 0.13, 0.294, 0.733, 0.359, 0.905, 0.771, 0.898, 0.621, 0.794, 0.614, 0.723, 0.926, 0.363, 0.77, 0.054, 0.338, 0.739, 0.794, 0.854, 0.546, 0.403, 0.909, 0.375, 0.797, 0.129, 0.3, 0.502, 0.65, 0.659, 0.25, 0.504, 0.247, 0.039, 0.784, 0.467, 0.042, 0.155, 0.538, 0.017, 0.154, 0.302, 0.277, 0.877, 0.48, 0.054, 0.817, 0.742, 0.524, 0.376, 0.814, 0.319, 0.157, 0.797, 0.445, 0.519, 0.95, 0.701, 0.522, 0.963, 0.08, 0.51, 0.822, 0.572, 0.209, 0.571, 0.684, 0.194, 0.584, 0.232, 0.562, 0.848, 0.738, 0.812, 0.157, 0.156, 0.945, 0.635, 0.425, 0.572, 0.84, 0.126, 0.479, 0.949, 0.66, 0.966, 0.887, 0.367, 0.75, 0.396, 0.25, 0.576, 0.43, 0.134, 0.125, 0.028, 0.385, 0.962, 0.073, 0.119, 0.397, 0.787, 0.951, 0.813, 0.201, 0.586, 0.87, 0.778, 0.683, 0.037, 0.259, 0.265, 0.873, 0.993, 0.577, 0.103, 0.803, 0.504, 0.894, 0.519, 0.076, 0.023, 0.823, 0.78, 0.847, 0.735, 0.62, 0.775, 0.863, 0.641, 0.523, 0.808, 0.521, 0.286, 0.586, 0.079, 0.469, 0.402, 0.407, 0.407, 0.935, 0.283, 0.481, 0.541, 0.445, 0.453, 0.945, 0.731, 0.324, 0.147, 0.632, 0.928, 0.443, 0.438, 0.338, 0.861, 0.012, 0.25, 0.666, 0.509, 0.246, 0.776, 0.056, 0.635, 0.93, 0.856, 0.465, 0.057, 0.044, 0.671, 0.102, 0.74, 0.64, 0.38, 0.217, 0.823, 0.62, 0.092, 0.462, 0.925, 0.359, 0.087, 0.064, 0.038, 0.168, 0.828, 0.343, 0.32, 0.06, 0.882, 0.952]
global q = [0.77, 0.882, 0.736, 0.977, 0.841, 0.408, 0.687, 0.609, 0.971, 0.53, 0.936, 0.997, 0.84, 0.575, 0.961, 0.911, 0.838, 0.998, 0.786, 0.543, 0.976, 0.93, 0.479, 0.925, 0.483, 0.924, 0.849, 0.586, 0.623, 0.744, 0.869, 0.422, 0.488, 0.995, 0.927, 0.747, 0.454, 0.547, 0.157, 0.56, 0.595, 0.935, 0.986, 0.953, 0.783, 0.078, 0.622, 0.885, 0.935, 0.811, 0.997, 0.611, 0.42, 0.993, 0.912, 0.879, 0.966, 0.969, 0.873, 0.794, 0.87, 0.69, 0.928, 0.935, 0.524, 0.96, 0.337, 0.545, 0.618, 0.812, 0.739, 0.647, 0.793, 0.207, 0.892, 0.832, 0.772, 0.997, 0.77, 0.677, 0.496, 0.576, 0.96, 0.551, 0.97, 0.606, 0.896, 0.601, 0.618, 0.411, 0.264, 0.818, 0.966, 0.882, 0.886, 0.931, 0.923, 0.841, 0.964, 0.48, 0.632, 0.945, 0.768, 0.878, 0.537, 0.831, 0.937, 0.367, 0.414, 0.882, 0.859, 0.734, 0.863, 0.566, 0.929, 0.619, 0.699, 0.893, 0.78, 0.793, 0.803, 0.979, 0.905, 0.923, 0.938, 0.827, 0.849, 0.866, 0.997, 0.837, 0.995, 0.594, 0.493, 0.948, 0.936, 0.993, 0.576, 0.616, 0.928, 0.594, 0.962, 0.758, 0.959, 0.558, 0.959, 0.675, 0.899, 0.88, 0.681, 0.756, 0.914, 0.555, 0.653, 0.709, 0.99, 0.044, 0.933, 0.78, 0.393, 0.972, 0.488, 0.852, 0.817, 0.816, 0.556, 0.387, 0.856, 0.948, 0.987, 0.943, 0.839, 0.727, 0.978, 0.906, 0.666, 0.984, 0.157, 0.742, 0.925, 0.782, 0.571, 0.745, 0.831, 0.339, 0.736, 0.563, 0.809, 0.854, 0.771, 0.988, 0.794, 0.621, 0.974, 0.888, 0.951, 0.777, 0.897, 0.147, 0.661, 0.953, 0.976, 0.996, 0.974, 0.8, 0.833, 0.967, 0.326, 0.835, 0.855, 0.43, 0.276, 0.324, 0.929, 0.974, 0.564, 0.559, 0.403, 0.872, 0.962, 0.873, 0.561, 0.875, 0.874, 0.823, 0.786, 0.937, 0.568, 0.381, 0.908, 0.996, 0.799, 0.951, 0.968, 0.573, 0.973, 0.554, 0.164, 0.475, 0.984, 0.902, 0.946, 0.788, 0.675, 0.871, 0.994, 0.8, 0.6, 0.978, 0.943, 0.317, 0.85, 0.334, 0.798, 0.993, 0.523, 0.747, 0.95, 0.603, 0.693, 0.69, 0.761, 0.486, 0.951, 0.777, 0.698, 0.894, 0.709, 0.975, 0.962, 0.885, 0.476, 0.925, 0.179, 0.434, 0.852, 0.77, 0.771, 0.946, 0.724, 0.863, 0.961, 0.869, 0.819, 0.361, 0.594, 0.712, 0.492, 0.785, 0.857, 0.459, 0.835, 0.897, 0.88, 0.883, 0.642, 0.998, 0.943, 0.811, 0.768, 0.592, 0.89, 0.961, 0.876, 0.473, 0.903, 0.944, 0.999]
global origin = 1
global destination = 60