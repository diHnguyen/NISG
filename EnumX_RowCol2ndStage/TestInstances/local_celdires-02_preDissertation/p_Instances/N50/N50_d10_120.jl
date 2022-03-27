global arcs = [1 2; 1 7; 1 8; 1 14; 1 23; 1 26; 1 31; 1 32; 2 4; 2 6; 2 29; 2 41; 2 50; 3 6; 3 40; 3 42; 4 29; 4 30; 4 50; 5 16; 5 30; 5 33; 5 34; 5 48; 5 50; 6 5; 6 8; 6 15; 6 20; 6 23; 6 33; 6 44; 6 50; 7 4; 7 9; 7 25; 7 33; 7 34; 7 45; 8 16; 9 10; 9 35; 9 37; 10 5; 10 6; 10 32; 10 38; 10 42; 11 12; 11 38; 11 46; 12 10; 12 17; 12 20; 12 23; 12 24; 12 27; 12 42; 12 45; 12 49; 13 14; 13 28; 13 30; 13 41; 14 30; 14 36; 14 46; 15 8; 15 30; 15 42; 16 14; 16 22; 16 27; 16 39; 16 42; 17 9; 17 11; 17 13; 17 37; 18 14; 18 15; 18 17; 18 30; 18 31; 19 14; 19 38; 20 40; 21 10; 21 12; 21 23; 21 34; 21 42; 21 43; 22 11; 22 18; 22 23; 22 33; 22 46; 22 47; 23 2; 23 14; 23 20; 23 30; 23 36; 23 50; 24 7; 24 10; 24 31; 25 11; 25 41; 25 44; 25 47; 26 2; 26 3; 26 7; 26 27; 26 28; 26 50; 27 22; 27 28; 27 37; 27 40; 27 46; 28 17; 28 19; 28 23; 28 26; 28 27; 28 38; 29 10; 29 19; 29 20; 29 22; 29 33; 29 38; 29 45; 29 47; 29 48; 29 49; 30 2; 30 20; 30 26; 30 28; 30 32; 30 43; 31 2; 31 7; 31 16; 31 29; 31 42; 31 47; 32 10; 32 22; 32 43; 33 4; 33 11; 33 12; 33 14; 33 20; 33 26; 33 48; 34 23; 34 43; 35 6; 35 17; 35 46; 36 8; 36 10; 36 29; 36 40; 37 28; 37 38; 38 6; 38 14; 38 16; 38 20; 38 44; 39 6; 39 10; 39 21; 40 4; 40 7; 40 13; 40 14; 40 24; 40 27; 40 43; 41 6; 41 7; 41 13; 41 16; 41 31; 41 35; 41 46; 42 17; 42 30; 42 33; 42 36; 42 37; 43 4; 43 13; 43 18; 43 27; 43 33; 43 40; 43 46; 43 48; 44 31; 44 39; 44 43; 44 47; 44 48; 45 14; 45 15; 45 20; 45 26; 45 28; 45 29; 45 34; 45 39; 45 44; 46 22; 47 15; 47 16; 47 20; 47 25; 47 26; 47 36; 48 2; 48 5; 48 39; 49 20; 49 25; 49 50]
global d_x = [7.0, 5.0, 7.0, 2.0, 1.0, 2.0, 7.0, 1.0, 7.0, 7.0, 10.0, 10.0, 8.0, 1.0, 9.0, 9.0, 8.0, 9.0, 5.0, 8.0, 1.0, 2.0, 5.0, 3.0, 5.0, 6.0, 4.0, 4.0, 4.0, 10.0, 7.0, 5.0, 7.0, 8.0, 2.0, 9.0, 9.0, 3.0, 2.0, 3.0, 7.0, 10.0, 4.0, 6.0, 6.0, 3.0, 9.0, 4.0, 7.0, 6.0, 9.0, 7.0, 10.0, 10.0, 9.0, 1.0, 6.0, 3.0, 3.0, 10.0, 10.0, 5.0, 7.0, 10.0, 8.0, 7.0, 2.0, 10.0, 9.0, 5.0, 9.0, 3.0, 2.0, 8.0, 2.0, 9.0, 1.0, 10.0, 10.0, 3.0, 8.0, 7.0, 2.0, 10.0, 7.0, 10.0, 9.0, 2.0, 6.0, 8.0, 5.0, 2.0, 4.0, 6.0, 3.0, 4.0, 6.0, 2.0, 7.0, 6.0, 5.0, 5.0, 2.0, 8.0, 1.0, 6.0, 7.0, 10.0, 2.0, 6.0, 8.0, 8.0, 1.0, 2.0, 10.0, 4.0, 6.0, 4.0, 10.0, 7.0, 1.0, 8.0, 9.0, 5.0, 7.0, 9.0, 2.0, 8.0, 3.0, 7.0, 2.0, 3.0, 10.0, 6.0, 5.0, 2.0, 3.0, 1.0, 9.0, 9.0, 10.0, 10.0, 2.0, 9.0, 6.0, 6.0, 6.0, 9.0, 8.0, 1.0, 10.0, 10.0, 5.0, 4.0, 2.0, 6.0, 1.0, 1.0, 6.0, 2.0, 6.0, 7.0, 2.0, 1.0, 3.0, 4.0, 7.0, 2.0, 5.0, 2.0, 6.0, 2.0, 6.0, 5.0, 4.0, 5.0, 9.0, 1.0, 10.0, 9.0, 9.0, 6.0, 8.0, 9.0, 7.0, 8.0, 7.0, 8.0, 1.0, 3.0, 1.0, 3.0, 5.0, 1.0, 6.0, 6.0, 10.0, 5.0, 9.0, 1.0, 9.0, 1.0, 7.0, 3.0, 9.0, 3.0, 2.0, 8.0, 8.0, 5.0, 9.0, 5.0, 3.0, 7.0, 10.0, 6.0, 2.0, 1.0, 10.0, 6.0, 4.0, 9.0, 6.0, 6.0, 5.0, 3.0, 7.0, 9.0, 2.0, 9.0, 5.0, 3.0, 1.0, 8.0]
global b_x = 5
global d_y = [8.0, 3.0, 1.0, 4.0, 7.0, 9.0, 5.0, 1.0, 5.0, 3.0, 5.0, 7.0, 7.0, 8.0, 1.0, 5.0, 1.0, 3.0, 7.0, 10.0, 5.0, 6.0, 6.0, 2.0, 8.0, 9.0, 7.0, 3.0, 4.0, 3.0, 4.0, 3.0, 1.0, 5.0, 3.0, 6.0, 7.0, 10.0, 4.0, 10.0, 3.0, 4.0, 6.0, 7.0, 5.0, 4.0, 10.0, 10.0, 2.0, 8.0, 2.0, 4.0, 6.0, 3.0, 10.0, 10.0, 7.0, 9.0, 8.0, 5.0, 3.0, 9.0, 9.0, 4.0, 9.0, 2.0, 1.0, 1.0, 4.0, 7.0, 5.0, 1.0, 8.0, 8.0, 7.0, 3.0, 6.0, 5.0, 7.0, 2.0, 1.0, 6.0, 3.0, 5.0, 9.0, 3.0, 6.0, 5.0, 9.0, 7.0, 7.0, 9.0, 1.0, 7.0, 4.0, 2.0, 2.0, 1.0, 7.0, 7.0, 10.0, 2.0, 6.0, 8.0, 1.0, 8.0, 6.0, 2.0, 2.0, 6.0, 3.0, 2.0, 6.0, 2.0, 2.0, 4.0, 7.0, 2.0, 4.0, 4.0, 8.0, 3.0, 8.0, 2.0, 4.0, 3.0, 4.0, 1.0, 8.0, 8.0, 9.0, 3.0, 4.0, 7.0, 4.0, 5.0, 10.0, 9.0, 8.0, 5.0, 10.0, 9.0, 10.0, 9.0, 10.0, 10.0, 4.0, 9.0, 1.0, 7.0, 4.0, 4.0, 4.0, 6.0, 4.0, 5.0, 6.0, 7.0, 6.0, 10.0, 2.0, 8.0, 10.0, 9.0, 10.0, 9.0, 4.0, 2.0, 4.0, 3.0, 9.0, 7.0, 7.0, 9.0, 7.0, 10.0, 8.0, 9.0, 7.0, 1.0, 7.0, 10.0, 8.0, 2.0, 9.0, 1.0, 2.0, 8.0, 1.0, 2.0, 10.0, 9.0, 5.0, 7.0, 1.0, 8.0, 1.0, 4.0, 4.0, 7.0, 4.0, 4.0, 7.0, 9.0, 7.0, 6.0, 7.0, 5.0, 9.0, 9.0, 1.0, 3.0, 3.0, 3.0, 5.0, 6.0, 7.0, 10.0, 6.0, 7.0, 6.0, 8.0, 4.0, 8.0, 7.0, 6.0, 2.0, 7.0, 4.0, 6.0, 8.0, 7.0, 1.0, 1.0]
global b_y = 10
global p = [0.639, 0.596, 0.035, 0.343, 0.456, 0.442, 0.498, 0.182, 0.835, 0.76, 0.548, 0.539, 0.393, 0.02, 0.952, 0.49, 0.142, 0.729, 0.286, 0.406, 0.107, 0.68, 0.381, 0.673, 0.225, 0.626, 0.897, 0.103, 0.912, 0.296, 0.562, 0.331, 0.226, 0.262, 0.249, 0.317, 0.58, 0.97, 0.078, 0.186, 0.261, 0.327, 0.218, 0.025, 0.632, 0.999, 0.98, 0.667, 0.44, 0.502, 0.206, 0.581, 0.684, 0.478, 0.672, 0.875, 0.264, 0.624, 0.101, 0.285, 0.923, 0.541, 0.64, 0.511, 0.682, 0.656, 0.28, 0.293, 0.117, 0.737, 0.742, 0.814, 0.27, 0.137, 0.835, 0.22, 0.097, 0.228, 0.011, 0.594, 0.085, 0.85, 0.393, 0.434, 0.131, 0.212, 0.164, 0.374, 0.284, 0.833, 0.672, 0.267, 0.928, 0.092, 0.64, 0.91, 0.473, 0.819, 0.374, 0.221, 0.297, 0.146, 0.742, 0.329, 0.8, 0.12, 0.491, 0.116, 0.823, 0.358, 0.883, 0.398, 0.112, 0.029, 0.899, 0.791, 0.865, 0.828, 0.278, 0.148, 0.076, 0.193, 0.032, 0.13, 0.45, 0.365, 0.917, 0.857, 0.411, 0.823, 0.73, 0.753, 0.075, 0.425, 0.714, 0.916, 0.105, 0.922, 0.526, 0.345, 0.17, 0.411, 0.741, 0.485, 0.034, 0.619, 0.573, 0.583, 0.385, 0.942, 0.034, 0.045, 0.687, 0.687, 0.369, 0.925, 0.095, 0.55, 0.339, 0.834, 0.701, 0.188, 0.434, 0.295, 0.856, 0.033, 0.933, 0.157, 0.466, 0.408, 0.851, 0.396, 0.605, 0.448, 0.322, 0.577, 0.972, 0.833, 0.797, 0.426, 0.183, 0.486, 0.163, 0.48, 0.439, 0.662, 0.039, 0.658, 0.522, 0.783, 0.73, 0.833, 0.105, 0.448, 0.387, 0.813, 0.071, 0.926, 0.149, 0.357, 0.315, 0.587, 0.733, 0.728, 0.456, 0.128, 0.938, 0.259, 0.317, 0.22, 0.565, 0.422, 0.187, 0.485, 0.7, 0.597, 0.717, 0.699, 0.763, 0.278, 0.093, 0.339, 0.91, 0.843, 0.646, 0.804, 0.276, 0.023, 0.117, 0.958, 0.743, 0.289, 0.211, 0.1]
global q = [0.966, 0.995, 0.128, 0.514, 0.832, 0.915, 0.948, 0.699, 0.93, 0.896, 0.972, 0.679, 0.741, 0.138, 0.956, 0.68, 0.552, 0.902, 0.3, 0.463, 0.702, 0.949, 0.741, 0.858, 0.625, 0.784, 0.987, 0.771, 0.98, 0.454, 0.964, 0.854, 0.333, 0.868, 0.378, 0.834, 0.832, 0.982, 0.202, 0.91, 0.915, 0.681, 0.478, 0.274, 0.717, 0.999, 0.982, 0.895, 0.836, 0.527, 0.222, 0.919, 0.92, 0.512, 0.677, 0.953, 0.477, 0.835, 0.224, 0.855, 0.995, 0.654, 0.793, 0.557, 0.827, 0.821, 0.304, 0.337, 0.249, 0.763, 0.778, 0.854, 0.95, 0.824, 0.855, 0.802, 0.522, 0.268, 0.411, 0.95, 0.436, 0.991, 0.936, 0.569, 0.939, 0.242, 0.321, 0.621, 0.966, 0.933, 0.684, 0.433, 0.963, 0.679, 0.663, 0.961, 0.572, 0.912, 0.474, 0.375, 0.71, 0.257, 0.864, 0.892, 0.889, 0.736, 0.868, 0.487, 0.85, 0.7, 0.943, 0.481, 0.457, 0.1, 0.932, 0.893, 0.901, 0.998, 0.766, 0.369, 0.596, 0.668, 0.91, 0.719, 0.901, 0.644, 0.962, 0.953, 0.792, 0.919, 0.806, 0.969, 0.573, 0.576, 0.836, 0.992, 0.765, 0.923, 0.576, 0.735, 0.724, 0.455, 0.795, 0.571, 0.934, 0.785, 0.948, 0.717, 0.887, 0.972, 0.592, 0.74, 0.791, 0.887, 0.658, 0.963, 0.945, 0.878, 0.382, 0.836, 0.714, 0.521, 0.949, 0.669, 0.901, 0.554, 0.945, 0.602, 0.997, 0.856, 0.946, 0.826, 0.875, 0.621, 0.543, 0.957, 0.982, 0.895, 0.935, 0.798, 0.539, 0.723, 0.784, 0.67, 0.888, 0.789, 0.743, 0.662, 0.623, 0.845, 0.98, 0.952, 0.86, 0.575, 0.998, 0.914, 0.636, 0.998, 0.672, 0.991, 0.616, 0.617, 0.85, 0.973, 0.717, 0.281, 0.993, 0.708, 0.657, 0.754, 0.84, 0.999, 0.275, 0.552, 0.941, 0.934, 0.918, 0.706, 0.954, 0.982, 0.331, 0.708, 0.969, 0.975, 0.98, 0.993, 0.713, 0.813, 0.977, 0.973, 0.889, 0.352, 0.538, 0.437]
global origin = 1
global destination = 50