global arcs = [1 2; 1 13; 1 27; 1 31; 1 37; 2 9; 2 22; 2 28; 2 29; 2 30; 2 35; 2 44; 2 46; 3 5; 3 23; 3 28; 3 29; 3 36; 3 44; 4 21; 5 9; 5 16; 5 32; 6 14; 6 25; 6 29; 6 30; 6 31; 6 41; 6 46; 6 49; 7 3; 7 9; 7 23; 7 26; 7 29; 7 39; 8 2; 8 11; 8 15; 8 20; 8 21; 8 27; 8 37; 8 44; 9 4; 9 11; 9 18; 9 25; 9 39; 9 41; 9 46; 9 50; 10 15; 10 23; 10 37; 11 5; 11 46; 12 37; 12 46; 12 48; 13 3; 13 16; 13 28; 13 30; 13 40; 14 13; 14 17; 14 18; 14 30; 14 49; 15 5; 15 9; 15 14; 15 22; 15 34; 15 45; 16 3; 16 29; 16 33; 16 37; 16 49; 17 2; 17 4; 17 11; 17 27; 17 30; 17 39; 17 40; 18 20; 19 13; 19 21; 19 29; 20 5; 20 27; 20 48; 21 12; 21 25; 21 33; 21 45; 22 14; 22 16; 22 37; 22 45; 22 47; 22 50; 23 3; 23 6; 23 14; 23 15; 23 27; 24 5; 24 14; 24 17; 24 29; 24 37; 24 43; 24 48; 24 49; 25 5; 25 15; 25 32; 25 47; 26 9; 26 19; 26 23; 27 20; 27 21; 27 33; 27 44; 27 47; 28 13; 28 15; 28 27; 28 39; 28 41; 29 3; 29 13; 29 16; 29 19; 29 22; 29 27; 29 28; 29 32; 29 41; 30 5; 30 6; 30 10; 30 14; 30 15; 30 36; 30 44; 31 2; 31 7; 31 23; 31 25; 31 45; 31 49; 32 9; 32 15; 32 27; 32 33; 32 37; 32 49; 32 50; 33 2; 33 4; 33 8; 33 16; 33 30; 33 37; 33 41; 33 46; 34 3; 34 17; 34 19; 34 20; 34 35; 34 47; 34 49; 35 28; 36 13; 36 37; 37 19; 37 22; 37 46; 38 16; 38 23; 38 25; 39 9; 39 11; 39 17; 39 41; 40 9; 40 16; 40 23; 40 24; 40 36; 41 3; 41 22; 41 27; 41 37; 42 10; 42 45; 43 3; 43 12; 43 38; 43 45; 44 2; 44 6; 44 15; 44 26; 44 46; 45 8; 45 13; 45 42; 45 48; 45 49; 46 19; 46 29; 46 35; 46 38; 47 16; 47 23; 47 24; 47 28; 47 29; 47 40; 48 8; 48 10; 48 16; 48 23; 48 30; 48 33; 48 37; 48 39; 48 42; 49 3; 49 13; 49 29; 49 30; 49 35]
global d_x = [4.0, 4.0, 1.0, 4.0, 2.0, 2.0, 1.0, 7.0, 8.0, 2.0, 8.0, 7.0, 1.0, 4.0, 9.0, 2.0, 5.0, 9.0, 5.0, 5.0, 2.0, 9.0, 3.0, 3.0, 2.0, 1.0, 4.0, 5.0, 9.0, 1.0, 1.0, 5.0, 6.0, 9.0, 2.0, 8.0, 7.0, 5.0, 5.0, 7.0, 2.0, 5.0, 5.0, 4.0, 9.0, 1.0, 9.0, 8.0, 3.0, 3.0, 6.0, 9.0, 7.0, 8.0, 5.0, 7.0, 5.0, 7.0, 10.0, 2.0, 2.0, 10.0, 3.0, 1.0, 6.0, 6.0, 6.0, 8.0, 2.0, 7.0, 7.0, 10.0, 1.0, 3.0, 10.0, 4.0, 6.0, 1.0, 6.0, 1.0, 9.0, 2.0, 6.0, 8.0, 10.0, 1.0, 9.0, 8.0, 1.0, 1.0, 1.0, 6.0, 9.0, 8.0, 9.0, 10.0, 1.0, 1.0, 7.0, 3.0, 10.0, 4.0, 1.0, 1.0, 4.0, 8.0, 10.0, 6.0, 5.0, 8.0, 3.0, 2.0, 1.0, 4.0, 7.0, 7.0, 8.0, 9.0, 3.0, 5.0, 10.0, 5.0, 10.0, 5.0, 1.0, 5.0, 7.0, 2.0, 5.0, 6.0, 7.0, 6.0, 9.0, 7.0, 5.0, 9.0, 6.0, 4.0, 6.0, 6.0, 5.0, 4.0, 2.0, 5.0, 7.0, 4.0, 10.0, 10.0, 10.0, 7.0, 6.0, 8.0, 1.0, 10.0, 1.0, 3.0, 8.0, 5.0, 10.0, 9.0, 8.0, 3.0, 4.0, 5.0, 10.0, 9.0, 9.0, 6.0, 3.0, 7.0, 6.0, 10.0, 4.0, 3.0, 10.0, 4.0, 6.0, 10.0, 7.0, 9.0, 5.0, 3.0, 9.0, 7.0, 5.0, 10.0, 3.0, 8.0, 1.0, 1.0, 3.0, 5.0, 10.0, 4.0, 2.0, 8.0, 7.0, 10.0, 2.0, 10.0, 7.0, 9.0, 1.0, 8.0, 8.0, 8.0, 3.0, 3.0, 8.0, 3.0, 7.0, 3.0, 8.0, 7.0, 2.0, 1.0, 10.0, 8.0, 10.0, 1.0, 10.0, 10.0, 8.0, 1.0, 2.0, 5.0, 6.0, 8.0, 10.0, 1.0, 9.0, 2.0, 1.0, 5.0, 10.0, 3.0, 4.0, 9.0, 4.0, 8.0, 5.0, 7.0]
global b_x = 5
global d_y = [7.0, 4.0, 5.0, 1.0, 10.0, 3.0, 8.0, 4.0, 9.0, 10.0, 5.0, 5.0, 7.0, 9.0, 2.0, 2.0, 4.0, 5.0, 2.0, 5.0, 6.0, 3.0, 6.0, 4.0, 9.0, 7.0, 5.0, 6.0, 4.0, 4.0, 7.0, 6.0, 8.0, 2.0, 8.0, 10.0, 7.0, 4.0, 1.0, 5.0, 8.0, 10.0, 6.0, 2.0, 4.0, 9.0, 10.0, 7.0, 1.0, 9.0, 8.0, 6.0, 7.0, 8.0, 2.0, 9.0, 2.0, 1.0, 6.0, 5.0, 8.0, 8.0, 9.0, 4.0, 1.0, 3.0, 6.0, 9.0, 3.0, 6.0, 8.0, 9.0, 4.0, 6.0, 5.0, 3.0, 10.0, 8.0, 8.0, 10.0, 3.0, 3.0, 2.0, 8.0, 10.0, 8.0, 7.0, 6.0, 5.0, 3.0, 5.0, 6.0, 1.0, 1.0, 9.0, 6.0, 3.0, 1.0, 3.0, 4.0, 4.0, 5.0, 2.0, 6.0, 4.0, 10.0, 8.0, 6.0, 1.0, 6.0, 5.0, 3.0, 5.0, 3.0, 10.0, 9.0, 4.0, 6.0, 6.0, 7.0, 9.0, 5.0, 2.0, 3.0, 7.0, 8.0, 4.0, 8.0, 8.0, 7.0, 5.0, 10.0, 10.0, 4.0, 8.0, 9.0, 4.0, 1.0, 4.0, 3.0, 5.0, 10.0, 1.0, 9.0, 4.0, 1.0, 10.0, 5.0, 4.0, 10.0, 10.0, 10.0, 5.0, 2.0, 10.0, 1.0, 6.0, 8.0, 2.0, 4.0, 4.0, 8.0, 10.0, 1.0, 8.0, 5.0, 4.0, 10.0, 1.0, 4.0, 5.0, 10.0, 6.0, 9.0, 4.0, 2.0, 6.0, 2.0, 5.0, 2.0, 8.0, 3.0, 1.0, 7.0, 6.0, 10.0, 2.0, 7.0, 1.0, 7.0, 10.0, 1.0, 3.0, 5.0, 2.0, 5.0, 7.0, 5.0, 1.0, 10.0, 1.0, 9.0, 5.0, 4.0, 6.0, 3.0, 2.0, 8.0, 2.0, 9.0, 3.0, 1.0, 1.0, 9.0, 1.0, 5.0, 8.0, 4.0, 10.0, 1.0, 6.0, 3.0, 2.0, 5.0, 6.0, 10.0, 1.0, 9.0, 2.0, 3.0, 3.0, 9.0, 6.0, 6.0, 4.0, 2.0, 9.0, 4.0, 8.0, 1.0, 9.0, 2.0]
global b_y = 10
global p = [0.673, 0.307, 0.696, 0.461, 0.371, 0.213, 0.421, 0.621, 0.652, 0.465, 0.718, 0.355, 0.823, 0.557, 0.931, 0.044, 0.06, 0.709, 0.936, 0.018, 0.032, 0.247, 0.647, 0.376, 0.247, 0.675, 0.504, 0.372, 0.519, 0.618, 0.171, 0.221, 0.474, 0.175, 0.23, 0.775, 0.838, 0.12, 0.946, 0.785, 0.115, 0.371, 0.831, 0.382, 0.938, 0.84, 0.455, 0.81, 0.987, 0.963, 0.039, 0.234, 0.733, 0.001, 0.282, 0.771, 0.995, 0.048, 0.581, 0.527, 0.456, 0.116, 0.937, 0.068, 0.408, 0.631, 0.353, 0.487, 0.24, 0.167, 0.185, 0.328, 0.625, 0.58, 0.634, 0.628, 0.059, 0.46, 0.439, 0.076, 0.666, 0.387, 0.933, 0.2, 0.206, 0.342, 0.443, 0.564, 0.899, 0.385, 0.84, 0.092, 0.563, 0.787, 0.86, 0.151, 0.509, 0.878, 0.806, 0.847, 0.279, 0.609, 0.378, 0.162, 0.261, 0.073, 0.731, 0.682, 0.279, 0.589, 0.725, 0.72, 0.39, 0.725, 0.111, 0.777, 0.959, 0.657, 0.115, 0.377, 0.136, 0.287, 0.661, 0.224, 0.801, 0.161, 0.446, 0.828, 0.137, 0.724, 0.854, 0.086, 0.47, 0.208, 0.395, 0.53, 0.447, 0.453, 0.85, 0.757, 0.086, 0.649, 0.232, 0.421, 0.352, 0.725, 0.648, 0.033, 0.491, 0.822, 0.043, 0.473, 0.277, 0.013, 0.55, 0.922, 0.61, 0.519, 0.305, 0.272, 0.961, 0.405, 0.255, 0.652, 0.606, 0.35, 0.613, 0.115, 0.301, 0.121, 0.959, 0.847, 0.831, 0.893, 0.856, 0.271, 0.506, 0.91, 0.135, 0.069, 0.662, 0.389, 0.55, 0.959, 0.335, 0.836, 0.873, 0.484, 0.937, 0.665, 0.933, 0.117, 0.251, 0.9, 0.673, 0.323, 0.053, 0.113, 0.455, 0.087, 0.271, 0.439, 0.253, 0.028, 0.154, 0.878, 0.834, 0.556, 0.718, 0.234, 0.374, 0.924, 0.19, 0.499, 0.07, 0.874, 0.308, 0.704, 0.68, 0.191, 0.464, 0.337, 0.175, 0.858, 0.475, 0.126, 0.056, 0.31, 0.405, 0.552, 0.966, 0.612, 0.571, 0.943, 0.537, 0.264, 0.149, 0.543, 0.586, 0.587, 0.134, 0.895]
global q = [0.91, 0.444, 0.942, 0.829, 0.761, 0.633, 0.913, 0.637, 0.772, 0.98, 0.833, 0.464, 0.955, 0.686, 0.975, 0.443, 0.909, 0.923, 0.949, 0.453, 0.209, 0.899, 0.713, 0.773, 0.788, 0.927, 0.874, 0.889, 0.8, 0.878, 0.68, 0.561, 0.931, 0.243, 0.889, 0.979, 0.937, 0.76, 0.948, 0.861, 0.215, 0.379, 0.95, 0.5, 0.966, 0.999, 0.759, 0.978, 0.996, 0.985, 0.494, 0.659, 0.934, 0.554, 0.537, 0.899, 0.999, 0.999, 0.781, 0.791, 0.623, 0.546, 0.982, 0.518, 0.531, 0.676, 0.589, 0.911, 0.289, 0.704, 0.896, 0.926, 0.78, 0.947, 0.641, 0.707, 0.88, 0.758, 0.886, 0.17, 0.916, 0.782, 0.947, 0.684, 0.223, 0.49, 0.866, 0.793, 0.906, 0.616, 0.875, 0.664, 0.872, 0.874, 0.982, 0.938, 0.743, 0.91, 0.814, 0.88, 0.694, 0.908, 0.72, 0.449, 0.628, 0.333, 0.984, 0.741, 0.312, 0.641, 0.911, 0.837, 0.575, 0.846, 0.669, 0.904, 0.965, 0.731, 0.873, 0.803, 0.289, 0.72, 0.852, 0.959, 0.965, 0.891, 0.463, 0.953, 0.979, 0.856, 0.977, 0.616, 0.668, 0.278, 0.605, 0.859, 0.785, 0.803, 0.905, 0.849, 0.757, 0.686, 0.83, 0.477, 0.553, 0.755, 0.973, 0.851, 0.787, 0.874, 0.89, 0.817, 0.685, 0.29, 0.575, 0.923, 0.808, 0.921, 0.471, 0.629, 0.968, 0.454, 0.622, 0.804, 0.609, 0.462, 0.655, 0.401, 0.901, 0.84, 0.995, 0.962, 0.907, 0.967, 0.887, 0.892, 0.687, 0.914, 0.674, 0.976, 0.845, 0.947, 0.813, 0.965, 0.69, 0.931, 0.88, 0.514, 0.997, 0.79, 0.946, 0.701, 0.525, 0.971, 0.684, 0.446, 0.427, 0.882, 0.676, 0.718, 0.782, 0.821, 0.683, 0.351, 0.268, 0.981, 0.951, 0.786, 0.875, 0.416, 0.674, 0.987, 0.223, 0.527, 0.414, 0.984, 0.995, 0.794, 0.681, 0.379, 0.701, 0.953, 0.51, 0.958, 0.638, 0.446, 0.295, 0.592, 0.652, 0.86, 0.995, 0.819, 0.935, 0.947, 0.666, 0.997, 0.539, 0.784, 0.736, 0.704, 0.585, 0.939]
global origin = 1
global destination = 50