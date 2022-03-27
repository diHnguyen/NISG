global arcs = [1 11; 1 13; 1 46; 2 8; 2 19; 2 20; 2 34; 2 36; 3 7; 3 9; 3 15; 3 18; 3 44; 3 45; 4 29; 4 46; 5 6; 5 24; 5 39; 5 40; 6 31; 6 42; 7 41; 7 46; 7 49; 8 16; 8 26; 8 28; 8 31; 8 32; 8 37; 8 39; 8 46; 9 3; 9 7; 9 10; 9 43; 10 8; 10 21; 10 28; 10 30; 11 6; 11 31; 11 36; 12 2; 12 8; 12 15; 12 16; 12 19; 12 40; 13 44; 14 33; 14 35; 14 36; 15 5; 15 6; 15 9; 16 7; 16 34; 16 39; 16 47; 17 3; 17 13; 17 14; 17 41; 18 24; 18 26; 18 28; 18 32; 18 46; 19 16; 19 18; 19 27; 19 29; 19 30; 19 38; 19 42; 20 12; 20 46; 21 2; 21 6; 21 8; 21 15; 21 30; 21 33; 21 39; 21 48; 22 16; 22 25; 22 32; 22 39; 22 43; 22 47; 23 25; 23 34; 23 39; 24 21; 24 36; 24 38; 25 2; 25 5; 25 16; 25 27; 25 35; 25 40; 26 32; 26 35; 26 37; 26 45; 27 14; 27 24; 27 44; 27 47; 28 2; 28 6; 28 14; 28 22; 28 26; 28 36; 28 39; 28 46; 29 7; 29 21; 29 33; 29 35; 29 44; 30 6; 30 10; 30 11; 30 15; 30 17; 30 20; 30 23; 31 23; 31 24; 31 28; 31 34; 31 35; 31 45; 31 47; 32 4; 32 11; 32 18; 32 22; 32 44; 32 49; 33 8; 33 28; 33 40; 33 41; 34 11; 34 13; 34 23; 34 24; 34 40; 35 8; 35 50; 36 3; 36 13; 36 15; 36 31; 36 41; 37 5; 37 10; 37 15; 37 22; 37 23; 37 44; 37 46; 38 2; 38 11; 38 29; 38 32; 38 35; 38 49; 39 4; 39 24; 39 31; 39 34; 39 35; 39 44; 39 48; 39 50; 40 8; 40 21; 40 27; 41 4; 41 9; 41 24; 41 39; 41 42; 41 49; 42 4; 42 6; 42 13; 42 25; 42 32; 42 33; 42 35; 42 49; 43 8; 43 12; 43 22; 43 38; 43 39; 43 49; 44 12; 44 28; 44 43; 45 6; 45 9; 45 14; 45 18; 45 19; 45 24; 45 25; 45 28; 45 31; 45 36; 46 16; 46 20; 46 25; 46 36; 47 13; 47 22; 47 25; 47 29; 48 6; 48 20; 48 39; 48 50; 49 2; 49 6; 49 20; 49 47]
global d_x = [3.0, 4.0, 8.0, 6.0, 1.0, 10.0, 3.0, 10.0, 3.0, 2.0, 1.0, 1.0, 4.0, 2.0, 4.0, 9.0, 6.0, 5.0, 3.0, 1.0, 2.0, 6.0, 8.0, 9.0, 5.0, 2.0, 7.0, 4.0, 10.0, 1.0, 5.0, 5.0, 8.0, 2.0, 3.0, 2.0, 9.0, 5.0, 4.0, 6.0, 6.0, 5.0, 1.0, 3.0, 7.0, 8.0, 10.0, 4.0, 4.0, 7.0, 4.0, 8.0, 6.0, 7.0, 5.0, 8.0, 6.0, 4.0, 4.0, 6.0, 4.0, 7.0, 4.0, 4.0, 1.0, 4.0, 10.0, 4.0, 8.0, 2.0, 5.0, 4.0, 7.0, 1.0, 3.0, 8.0, 6.0, 3.0, 9.0, 5.0, 4.0, 1.0, 4.0, 9.0, 5.0, 4.0, 10.0, 3.0, 9.0, 7.0, 5.0, 2.0, 5.0, 6.0, 2.0, 7.0, 2.0, 8.0, 9.0, 8.0, 10.0, 2.0, 9.0, 10.0, 6.0, 4.0, 1.0, 3.0, 8.0, 4.0, 2.0, 5.0, 2.0, 10.0, 3.0, 7.0, 9.0, 10.0, 7.0, 10.0, 4.0, 10.0, 8.0, 3.0, 10.0, 9.0, 9.0, 9.0, 6.0, 2.0, 5.0, 3.0, 5.0, 8.0, 2.0, 4.0, 5.0, 8.0, 4.0, 8.0, 3.0, 10.0, 2.0, 1.0, 10.0, 9.0, 7.0, 10.0, 4.0, 7.0, 8.0, 5.0, 8.0, 4.0, 6.0, 3.0, 2.0, 3.0, 3.0, 10.0, 2.0, 8.0, 8.0, 5.0, 7.0, 8.0, 8.0, 8.0, 9.0, 10.0, 7.0, 3.0, 8.0, 8.0, 7.0, 6.0, 9.0, 1.0, 9.0, 7.0, 10.0, 1.0, 7.0, 1.0, 8.0, 10.0, 6.0, 1.0, 10.0, 5.0, 7.0, 4.0, 3.0, 1.0, 8.0, 1.0, 2.0, 9.0, 10.0, 4.0, 1.0, 9.0, 10.0, 2.0, 5.0, 1.0, 1.0, 3.0, 10.0, 1.0, 2.0, 3.0, 5.0, 2.0, 1.0, 4.0, 10.0, 4.0, 5.0, 8.0, 1.0, 8.0, 2.0, 9.0, 3.0, 3.0, 3.0, 3.0, 8.0, 4.0, 5.0, 10.0, 3.0, 9.0, 10.0]
global b_x = 5
global d_y = [5.0, 9.0, 8.0, 6.0, 1.0, 10.0, 6.0, 9.0, 7.0, 10.0, 3.0, 5.0, 5.0, 8.0, 3.0, 7.0, 6.0, 5.0, 5.0, 5.0, 8.0, 7.0, 4.0, 2.0, 7.0, 4.0, 1.0, 4.0, 10.0, 3.0, 2.0, 6.0, 4.0, 5.0, 7.0, 4.0, 6.0, 8.0, 1.0, 4.0, 5.0, 6.0, 6.0, 6.0, 10.0, 5.0, 4.0, 6.0, 1.0, 2.0, 9.0, 5.0, 1.0, 6.0, 7.0, 8.0, 2.0, 6.0, 8.0, 3.0, 4.0, 1.0, 7.0, 5.0, 4.0, 4.0, 1.0, 10.0, 2.0, 7.0, 10.0, 10.0, 10.0, 4.0, 5.0, 4.0, 1.0, 7.0, 4.0, 6.0, 1.0, 2.0, 8.0, 4.0, 10.0, 8.0, 5.0, 10.0, 5.0, 9.0, 1.0, 8.0, 8.0, 7.0, 1.0, 9.0, 4.0, 3.0, 7.0, 2.0, 5.0, 6.0, 5.0, 8.0, 4.0, 3.0, 2.0, 1.0, 5.0, 10.0, 4.0, 1.0, 9.0, 6.0, 1.0, 10.0, 6.0, 10.0, 2.0, 8.0, 5.0, 3.0, 4.0, 2.0, 2.0, 7.0, 5.0, 5.0, 4.0, 8.0, 1.0, 9.0, 6.0, 3.0, 2.0, 6.0, 3.0, 9.0, 6.0, 5.0, 3.0, 2.0, 1.0, 2.0, 6.0, 1.0, 7.0, 7.0, 3.0, 5.0, 4.0, 7.0, 4.0, 8.0, 5.0, 9.0, 1.0, 3.0, 1.0, 7.0, 3.0, 2.0, 4.0, 4.0, 5.0, 8.0, 3.0, 1.0, 2.0, 7.0, 2.0, 2.0, 2.0, 7.0, 9.0, 5.0, 6.0, 1.0, 6.0, 2.0, 9.0, 3.0, 1.0, 6.0, 9.0, 2.0, 1.0, 4.0, 4.0, 3.0, 1.0, 1.0, 9.0, 8.0, 8.0, 7.0, 7.0, 6.0, 10.0, 8.0, 1.0, 1.0, 1.0, 3.0, 1.0, 10.0, 10.0, 5.0, 8.0, 10.0, 3.0, 4.0, 4.0, 6.0, 4.0, 4.0, 3.0, 3.0, 2.0, 4.0, 5.0, 5.0, 4.0, 8.0, 6.0, 2.0, 8.0, 3.0, 9.0, 2.0, 6.0, 7.0, 10.0, 8.0, 1.0]
global b_y = 10
global p = [0.028, 0.702, 0.801, 0.67, 0.76, 0.006, 0.774, 0.278, 0.686, 0.896, 0.402, 0.218, 0.393, 0.505, 0.415, 0.526, 0.29, 0.621, 0.941, 0.977, 0.371, 0.717, 0.226, 0.225, 0.597, 0.056, 0.06, 0.83, 0.833, 0.835, 0.333, 0.148, 0.9, 0.337, 0.729, 0.705, 0.702, 0.513, 0.111, 0.799, 0.613, 0.036, 0.757, 0.695, 0.32, 0.851, 0.956, 0.7, 0.598, 0.129, 0.781, 0.042, 0.437, 0.624, 0.477, 0.713, 0.026, 0.128, 0.608, 0.523, 0.705, 0.099, 0.958, 0.638, 0.694, 0.302, 0.64, 0.755, 0.875, 0.674, 0.01, 0.788, 0.913, 0.718, 0.894, 0.09, 0.218, 0.601, 0.168, 0.473, 0.528, 0.551, 0.401, 0.556, 0.295, 0.788, 0.451, 0.514, 0.494, 0.165, 0.587, 0.938, 0.462, 0.336, 0.267, 0.114, 0.486, 0.197, 0.615, 0.831, 0.04, 0.588, 0.669, 0.789, 0.484, 0.074, 0.519, 0.671, 0.285, 0.491, 0.327, 0.359, 0.859, 0.531, 0.967, 0.272, 0.023, 0.821, 0.342, 0.506, 0.099, 0.155, 0.244, 0.303, 0.294, 0.493, 0.984, 0.293, 0.467, 0.67, 0.491, 0.167, 0.192, 0.544, 0.576, 0.037, 0.018, 0.498, 0.147, 0.016, 0.899, 0.131, 0.804, 0.735, 0.715, 0.917, 0.12, 0.852, 0.046, 0.551, 0.857, 0.635, 0.17, 0.851, 0.72, 0.274, 0.954, 0.173, 0.761, 0.478, 0.274, 0.992, 0.705, 0.268, 0.242, 0.476, 0.803, 0.562, 0.277, 0.935, 0.51, 0.16, 0.221, 0.06, 0.755, 0.403, 0.161, 0.902, 0.077, 0.383, 0.361, 0.558, 0.832, 0.672, 0.558, 0.907, 0.419, 0.339, 0.562, 0.676, 0.665, 0.413, 0.217, 0.329, 0.652, 0.665, 0.901, 0.084, 0.194, 0.711, 0.729, 0.704, 0.952, 0.885, 0.431, 0.518, 0.227, 0.203, 0.799, 0.342, 0.852, 0.274, 0.196, 0.37, 0.531, 0.751, 0.231, 0.281, 0.031, 0.665, 0.816, 0.538, 0.119, 0.992, 0.516, 0.882, 0.229, 0.271, 0.241, 0.885, 0.694, 0.637, 0.125, 0.281, 0.887]
global q = [0.418, 0.792, 0.968, 0.901, 0.881, 0.47, 0.978, 0.633, 0.838, 0.958, 0.829, 0.24, 0.876, 0.689, 0.723, 0.649, 0.367, 0.84, 0.964, 0.997, 0.745, 0.812, 0.646, 0.823, 0.634, 0.637, 0.439, 0.852, 0.906, 0.954, 0.997, 0.28, 0.922, 0.639, 0.894, 0.729, 0.813, 0.778, 0.309, 0.95, 0.967, 0.645, 0.992, 0.72, 0.587, 0.854, 0.966, 0.758, 0.871, 0.765, 0.949, 0.552, 0.63, 0.843, 0.664, 0.959, 0.621, 0.248, 0.735, 0.606, 0.951, 0.302, 0.962, 0.821, 0.788, 0.584, 0.915, 0.895, 0.922, 0.843, 0.357, 0.85, 0.941, 0.89, 0.903, 0.147, 0.919, 0.799, 0.256, 0.573, 0.744, 0.757, 0.807, 0.71, 0.944, 0.972, 0.681, 0.572, 0.827, 0.848, 0.989, 0.943, 0.814, 0.461, 0.442, 0.493, 0.935, 0.674, 0.876, 0.871, 0.963, 0.672, 0.955, 0.935, 0.883, 0.95, 0.677, 0.934, 0.83, 0.644, 0.771, 0.969, 0.948, 0.949, 0.989, 0.299, 0.408, 0.968, 0.448, 0.735, 0.957, 0.977, 0.385, 0.852, 0.818, 0.517, 0.985, 0.799, 0.525, 0.862, 0.776, 0.275, 0.416, 0.567, 0.652, 0.677, 0.852, 0.878, 0.911, 0.557, 0.908, 0.717, 0.81, 0.896, 0.979, 0.956, 0.143, 0.882, 0.968, 0.85, 0.951, 0.904, 0.266, 0.892, 0.877, 0.457, 0.97, 0.72, 0.89, 0.584, 0.948, 0.992, 0.821, 0.555, 0.964, 0.66, 0.806, 0.74, 0.582, 0.948, 0.833, 0.685, 0.274, 0.525, 0.97, 0.8, 0.891, 0.935, 0.698, 0.927, 0.639, 0.79, 0.985, 0.861, 0.888, 0.948, 0.526, 0.76, 0.995, 0.994, 0.868, 0.531, 0.541, 0.531, 0.969, 0.978, 0.91, 0.865, 0.854, 0.762, 0.973, 0.836, 0.995, 0.92, 0.543, 0.859, 0.649, 0.336, 0.955, 0.93, 0.861, 0.393, 0.849, 0.449, 0.554, 0.752, 0.695, 0.48, 0.101, 0.841, 0.983, 0.551, 0.779, 0.994, 0.97, 0.892, 0.74, 0.348, 0.851, 0.93, 0.866, 0.742, 0.882, 0.293, 0.907]
global origin = 1
global destination = 50