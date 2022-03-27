global arcs = [1 5; 1 11; 1 12; 1 27; 1 28; 1 33; 1 45; 2 12; 2 16; 2 22; 2 36; 2 44; 2 47; 2 48; 2 54; 2 56; 3 5; 3 9; 3 26; 3 29; 3 54; 3 56; 3 57; 4 3; 4 19; 4 21; 4 28; 4 44; 4 52; 5 3; 5 15; 5 34; 5 42; 5 56; 6 9; 6 10; 6 18; 6 33; 6 35; 6 57; 7 2; 7 5; 7 10; 7 15; 7 19; 7 28; 7 32; 7 46; 7 49; 7 53; 8 14; 8 18; 8 33; 8 54; 9 8; 9 12; 9 13; 9 30; 9 44; 9 57; 10 11; 10 21; 10 55; 10 57; 11 6; 11 12; 11 38; 11 47; 11 53; 12 10; 12 11; 12 15; 12 16; 12 18; 12 28; 12 40; 13 5; 13 21; 13 30; 13 32; 13 48; 14 4; 14 10; 14 13; 14 23; 14 34; 14 40; 14 42; 14 49; 15 22; 15 23; 15 30; 16 11; 16 35; 16 44; 16 46; 17 12; 17 16; 17 23; 17 33; 17 40; 18 4; 18 32; 18 45; 18 51; 19 14; 19 20; 19 32; 19 42; 19 56; 20 6; 20 10; 20 49; 20 52; 20 55; 21 19; 21 26; 21 35; 21 46; 21 53; 22 9; 22 11; 22 40; 22 46; 23 5; 23 24; 23 34; 23 51; 23 54; 23 60; 24 12; 24 23; 24 38; 24 54; 24 59; 25 3; 25 7; 25 26; 25 50; 25 52; 26 57; 26 60; 27 2; 27 12; 27 17; 27 55; 28 3; 28 26; 28 38; 28 40; 28 46; 28 51; 29 7; 29 10; 29 16; 29 24; 29 27; 29 48; 29 54; 30 8; 30 9; 30 11; 30 24; 30 33; 30 36; 30 42; 30 50; 30 58; 30 59; 31 14; 31 17; 31 50; 32 18; 32 29; 32 43; 32 53; 32 59; 33 3; 33 5; 33 16; 33 40; 33 56; 34 2; 34 27; 34 30; 34 39; 34 55; 34 56; 34 59; 34 60; 35 32; 35 45; 35 58; 35 60; 36 17; 36 20; 36 28; 36 37; 36 45; 36 54; 37 6; 37 8; 37 9; 37 43; 38 8; 38 14; 38 21; 38 24; 38 25; 38 40; 38 42; 38 60; 39 2; 39 14; 39 15; 39 30; 39 33; 39 47; 39 57; 40 4; 40 18; 40 31; 40 34; 40 37; 40 41; 40 42; 40 58; 41 5; 41 12; 41 17; 41 18; 41 28; 41 40; 41 43; 41 46; 42 9; 42 24; 42 30; 42 36; 42 39; 42 41; 42 48; 42 56; 43 4; 43 42; 43 48; 43 51; 44 22; 45 10; 45 27; 45 38; 46 8; 46 12; 46 15; 46 16; 46 42; 46 45; 46 58; 47 14; 47 23; 47 33; 47 37; 47 43; 47 52; 48 12; 48 13; 48 30; 48 40; 48 43; 48 46; 48 49; 49 8; 49 10; 49 19; 49 36; 49 38; 49 43; 49 57; 50 2; 50 16; 50 38; 50 41; 50 43; 50 46; 51 14; 51 22; 51 34; 51 50; 51 57; 52 7; 52 11; 52 26; 52 44; 52 55; 52 57; 52 58; 53 19; 53 21; 53 23; 53 25; 53 46; 53 48; 54 5; 54 12; 54 26; 54 37; 54 57; 55 10; 55 18; 55 27; 55 38; 55 42; 55 53; 55 60; 56 11; 56 48; 56 59; 57 6; 57 7; 57 13; 57 23; 57 24; 57 29; 57 36; 57 52; 58 7; 58 9; 58 19; 58 28; 58 52; 58 55; 58 59; 59 4; 59 8; 59 9; 59 14; 59 16; 59 19; 59 29; 59 35; 59 50]
global d_x = [4.0, 2.0, 2.0, 4.0, 6.0, 3.0, 5.0, 8.0, 2.0, 8.0, 7.0, 2.0, 10.0, 7.0, 2.0, 6.0, 10.0, 4.0, 2.0, 1.0, 8.0, 7.0, 5.0, 9.0, 10.0, 10.0, 7.0, 3.0, 4.0, 6.0, 9.0, 3.0, 8.0, 9.0, 9.0, 3.0, 7.0, 4.0, 2.0, 8.0, 7.0, 7.0, 6.0, 2.0, 4.0, 9.0, 5.0, 8.0, 9.0, 3.0, 10.0, 2.0, 4.0, 3.0, 1.0, 1.0, 3.0, 9.0, 7.0, 4.0, 1.0, 10.0, 8.0, 7.0, 3.0, 3.0, 9.0, 5.0, 8.0, 6.0, 6.0, 1.0, 1.0, 3.0, 9.0, 9.0, 9.0, 8.0, 1.0, 5.0, 3.0, 7.0, 7.0, 3.0, 6.0, 6.0, 10.0, 1.0, 4.0, 5.0, 9.0, 9.0, 10.0, 5.0, 5.0, 5.0, 3.0, 9.0, 4.0, 3.0, 5.0, 8.0, 5.0, 4.0, 6.0, 6.0, 4.0, 5.0, 6.0, 7.0, 5.0, 8.0, 8.0, 9.0, 8.0, 6.0, 1.0, 5.0, 7.0, 6.0, 8.0, 3.0, 6.0, 6.0, 8.0, 6.0, 6.0, 2.0, 6.0, 7.0, 6.0, 2.0, 9.0, 8.0, 3.0, 1.0, 7.0, 4.0, 9.0, 7.0, 4.0, 8.0, 4.0, 10.0, 4.0, 4.0, 5.0, 6.0, 7.0, 10.0, 4.0, 8.0, 1.0, 9.0, 1.0, 9.0, 6.0, 5.0, 5.0, 3.0, 7.0, 5.0, 6.0, 2.0, 4.0, 2.0, 7.0, 10.0, 8.0, 7.0, 9.0, 8.0, 6.0, 5.0, 9.0, 8.0, 5.0, 7.0, 8.0, 3.0, 2.0, 4.0, 7.0, 4.0, 3.0, 6.0, 9.0, 8.0, 6.0, 6.0, 5.0, 4.0, 7.0, 10.0, 5.0, 5.0, 7.0, 4.0, 3.0, 7.0, 4.0, 10.0, 4.0, 3.0, 8.0, 2.0, 8.0, 5.0, 3.0, 10.0, 1.0, 9.0, 8.0, 2.0, 8.0, 7.0, 9.0, 4.0, 4.0, 10.0, 9.0, 4.0, 8.0, 4.0, 6.0, 8.0, 8.0, 7.0, 1.0, 7.0, 6.0, 3.0, 3.0, 7.0, 6.0, 2.0, 9.0, 9.0, 1.0, 1.0, 1.0, 4.0, 6.0, 2.0, 4.0, 4.0, 4.0, 6.0, 4.0, 9.0, 2.0, 2.0, 6.0, 9.0, 10.0, 1.0, 4.0, 3.0, 4.0, 1.0, 7.0, 5.0, 6.0, 5.0, 6.0, 9.0, 2.0, 7.0, 6.0, 10.0, 8.0, 7.0, 4.0, 1.0, 4.0, 7.0, 4.0, 2.0, 5.0, 4.0, 7.0, 2.0, 1.0, 1.0, 7.0, 4.0, 7.0, 7.0, 9.0, 4.0, 5.0, 1.0, 2.0, 8.0, 1.0, 1.0, 2.0, 1.0, 5.0, 8.0, 10.0, 1.0, 7.0, 1.0, 3.0, 4.0, 4.0, 1.0, 5.0, 4.0, 3.0, 4.0, 7.0, 5.0, 8.0, 5.0, 1.0, 9.0, 6.0, 9.0, 2.0, 5.0, 10.0, 4.0, 1.0, 9.0, 3.0, 8.0, 9.0, 3.0, 9.0, 7.0, 7.0, 3.0, 1.0, 1.0, 4.0, 5.0, 6.0, 8.0, 2.0]
global b_x = 5
global d_y = [4.0, 10.0, 1.0, 5.0, 6.0, 7.0, 10.0, 9.0, 10.0, 7.0, 5.0, 7.0, 7.0, 7.0, 8.0, 3.0, 4.0, 1.0, 6.0, 9.0, 3.0, 6.0, 5.0, 4.0, 5.0, 1.0, 10.0, 4.0, 4.0, 5.0, 2.0, 10.0, 5.0, 2.0, 10.0, 4.0, 9.0, 9.0, 1.0, 7.0, 3.0, 8.0, 6.0, 6.0, 6.0, 1.0, 1.0, 5.0, 7.0, 3.0, 8.0, 6.0, 6.0, 6.0, 6.0, 3.0, 4.0, 5.0, 2.0, 3.0, 1.0, 3.0, 3.0, 6.0, 6.0, 1.0, 6.0, 5.0, 8.0, 3.0, 10.0, 1.0, 3.0, 1.0, 10.0, 4.0, 7.0, 4.0, 4.0, 4.0, 10.0, 10.0, 9.0, 10.0, 7.0, 2.0, 7.0, 2.0, 1.0, 3.0, 2.0, 4.0, 3.0, 9.0, 1.0, 8.0, 5.0, 8.0, 7.0, 1.0, 3.0, 2.0, 5.0, 2.0, 8.0, 7.0, 7.0, 8.0, 2.0, 7.0, 3.0, 2.0, 1.0, 5.0, 4.0, 9.0, 8.0, 7.0, 5.0, 10.0, 5.0, 1.0, 7.0, 3.0, 3.0, 1.0, 6.0, 9.0, 10.0, 2.0, 3.0, 3.0, 10.0, 9.0, 4.0, 8.0, 1.0, 4.0, 7.0, 6.0, 7.0, 8.0, 10.0, 5.0, 4.0, 1.0, 7.0, 1.0, 2.0, 6.0, 2.0, 9.0, 7.0, 5.0, 10.0, 4.0, 5.0, 2.0, 7.0, 1.0, 5.0, 4.0, 4.0, 7.0, 8.0, 5.0, 6.0, 6.0, 4.0, 3.0, 2.0, 3.0, 2.0, 4.0, 6.0, 8.0, 6.0, 7.0, 2.0, 6.0, 3.0, 4.0, 8.0, 8.0, 5.0, 6.0, 9.0, 3.0, 6.0, 7.0, 10.0, 7.0, 4.0, 10.0, 4.0, 9.0, 9.0, 3.0, 7.0, 6.0, 8.0, 3.0, 9.0, 3.0, 8.0, 2.0, 6.0, 5.0, 2.0, 10.0, 6.0, 10.0, 8.0, 10.0, 4.0, 1.0, 7.0, 3.0, 3.0, 3.0, 3.0, 5.0, 10.0, 9.0, 7.0, 10.0, 7.0, 10.0, 6.0, 8.0, 7.0, 5.0, 9.0, 5.0, 10.0, 1.0, 1.0, 5.0, 1.0, 8.0, 8.0, 10.0, 3.0, 3.0, 2.0, 5.0, 1.0, 7.0, 9.0, 3.0, 10.0, 9.0, 10.0, 5.0, 2.0, 1.0, 6.0, 4.0, 8.0, 6.0, 10.0, 2.0, 10.0, 5.0, 9.0, 10.0, 4.0, 2.0, 5.0, 9.0, 2.0, 10.0, 9.0, 10.0, 7.0, 2.0, 8.0, 8.0, 4.0, 2.0, 6.0, 8.0, 7.0, 6.0, 5.0, 1.0, 5.0, 1.0, 9.0, 10.0, 3.0, 5.0, 4.0, 4.0, 9.0, 7.0, 9.0, 5.0, 9.0, 5.0, 9.0, 2.0, 1.0, 5.0, 5.0, 1.0, 5.0, 2.0, 10.0, 3.0, 6.0, 5.0, 4.0, 4.0, 7.0, 9.0, 8.0, 9.0, 8.0, 5.0, 8.0, 5.0, 9.0, 7.0, 6.0, 1.0, 5.0, 3.0, 1.0, 10.0, 2.0, 2.0, 3.0, 2.0, 1.0, 1.0, 2.0, 5.0, 3.0, 9.0, 4.0]
global b_y = 10
global p = [0.97, 0.891, 0.054, 0.228, 0.685, 0.967, 0.724, 0.157, 0.749, 0.923, 0.359, 0.362, 0.268, 0.607, 0.156, 0.844, 0.493, 0.913, 0.027, 0.308, 0.713, 0.178, 0.648, 0.866, 0.768, 0.955, 0.571, 0.716, 0.485, 0.466, 0.788, 0.829, 0.231, 0.776, 0.79, 0.666, 0.263, 0.004, 0.026, 0.344, 0.905, 0.655, 0.542, 0.718, 0.382, 0.094, 0.142, 0.504, 0.088, 0.832, 0.911, 0.225, 0.944, 0.371, 0.312, 0.221, 0.771, 0.116, 0.14, 0.617, 0.663, 0.449, 0.782, 0.745, 0.277, 0.292, 0.325, 0.455, 0.142, 0.124, 0.413, 0.039, 0.711, 0.852, 0.366, 0.915, 0.644, 0.288, 0.899, 0.247, 0.466, 0.046, 0.962, 0.628, 0.239, 0.438, 0.665, 0.309, 0.065, 0.877, 0.425, 0.849, 0.795, 0.626, 0.219, 0.909, 0.845, 0.54, 0.472, 0.519, 0.427, 0.116, 0.256, 0.198, 0.321, 0.124, 0.752, 0.584, 0.911, 0.602, 0.854, 0.081, 0.376, 0.197, 0.072, 0.852, 0.173, 0.176, 0.684, 0.477, 0.262, 0.856, 0.277, 0.862, 0.594, 0.271, 0.153, 0.64, 0.847, 0.442, 0.069, 0.256, 0.114, 0.073, 0.545, 0.15, 0.939, 0.24, 0.069, 0.829, 0.802, 0.514, 0.015, 0.05, 0.678, 0.671, 0.573, 0.375, 0.6, 0.412, 0.16, 0.782, 0.477, 0.744, 0.51, 0.922, 0.62, 0.319, 0.169, 0.535, 0.226, 0.688, 0.793, 0.164, 0.416, 0.313, 0.623, 0.821, 0.517, 0.76, 0.195, 0.587, 0.806, 0.573, 0.426, 0.675, 0.222, 0.723, 0.792, 0.863, 0.788, 0.193, 0.134, 0.197, 0.275, 0.253, 0.561, 0.194, 0.491, 0.67, 0.342, 0.99, 0.818, 0.11, 0.624, 0.383, 0.465, 0.981, 0.216, 0.789, 0.801, 0.584, 0.628, 0.302, 0.947, 0.592, 0.866, 0.448, 0.451, 0.672, 0.153, 0.584, 0.814, 0.807, 0.136, 0.823, 0.976, 0.61, 0.504, 0.784, 0.224, 0.179, 0.418, 0.82, 0.829, 0.78, 0.156, 0.034, 0.236, 0.495, 0.368, 0.736, 0.373, 0.089, 0.364, 0.68, 0.228, 0.375, 0.37, 0.891, 0.277, 0.165, 0.688, 0.975, 0.86, 0.136, 0.752, 0.411, 0.624, 0.293, 0.777, 0.903, 0.657, 0.531, 0.986, 0.99, 0.073, 0.693, 0.099, 0.695, 0.37, 0.131, 0.537, 0.256, 0.37, 0.151, 0.639, 0.265, 0.796, 0.17, 0.708, 0.238, 0.7, 0.04, 0.109, 0.644, 0.954, 0.128, 0.525, 0.289, 0.112, 0.331, 0.021, 0.59, 0.317, 0.52, 0.921, 0.945, 0.599, 0.219, 0.635, 0.541, 0.775, 0.739, 0.496, 0.334, 0.458, 0.15, 0.587, 0.191, 0.001, 0.913, 0.023, 0.553, 0.249, 0.047, 0.209, 0.35, 0.434, 0.645, 0.697, 0.669, 0.158, 0.71, 0.309, 0.898, 0.739, 0.78, 0.259, 0.377, 0.933, 0.884, 0.627, 0.276, 0.899, 0.446, 0.166, 0.291, 0.031, 0.778, 0.17, 0.326, 0.574, 0.428, 0.302, 0.266, 0.788, 0.249, 0.04, 0.175, 0.955]
global q = [0.987, 0.984, 0.298, 0.576, 0.745, 0.995, 0.861, 0.994, 0.908, 0.924, 0.552, 0.853, 0.7, 0.943, 0.828, 0.868, 0.545, 0.963, 0.962, 0.61, 0.83, 0.842, 0.721, 0.902, 0.881, 0.968, 0.678, 0.99, 0.889, 0.567, 0.801, 0.871, 0.862, 0.79, 0.928, 0.738, 0.732, 0.714, 0.965, 0.527, 0.905, 0.769, 0.681, 0.912, 0.418, 0.568, 0.212, 0.937, 0.72, 0.844, 0.928, 0.917, 0.983, 0.918, 0.772, 0.644, 0.966, 0.856, 0.423, 0.826, 0.677, 0.981, 0.971, 0.834, 0.583, 0.622, 0.479, 0.574, 0.935, 0.7, 0.961, 0.855, 0.913, 0.9, 0.58, 0.92, 0.866, 0.647, 0.982, 0.323, 0.992, 0.769, 0.977, 0.707, 0.393, 0.742, 0.72, 0.544, 0.683, 0.906, 0.71, 0.92, 0.965, 0.877, 0.427, 0.94, 0.916, 0.879, 0.477, 0.676, 0.514, 0.366, 0.38, 0.22, 0.718, 0.223, 0.88, 0.597, 0.911, 0.808, 0.932, 0.154, 0.952, 0.555, 0.528, 0.999, 0.822, 0.546, 0.892, 0.989, 0.685, 0.947, 0.339, 0.914, 0.946, 0.546, 0.601, 0.707, 0.993, 0.903, 0.422, 0.801, 0.479, 0.694, 0.712, 0.555, 0.959, 0.7, 0.972, 0.87, 0.835, 0.637, 0.523, 0.546, 0.881, 0.85, 0.583, 0.541, 0.75, 0.818, 0.738, 0.993, 0.56, 0.767, 0.877, 0.954, 0.638, 0.392, 0.466, 0.562, 0.566, 0.866, 0.961, 0.844, 0.984, 0.582, 0.695, 0.939, 0.942, 0.81, 0.776, 0.846, 0.887, 0.681, 0.623, 0.754, 0.414, 0.738, 0.84, 0.872, 0.832, 0.521, 0.844, 0.686, 0.732, 0.699, 0.756, 0.95, 0.494, 0.807, 0.631, 0.995, 0.906, 0.763, 0.893, 0.719, 0.756, 0.984, 0.837, 0.899, 0.986, 0.896, 0.939, 0.754, 0.994, 0.659, 0.884, 0.776, 0.671, 0.984, 0.516, 0.867, 0.893, 0.974, 0.964, 0.842, 0.981, 0.683, 0.656, 0.949, 0.304, 0.233, 0.63, 0.837, 0.911, 0.883, 0.314, 0.283, 0.339, 0.748, 0.412, 0.807, 0.437, 0.419, 0.904, 0.68, 0.9, 0.507, 0.584, 0.939, 0.885, 0.788, 0.833, 0.985, 0.997, 0.608, 0.947, 0.861, 0.748, 0.885, 0.819, 0.924, 0.957, 0.883, 0.99, 0.993, 0.215, 0.775, 0.805, 0.902, 0.896, 0.303, 0.588, 0.445, 0.56, 0.596, 0.77, 0.628, 0.895, 0.938, 0.755, 0.449, 0.92, 0.417, 0.477, 0.838, 0.985, 0.727, 0.935, 0.961, 0.424, 0.724, 0.081, 0.645, 0.562, 0.893, 0.936, 0.945, 0.646, 0.636, 0.923, 0.857, 0.908, 0.913, 0.567, 0.799, 0.587, 0.641, 0.844, 0.81, 0.651, 0.963, 0.455, 0.765, 0.913, 0.763, 0.711, 0.627, 0.939, 0.781, 0.739, 0.676, 0.691, 0.964, 0.676, 0.926, 0.851, 0.86, 0.758, 0.765, 0.941, 0.995, 0.881, 0.853, 0.931, 0.781, 0.751, 0.442, 0.41, 0.998, 0.174, 0.399, 0.605, 0.584, 0.556, 0.831, 0.79, 0.471, 0.893, 0.741, 0.999]
global origin = 1
global destination = 60