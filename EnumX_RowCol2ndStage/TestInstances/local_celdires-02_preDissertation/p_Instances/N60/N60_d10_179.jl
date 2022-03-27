global arcs = [1 2; 1 6; 1 7; 1 8; 1 21; 1 25; 1 29; 1 31; 1 49; 1 58; 1 59; 2 4; 2 17; 2 19; 2 30; 2 31; 2 40; 2 54; 3 31; 3 37; 3 40; 3 56; 3 57; 4 17; 4 33; 4 34; 4 59; 5 10; 5 11; 5 17; 5 34; 5 41; 6 3; 6 4; 6 17; 6 27; 6 39; 7 2; 7 10; 7 21; 7 27; 7 35; 7 40; 7 47; 7 53; 8 4; 8 10; 8 13; 8 19; 8 21; 8 22; 8 26; 8 36; 8 37; 8 51; 9 50; 9 59; 10 2; 10 3; 10 11; 10 27; 10 37; 10 55; 11 2; 11 6; 11 25; 11 31; 11 34; 11 39; 12 20; 12 21; 12 22; 12 31; 12 33; 12 39; 12 49; 12 54; 13 17; 13 34; 13 36; 13 42; 13 49; 13 51; 13 52; 13 54; 13 55; 13 56; 14 8; 14 20; 14 25; 14 36; 14 46; 14 54; 15 2; 15 5; 15 28; 15 49; 15 60; 16 45; 16 56; 16 60; 17 9; 17 12; 17 13; 17 20; 17 34; 17 38; 18 4; 18 5; 18 31; 18 36; 18 42; 18 50; 18 51; 18 53; 19 13; 19 17; 20 32; 20 38; 20 49; 21 27; 21 35; 21 54; 21 56; 22 3; 22 6; 22 7; 22 10; 22 26; 23 16; 23 17; 23 19; 23 22; 23 33; 23 34; 23 38; 23 55; 24 17; 24 18; 24 26; 24 49; 24 52; 25 4; 25 7; 25 37; 25 45; 25 55; 26 2; 26 14; 26 21; 26 25; 26 30; 26 47; 26 57; 27 8; 27 13; 27 23; 27 24; 27 26; 27 29; 27 51; 27 53; 27 60; 28 12; 28 38; 28 46; 28 48; 28 54; 29 6; 29 26; 29 44; 29 56; 30 8; 30 9; 30 16; 30 33; 30 37; 30 42; 30 49; 30 52; 30 58; 31 8; 31 20; 31 26; 31 60; 32 10; 32 11; 32 17; 32 22; 32 33; 32 34; 32 41; 32 54; 33 16; 33 19; 33 22; 33 23; 33 27; 33 51; 33 55; 33 60; 34 9; 34 12; 34 37; 34 44; 34 57; 34 59; 35 6; 35 15; 35 25; 35 47; 35 48; 36 16; 36 33; 36 42; 36 48; 36 58; 36 60; 37 16; 37 19; 37 26; 37 40; 38 9; 38 13; 38 21; 38 35; 38 47; 38 49; 38 52; 38 59; 39 19; 39 41; 39 43; 40 22; 40 33; 40 34; 41 4; 41 20; 41 23; 41 40; 41 42; 41 60; 42 16; 42 18; 42 25; 42 28; 42 36; 43 5; 43 29; 43 35; 43 46; 43 47; 43 53; 43 56; 44 11; 44 24; 44 42; 44 59; 45 37; 46 7; 46 24; 46 28; 46 32; 46 33; 46 50; 47 18; 47 43; 48 8; 48 11; 48 12; 48 22; 48 23; 48 24; 48 31; 48 42; 48 44; 48 57; 49 9; 49 10; 49 11; 49 30; 49 36; 49 40; 49 53; 49 54; 50 3; 50 19; 50 24; 50 27; 50 35; 50 48; 50 55; 50 58; 51 11; 51 28; 51 54; 51 57; 51 60; 52 14; 52 24; 52 25; 52 36; 52 56; 53 4; 53 14; 53 44; 53 60; 54 3; 54 31; 54 37; 54 41; 54 47; 54 50; 55 12; 55 21; 55 42; 55 50; 55 53; 55 59; 56 16; 56 36; 56 52; 57 16; 57 21; 57 39; 58 4; 58 14; 58 29; 58 52; 58 57; 59 8; 59 9; 59 20; 59 32; 59 49; 59 53]
global d_x = [2.0, 7.0, 7.0, 8.0, 6.0, 3.0, 10.0, 8.0, 9.0, 2.0, 10.0, 4.0, 6.0, 1.0, 5.0, 10.0, 10.0, 3.0, 3.0, 2.0, 3.0, 1.0, 6.0, 4.0, 6.0, 2.0, 3.0, 10.0, 1.0, 2.0, 5.0, 8.0, 7.0, 1.0, 8.0, 9.0, 6.0, 2.0, 7.0, 4.0, 6.0, 3.0, 9.0, 7.0, 6.0, 9.0, 9.0, 8.0, 10.0, 6.0, 4.0, 5.0, 10.0, 8.0, 5.0, 10.0, 3.0, 5.0, 5.0, 9.0, 8.0, 3.0, 6.0, 8.0, 10.0, 10.0, 3.0, 5.0, 2.0, 6.0, 1.0, 1.0, 8.0, 6.0, 6.0, 5.0, 4.0, 6.0, 5.0, 4.0, 1.0, 10.0, 7.0, 3.0, 9.0, 9.0, 1.0, 6.0, 3.0, 9.0, 6.0, 9.0, 8.0, 1.0, 3.0, 8.0, 5.0, 4.0, 6.0, 3.0, 3.0, 1.0, 3.0, 9.0, 2.0, 9.0, 1.0, 1.0, 6.0, 9.0, 1.0, 9.0, 9.0, 10.0, 7.0, 6.0, 6.0, 5.0, 3.0, 5.0, 5.0, 1.0, 9.0, 7.0, 1.0, 7.0, 9.0, 10.0, 2.0, 5.0, 9.0, 10.0, 3.0, 6.0, 9.0, 6.0, 1.0, 2.0, 8.0, 3.0, 3.0, 4.0, 7.0, 6.0, 1.0, 7.0, 8.0, 10.0, 1.0, 7.0, 4.0, 9.0, 10.0, 4.0, 3.0, 5.0, 2.0, 10.0, 6.0, 1.0, 4.0, 5.0, 8.0, 5.0, 1.0, 7.0, 10.0, 7.0, 1.0, 3.0, 2.0, 4.0, 6.0, 4.0, 10.0, 8.0, 9.0, 1.0, 4.0, 1.0, 3.0, 5.0, 4.0, 3.0, 10.0, 6.0, 5.0, 7.0, 1.0, 1.0, 9.0, 8.0, 1.0, 2.0, 4.0, 8.0, 9.0, 8.0, 6.0, 4.0, 9.0, 6.0, 1.0, 2.0, 4.0, 4.0, 4.0, 2.0, 8.0, 3.0, 5.0, 6.0, 2.0, 6.0, 6.0, 8.0, 6.0, 7.0, 1.0, 9.0, 8.0, 5.0, 2.0, 8.0, 5.0, 1.0, 3.0, 8.0, 2.0, 4.0, 6.0, 6.0, 7.0, 5.0, 6.0, 4.0, 10.0, 8.0, 5.0, 10.0, 1.0, 10.0, 5.0, 2.0, 10.0, 8.0, 4.0, 2.0, 3.0, 2.0, 6.0, 1.0, 8.0, 10.0, 5.0, 1.0, 8.0, 9.0, 8.0, 8.0, 1.0, 9.0, 7.0, 6.0, 6.0, 9.0, 1.0, 4.0, 3.0, 9.0, 1.0, 4.0, 4.0, 5.0, 9.0, 4.0, 2.0, 6.0, 8.0, 2.0, 1.0, 5.0, 7.0, 6.0, 1.0, 3.0, 9.0, 7.0, 9.0, 3.0, 4.0, 8.0, 6.0, 6.0, 10.0, 3.0, 1.0, 9.0, 3.0, 9.0, 4.0, 9.0, 5.0, 10.0, 10.0, 4.0, 6.0, 10.0, 3.0, 2.0, 10.0, 9.0, 6.0, 10.0, 6.0, 4.0, 10.0, 10.0, 1.0, 1.0, 4.0, 4.0, 6.0, 8.0, 10.0, 5.0, 9.0, 5.0, 2.0, 7.0, 10.0, 3.0, 1.0, 2.0, 2.0, 8.0]
global b_x = 5
global d_y = [7.0, 5.0, 9.0, 5.0, 9.0, 1.0, 9.0, 4.0, 7.0, 10.0, 8.0, 5.0, 1.0, 7.0, 10.0, 10.0, 2.0, 3.0, 2.0, 2.0, 5.0, 2.0, 4.0, 2.0, 3.0, 9.0, 7.0, 8.0, 2.0, 1.0, 4.0, 3.0, 3.0, 7.0, 4.0, 5.0, 2.0, 9.0, 8.0, 6.0, 7.0, 6.0, 5.0, 1.0, 1.0, 9.0, 10.0, 3.0, 8.0, 4.0, 3.0, 4.0, 1.0, 4.0, 2.0, 10.0, 5.0, 2.0, 4.0, 8.0, 5.0, 9.0, 9.0, 8.0, 3.0, 4.0, 9.0, 9.0, 8.0, 1.0, 5.0, 8.0, 7.0, 9.0, 6.0, 7.0, 2.0, 9.0, 10.0, 3.0, 4.0, 3.0, 5.0, 7.0, 7.0, 6.0, 3.0, 1.0, 9.0, 2.0, 2.0, 1.0, 2.0, 3.0, 10.0, 5.0, 4.0, 10.0, 1.0, 6.0, 1.0, 6.0, 1.0, 5.0, 2.0, 10.0, 6.0, 7.0, 6.0, 8.0, 5.0, 7.0, 2.0, 8.0, 6.0, 9.0, 8.0, 10.0, 2.0, 6.0, 6.0, 3.0, 2.0, 10.0, 8.0, 4.0, 8.0, 7.0, 4.0, 5.0, 1.0, 7.0, 9.0, 2.0, 7.0, 4.0, 1.0, 4.0, 3.0, 4.0, 5.0, 6.0, 8.0, 9.0, 7.0, 6.0, 3.0, 10.0, 10.0, 1.0, 1.0, 1.0, 5.0, 10.0, 1.0, 6.0, 8.0, 2.0, 2.0, 9.0, 9.0, 8.0, 1.0, 8.0, 7.0, 6.0, 9.0, 5.0, 8.0, 8.0, 5.0, 1.0, 4.0, 1.0, 10.0, 8.0, 5.0, 3.0, 9.0, 10.0, 10.0, 10.0, 3.0, 3.0, 5.0, 3.0, 4.0, 2.0, 2.0, 7.0, 4.0, 9.0, 2.0, 10.0, 5.0, 10.0, 2.0, 1.0, 3.0, 3.0, 5.0, 10.0, 1.0, 2.0, 2.0, 9.0, 9.0, 10.0, 5.0, 4.0, 9.0, 8.0, 6.0, 8.0, 4.0, 10.0, 6.0, 6.0, 2.0, 10.0, 9.0, 6.0, 7.0, 8.0, 1.0, 8.0, 1.0, 2.0, 7.0, 9.0, 9.0, 10.0, 7.0, 8.0, 1.0, 9.0, 3.0, 4.0, 9.0, 1.0, 4.0, 9.0, 3.0, 9.0, 7.0, 7.0, 5.0, 8.0, 2.0, 7.0, 5.0, 7.0, 9.0, 10.0, 7.0, 7.0, 8.0, 9.0, 9.0, 4.0, 6.0, 8.0, 2.0, 5.0, 10.0, 7.0, 10.0, 4.0, 5.0, 8.0, 2.0, 3.0, 3.0, 5.0, 9.0, 4.0, 5.0, 1.0, 3.0, 1.0, 6.0, 4.0, 6.0, 3.0, 7.0, 9.0, 9.0, 8.0, 9.0, 1.0, 10.0, 9.0, 5.0, 8.0, 6.0, 6.0, 8.0, 5.0, 6.0, 2.0, 7.0, 4.0, 7.0, 10.0, 1.0, 8.0, 8.0, 9.0, 3.0, 7.0, 3.0, 5.0, 2.0, 8.0, 7.0, 3.0, 4.0, 4.0, 9.0, 5.0, 2.0, 9.0, 8.0, 5.0, 7.0, 9.0, 5.0, 1.0, 2.0, 10.0, 2.0, 7.0, 1.0, 5.0, 9.0, 7.0]
global b_y = 10
global p = [0.374, 0.038, 0.574, 0.629, 0.31, 0.292, 0.026, 0.129, 0.169, 0.805, 0.046, 0.372, 0.185, 0.67, 0.03, 0.491, 0.312, 0.934, 0.003, 0.123, 0.925, 0.153, 0.242, 0.567, 0.553, 0.373, 0.17, 0.264, 0.603, 0.137, 0.31, 0.812, 0.128, 0.114, 0.474, 0.386, 0.335, 0.987, 0.476, 0.971, 0.746, 0.013, 0.065, 0.239, 0.172, 0.503, 0.227, 0.818, 0.261, 0.647, 0.653, 0.486, 0.824, 0.496, 0.468, 0.786, 0.351, 0.972, 0.025, 0.633, 0.97, 0.841, 0.819, 0.112, 0.392, 0.93, 0.764, 0.275, 0.259, 0.168, 0.619, 0.332, 0.97, 0.652, 0.653, 0.725, 0.81, 0.557, 0.605, 0.952, 0.245, 0.567, 0.279, 0.774, 0.518, 0.348, 0.784, 0.88, 0.573, 0.487, 0.442, 0.797, 0.035, 0.328, 0.378, 0.539, 0.344, 0.208, 0.908, 0.104, 0.609, 0.886, 0.508, 0.567, 0.453, 0.089, 0.997, 0.077, 0.238, 0.611, 0.228, 0.407, 0.859, 0.394, 0.541, 0.679, 0.18, 0.59, 0.201, 0.022, 0.133, 0.953, 0.445, 0.167, 0.518, 0.696, 0.514, 0.396, 0.307, 0.522, 0.507, 0.863, 0.574, 0.579, 0.343, 0.948, 0.25, 0.747, 0.283, 0.413, 0.219, 0.529, 0.538, 0.956, 0.347, 0.533, 0.369, 0.747, 0.753, 0.901, 0.314, 0.629, 0.335, 0.101, 0.253, 0.484, 0.179, 0.552, 0.151, 0.45, 0.572, 0.417, 0.412, 0.635, 0.547, 0.879, 0.494, 0.457, 0.895, 0.524, 0.459, 0.735, 0.018, 0.696, 0.82, 0.808, 0.178, 0.472, 0.087, 0.94, 0.468, 0.583, 0.313, 0.664, 0.223, 0.703, 0.806, 0.217, 0.235, 0.705, 0.32, 0.978, 0.081, 0.974, 0.789, 0.925, 0.694, 0.045, 0.005, 0.833, 0.871, 0.914, 0.712, 0.521, 0.899, 0.623, 0.268, 0.911, 0.783, 0.851, 0.084, 0.942, 0.437, 0.826, 0.66, 0.131, 0.119, 0.134, 0.897, 0.797, 0.001, 0.064, 0.184, 0.738, 0.182, 0.14, 0.902, 0.022, 0.563, 0.344, 0.339, 0.512, 0.249, 0.494, 0.874, 0.017, 0.167, 0.702, 0.114, 0.096, 0.634, 0.783, 0.76, 0.702, 0.496, 0.421, 0.541, 0.847, 0.79, 0.755, 0.48, 0.555, 0.437, 0.921, 0.74, 0.902, 0.662, 0.726, 0.903, 0.233, 0.63, 0.904, 0.001, 0.439, 0.339, 0.5, 0.052, 0.65, 0.935, 0.525, 0.155, 0.748, 0.729, 0.506, 0.686, 0.257, 0.306, 0.271, 0.978, 0.872, 0.925, 0.273, 0.019, 0.22, 0.112, 0.612, 0.552, 0.692, 0.859, 0.375, 0.176, 0.772, 0.962, 0.126, 0.029, 0.955, 0.319, 0.423, 0.236, 0.119, 0.184, 0.275, 0.235, 0.34, 0.906, 0.261, 0.421, 0.548, 0.441, 0.576, 0.703, 0.041, 0.239, 0.212, 0.7, 0.008, 0.756, 0.969, 0.666, 0.384, 0.239, 0.921, 0.303, 0.222, 0.866, 0.008, 0.59, 0.578, 0.088, 0.277, 0.034, 0.045, 0.876, 0.661, 0.691, 0.811]
global q = [0.383, 0.075, 0.994, 0.837, 0.328, 0.404, 0.759, 0.32, 0.384, 0.969, 0.446, 0.661, 0.741, 0.725, 0.552, 0.684, 0.429, 0.999, 0.804, 0.304, 0.972, 0.928, 0.918, 0.95, 0.747, 0.442, 0.569, 0.399, 0.753, 0.513, 0.912, 0.886, 0.205, 0.789, 0.558, 0.653, 0.902, 0.993, 0.917, 0.999, 0.869, 0.712, 0.285, 0.272, 0.329, 0.962, 0.531, 0.95, 0.343, 0.88, 0.96, 0.791, 0.903, 0.63, 0.534, 0.858, 0.482, 0.972, 0.316, 0.974, 0.989, 0.916, 0.838, 0.504, 0.948, 0.995, 0.777, 0.939, 0.365, 0.259, 0.862, 0.959, 0.992, 0.93, 0.799, 0.961, 0.886, 0.976, 0.922, 0.988, 0.519, 0.832, 0.692, 0.792, 0.915, 0.543, 0.818, 0.901, 0.913, 0.647, 0.718, 0.866, 0.945, 0.716, 0.45, 0.623, 0.443, 0.423, 0.957, 0.218, 0.828, 0.958, 0.869, 0.591, 0.731, 0.498, 0.997, 0.673, 0.402, 0.859, 0.257, 0.993, 0.902, 0.952, 0.815, 0.873, 0.187, 0.633, 0.536, 0.637, 0.164, 0.98, 0.985, 0.637, 0.857, 0.972, 0.56, 0.717, 0.672, 0.589, 0.671, 0.985, 0.684, 0.835, 0.978, 0.982, 0.253, 0.908, 0.352, 0.442, 0.342, 0.893, 0.671, 0.973, 0.842, 0.855, 0.503, 0.983, 0.973, 0.922, 0.374, 0.915, 0.785, 0.376, 0.644, 0.831, 0.84, 0.899, 0.422, 0.772, 0.986, 0.47, 0.463, 0.657, 0.94, 0.984, 0.591, 0.572, 0.903, 0.737, 0.865, 0.955, 0.179, 0.772, 0.915, 0.993, 0.805, 0.545, 0.422, 0.986, 0.942, 0.838, 0.728, 0.944, 0.411, 0.812, 0.848, 0.382, 0.297, 0.803, 0.524, 0.999, 0.198, 0.993, 0.95, 0.978, 0.932, 0.784, 0.641, 0.949, 0.976, 0.956, 0.735, 0.779, 0.986, 0.689, 0.379, 0.945, 0.867, 0.951, 0.298, 0.995, 0.971, 0.843, 0.753, 0.494, 0.609, 0.19, 0.938, 0.885, 0.067, 0.105, 0.951, 0.994, 0.645, 0.319, 0.915, 0.39, 0.655, 0.921, 0.517, 0.6, 0.325, 0.867, 0.957, 0.72, 0.515, 0.738, 0.931, 0.347, 0.976, 0.979, 0.997, 0.722, 0.504, 0.472, 0.603, 0.989, 0.804, 0.907, 0.617, 0.784, 0.743, 0.927, 0.902, 0.962, 0.736, 0.93, 0.945, 0.396, 0.864, 0.928, 0.891, 0.973, 0.991, 0.939, 0.401, 0.909, 0.946, 0.602, 0.157, 0.812, 0.955, 0.741, 0.752, 0.92, 0.766, 0.285, 0.987, 0.94, 0.943, 0.746, 0.569, 0.964, 0.893, 0.82, 0.837, 0.981, 0.878, 0.897, 0.585, 0.781, 0.98, 0.666, 0.99, 0.968, 0.845, 0.552, 0.461, 0.931, 0.64, 0.831, 0.639, 0.671, 0.927, 0.634, 0.546, 0.942, 0.615, 0.591, 0.919, 0.867, 0.355, 0.779, 0.829, 0.218, 0.845, 0.976, 0.791, 0.632, 0.411, 0.939, 0.635, 0.795, 0.869, 0.729, 0.697, 0.587, 0.409, 0.632, 0.386, 0.346, 0.967, 0.825, 0.826, 0.901]
global origin = 1
global destination = 60