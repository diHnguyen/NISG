global arcs = [1 12; 1 16; 1 17; 1 28; 2 3; 2 20; 2 48; 3 2; 3 4; 3 5; 3 10; 3 11; 3 13; 3 27; 3 45; 3 46; 3 48; 4 21; 5 24; 6 10; 6 27; 6 28; 6 30; 6 34; 7 4; 7 10; 7 17; 7 28; 7 36; 7 38; 7 40; 7 47; 8 4; 8 27; 8 32; 8 34; 8 36; 9 11; 9 14; 9 17; 9 21; 9 23; 9 27; 9 38; 9 46; 10 6; 10 7; 10 29; 10 34; 10 38; 11 6; 11 18; 11 21; 11 45; 11 50; 12 14; 12 21; 12 33; 13 9; 13 14; 13 20; 13 30; 13 33; 13 34; 13 46; 14 10; 15 6; 15 12; 15 23; 15 32; 15 39; 15 45; 15 48; 16 11; 16 41; 17 4; 17 8; 17 16; 17 18; 17 32; 17 38; 18 4; 18 13; 18 22; 18 28; 18 37; 18 39; 18 46; 19 8; 19 30; 19 34; 19 41; 20 3; 20 8; 20 11; 20 25; 20 26; 20 27; 20 40; 20 49; 20 50; 21 5; 21 8; 21 14; 21 43; 22 21; 22 25; 22 30; 22 32; 22 44; 23 2; 23 10; 23 18; 23 32; 23 33; 23 46; 23 50; 24 2; 24 22; 24 29; 24 42; 24 47; 24 48; 25 3; 25 14; 25 23; 25 34; 25 45; 25 47; 26 4; 26 43; 26 49; 27 16; 27 17; 27 18; 27 35; 27 44; 27 48; 28 13; 28 22; 28 25; 28 30; 28 35; 28 37; 29 39; 29 49; 29 50; 30 7; 30 12; 30 15; 30 18; 30 28; 31 3; 31 7; 31 19; 31 24; 31 35; 31 40; 31 43; 31 49; 31 50; 32 11; 32 14; 32 17; 32 41; 32 45; 33 9; 33 10; 34 2; 34 3; 34 17; 34 44; 34 46; 35 31; 35 37; 35 50; 36 11; 36 15; 36 25; 36 37; 36 43; 36 50; 37 5; 37 13; 37 14; 37 18; 37 21; 37 23; 37 44; 38 3; 38 7; 38 12; 38 13; 38 19; 38 30; 38 48; 39 6; 39 7; 39 14; 39 33; 39 43; 39 47; 40 50; 41 27; 41 28; 42 35; 42 40; 43 3; 43 8; 43 13; 43 23; 43 28; 43 35; 43 40; 43 48; 44 22; 44 25; 44 28; 44 34; 45 23; 45 25; 45 28; 45 33; 45 36; 45 50; 46 2; 46 13; 46 43; 46 47; 46 50; 47 12; 47 26; 47 37; 47 39; 47 41; 47 44; 48 12; 48 22; 48 28; 48 29; 48 36; 48 47; 49 20; 49 34]
global d_x = [8.0, 7.0, 1.0, 6.0, 4.0, 7.0, 3.0, 6.0, 10.0, 9.0, 8.0, 7.0, 1.0, 5.0, 1.0, 2.0, 6.0, 1.0, 9.0, 6.0, 9.0, 1.0, 3.0, 3.0, 9.0, 7.0, 1.0, 7.0, 4.0, 6.0, 10.0, 1.0, 5.0, 9.0, 3.0, 8.0, 4.0, 1.0, 9.0, 5.0, 4.0, 4.0, 9.0, 6.0, 5.0, 8.0, 10.0, 1.0, 9.0, 10.0, 3.0, 7.0, 9.0, 1.0, 4.0, 10.0, 7.0, 9.0, 7.0, 10.0, 10.0, 4.0, 10.0, 7.0, 6.0, 2.0, 6.0, 6.0, 9.0, 3.0, 5.0, 10.0, 9.0, 7.0, 9.0, 1.0, 4.0, 7.0, 5.0, 3.0, 1.0, 7.0, 5.0, 6.0, 2.0, 3.0, 9.0, 1.0, 6.0, 1.0, 9.0, 4.0, 8.0, 7.0, 4.0, 7.0, 4.0, 5.0, 6.0, 4.0, 9.0, 9.0, 8.0, 8.0, 1.0, 7.0, 5.0, 7.0, 5.0, 2.0, 4.0, 1.0, 2.0, 1.0, 1.0, 7.0, 7.0, 5.0, 6.0, 8.0, 3.0, 1.0, 8.0, 1.0, 5.0, 7.0, 10.0, 6.0, 10.0, 1.0, 10.0, 6.0, 9.0, 7.0, 3.0, 7.0, 10.0, 1.0, 4.0, 4.0, 9.0, 2.0, 10.0, 9.0, 8.0, 5.0, 3.0, 4.0, 1.0, 3.0, 5.0, 5.0, 8.0, 9.0, 5.0, 7.0, 1.0, 8.0, 9.0, 6.0, 6.0, 7.0, 1.0, 8.0, 10.0, 8.0, 9.0, 7.0, 6.0, 6.0, 2.0, 6.0, 3.0, 10.0, 9.0, 6.0, 9.0, 5.0, 3.0, 9.0, 10.0, 8.0, 10.0, 10.0, 2.0, 3.0, 5.0, 2.0, 3.0, 9.0, 9.0, 2.0, 5.0, 5.0, 3.0, 3.0, 3.0, 7.0, 2.0, 1.0, 10.0, 5.0, 4.0, 2.0, 10.0, 3.0, 6.0, 10.0, 9.0, 6.0, 3.0, 7.0, 2.0, 7.0, 9.0, 8.0, 2.0, 6.0, 1.0, 5.0, 1.0, 1.0, 8.0, 4.0, 10.0, 3.0, 7.0, 5.0, 7.0, 2.0, 8.0, 9.0, 4.0, 2.0, 7.0, 4.0, 7.0, 7.0, 4.0, 4.0, 3.0, 7.0, 10.0, 10.0]
global b_x = 5
global d_y = [4.0, 6.0, 6.0, 2.0, 4.0, 7.0, 5.0, 4.0, 3.0, 6.0, 5.0, 6.0, 4.0, 5.0, 7.0, 8.0, 10.0, 8.0, 6.0, 5.0, 4.0, 1.0, 8.0, 4.0, 6.0, 9.0, 3.0, 10.0, 5.0, 3.0, 1.0, 3.0, 9.0, 7.0, 2.0, 8.0, 9.0, 9.0, 6.0, 7.0, 8.0, 5.0, 5.0, 2.0, 8.0, 2.0, 7.0, 3.0, 5.0, 10.0, 3.0, 6.0, 4.0, 3.0, 5.0, 7.0, 6.0, 7.0, 6.0, 4.0, 3.0, 6.0, 9.0, 5.0, 2.0, 5.0, 9.0, 3.0, 10.0, 8.0, 2.0, 9.0, 4.0, 2.0, 6.0, 6.0, 6.0, 1.0, 5.0, 10.0, 9.0, 3.0, 5.0, 4.0, 10.0, 8.0, 1.0, 5.0, 7.0, 6.0, 5.0, 2.0, 1.0, 10.0, 4.0, 1.0, 7.0, 10.0, 4.0, 3.0, 9.0, 10.0, 3.0, 7.0, 2.0, 10.0, 2.0, 3.0, 3.0, 1.0, 4.0, 7.0, 4.0, 5.0, 1.0, 7.0, 3.0, 2.0, 8.0, 1.0, 8.0, 2.0, 10.0, 6.0, 8.0, 3.0, 3.0, 6.0, 6.0, 6.0, 5.0, 10.0, 10.0, 3.0, 9.0, 7.0, 9.0, 1.0, 1.0, 10.0, 8.0, 7.0, 3.0, 10.0, 10.0, 4.0, 8.0, 2.0, 3.0, 8.0, 8.0, 8.0, 3.0, 4.0, 9.0, 3.0, 5.0, 2.0, 1.0, 8.0, 4.0, 5.0, 1.0, 5.0, 7.0, 9.0, 7.0, 3.0, 10.0, 3.0, 5.0, 8.0, 1.0, 10.0, 6.0, 6.0, 2.0, 9.0, 10.0, 6.0, 6.0, 7.0, 4.0, 5.0, 10.0, 2.0, 10.0, 3.0, 7.0, 1.0, 7.0, 10.0, 7.0, 2.0, 10.0, 4.0, 6.0, 8.0, 8.0, 10.0, 2.0, 8.0, 5.0, 8.0, 1.0, 10.0, 6.0, 7.0, 8.0, 5.0, 8.0, 5.0, 4.0, 3.0, 9.0, 4.0, 2.0, 8.0, 1.0, 5.0, 9.0, 6.0, 3.0, 2.0, 10.0, 2.0, 5.0, 4.0, 9.0, 1.0, 6.0, 9.0, 2.0, 2.0, 7.0, 6.0, 2.0, 6.0, 5.0, 1.0, 10.0, 9.0, 7.0, 6.0]
global b_y = 10
global p = [0.236, 0.657, 0.024, 0.131, 0.065, 0.888, 0.981, 0.781, 0.712, 0.555, 0.44, 0.768, 0.045, 0.311, 0.15, 0.576, 0.747, 0.088, 0.054, 0.259, 0.212, 0.289, 0.467, 0.962, 0.656, 0.087, 0.038, 0.492, 0.357, 0.537, 0.809, 0.724, 0.291, 0.725, 0.874, 0.519, 0.335, 0.684, 0.788, 0.583, 0.964, 0.027, 0.744, 0.29, 0.46, 0.236, 0.256, 0.504, 0.85, 0.953, 0.725, 0.482, 0.903, 0.12, 0.048, 0.053, 0.215, 0.595, 0.273, 0.526, 0.346, 0.317, 0.22, 0.173, 0.037, 0.608, 0.127, 0.789, 0.057, 0.41, 0.427, 0.312, 0.747, 0.069, 0.899, 0.416, 0.143, 0.499, 0.108, 0.82, 0.863, 0.313, 0.401, 0.104, 0.088, 0.567, 0.009, 0.422, 0.647, 0.674, 0.985, 0.325, 0.974, 0.423, 0.25, 0.224, 0.727, 0.529, 0.222, 0.029, 0.409, 0.409, 0.607, 0.49, 0.522, 0.152, 0.514, 0.336, 0.803, 0.304, 0.231, 0.648, 0.395, 0.667, 0.654, 0.893, 0.218, 0.495, 0.023, 0.876, 0.287, 0.144, 0.049, 0.363, 0.208, 0.281, 0.307, 0.508, 0.744, 0.521, 0.558, 0.222, 0.551, 0.491, 0.674, 0.415, 0.358, 0.434, 0.876, 0.972, 0.098, 0.954, 0.206, 0.497, 0.482, 0.127, 0.995, 0.152, 0.485, 0.435, 0.046, 0.497, 0.255, 0.573, 0.8, 0.144, 0.618, 0.456, 0.863, 0.142, 0.792, 0.915, 0.259, 0.194, 0.759, 0.08, 0.919, 0.995, 0.027, 0.068, 0.115, 0.778, 0.157, 0.511, 0.029, 0.367, 0.094, 0.637, 0.426, 0.55, 0.705, 0.725, 0.624, 0.726, 0.8, 0.008, 0.279, 0.652, 0.663, 0.973, 0.869, 0.588, 0.316, 0.44, 0.522, 0.428, 0.243, 0.059, 0.08, 0.428, 0.887, 0.181, 0.05, 0.18, 0.027, 0.789, 0.169, 0.158, 0.892, 0.096, 0.17, 0.534, 0.716, 0.866, 0.161, 0.871, 0.221, 0.025, 0.721, 0.137, 0.093, 0.527, 0.831, 0.345, 0.431, 0.662, 0.154, 0.408, 0.826, 0.981, 0.65, 0.944, 0.939, 0.739, 0.752, 0.391, 0.154, 0.741, 0.164, 0.339, 0.098, 0.853, 0.742, 0.987]
global q = [0.61, 0.947, 0.343, 0.234, 0.616, 0.957, 0.982, 0.879, 0.93, 0.719, 0.995, 0.971, 0.472, 0.886, 0.876, 0.805, 0.859, 0.48, 0.211, 0.313, 0.995, 0.735, 0.896, 0.962, 0.839, 0.601, 0.113, 0.644, 0.358, 0.877, 0.876, 0.793, 0.394, 0.752, 0.905, 0.845, 0.553, 0.773, 0.936, 0.649, 0.974, 0.976, 0.825, 0.722, 0.65, 0.537, 0.299, 0.685, 0.975, 0.978, 0.962, 0.945, 0.992, 0.356, 0.736, 0.119, 0.628, 0.91, 0.544, 0.812, 0.447, 0.68, 0.739, 0.632, 0.6, 0.616, 0.842, 0.881, 0.645, 0.497, 0.556, 0.529, 0.923, 0.361, 0.967, 0.852, 0.4, 0.625, 0.767, 0.915, 0.891, 0.345, 0.772, 0.779, 0.315, 0.63, 0.56, 0.959, 0.961, 0.967, 0.993, 0.872, 0.974, 0.982, 0.684, 0.556, 0.916, 0.574, 0.329, 0.469, 0.61, 0.437, 0.87, 0.771, 0.775, 0.436, 0.6, 0.35, 0.957, 0.343, 0.81, 0.975, 0.659, 0.886, 0.772, 0.904, 0.791, 0.901, 0.628, 0.89, 0.775, 0.817, 0.966, 0.622, 0.746, 0.343, 0.997, 0.78, 0.913, 0.817, 0.994, 0.536, 0.841, 0.725, 0.963, 0.89, 0.561, 0.488, 0.882, 0.983, 0.995, 0.954, 0.936, 0.66, 0.628, 0.868, 0.999, 0.76, 0.686, 0.621, 0.557, 0.787, 0.462, 0.687, 0.952, 0.194, 0.785, 0.647, 0.937, 0.692, 0.866, 0.96, 0.587, 0.66, 0.831, 0.29, 0.924, 0.995, 0.194, 0.517, 0.506, 0.961, 0.708, 0.91, 0.816, 0.415, 0.106, 0.816, 0.967, 0.813, 0.893, 0.895, 0.785, 0.872, 0.862, 0.559, 0.602, 0.871, 0.734, 0.981, 0.946, 0.78, 0.448, 0.896, 0.704, 0.654, 0.81, 0.343, 0.371, 0.776, 0.951, 0.839, 0.104, 0.512, 0.894, 0.947, 0.454, 0.302, 0.992, 0.408, 0.365, 0.727, 0.867, 0.919, 0.493, 0.976, 0.732, 0.718, 0.792, 0.195, 0.33, 0.579, 0.952, 0.824, 0.867, 0.944, 0.83, 0.504, 0.934, 0.994, 0.797, 0.984, 0.977, 0.965, 0.875, 0.638, 0.881, 0.805, 0.432, 0.654, 0.366, 0.916, 0.753, 0.998]
global origin = 1
global destination = 50