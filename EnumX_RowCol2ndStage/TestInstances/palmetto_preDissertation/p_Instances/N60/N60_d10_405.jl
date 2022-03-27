global arcs = [1 4; 1 5; 1 11; 1 21; 1 25; 1 39; 1 42; 1 45; 2 12; 2 13; 2 17; 2 22; 2 29; 2 37; 2 39; 2 42; 2 56; 2 57; 3 6; 3 9; 3 19; 3 22; 3 26; 3 39; 3 42; 3 43; 3 45; 3 46; 3 47; 3 48; 3 51; 3 60; 4 10; 4 15; 4 20; 4 23; 4 29; 4 33; 4 51; 5 50; 6 23; 6 30; 6 31; 6 60; 7 5; 7 12; 7 14; 7 17; 7 26; 7 28; 7 32; 7 41; 7 47; 8 10; 8 11; 8 18; 8 19; 8 21; 8 45; 8 48; 8 57; 9 4; 9 7; 9 33; 9 50; 10 7; 10 17; 10 28; 10 50; 10 51; 10 54; 10 58; 11 5; 11 7; 11 10; 11 27; 11 43; 11 46; 11 57; 12 5; 12 6; 12 9; 12 26; 12 41; 12 44; 12 50; 12 52; 13 11; 13 32; 13 37; 13 40; 13 48; 13 49; 13 54; 13 57; 14 27; 14 39; 14 42; 14 44; 14 48; 15 5; 15 9; 15 12; 15 13; 15 26; 15 33; 15 40; 15 41; 16 4; 16 8; 16 10; 16 25; 16 53; 17 6; 17 9; 17 10; 17 11; 17 23; 17 28; 17 31; 17 44; 17 48; 17 51; 18 14; 18 21; 18 22; 18 33; 18 41; 19 7; 19 9; 19 13; 19 15; 19 23; 19 30; 19 31; 19 32; 19 39; 19 42; 19 51; 19 59; 19 60; 20 2; 20 25; 20 30; 20 37; 20 58; 21 20; 21 23; 21 27; 21 29; 21 40; 21 46; 22 6; 22 14; 22 39; 23 3; 23 5; 23 6; 23 20; 23 34; 23 55; 24 4; 24 16; 24 21; 24 23; 24 36; 24 39; 24 48; 24 53; 24 60; 25 6; 25 10; 25 11; 25 23; 25 24; 26 10; 26 19; 26 49; 27 12; 27 15; 27 21; 27 32; 27 43; 28 18; 28 34; 28 43; 28 47; 29 13; 29 16; 29 24; 29 34; 30 10; 30 13; 30 20; 30 29; 30 55; 31 14; 31 38; 31 58; 32 4; 32 5; 32 29; 32 31; 32 40; 32 51; 32 56; 33 2; 33 4; 33 21; 33 35; 33 59; 34 4; 34 27; 34 41; 34 44; 34 46; 34 60; 35 12; 35 37; 35 39; 35 40; 35 44; 35 45; 35 47; 35 51; 35 58; 36 16; 36 20; 36 23; 36 31; 36 43; 36 54; 37 20; 37 23; 37 24; 37 29; 37 31; 37 33; 37 54; 37 55; 37 56; 38 5; 38 15; 38 22; 38 35; 38 49; 38 54; 39 3; 39 15; 39 19; 39 21; 39 33; 39 34; 39 42; 39 54; 39 59; 40 27; 40 28; 40 34; 40 44; 40 50; 40 58; 41 2; 41 7; 41 24; 41 43; 41 53; 42 13; 42 25; 42 39; 42 48; 42 52; 42 57; 43 15; 43 24; 43 31; 43 46; 44 2; 44 7; 44 9; 44 12; 44 14; 44 16; 44 21; 44 35; 44 48; 45 18; 45 28; 45 47; 46 10; 46 15; 46 29; 46 36; 46 56; 47 2; 47 16; 47 46; 48 4; 48 7; 48 12; 48 15; 48 46; 48 56; 49 6; 49 7; 49 19; 49 25; 49 33; 49 47; 49 56; 50 12; 50 18; 50 25; 50 37; 50 47; 51 9; 51 14; 51 16; 51 32; 51 39; 51 56; 51 58; 51 59; 52 9; 52 32; 52 41; 53 17; 53 25; 54 14; 54 26; 54 29; 54 51; 54 59; 55 2; 55 5; 55 13; 55 16; 55 44; 56 12; 56 19; 56 22; 56 24; 56 40; 57 11; 57 17; 57 19; 57 30; 57 44; 57 58; 58 9; 58 19; 58 28; 59 2; 59 4; 59 6; 59 11; 59 21; 59 27; 59 30; 59 32; 59 58]
global d_x = [8.0, 7.0, 7.0, 5.0, 10.0, 4.0, 9.0, 7.0, 9.0, 7.0, 7.0, 3.0, 1.0, 2.0, 9.0, 9.0, 4.0, 1.0, 2.0, 7.0, 8.0, 7.0, 7.0, 10.0, 10.0, 10.0, 10.0, 2.0, 4.0, 5.0, 9.0, 9.0, 9.0, 3.0, 4.0, 5.0, 8.0, 3.0, 9.0, 5.0, 1.0, 3.0, 2.0, 3.0, 4.0, 7.0, 1.0, 8.0, 6.0, 10.0, 4.0, 9.0, 4.0, 8.0, 6.0, 10.0, 5.0, 4.0, 2.0, 5.0, 5.0, 2.0, 2.0, 8.0, 4.0, 6.0, 6.0, 8.0, 2.0, 8.0, 7.0, 4.0, 7.0, 5.0, 2.0, 2.0, 6.0, 8.0, 10.0, 4.0, 7.0, 9.0, 9.0, 6.0, 2.0, 4.0, 4.0, 3.0, 3.0, 2.0, 5.0, 4.0, 1.0, 10.0, 4.0, 10.0, 2.0, 4.0, 6.0, 9.0, 5.0, 8.0, 4.0, 5.0, 7.0, 2.0, 7.0, 7.0, 6.0, 1.0, 8.0, 6.0, 10.0, 9.0, 1.0, 1.0, 5.0, 3.0, 10.0, 5.0, 4.0, 7.0, 9.0, 2.0, 6.0, 10.0, 2.0, 3.0, 5.0, 8.0, 7.0, 3.0, 2.0, 4.0, 4.0, 7.0, 1.0, 8.0, 8.0, 8.0, 4.0, 4.0, 9.0, 1.0, 7.0, 8.0, 3.0, 7.0, 3.0, 5.0, 1.0, 2.0, 1.0, 3.0, 9.0, 8.0, 4.0, 8.0, 1.0, 9.0, 6.0, 5.0, 7.0, 2.0, 1.0, 5.0, 3.0, 4.0, 2.0, 9.0, 2.0, 2.0, 5.0, 1.0, 3.0, 3.0, 7.0, 3.0, 8.0, 1.0, 10.0, 4.0, 7.0, 8.0, 4.0, 3.0, 5.0, 2.0, 4.0, 1.0, 6.0, 1.0, 9.0, 6.0, 6.0, 5.0, 6.0, 4.0, 3.0, 10.0, 7.0, 9.0, 5.0, 9.0, 4.0, 7.0, 7.0, 10.0, 2.0, 2.0, 6.0, 1.0, 10.0, 4.0, 7.0, 9.0, 4.0, 9.0, 2.0, 10.0, 10.0, 8.0, 9.0, 3.0, 3.0, 1.0, 8.0, 8.0, 7.0, 7.0, 8.0, 1.0, 10.0, 10.0, 6.0, 1.0, 5.0, 4.0, 7.0, 3.0, 2.0, 3.0, 8.0, 2.0, 9.0, 3.0, 5.0, 7.0, 1.0, 7.0, 1.0, 3.0, 6.0, 9.0, 9.0, 3.0, 5.0, 7.0, 10.0, 6.0, 8.0, 3.0, 2.0, 3.0, 8.0, 1.0, 4.0, 1.0, 4.0, 3.0, 5.0, 3.0, 7.0, 5.0, 7.0, 9.0, 1.0, 7.0, 2.0, 9.0, 10.0, 4.0, 6.0, 7.0, 5.0, 4.0, 4.0, 7.0, 7.0, 4.0, 9.0, 3.0, 1.0, 2.0, 4.0, 3.0, 1.0, 4.0, 1.0, 10.0, 9.0, 4.0, 6.0, 7.0, 3.0, 1.0, 5.0, 5.0, 9.0, 6.0, 7.0, 2.0, 5.0, 9.0, 10.0, 10.0, 6.0, 8.0, 10.0, 2.0, 5.0, 3.0, 10.0, 3.0, 5.0, 7.0, 6.0, 3.0, 4.0, 7.0, 3.0, 2.0, 1.0, 10.0, 1.0, 8.0, 6.0, 10.0, 6.0, 7.0, 2.0, 4.0, 2.0, 4.0, 9.0, 1.0, 8.0, 3.0, 1.0, 4.0, 1.0, 7.0, 7.0, 2.0, 2.0, 8.0, 8.0, 5.0, 8.0, 9.0, 8.0]
global b_x = 5
global d_y = [6.0, 4.0, 4.0, 6.0, 5.0, 8.0, 1.0, 9.0, 7.0, 3.0, 3.0, 4.0, 9.0, 7.0, 1.0, 3.0, 5.0, 8.0, 2.0, 3.0, 3.0, 7.0, 2.0, 9.0, 6.0, 10.0, 7.0, 4.0, 6.0, 1.0, 9.0, 8.0, 10.0, 9.0, 6.0, 7.0, 7.0, 1.0, 9.0, 7.0, 1.0, 10.0, 4.0, 2.0, 4.0, 8.0, 9.0, 6.0, 5.0, 5.0, 7.0, 4.0, 7.0, 8.0, 2.0, 3.0, 6.0, 9.0, 6.0, 8.0, 9.0, 9.0, 4.0, 1.0, 7.0, 1.0, 8.0, 5.0, 7.0, 10.0, 10.0, 4.0, 9.0, 9.0, 4.0, 9.0, 7.0, 4.0, 5.0, 7.0, 3.0, 10.0, 6.0, 10.0, 10.0, 4.0, 4.0, 1.0, 9.0, 9.0, 6.0, 4.0, 6.0, 10.0, 7.0, 1.0, 4.0, 6.0, 1.0, 5.0, 1.0, 3.0, 4.0, 4.0, 4.0, 1.0, 1.0, 7.0, 5.0, 10.0, 5.0, 9.0, 5.0, 8.0, 1.0, 5.0, 3.0, 3.0, 5.0, 6.0, 8.0, 1.0, 9.0, 6.0, 8.0, 3.0, 9.0, 5.0, 1.0, 8.0, 1.0, 8.0, 10.0, 8.0, 9.0, 10.0, 5.0, 3.0, 5.0, 6.0, 2.0, 4.0, 10.0, 1.0, 3.0, 6.0, 4.0, 1.0, 9.0, 9.0, 8.0, 10.0, 10.0, 4.0, 10.0, 4.0, 3.0, 9.0, 1.0, 7.0, 5.0, 10.0, 6.0, 7.0, 1.0, 9.0, 3.0, 7.0, 1.0, 7.0, 5.0, 6.0, 2.0, 2.0, 5.0, 1.0, 1.0, 10.0, 10.0, 8.0, 9.0, 2.0, 10.0, 5.0, 1.0, 1.0, 2.0, 4.0, 7.0, 6.0, 2.0, 10.0, 5.0, 2.0, 2.0, 8.0, 6.0, 1.0, 6.0, 4.0, 1.0, 10.0, 6.0, 3.0, 1.0, 10.0, 8.0, 10.0, 6.0, 6.0, 3.0, 2.0, 3.0, 6.0, 8.0, 10.0, 6.0, 1.0, 2.0, 6.0, 9.0, 10.0, 3.0, 8.0, 5.0, 4.0, 3.0, 7.0, 2.0, 8.0, 8.0, 2.0, 10.0, 4.0, 5.0, 4.0, 2.0, 3.0, 9.0, 3.0, 8.0, 3.0, 5.0, 7.0, 6.0, 6.0, 1.0, 7.0, 9.0, 5.0, 8.0, 3.0, 8.0, 5.0, 2.0, 10.0, 8.0, 1.0, 5.0, 2.0, 7.0, 7.0, 5.0, 4.0, 5.0, 5.0, 2.0, 8.0, 9.0, 3.0, 8.0, 2.0, 4.0, 9.0, 4.0, 1.0, 7.0, 4.0, 8.0, 10.0, 2.0, 4.0, 5.0, 6.0, 5.0, 3.0, 10.0, 6.0, 10.0, 3.0, 4.0, 4.0, 2.0, 8.0, 6.0, 8.0, 7.0, 5.0, 7.0, 2.0, 7.0, 1.0, 2.0, 8.0, 1.0, 2.0, 9.0, 4.0, 7.0, 9.0, 9.0, 10.0, 2.0, 8.0, 9.0, 2.0, 5.0, 9.0, 1.0, 8.0, 9.0, 2.0, 9.0, 8.0, 7.0, 7.0, 2.0, 2.0, 3.0, 10.0, 3.0, 1.0, 6.0, 8.0, 7.0, 7.0, 2.0, 1.0, 10.0, 10.0, 10.0, 1.0, 4.0, 10.0, 6.0, 9.0, 5.0, 9.0, 5.0, 9.0, 7.0, 7.0, 7.0, 4.0, 3.0, 6.0, 9.0, 5.0, 2.0, 10.0, 2.0]
global b_y = 10
global p = [0.838, 0.792, 0.483, 0.268, 0.918, 0.026, 0.921, 0.181, 0.825, 0.15, 0.339, 0.537, 0.224, 0.598, 0.043, 0.289, 0.903, 0.241, 0.943, 0.631, 0.292, 0.677, 0.77, 0.971, 0.384, 0.029, 0.349, 0.14, 0.677, 0.27, 0.435, 0.299, 0.431, 0.348, 0.767, 0.306, 0.907, 0.144, 0.83, 0.084, 0.285, 0.529, 0.697, 0.663, 0.705, 0.819, 0.949, 0.004, 0.636, 0.398, 0.783, 0.111, 0.138, 0.35, 0.795, 0.138, 0.596, 0.525, 0.594, 0.301, 0.093, 0.538, 0.19, 0.482, 0.104, 0.999, 0.988, 0.541, 0.334, 0.829, 0.628, 0.676, 0.196, 0.994, 0.58, 0.838, 0.191, 0.497, 0.136, 0.625, 0.773, 0.356, 0.256, 0.53, 0.938, 0.859, 0.552, 0.946, 0.219, 0.747, 0.349, 0.237, 0.976, 0.309, 0.115, 0.081, 0.307, 0.982, 0.947, 0.772, 0.383, 0.076, 0.518, 0.27, 0.287, 0.678, 0.63, 0.73, 0.537, 0.713, 0.313, 0.765, 0.067, 0.889, 0.242, 0.722, 0.58, 0.505, 0.25, 0.975, 0.461, 0.972, 0.191, 0.3, 0.247, 0.757, 0.304, 0.015, 0.382, 0.004, 0.039, 0.63, 0.506, 0.444, 0.068, 0.973, 0.87, 0.024, 0.775, 0.091, 0.254, 0.376, 0.757, 0.426, 0.078, 0.334, 0.419, 0.599, 0.057, 0.503, 0.277, 0.622, 0.547, 0.597, 0.936, 0.335, 0.486, 0.711, 0.304, 0.474, 0.404, 0.107, 0.048, 0.146, 0.705, 0.377, 0.87, 0.922, 0.866, 0.582, 0.258, 0.996, 0.763, 0.403, 0.934, 0.594, 0.195, 0.463, 0.712, 0.34, 0.789, 0.557, 0.3, 0.138, 0.371, 0.287, 0.23, 0.858, 0.299, 0.389, 0.058, 0.979, 0.3, 0.005, 0.16, 0.368, 0.505, 0.21, 0.507, 0.61, 0.248, 0.864, 0.147, 0.454, 0.897, 0.065, 0.353, 0.389, 0.981, 0.328, 0.901, 0.377, 0.209, 0.65, 0.301, 0.413, 0.121, 0.048, 0.109, 0.897, 0.142, 0.458, 0.575, 0.132, 0.392, 0.203, 0.045, 0.336, 0.724, 0.703, 0.694, 0.595, 0.627, 0.077, 0.673, 0.624, 0.16, 0.409, 0.211, 0.757, 0.604, 0.075, 0.207, 0.466, 0.575, 0.488, 0.315, 0.944, 0.352, 0.172, 0.026, 0.058, 0.942, 0.989, 0.94, 0.733, 0.115, 0.563, 0.124, 0.592, 0.008, 0.743, 0.342, 0.702, 0.86, 0.231, 0.879, 0.115, 0.042, 0.83, 0.62, 0.582, 0.595, 0.515, 0.478, 0.6, 0.174, 0.434, 0.687, 0.98, 0.847, 0.454, 0.701, 0.114, 0.649, 0.514, 0.722, 0.177, 0.344, 0.342, 0.956, 0.098, 0.547, 0.043, 0.498, 0.66, 0.91, 0.567, 0.502, 0.141, 0.902, 0.328, 0.859, 0.671, 0.362, 0.922, 0.82, 0.205, 0.813, 0.828, 0.965, 0.988, 0.531, 0.71, 0.991, 0.926, 0.337, 0.889, 0.059, 0.897, 0.397, 0.044, 0.698, 0.723, 0.68, 0.973, 0.239, 0.453, 0.22, 0.083, 0.766, 0.797, 0.853, 0.019, 0.079, 0.563, 0.009, 0.818, 0.448, 0.137, 0.944, 0.208, 0.611, 0.318, 0.173, 0.708, 0.699, 0.889, 0.17, 0.506, 0.15, 0.763, 0.474, 0.127, 0.459, 0.895, 0.738, 0.546, 0.991, 0.834, 0.445]
global q = [0.925, 0.975, 0.632, 0.932, 0.987, 0.622, 0.989, 0.764, 0.982, 0.947, 0.985, 0.79, 0.64, 0.658, 0.547, 0.323, 0.96, 0.962, 0.959, 0.643, 0.474, 0.754, 0.94, 0.992, 0.746, 0.176, 0.93, 0.717, 0.714, 0.318, 0.555, 0.633, 0.904, 0.549, 0.996, 0.723, 0.987, 0.804, 0.971, 0.092, 0.589, 0.591, 0.994, 0.715, 0.996, 0.985, 0.959, 0.869, 0.964, 0.595, 0.881, 0.176, 0.326, 0.388, 0.854, 0.183, 0.912, 0.724, 0.884, 0.727, 0.289, 0.983, 0.396, 0.606, 0.654, 0.999, 0.997, 0.89, 0.665, 0.866, 0.86, 0.949, 0.38, 0.995, 0.734, 0.884, 0.731, 0.669, 0.596, 0.865, 0.993, 0.56, 0.522, 0.681, 0.985, 0.927, 0.581, 0.966, 0.675, 0.907, 0.579, 0.866, 0.996, 0.694, 0.271, 0.366, 0.504, 0.983, 0.96, 0.901, 0.999, 0.44, 0.73, 0.299, 0.931, 0.978, 0.772, 0.976, 0.903, 0.989, 0.503, 0.928, 0.246, 0.901, 0.632, 0.948, 0.672, 0.784, 0.518, 0.992, 0.643, 0.981, 0.206, 0.673, 0.388, 0.889, 0.449, 0.919, 0.882, 0.055, 0.822, 0.807, 0.89, 0.93, 0.353, 0.977, 0.992, 0.954, 0.957, 0.877, 0.995, 0.991, 0.848, 0.806, 0.818, 0.444, 0.788, 0.784, 0.423, 0.77, 0.968, 0.674, 0.741, 0.685, 0.952, 0.806, 0.88, 0.982, 0.925, 0.669, 0.47, 0.951, 0.493, 0.656, 0.971, 0.537, 0.904, 0.924, 0.998, 0.81, 0.481, 0.998, 0.907, 0.423, 0.978, 0.749, 0.707, 0.859, 0.812, 0.36, 0.9, 0.61, 0.549, 0.571, 0.933, 0.547, 0.356, 0.909, 0.345, 0.957, 0.952, 0.991, 0.851, 0.717, 0.79, 0.609, 0.699, 0.91, 0.72, 0.74, 0.396, 0.996, 0.97, 0.715, 0.918, 0.367, 0.374, 0.49, 0.994, 0.401, 0.937, 0.906, 0.803, 0.861, 0.513, 0.572, 0.577, 0.791, 0.813, 0.979, 0.383, 0.855, 0.895, 0.169, 0.711, 0.859, 0.994, 0.696, 0.791, 0.991, 0.7, 0.606, 0.876, 0.167, 0.72, 0.923, 0.465, 0.729, 0.222, 0.915, 0.613, 0.251, 0.307, 0.997, 0.659, 0.786, 0.447, 0.961, 0.584, 0.781, 0.822, 0.092, 0.998, 0.998, 0.982, 0.827, 0.688, 0.643, 0.134, 0.752, 0.094, 0.997, 0.634, 0.9, 0.975, 0.77, 0.978, 0.788, 0.309, 0.865, 0.768, 0.665, 0.602, 0.882, 0.729, 0.791, 0.322, 0.78, 0.83, 0.988, 0.876, 0.773, 0.946, 0.289, 0.989, 0.808, 0.93, 0.511, 0.799, 0.936, 0.969, 0.849, 0.599, 0.215, 0.569, 0.709, 0.976, 0.776, 0.564, 0.235, 0.913, 0.919, 0.915, 0.798, 0.392, 0.925, 0.851, 0.477, 0.87, 0.841, 0.985, 0.995, 0.616, 0.803, 0.993, 0.953, 0.904, 0.894, 0.624, 0.968, 0.777, 0.931, 0.725, 0.978, 0.713, 0.994, 0.748, 0.879, 0.356, 0.339, 0.808, 0.879, 0.899, 0.473, 0.706, 0.57, 0.022, 0.824, 0.599, 0.718, 0.989, 0.865, 0.937, 0.524, 0.609, 0.971, 0.832, 0.924, 0.403, 0.911, 0.505, 0.951, 0.784, 0.866, 0.996, 0.951, 0.9, 0.935, 0.994, 0.892, 0.772]
global origin = 1
global destination = 60