global arcs = [1 8; 1 11; 1 19; 1 24; 1 26; 1 30; 1 35; 1 39; 1 42; 2 19; 2 24; 2 31; 2 38; 3 5; 3 8; 3 11; 3 12; 3 14; 3 15; 3 37; 4 10; 4 38; 4 42; 4 46; 5 34; 6 8; 6 11; 6 21; 6 25; 6 45; 6 46; 6 47; 7 4; 7 14; 7 41; 8 13; 8 26; 8 34; 8 48; 9 20; 9 42; 9 50; 10 9; 10 20; 10 22; 10 29; 10 30; 10 34; 10 39; 10 47; 11 37; 12 4; 12 5; 12 8; 12 24; 12 25; 13 9; 13 15; 13 42; 14 12; 14 28; 14 45; 15 2; 15 4; 15 11; 15 19; 15 24; 15 28; 15 33; 15 41; 15 42; 16 11; 16 13; 16 17; 16 24; 16 37; 16 41; 16 42; 17 8; 17 16; 17 23; 17 31; 17 46; 18 6; 18 15; 18 22; 18 28; 18 37; 19 15; 19 30; 19 45; 20 3; 20 14; 20 19; 20 42; 20 43; 20 50; 21 17; 21 20; 21 23; 21 27; 21 44; 22 6; 22 13; 22 14; 22 15; 22 18; 22 19; 22 25; 22 35; 22 37; 23 5; 23 31; 23 36; 23 42; 24 2; 24 4; 24 14; 24 18; 24 38; 24 44; 24 48; 25 4; 25 8; 25 28; 25 30; 25 35; 25 37; 26 6; 26 20; 26 28; 26 34; 26 50; 27 8; 27 14; 27 28; 27 35; 27 40; 27 45; 27 47; 28 6; 28 16; 28 20; 28 23; 28 27; 28 31; 28 39; 28 42; 28 44; 28 47; 28 50; 29 27; 29 40; 29 45; 30 3; 30 4; 30 20; 30 23; 30 35; 30 45; 30 50; 31 12; 31 26; 31 28; 31 50; 32 26; 32 41; 32 45; 32 48; 33 12; 33 17; 33 20; 33 30; 33 35; 33 41; 34 9; 34 18; 34 21; 34 26; 34 28; 34 44; 35 18; 35 29; 35 32; 35 47; 36 9; 36 33; 37 2; 37 33; 37 36; 37 40; 37 41; 38 24; 38 26; 38 34; 38 39; 38 44; 39 11; 39 15; 39 16; 39 23; 39 40; 39 46; 40 4; 40 9; 40 18; 40 19; 40 31; 40 32; 40 48; 41 22; 41 30; 41 47; 41 49; 42 9; 42 10; 42 11; 42 16; 42 35; 42 44; 42 46; 43 14; 43 22; 43 36; 43 46; 44 2; 44 6; 44 7; 44 25; 44 31; 44 40; 44 41; 44 45; 45 9; 45 18; 45 21; 45 25; 45 42; 46 7; 46 19; 46 31; 46 48; 47 3; 47 4; 47 10; 47 16; 47 49; 48 43; 49 10; 49 13; 49 22; 49 32; 49 40]
global d_x = [10.0, 1.0, 10.0, 5.0, 1.0, 1.0, 9.0, 8.0, 5.0, 3.0, 9.0, 8.0, 4.0, 9.0, 8.0, 9.0, 2.0, 1.0, 6.0, 1.0, 6.0, 1.0, 4.0, 8.0, 4.0, 8.0, 9.0, 9.0, 2.0, 10.0, 7.0, 3.0, 5.0, 8.0, 7.0, 10.0, 9.0, 2.0, 4.0, 8.0, 2.0, 5.0, 4.0, 10.0, 4.0, 9.0, 7.0, 6.0, 2.0, 4.0, 9.0, 7.0, 6.0, 3.0, 8.0, 1.0, 5.0, 2.0, 3.0, 4.0, 8.0, 2.0, 2.0, 7.0, 9.0, 1.0, 7.0, 10.0, 5.0, 6.0, 10.0, 2.0, 2.0, 1.0, 1.0, 10.0, 2.0, 8.0, 6.0, 8.0, 8.0, 2.0, 5.0, 10.0, 6.0, 1.0, 9.0, 7.0, 4.0, 1.0, 9.0, 5.0, 9.0, 9.0, 4.0, 9.0, 2.0, 7.0, 10.0, 10.0, 3.0, 8.0, 10.0, 8.0, 4.0, 10.0, 3.0, 1.0, 5.0, 8.0, 8.0, 6.0, 3.0, 2.0, 2.0, 9.0, 1.0, 10.0, 10.0, 9.0, 3.0, 7.0, 4.0, 4.0, 5.0, 4.0, 4.0, 2.0, 6.0, 6.0, 10.0, 6.0, 3.0, 10.0, 7.0, 7.0, 6.0, 7.0, 7.0, 5.0, 4.0, 4.0, 5.0, 4.0, 10.0, 1.0, 5.0, 10.0, 3.0, 10.0, 3.0, 6.0, 3.0, 7.0, 4.0, 4.0, 4.0, 3.0, 9.0, 2.0, 9.0, 9.0, 5.0, 1.0, 4.0, 5.0, 4.0, 8.0, 9.0, 7.0, 2.0, 8.0, 4.0, 7.0, 2.0, 5.0, 7.0, 9.0, 3.0, 5.0, 10.0, 9.0, 2.0, 7.0, 9.0, 5.0, 5.0, 6.0, 4.0, 7.0, 9.0, 9.0, 6.0, 10.0, 4.0, 1.0, 3.0, 4.0, 5.0, 7.0, 9.0, 7.0, 2.0, 8.0, 9.0, 10.0, 6.0, 2.0, 10.0, 2.0, 10.0, 5.0, 5.0, 10.0, 2.0, 4.0, 3.0, 5.0, 4.0, 1.0, 3.0, 4.0, 1.0, 4.0, 8.0, 3.0, 7.0, 5.0, 9.0, 8.0, 10.0, 3.0, 2.0, 1.0, 7.0, 9.0, 4.0, 9.0, 10.0, 3.0, 4.0, 5.0, 2.0, 3.0, 7.0, 7.0, 2.0, 7.0, 5.0, 1.0, 2.0, 5.0, 5.0]
global b_x = 5
global d_y = [7.0, 9.0, 8.0, 9.0, 10.0, 9.0, 9.0, 8.0, 1.0, 6.0, 2.0, 3.0, 9.0, 10.0, 4.0, 9.0, 1.0, 4.0, 6.0, 6.0, 3.0, 5.0, 3.0, 9.0, 2.0, 2.0, 4.0, 4.0, 1.0, 4.0, 7.0, 9.0, 3.0, 1.0, 10.0, 3.0, 5.0, 1.0, 10.0, 3.0, 6.0, 10.0, 1.0, 6.0, 6.0, 6.0, 4.0, 5.0, 1.0, 7.0, 5.0, 8.0, 5.0, 7.0, 7.0, 9.0, 2.0, 7.0, 6.0, 7.0, 5.0, 2.0, 6.0, 3.0, 2.0, 1.0, 1.0, 4.0, 4.0, 6.0, 10.0, 9.0, 9.0, 2.0, 9.0, 4.0, 7.0, 2.0, 9.0, 2.0, 10.0, 7.0, 8.0, 9.0, 9.0, 5.0, 2.0, 7.0, 2.0, 5.0, 8.0, 8.0, 7.0, 6.0, 1.0, 10.0, 9.0, 5.0, 4.0, 5.0, 8.0, 6.0, 5.0, 1.0, 3.0, 5.0, 10.0, 3.0, 5.0, 5.0, 10.0, 7.0, 1.0, 4.0, 3.0, 10.0, 8.0, 9.0, 6.0, 4.0, 5.0, 5.0, 6.0, 2.0, 4.0, 4.0, 9.0, 1.0, 8.0, 7.0, 5.0, 9.0, 3.0, 3.0, 9.0, 1.0, 9.0, 8.0, 2.0, 3.0, 6.0, 1.0, 1.0, 5.0, 7.0, 5.0, 6.0, 10.0, 10.0, 4.0, 5.0, 1.0, 7.0, 4.0, 2.0, 10.0, 5.0, 7.0, 7.0, 8.0, 10.0, 4.0, 6.0, 1.0, 6.0, 10.0, 6.0, 10.0, 8.0, 6.0, 5.0, 10.0, 6.0, 6.0, 9.0, 8.0, 2.0, 10.0, 8.0, 4.0, 4.0, 2.0, 2.0, 5.0, 7.0, 2.0, 7.0, 1.0, 1.0, 4.0, 2.0, 5.0, 6.0, 3.0, 1.0, 8.0, 9.0, 8.0, 1.0, 1.0, 3.0, 10.0, 3.0, 10.0, 3.0, 9.0, 1.0, 5.0, 5.0, 8.0, 4.0, 6.0, 1.0, 1.0, 5.0, 3.0, 4.0, 3.0, 9.0, 3.0, 8.0, 8.0, 4.0, 6.0, 4.0, 3.0, 2.0, 4.0, 9.0, 1.0, 2.0, 1.0, 7.0, 8.0, 4.0, 4.0, 7.0, 5.0, 5.0, 4.0, 9.0, 5.0, 6.0, 1.0, 10.0, 2.0, 6.0, 4.0, 9.0, 3.0, 4.0, 4.0, 5.0]
global b_y = 10
global p = [0.457, 0.796, 0.591, 0.215, 0.435, 0.626, 0.105, 0.777, 0.714, 0.528, 0.013, 0.287, 0.902, 0.836, 0.082, 0.935, 0.3, 0.723, 0.188, 0.8, 0.871, 0.053, 0.453, 0.004, 0.954, 0.551, 0.716, 0.715, 0.15, 0.09, 0.27, 0.318, 0.313, 0.099, 0.164, 0.059, 0.008, 0.825, 0.93, 0.545, 0.887, 0.848, 0.318, 0.612, 0.918, 0.204, 0.774, 0.534, 0.625, 0.348, 0.024, 0.755, 0.967, 0.478, 0.713, 0.217, 0.414, 0.001, 0.82, 0.195, 0.387, 0.786, 0.834, 0.268, 0.035, 0.73, 0.641, 0.962, 0.401, 0.569, 0.212, 0.793, 0.614, 0.266, 0.698, 0.721, 0.735, 0.562, 0.863, 0.88, 0.269, 0.462, 0.394, 0.429, 0.611, 0.055, 0.402, 0.859, 0.967, 0.938, 0.583, 0.569, 0.507, 0.717, 0.395, 0.726, 0.179, 0.818, 0.035, 0.791, 0.737, 0.485, 0.732, 0.679, 0.002, 0.106, 0.294, 0.908, 0.328, 0.737, 0.043, 0.759, 0.028, 0.267, 0.616, 0.131, 0.319, 0.01, 0.465, 0.347, 0.672, 0.449, 0.108, 0.024, 0.273, 0.823, 0.361, 0.066, 0.574, 0.526, 0.776, 0.82, 0.781, 0.824, 0.103, 0.425, 0.083, 0.508, 0.689, 0.706, 0.56, 0.945, 0.16, 0.773, 0.462, 0.363, 0.941, 0.16, 0.141, 0.005, 0.288, 0.442, 0.228, 0.953, 0.659, 0.093, 0.207, 0.251, 0.814, 0.04, 0.408, 0.789, 0.536, 0.713, 0.393, 0.268, 0.415, 0.662, 0.291, 0.301, 0.196, 0.557, 0.894, 0.333, 0.925, 0.627, 0.744, 0.488, 0.645, 0.653, 0.328, 0.404, 0.823, 0.721, 0.483, 0.789, 0.398, 0.034, 0.403, 0.342, 0.396, 0.811, 0.933, 0.078, 0.574, 0.453, 0.263, 0.439, 0.612, 0.84, 0.441, 0.9, 0.838, 0.203, 0.855, 0.35, 0.178, 0.677, 0.438, 0.482, 0.106, 0.614, 0.936, 0.641, 0.993, 0.581, 0.653, 0.934, 0.944, 0.849, 0.319, 0.352, 0.109, 0.436, 0.583, 0.889, 0.734, 0.747, 0.065, 0.143, 0.981, 0.55, 0.973, 0.679, 0.439, 0.949, 0.772, 0.015, 0.413, 0.435, 0.205, 0.664, 0.368, 0.927, 0.708, 0.352, 0.438, 0.714, 0.685, 0.119, 0.543, 0.667, 0.68]
global q = [0.471, 0.884, 0.765, 0.521, 0.671, 0.987, 0.9, 0.934, 0.754, 0.705, 0.511, 0.333, 0.993, 0.961, 0.953, 0.966, 0.955, 0.881, 0.662, 0.998, 0.948, 0.126, 0.841, 0.763, 0.981, 0.875, 0.916, 0.861, 0.279, 0.168, 0.875, 0.713, 0.696, 0.814, 0.276, 0.399, 0.818, 0.986, 0.99, 0.583, 0.996, 0.864, 0.65, 0.961, 0.984, 0.575, 0.856, 0.603, 0.89, 0.523, 0.37, 0.843, 0.971, 0.861, 0.932, 0.416, 0.806, 0.376, 0.862, 0.933, 0.481, 0.836, 0.84, 0.28, 0.707, 0.906, 0.75, 0.979, 0.99, 0.782, 0.647, 0.844, 0.695, 0.508, 0.804, 0.931, 0.745, 0.699, 0.968, 0.893, 0.453, 0.833, 0.545, 0.512, 0.955, 0.471, 0.698, 0.96, 0.998, 0.942, 0.722, 0.782, 0.578, 0.823, 0.742, 0.906, 0.887, 0.886, 0.242, 0.889, 0.842, 0.903, 0.789, 0.895, 0.593, 0.997, 0.882, 0.986, 0.528, 0.925, 0.134, 0.915, 0.106, 0.689, 0.757, 0.619, 0.415, 0.049, 0.612, 0.379, 0.836, 0.458, 0.891, 0.313, 0.683, 0.901, 0.935, 0.401, 0.969, 0.987, 0.778, 0.969, 0.979, 0.901, 0.579, 0.497, 0.324, 0.669, 0.712, 0.894, 0.759, 0.954, 0.498, 0.899, 0.555, 0.659, 0.985, 0.429, 0.19, 0.84, 0.363, 0.561, 0.926, 0.982, 0.727, 0.105, 0.924, 0.302, 0.864, 0.993, 0.919, 0.934, 0.95, 0.777, 0.85, 0.927, 0.961, 0.945, 0.964, 0.544, 0.593, 0.57, 0.986, 0.463, 0.956, 0.857, 0.772, 0.873, 0.871, 0.963, 0.491, 0.437, 0.915, 0.762, 0.801, 0.873, 0.486, 0.341, 0.578, 0.432, 0.788, 0.916, 0.971, 0.102, 0.793, 0.685, 0.436, 0.935, 0.834, 0.964, 0.56, 0.988, 0.885, 0.441, 0.92, 0.519, 0.72, 0.772, 0.658, 0.549, 0.285, 0.775, 0.976, 0.646, 0.997, 0.966, 0.92, 0.962, 0.955, 0.889, 0.517, 0.452, 0.849, 0.625, 0.695, 0.951, 0.997, 0.963, 0.776, 0.527, 0.989, 0.916, 0.982, 0.993, 0.592, 0.963, 0.789, 0.533, 0.992, 0.928, 0.483, 0.855, 0.489, 0.934, 0.778, 0.692, 0.871, 0.961, 0.771, 0.823, 0.723, 0.678, 0.824]
global origin = 1
global destination = 50