global arcs = [1 33; 1 34; 1 37; 2 17; 2 28; 3 11; 3 16; 3 35; 3 36; 4 2; 4 6; 4 36; 5 14; 5 15; 5 40; 6 11; 6 15; 6 28; 6 29; 6 41; 6 50; 7 3; 7 4; 7 11; 7 12; 7 19; 7 35; 7 40; 7 42; 7 46; 7 50; 8 7; 8 10; 8 24; 8 37; 8 39; 9 15; 9 18; 9 22; 9 26; 9 27; 9 46; 9 47; 10 8; 10 28; 10 47; 11 29; 11 40; 11 45; 11 48; 12 37; 12 41; 12 44; 13 2; 13 34; 13 37; 13 39; 13 40; 13 41; 14 20; 14 29; 14 32; 14 40; 15 14; 15 17; 15 32; 16 2; 16 13; 16 18; 16 37; 17 6; 17 10; 17 13; 17 19; 17 35; 18 10; 18 11; 18 30; 18 34; 18 38; 18 43; 18 48; 19 25; 19 38; 19 40; 20 3; 20 27; 20 32; 20 47; 21 19; 21 22; 21 40; 21 46; 22 9; 22 11; 22 23; 22 27; 22 49; 23 5; 23 50; 24 29; 24 32; 24 50; 25 3; 25 12; 25 17; 25 36; 26 3; 26 9; 26 13; 26 19; 26 20; 26 44; 27 2; 27 4; 27 11; 27 39; 28 13; 28 15; 28 20; 28 35; 29 5; 29 10; 29 13; 29 22; 29 32; 29 34; 29 38; 29 41; 30 18; 30 35; 30 37; 30 38; 30 44; 30 45; 30 49; 31 2; 31 5; 31 6; 31 7; 31 9; 31 13; 31 26; 31 39; 31 41; 31 44; 32 8; 32 13; 32 22; 32 26; 32 40; 33 7; 33 27; 33 28; 33 47; 34 22; 34 24; 34 40; 34 49; 34 50; 35 8; 35 27; 35 40; 36 10; 36 14; 36 16; 36 23; 36 28; 36 33; 36 42; 36 43; 36 46; 36 48; 37 11; 37 21; 37 42; 38 50; 39 8; 39 13; 39 40; 40 9; 40 42; 40 47; 41 10; 41 11; 41 17; 41 19; 41 22; 41 25; 41 31; 42 7; 42 27; 42 37; 43 10; 43 15; 43 17; 43 21; 43 31; 43 46; 44 28; 45 10; 45 14; 45 18; 45 19; 45 25; 45 30; 45 48; 46 4; 46 6; 46 16; 46 17; 46 38; 46 39; 47 4; 47 34; 47 39; 47 45; 47 48; 47 50; 48 11; 48 21; 48 36; 48 37; 48 42; 48 47; 49 2; 49 16; 49 22; 49 26; 49 35; 49 38]
global d_x = [5.0, 9.0, 10.0, 5.0, 7.0, 6.0, 2.0, 5.0, 9.0, 1.0, 1.0, 7.0, 7.0, 3.0, 2.0, 2.0, 7.0, 6.0, 5.0, 1.0, 10.0, 2.0, 6.0, 6.0, 7.0, 7.0, 2.0, 1.0, 7.0, 3.0, 5.0, 6.0, 3.0, 7.0, 10.0, 6.0, 9.0, 9.0, 9.0, 10.0, 4.0, 5.0, 1.0, 5.0, 1.0, 2.0, 2.0, 7.0, 5.0, 2.0, 3.0, 4.0, 3.0, 9.0, 6.0, 1.0, 8.0, 2.0, 3.0, 10.0, 8.0, 7.0, 3.0, 8.0, 10.0, 8.0, 2.0, 2.0, 10.0, 1.0, 1.0, 2.0, 10.0, 9.0, 2.0, 4.0, 6.0, 1.0, 2.0, 10.0, 1.0, 7.0, 10.0, 6.0, 9.0, 3.0, 7.0, 8.0, 10.0, 4.0, 6.0, 6.0, 6.0, 4.0, 8.0, 6.0, 6.0, 7.0, 8.0, 4.0, 7.0, 10.0, 5.0, 9.0, 5.0, 6.0, 8.0, 5.0, 1.0, 7.0, 1.0, 9.0, 8.0, 9.0, 8.0, 3.0, 10.0, 6.0, 10.0, 1.0, 7.0, 6.0, 9.0, 10.0, 7.0, 2.0, 5.0, 6.0, 4.0, 9.0, 10.0, 10.0, 4.0, 9.0, 6.0, 5.0, 1.0, 2.0, 5.0, 10.0, 6.0, 6.0, 6.0, 9.0, 6.0, 1.0, 2.0, 6.0, 3.0, 7.0, 10.0, 4.0, 2.0, 9.0, 8.0, 1.0, 3.0, 9.0, 8.0, 6.0, 6.0, 4.0, 10.0, 10.0, 1.0, 4.0, 2.0, 10.0, 8.0, 2.0, 9.0, 7.0, 6.0, 1.0, 2.0, 3.0, 1.0, 8.0, 3.0, 1.0, 5.0, 6.0, 4.0, 1.0, 9.0, 6.0, 5.0, 4.0, 6.0, 1.0, 9.0, 1.0, 9.0, 7.0, 2.0, 8.0, 6.0, 8.0, 8.0, 4.0, 9.0, 7.0, 7.0, 8.0, 3.0, 9.0, 8.0, 2.0, 10.0, 2.0, 10.0, 3.0, 8.0, 4.0, 6.0, 1.0, 2.0, 6.0, 9.0, 8.0, 2.0, 4.0, 10.0, 4.0, 7.0, 5.0, 6.0, 10.0, 9.0, 3.0, 9.0]
global b_x = 5
global d_y = [4.0, 6.0, 9.0, 6.0, 6.0, 2.0, 3.0, 2.0, 5.0, 2.0, 4.0, 4.0, 7.0, 5.0, 7.0, 2.0, 5.0, 1.0, 1.0, 10.0, 6.0, 2.0, 9.0, 1.0, 5.0, 3.0, 8.0, 8.0, 1.0, 5.0, 3.0, 9.0, 4.0, 8.0, 4.0, 4.0, 8.0, 10.0, 10.0, 7.0, 7.0, 4.0, 3.0, 5.0, 10.0, 4.0, 5.0, 7.0, 10.0, 9.0, 5.0, 2.0, 3.0, 7.0, 1.0, 9.0, 1.0, 3.0, 9.0, 9.0, 7.0, 1.0, 2.0, 3.0, 8.0, 1.0, 10.0, 5.0, 2.0, 3.0, 2.0, 4.0, 2.0, 7.0, 5.0, 1.0, 5.0, 9.0, 3.0, 6.0, 8.0, 4.0, 8.0, 3.0, 5.0, 5.0, 5.0, 9.0, 4.0, 6.0, 6.0, 10.0, 5.0, 5.0, 7.0, 10.0, 2.0, 10.0, 6.0, 5.0, 6.0, 7.0, 6.0, 8.0, 3.0, 5.0, 9.0, 6.0, 10.0, 1.0, 9.0, 10.0, 2.0, 8.0, 7.0, 2.0, 9.0, 9.0, 3.0, 8.0, 9.0, 5.0, 2.0, 3.0, 7.0, 5.0, 9.0, 3.0, 7.0, 10.0, 3.0, 9.0, 4.0, 3.0, 6.0, 3.0, 7.0, 7.0, 2.0, 6.0, 10.0, 3.0, 10.0, 1.0, 10.0, 8.0, 8.0, 7.0, 7.0, 1.0, 6.0, 7.0, 8.0, 3.0, 8.0, 7.0, 7.0, 1.0, 7.0, 7.0, 2.0, 9.0, 3.0, 2.0, 6.0, 10.0, 5.0, 10.0, 10.0, 3.0, 6.0, 6.0, 5.0, 10.0, 2.0, 6.0, 6.0, 9.0, 9.0, 8.0, 7.0, 2.0, 3.0, 8.0, 1.0, 10.0, 8.0, 7.0, 3.0, 10.0, 4.0, 2.0, 2.0, 2.0, 4.0, 8.0, 4.0, 3.0, 4.0, 10.0, 8.0, 5.0, 1.0, 2.0, 6.0, 4.0, 3.0, 3.0, 5.0, 2.0, 8.0, 1.0, 5.0, 9.0, 5.0, 8.0, 4.0, 9.0, 8.0, 10.0, 8.0, 6.0, 4.0, 5.0, 2.0, 8.0, 5.0, 1.0, 9.0, 7.0, 4.0]
global b_y = 10
global p = [0.15, 0.368, 0.514, 0.228, 0.547, 0.904, 0.794, 0.61, 0.878, 0.464, 0.515, 0.325, 0.368, 0.564, 0.876, 0.857, 0.973, 0.825, 0.311, 0.333, 0.103, 0.622, 0.061, 0.848, 0.549, 0.314, 0.756, 0.179, 0.149, 0.895, 0.91, 0.653, 0.578, 0.624, 0.938, 0.277, 0.531, 0.72, 0.184, 0.963, 0.823, 0.38, 0.948, 0.283, 0.314, 0.579, 0.71, 0.469, 0.323, 0.084, 0.976, 0.085, 0.431, 0.138, 0.221, 0.256, 0.857, 0.279, 0.622, 0.197, 0.283, 0.342, 0.776, 0.331, 0.96, 0.567, 0.621, 0.464, 0.474, 0.212, 0.812, 0.673, 0.173, 0.93, 0.8, 0.54, 0.327, 0.849, 0.035, 0.79, 0.615, 0.571, 0.415, 0.326, 0.258, 0.167, 0.326, 0.937, 0.08, 0.212, 0.249, 0.433, 0.003, 0.113, 0.049, 0.827, 0.237, 0.898, 0.345, 0.298, 0.631, 0.094, 0.516, 0.24, 0.537, 0.289, 0.754, 0.2, 0.728, 0.129, 0.787, 0.082, 0.037, 0.754, 0.885, 0.124, 0.449, 0.954, 0.896, 0.187, 0.549, 0.64, 0.028, 0.224, 0.116, 0.339, 0.694, 0.541, 0.447, 0.013, 0.277, 0.732, 0.006, 0.107, 0.378, 0.573, 0.706, 0.112, 0.79, 0.841, 0.296, 0.765, 0.363, 0.644, 0.346, 0.823, 0.786, 0.95, 0.177, 0.24, 0.561, 0.995, 0.116, 0.67, 0.418, 0.972, 0.432, 0.232, 0.516, 0.348, 0.118, 0.779, 0.078, 0.874, 0.941, 0.892, 0.394, 0.257, 0.401, 0.138, 0.982, 0.556, 0.156, 0.922, 0.056, 0.614, 0.686, 0.486, 0.403, 0.981, 0.827, 0.763, 0.489, 0.254, 0.948, 0.392, 0.535, 0.654, 0.356, 0.721, 0.747, 0.979, 0.864, 0.611, 0.506, 0.708, 0.627, 0.202, 0.645, 0.556, 0.588, 0.741, 0.821, 0.793, 0.714, 0.258, 0.981, 0.674, 0.914, 0.485, 0.663, 0.221, 0.005, 0.984, 0.782, 0.562, 0.008, 0.617, 0.36, 0.627, 0.36, 0.076, 0.664, 0.036, 0.105, 0.891, 0.394, 0.894, 0.302, 0.999, 0.473]
global q = [0.317, 0.908, 0.757, 0.881, 0.818, 0.944, 0.999, 0.913, 0.992, 0.935, 0.853, 0.504, 0.494, 0.704, 0.944, 0.916, 0.984, 0.922, 0.329, 0.477, 0.508, 0.991, 0.284, 0.956, 0.72, 0.91, 0.782, 0.466, 0.504, 0.945, 0.944, 0.946, 0.73, 0.742, 0.971, 0.438, 0.87, 0.764, 0.586, 0.982, 0.976, 0.479, 0.957, 0.302, 0.487, 0.918, 0.877, 0.607, 0.461, 0.727, 0.985, 0.38, 0.595, 0.958, 0.882, 0.853, 0.99, 0.591, 0.735, 0.904, 0.354, 0.838, 0.817, 0.702, 0.973, 0.654, 0.716, 0.709, 0.595, 0.632, 0.932, 0.791, 0.28, 0.932, 0.907, 0.818, 0.531, 0.87, 0.851, 0.995, 0.743, 0.905, 0.936, 0.62, 0.258, 0.501, 0.862, 0.983, 0.657, 0.586, 0.392, 0.436, 0.413, 0.791, 0.912, 0.838, 0.299, 0.941, 0.36, 0.652, 0.707, 0.274, 0.64, 0.942, 0.884, 0.73, 0.855, 0.353, 0.814, 0.439, 0.873, 0.091, 0.276, 0.76, 0.885, 0.726, 0.928, 0.98, 0.899, 0.61, 0.588, 0.932, 0.594, 0.753, 0.644, 0.97, 0.7, 0.958, 0.507, 0.562, 0.883, 0.872, 0.363, 0.666, 0.663, 0.748, 0.976, 0.149, 0.931, 0.927, 0.396, 0.924, 0.704, 0.935, 0.65, 0.96, 0.829, 0.998, 0.954, 0.888, 0.827, 0.995, 0.347, 0.72, 0.912, 0.994, 0.93, 0.627, 0.679, 0.937, 0.196, 0.827, 0.686, 0.906, 0.946, 0.948, 0.922, 0.263, 0.445, 0.34, 0.995, 0.703, 0.209, 0.95, 0.743, 0.682, 0.884, 0.763, 0.529, 0.983, 0.856, 0.968, 0.884, 0.508, 0.997, 0.927, 0.793, 0.679, 0.775, 0.763, 0.751, 0.988, 0.953, 0.654, 0.645, 0.738, 0.667, 0.763, 0.815, 0.724, 0.776, 0.992, 0.918, 0.82, 0.836, 0.734, 0.981, 0.714, 0.986, 0.58, 0.716, 0.982, 0.255, 0.984, 0.842, 0.624, 0.602, 0.828, 0.552, 0.731, 0.579, 0.383, 0.772, 0.186, 0.318, 0.896, 0.455, 0.926, 0.694, 0.999, 0.829]
global origin = 1
global destination = 50