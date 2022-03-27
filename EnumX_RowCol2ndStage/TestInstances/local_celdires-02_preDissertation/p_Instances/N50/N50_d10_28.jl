global arcs = [1 14; 1 17; 1 27; 1 34; 1 36; 1 44; 1 50; 2 4; 2 6; 2 8; 2 13; 2 23; 2 44; 2 48; 3 10; 3 34; 3 35; 3 39; 3 44; 4 27; 4 44; 4 45; 4 49; 5 2; 5 19; 5 34; 5 36; 5 40; 6 17; 7 12; 7 18; 7 24; 7 28; 7 32; 7 35; 7 36; 7 50; 8 3; 8 37; 8 46; 9 2; 9 4; 9 5; 9 8; 9 20; 9 22; 9 26; 9 34; 9 41; 9 45; 9 49; 10 9; 10 18; 10 23; 10 25; 10 26; 10 29; 10 31; 10 33; 11 23; 11 28; 11 44; 12 11; 12 15; 12 21; 12 24; 12 49; 13 2; 13 14; 13 22; 13 47; 14 19; 14 26; 14 45; 15 12; 15 19; 15 31; 16 30; 16 35; 16 45; 16 46; 17 18; 17 33; 18 3; 18 11; 18 14; 18 22; 18 24; 19 18; 19 30; 19 35; 19 37; 19 42; 20 23; 20 24; 20 29; 20 32; 20 35; 20 36; 20 37; 20 38; 20 44; 21 5; 21 13; 21 15; 21 31; 21 36; 21 38; 21 42; 21 45; 22 2; 22 38; 22 50; 23 11; 23 12; 23 18; 23 19; 23 24; 23 44; 24 18; 25 2; 25 3; 25 15; 25 20; 25 35; 25 39; 26 5; 26 9; 26 10; 26 14; 26 20; 26 25; 26 48; 27 14; 27 16; 27 25; 27 28; 27 41; 27 45; 27 50; 28 3; 28 6; 28 44; 28 45; 28 50; 29 21; 30 8; 30 10; 30 12; 30 28; 30 50; 31 17; 31 18; 31 24; 31 27; 31 35; 31 37; 31 42; 32 20; 32 23; 32 29; 32 37; 32 39; 33 5; 33 9; 33 11; 33 15; 33 22; 33 32; 33 42; 33 50; 34 11; 34 12; 34 43; 34 47; 35 13; 35 42; 36 9; 36 16; 36 20; 36 39; 36 41; 37 5; 37 21; 37 28; 37 47; 38 16; 38 23; 38 24; 38 36; 38 39; 38 48; 39 17; 39 24; 39 37; 40 2; 40 10; 40 21; 40 28; 40 35; 40 36; 40 46; 41 7; 41 13; 41 16; 41 34; 41 36; 41 42; 41 44; 42 2; 42 14; 42 18; 42 27; 42 29; 42 33; 42 37; 42 38; 42 46; 43 30; 43 39; 44 18; 44 23; 44 29; 44 35; 44 39; 44 47; 45 22; 45 31; 46 3; 46 4; 46 7; 46 15; 46 29; 46 33; 47 5; 47 11; 47 25; 47 39; 47 40; 48 33; 48 41; 48 43; 48 45; 49 18; 49 26; 49 28; 49 41]
global d_x = [10.0, 3.0, 10.0, 6.0, 5.0, 7.0, 1.0, 3.0, 10.0, 4.0, 8.0, 2.0, 1.0, 10.0, 9.0, 5.0, 9.0, 10.0, 7.0, 1.0, 3.0, 8.0, 8.0, 10.0, 10.0, 5.0, 2.0, 8.0, 7.0, 2.0, 6.0, 10.0, 9.0, 1.0, 10.0, 2.0, 8.0, 2.0, 8.0, 2.0, 1.0, 1.0, 10.0, 1.0, 2.0, 9.0, 6.0, 2.0, 7.0, 8.0, 9.0, 7.0, 1.0, 9.0, 10.0, 10.0, 7.0, 10.0, 3.0, 4.0, 9.0, 3.0, 1.0, 7.0, 2.0, 9.0, 10.0, 1.0, 6.0, 9.0, 1.0, 10.0, 10.0, 4.0, 8.0, 7.0, 5.0, 6.0, 5.0, 10.0, 1.0, 4.0, 1.0, 2.0, 1.0, 7.0, 1.0, 2.0, 3.0, 6.0, 6.0, 10.0, 6.0, 10.0, 7.0, 2.0, 5.0, 10.0, 10.0, 4.0, 7.0, 3.0, 4.0, 9.0, 9.0, 3.0, 2.0, 4.0, 6.0, 4.0, 3.0, 7.0, 4.0, 1.0, 3.0, 10.0, 6.0, 5.0, 5.0, 9.0, 10.0, 9.0, 4.0, 5.0, 5.0, 10.0, 4.0, 1.0, 1.0, 4.0, 8.0, 3.0, 7.0, 2.0, 10.0, 6.0, 3.0, 4.0, 3.0, 1.0, 4.0, 7.0, 2.0, 5.0, 1.0, 4.0, 5.0, 2.0, 7.0, 6.0, 10.0, 4.0, 8.0, 1.0, 8.0, 6.0, 5.0, 3.0, 7.0, 2.0, 8.0, 9.0, 7.0, 4.0, 3.0, 3.0, 8.0, 4.0, 8.0, 6.0, 8.0, 8.0, 6.0, 7.0, 3.0, 3.0, 5.0, 5.0, 2.0, 7.0, 10.0, 3.0, 3.0, 2.0, 10.0, 4.0, 2.0, 1.0, 6.0, 7.0, 10.0, 9.0, 8.0, 4.0, 6.0, 3.0, 6.0, 6.0, 1.0, 4.0, 4.0, 5.0, 1.0, 3.0, 4.0, 9.0, 2.0, 9.0, 2.0, 9.0, 1.0, 10.0, 1.0, 10.0, 2.0, 3.0, 1.0, 3.0, 2.0, 4.0, 10.0, 7.0, 6.0, 4.0, 7.0, 10.0, 1.0, 1.0, 9.0, 2.0, 1.0, 4.0, 6.0, 4.0, 3.0, 3.0, 2.0, 1.0, 2.0, 9.0, 6.0, 2.0, 4.0, 5.0, 1.0, 10.0, 1.0]
global b_x = 5
global d_y = [5.0, 5.0, 4.0, 10.0, 4.0, 8.0, 7.0, 7.0, 6.0, 5.0, 6.0, 6.0, 7.0, 6.0, 9.0, 4.0, 5.0, 4.0, 3.0, 1.0, 6.0, 9.0, 3.0, 7.0, 10.0, 6.0, 10.0, 3.0, 2.0, 9.0, 6.0, 6.0, 7.0, 8.0, 3.0, 7.0, 5.0, 10.0, 9.0, 8.0, 8.0, 4.0, 2.0, 5.0, 9.0, 1.0, 3.0, 3.0, 6.0, 10.0, 4.0, 10.0, 7.0, 7.0, 2.0, 1.0, 5.0, 2.0, 4.0, 1.0, 1.0, 3.0, 9.0, 4.0, 3.0, 5.0, 8.0, 10.0, 7.0, 9.0, 6.0, 7.0, 4.0, 7.0, 4.0, 6.0, 7.0, 5.0, 9.0, 4.0, 6.0, 4.0, 2.0, 6.0, 8.0, 5.0, 10.0, 10.0, 4.0, 1.0, 5.0, 3.0, 10.0, 1.0, 5.0, 7.0, 3.0, 4.0, 3.0, 8.0, 6.0, 9.0, 10.0, 2.0, 3.0, 1.0, 10.0, 1.0, 8.0, 6.0, 6.0, 3.0, 8.0, 7.0, 7.0, 8.0, 7.0, 8.0, 9.0, 8.0, 10.0, 4.0, 2.0, 10.0, 7.0, 3.0, 6.0, 1.0, 7.0, 10.0, 7.0, 3.0, 2.0, 1.0, 1.0, 7.0, 9.0, 4.0, 10.0, 8.0, 1.0, 10.0, 6.0, 5.0, 2.0, 10.0, 8.0, 5.0, 7.0, 5.0, 1.0, 2.0, 5.0, 4.0, 3.0, 10.0, 1.0, 8.0, 4.0, 4.0, 7.0, 9.0, 9.0, 1.0, 6.0, 8.0, 5.0, 5.0, 6.0, 1.0, 9.0, 10.0, 2.0, 7.0, 4.0, 10.0, 4.0, 1.0, 7.0, 2.0, 1.0, 8.0, 7.0, 3.0, 7.0, 1.0, 2.0, 5.0, 4.0, 8.0, 10.0, 5.0, 5.0, 7.0, 7.0, 6.0, 10.0, 8.0, 5.0, 8.0, 8.0, 5.0, 1.0, 8.0, 9.0, 2.0, 3.0, 5.0, 9.0, 6.0, 10.0, 6.0, 3.0, 2.0, 7.0, 10.0, 5.0, 1.0, 5.0, 2.0, 3.0, 10.0, 9.0, 7.0, 9.0, 6.0, 6.0, 8.0, 5.0, 7.0, 4.0, 1.0, 6.0, 4.0, 7.0, 5.0, 6.0, 2.0, 1.0, 1.0, 7.0, 6.0, 9.0, 8.0, 8.0, 9.0, 7.0]
global b_y = 10
global p = [0.903, 0.721, 0.158, 0.719, 0.597, 0.394, 0.986, 0.619, 0.061, 0.205, 0.284, 0.631, 0.232, 0.397, 0.301, 0.399, 0.096, 0.707, 0.768, 0.808, 0.053, 0.812, 0.719, 0.822, 0.035, 0.157, 0.851, 0.416, 0.57, 0.757, 0.403, 0.92, 0.514, 0.739, 0.057, 0.457, 0.957, 0.449, 0.236, 0.597, 0.322, 0.296, 0.885, 0.017, 0.073, 0.59, 0.567, 0.863, 0.804, 0.206, 0.616, 0.616, 0.381, 0.522, 0.74, 0.28, 0.406, 0.022, 0.29, 0.071, 0.704, 0.538, 0.388, 0.828, 0.601, 0.852, 0.777, 0.848, 0.337, 0.067, 0.683, 0.655, 0.762, 0.545, 0.719, 0.708, 0.794, 0.298, 0.858, 0.258, 0.281, 0.512, 0.4, 0.184, 0.056, 0.583, 0.808, 0.311, 0.537, 0.574, 0.616, 0.279, 0.347, 0.831, 0.172, 0.955, 0.697, 0.556, 0.199, 0.254, 0.141, 0.978, 0.649, 0.638, 0.307, 0.048, 0.678, 0.953, 0.639, 0.497, 0.141, 0.699, 0.624, 0.364, 0.489, 0.561, 0.018, 0.311, 0.359, 0.539, 0.108, 0.563, 0.822, 0.661, 0.723, 0.404, 0.775, 0.695, 0.38, 0.152, 0.357, 0.448, 0.518, 0.666, 0.931, 0.846, 0.486, 0.901, 0.125, 0.263, 0.232, 0.549, 0.982, 0.838, 0.65, 0.036, 0.515, 0.19, 0.813, 0.683, 0.566, 0.535, 0.045, 0.355, 0.376, 0.3, 0.058, 0.883, 0.799, 0.501, 0.866, 0.081, 0.198, 0.274, 0.56, 0.138, 0.663, 0.253, 0.901, 0.149, 0.567, 0.26, 0.203, 0.602, 0.28, 0.707, 0.915, 0.739, 0.342, 0.836, 0.769, 0.8, 0.189, 0.879, 0.01, 0.993, 0.353, 0.686, 0.708, 0.574, 0.22, 0.872, 0.031, 0.135, 0.169, 0.977, 0.447, 0.234, 0.122, 0.52, 0.143, 0.464, 0.52, 0.339, 0.612, 0.862, 0.403, 0.275, 0.864, 0.086, 0.272, 0.487, 0.786, 0.275, 0.907, 0.306, 0.46, 0.697, 0.454, 0.493, 0.406, 0.247, 0.894, 0.941, 0.34, 0.645, 0.93, 0.696, 0.364, 0.017, 0.07, 0.316, 0.54, 0.347, 0.201, 0.001, 0.965, 0.935, 0.345, 0.055, 0.976, 0.132, 0.776, 0.836, 0.91, 0.854, 0.685]
global q = [0.937, 0.768, 0.761, 0.759, 0.652, 0.494, 0.999, 0.651, 0.507, 0.51, 0.49, 0.642, 0.545, 0.762, 0.714, 0.561, 0.905, 0.851, 0.96, 0.825, 0.95, 0.835, 0.998, 0.927, 0.394, 0.456, 0.967, 0.495, 0.876, 0.841, 0.975, 0.929, 0.636, 0.882, 0.38, 0.998, 0.961, 0.809, 0.468, 0.954, 0.63, 0.999, 0.96, 0.392, 0.693, 0.86, 0.636, 0.935, 0.972, 0.584, 0.721, 0.819, 0.436, 0.531, 0.972, 0.726, 0.904, 0.61, 0.664, 0.834, 0.808, 0.936, 0.484, 0.874, 0.948, 0.897, 0.951, 0.997, 0.872, 0.913, 0.992, 0.932, 0.854, 0.91, 0.939, 0.787, 0.902, 0.584, 0.978, 0.502, 0.54, 0.916, 0.466, 0.209, 0.659, 0.979, 0.86, 0.399, 0.763, 0.974, 0.751, 0.358, 0.64, 0.917, 0.928, 0.99, 0.934, 0.984, 0.289, 0.379, 0.408, 0.999, 0.748, 0.931, 0.549, 0.379, 0.978, 0.985, 0.913, 0.935, 0.221, 0.885, 0.874, 0.891, 0.581, 0.868, 0.566, 0.431, 0.641, 0.656, 0.147, 0.912, 0.889, 0.759, 0.78, 0.83, 0.854, 0.969, 0.485, 0.562, 0.967, 0.786, 0.658, 0.816, 0.968, 0.854, 0.631, 0.997, 0.206, 0.764, 0.47, 0.715, 0.997, 0.969, 0.871, 0.288, 0.898, 0.679, 0.928, 0.728, 0.795, 0.704, 0.673, 0.692, 0.83, 0.561, 0.274, 0.884, 0.879, 0.569, 0.958, 0.89, 0.764, 0.675, 0.758, 0.286, 0.666, 0.297, 0.94, 0.891, 0.613, 0.843, 0.957, 0.997, 0.753, 0.898, 0.96, 0.957, 0.652, 0.916, 0.921, 0.809, 0.724, 0.982, 0.45, 0.998, 0.378, 0.792, 0.728, 0.819, 0.675, 0.975, 0.559, 0.46, 0.808, 0.995, 0.911, 0.805, 0.356, 0.585, 0.598, 0.846, 0.619, 0.407, 0.887, 0.902, 0.856, 0.765, 0.865, 0.13, 0.509, 0.942, 0.891, 0.414, 0.98, 0.629, 0.737, 0.985, 0.795, 0.51, 0.801, 0.698, 0.955, 0.983, 0.88, 0.835, 0.935, 0.698, 0.736, 0.588, 0.183, 0.93, 0.849, 0.571, 0.49, 0.869, 0.989, 0.952, 0.972, 0.831, 0.989, 0.442, 0.785, 0.927, 0.932, 0.883, 0.957]
global origin = 1
global destination = 50