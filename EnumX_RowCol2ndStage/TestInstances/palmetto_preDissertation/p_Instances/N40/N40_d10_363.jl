global arcs = [1 2; 1 6; 1 14; 1 17; 1 22; 1 27; 1 34; 2 6; 2 7; 2 12; 2 14; 2 20; 2 33; 2 35; 3 4; 3 24; 3 27; 3 28; 4 2; 4 3; 4 11; 4 13; 4 23; 4 35; 4 36; 4 38; 5 22; 5 34; 6 13; 6 14; 6 35; 7 2; 7 8; 7 15; 7 21; 7 23; 7 31; 7 32; 7 34; 7 35; 8 7; 8 35; 9 4; 9 8; 9 12; 9 16; 9 30; 10 25; 10 29; 10 38; 11 3; 11 27; 11 37; 11 38; 12 8; 12 24; 12 25; 12 27; 12 34; 12 35; 13 14; 14 12; 14 30; 14 35; 14 37; 14 38; 15 18; 15 26; 15 37; 15 40; 16 26; 16 38; 17 6; 17 40; 18 4; 18 15; 18 19; 18 20; 18 24; 18 38; 19 6; 19 15; 19 24; 19 30; 19 31; 19 34; 20 2; 20 8; 20 11; 20 27; 21 6; 21 8; 21 11; 21 15; 21 17; 21 26; 22 6; 22 8; 23 2; 23 5; 23 15; 23 19; 23 32; 23 38; 24 2; 24 4; 24 11; 24 26; 25 18; 25 30; 26 2; 26 10; 26 19; 26 32; 27 8; 27 25; 27 32; 28 4; 28 22; 29 2; 29 13; 29 26; 29 31; 29 37; 29 40; 30 17; 30 19; 30 20; 30 27; 30 29; 31 17; 31 34; 31 35; 31 39; 31 40; 32 7; 33 6; 33 8; 33 27; 34 30; 34 37; 35 13; 35 15; 35 20; 35 30; 35 32; 35 36; 35 37; 35 38; 35 40; 36 8; 36 9; 37 3; 37 5; 37 19; 37 24; 37 29; 38 8; 38 13; 38 28; 39 8; 39 12; 39 33; 39 35]
global d_x = [6.0, 7.0, 9.0, 7.0, 6.0, 3.0, 6.0, 2.0, 10.0, 5.0, 9.0, 9.0, 7.0, 1.0, 4.0, 3.0, 1.0, 2.0, 5.0, 9.0, 9.0, 7.0, 8.0, 7.0, 5.0, 10.0, 6.0, 4.0, 2.0, 10.0, 4.0, 5.0, 8.0, 1.0, 7.0, 6.0, 3.0, 7.0, 9.0, 4.0, 10.0, 3.0, 5.0, 8.0, 2.0, 4.0, 4.0, 3.0, 7.0, 4.0, 1.0, 10.0, 4.0, 3.0, 2.0, 5.0, 8.0, 9.0, 1.0, 8.0, 4.0, 1.0, 7.0, 10.0, 2.0, 2.0, 7.0, 8.0, 6.0, 5.0, 3.0, 7.0, 4.0, 2.0, 8.0, 4.0, 7.0, 3.0, 9.0, 8.0, 1.0, 7.0, 8.0, 7.0, 10.0, 8.0, 9.0, 4.0, 8.0, 1.0, 6.0, 7.0, 8.0, 3.0, 2.0, 9.0, 1.0, 4.0, 10.0, 9.0, 7.0, 4.0, 10.0, 9.0, 3.0, 7.0, 8.0, 9.0, 3.0, 9.0, 6.0, 5.0, 2.0, 1.0, 9.0, 2.0, 10.0, 8.0, 7.0, 7.0, 9.0, 10.0, 7.0, 6.0, 3.0, 10.0, 2.0, 1.0, 3.0, 1.0, 3.0, 8.0, 5.0, 2.0, 7.0, 4.0, 2.0, 4.0, 4.0, 3.0, 3.0, 6.0, 8.0, 3.0, 1.0, 5.0, 3.0, 10.0, 3.0, 1.0, 4.0, 7.0, 5.0, 4.0, 3.0, 6.0, 7.0, 6.0, 1.0, 10.0, 7.0, 1.0, 7.0, 3.0]
global b_x = 5
global d_y = [10.0, 6.0, 9.0, 2.0, 5.0, 8.0, 4.0, 9.0, 4.0, 10.0, 1.0, 5.0, 7.0, 8.0, 5.0, 6.0, 4.0, 7.0, 4.0, 3.0, 8.0, 8.0, 6.0, 4.0, 6.0, 8.0, 3.0, 3.0, 7.0, 8.0, 8.0, 10.0, 9.0, 9.0, 4.0, 6.0, 4.0, 10.0, 3.0, 2.0, 9.0, 4.0, 4.0, 10.0, 9.0, 7.0, 8.0, 4.0, 4.0, 1.0, 10.0, 6.0, 2.0, 9.0, 6.0, 1.0, 5.0, 3.0, 10.0, 8.0, 1.0, 2.0, 3.0, 4.0, 3.0, 6.0, 3.0, 9.0, 2.0, 6.0, 5.0, 10.0, 9.0, 10.0, 5.0, 7.0, 6.0, 3.0, 6.0, 7.0, 2.0, 10.0, 3.0, 3.0, 5.0, 2.0, 2.0, 5.0, 2.0, 2.0, 10.0, 9.0, 4.0, 10.0, 8.0, 2.0, 5.0, 5.0, 7.0, 4.0, 4.0, 5.0, 9.0, 7.0, 6.0, 5.0, 4.0, 2.0, 9.0, 1.0, 2.0, 5.0, 6.0, 3.0, 4.0, 8.0, 3.0, 3.0, 5.0, 7.0, 4.0, 7.0, 2.0, 1.0, 1.0, 5.0, 1.0, 9.0, 5.0, 5.0, 10.0, 7.0, 3.0, 1.0, 2.0, 2.0, 9.0, 8.0, 2.0, 2.0, 5.0, 8.0, 3.0, 7.0, 5.0, 10.0, 6.0, 2.0, 10.0, 6.0, 10.0, 7.0, 8.0, 5.0, 9.0, 1.0, 3.0, 3.0, 4.0, 1.0, 3.0, 3.0, 6.0, 2.0]
global b_y = 10
global p = [0.926, 0.093, 0.167, 0.453, 0.854, 0.564, 0.026, 0.273, 0.264, 0.28, 0.257, 0.073, 0.546, 0.927, 0.825, 0.111, 0.772, 0.408, 0.204, 0.715, 0.346, 0.127, 0.054, 0.033, 0.741, 0.413, 0.71, 0.07, 0.121, 0.154, 0.709, 0.056, 0.399, 0.384, 0.91, 0.796, 0.651, 0.045, 0.327, 0.795, 0.636, 0.388, 0.98, 0.444, 0.24, 0.268, 0.196, 0.829, 0.178, 0.989, 0.201, 0.882, 0.358, 0.276, 0.601, 0.553, 0.192, 0.066, 0.73, 0.923, 0.82, 0.227, 0.089, 0.883, 0.524, 0.595, 0.092, 0.295, 0.043, 0.917, 0.184, 0.575, 0.161, 0.645, 0.253, 0.318, 0.852, 0.31, 0.754, 0.178, 0.498, 0.204, 0.916, 0.441, 0.698, 0.34, 0.245, 0.331, 0.957, 0.03, 0.999, 0.957, 0.59, 0.638, 0.036, 0.854, 0.914, 0.346, 0.522, 0.122, 0.541, 0.403, 0.663, 0.803, 0.682, 0.709, 0.823, 0.742, 0.542, 0.242, 0.221, 0.158, 0.878, 0.421, 0.631, 0.827, 0.9, 0.352, 0.374, 0.02, 0.429, 0.736, 0.497, 0.741, 0.36, 0.799, 0.739, 0.525, 0.545, 0.011, 0.788, 0.233, 0.232, 0.757, 0.84, 0.58, 0.516, 0.477, 0.946, 0.15, 0.313, 0.838, 0.566, 0.732, 0.943, 0.028, 0.339, 0.632, 0.404, 0.384, 0.917, 0.579, 0.566, 0.11, 0.9, 0.124, 0.869, 0.497, 0.388, 0.105, 0.92, 0.23, 0.12, 0.443]
global q = [0.985, 0.857, 0.698, 0.587, 0.871, 0.589, 0.475, 0.483, 0.848, 0.467, 0.995, 0.263, 0.622, 0.963, 0.876, 0.321, 0.825, 0.994, 0.617, 0.791, 0.604, 0.316, 0.922, 0.653, 0.99, 0.938, 0.935, 0.435, 0.65, 0.488, 0.904, 0.775, 0.898, 0.738, 0.995, 0.847, 0.753, 0.126, 0.776, 0.979, 0.945, 0.581, 0.981, 0.511, 0.48, 0.35, 0.713, 0.944, 0.782, 0.992, 0.639, 0.96, 0.764, 0.31, 0.897, 0.693, 0.818, 0.656, 0.988, 0.994, 0.882, 0.74, 0.632, 0.961, 0.845, 0.73, 0.867, 0.468, 0.979, 0.958, 0.863, 0.828, 0.957, 0.772, 0.898, 0.91, 0.988, 0.725, 0.853, 0.286, 0.999, 0.688, 0.93, 0.77, 0.959, 0.634, 0.976, 0.386, 0.991, 0.959, 0.999, 0.985, 0.628, 0.666, 0.288, 0.943, 0.928, 0.466, 0.892, 0.38, 0.78, 0.42, 0.89, 0.819, 0.844, 0.727, 0.986, 0.813, 0.699, 0.6, 0.431, 0.667, 0.892, 0.693, 0.722, 0.986, 0.932, 0.976, 0.61, 0.438, 0.555, 0.737, 0.885, 0.778, 0.822, 0.958, 0.786, 0.944, 0.619, 0.041, 0.847, 0.788, 0.676, 0.891, 0.957, 0.635, 0.607, 0.724, 0.998, 0.686, 0.828, 0.977, 0.798, 0.81, 0.966, 0.408, 0.452, 0.812, 0.672, 0.957, 0.96, 0.693, 0.97, 0.228, 0.995, 0.382, 0.94, 0.917, 0.463, 0.741, 0.981, 0.324, 0.981, 0.558]
global origin = 1
global destination = 40