global arcs = [1 2; 1 17; 1 33; 1 37; 2 31; 2 35; 2 40; 3 23; 3 31; 3 39; 4 5; 4 11; 4 14; 4 38; 5 2; 5 13; 5 15; 5 22; 5 36; 5 37; 5 39; 6 2; 6 5; 6 19; 6 25; 6 28; 6 40; 7 9; 7 10; 7 11; 8 18; 8 31; 8 36; 9 12; 9 14; 9 18; 9 39; 10 6; 10 15; 10 34; 10 37; 11 2; 11 8; 11 15; 11 33; 12 3; 12 23; 12 26; 12 31; 12 33; 12 34; 13 14; 13 17; 14 5; 14 9; 14 30; 15 3; 15 22; 15 27; 15 35; 15 39; 16 3; 16 21; 17 38; 18 5; 18 19; 18 20; 18 27; 18 31; 18 32; 18 36; 19 8; 19 17; 19 24; 19 33; 20 8; 20 9; 20 11; 20 21; 20 27; 21 20; 22 4; 22 11; 22 14; 22 15; 22 40; 23 4; 23 11; 23 20; 23 25; 23 39; 24 13; 24 21; 24 26; 24 39; 25 14; 25 19; 25 27; 26 31; 26 32; 26 33; 26 37; 27 40; 28 6; 28 16; 28 21; 28 35; 28 40; 29 11; 29 23; 29 24; 29 33; 29 36; 30 11; 30 16; 30 26; 30 29; 30 32; 30 40; 31 2; 31 19; 31 27; 31 29; 32 9; 32 19; 32 33; 32 34; 32 36; 32 37; 32 39; 33 5; 33 7; 33 13; 33 26; 33 32; 33 38; 33 39; 34 2; 34 9; 34 13; 34 15; 34 23; 35 29; 35 32; 36 12; 36 25; 36 29; 36 34; 36 38; 37 3; 37 5; 37 8; 37 9; 37 13; 37 29; 37 34; 37 35; 37 36; 38 8; 38 11; 38 19; 38 20; 38 25; 38 31; 39 4; 39 13; 39 17; 39 18; 39 29; 39 31; 39 32; 39 35]
global d_x = [3.0, 6.0, 1.0, 9.0, 2.0, 3.0, 10.0, 10.0, 8.0, 7.0, 8.0, 8.0, 5.0, 3.0, 8.0, 9.0, 10.0, 9.0, 9.0, 3.0, 1.0, 6.0, 2.0, 6.0, 6.0, 7.0, 8.0, 10.0, 2.0, 3.0, 5.0, 9.0, 9.0, 4.0, 9.0, 4.0, 1.0, 8.0, 2.0, 8.0, 6.0, 2.0, 9.0, 10.0, 6.0, 7.0, 9.0, 1.0, 1.0, 3.0, 3.0, 6.0, 4.0, 1.0, 8.0, 5.0, 8.0, 5.0, 3.0, 1.0, 10.0, 9.0, 5.0, 5.0, 3.0, 9.0, 10.0, 1.0, 1.0, 9.0, 9.0, 10.0, 9.0, 1.0, 5.0, 6.0, 3.0, 1.0, 4.0, 5.0, 4.0, 1.0, 8.0, 8.0, 9.0, 8.0, 9.0, 8.0, 7.0, 4.0, 4.0, 3.0, 1.0, 6.0, 10.0, 7.0, 4.0, 2.0, 1.0, 7.0, 5.0, 7.0, 8.0, 4.0, 2.0, 4.0, 6.0, 7.0, 1.0, 6.0, 8.0, 6.0, 8.0, 2.0, 10.0, 9.0, 7.0, 1.0, 1.0, 7.0, 7.0, 9.0, 6.0, 4.0, 7.0, 6.0, 3.0, 2.0, 10.0, 2.0, 5.0, 6.0, 4.0, 8.0, 8.0, 7.0, 2.0, 10.0, 2.0, 1.0, 7.0, 9.0, 8.0, 4.0, 6.0, 9.0, 1.0, 4.0, 7.0, 2.0, 7.0, 4.0, 5.0, 4.0, 4.0, 4.0, 8.0, 2.0, 1.0, 5.0, 3.0, 1.0, 1.0, 5.0, 9.0, 9.0, 7.0, 9.0, 3.0, 6.0, 6.0, 6.0]
global b_x = 5
global d_y = [6.0, 7.0, 5.0, 9.0, 10.0, 10.0, 1.0, 10.0, 1.0, 2.0, 8.0, 4.0, 4.0, 3.0, 7.0, 9.0, 10.0, 7.0, 6.0, 4.0, 4.0, 2.0, 6.0, 3.0, 5.0, 10.0, 9.0, 10.0, 10.0, 8.0, 2.0, 10.0, 7.0, 5.0, 9.0, 4.0, 8.0, 3.0, 8.0, 3.0, 9.0, 9.0, 5.0, 1.0, 9.0, 7.0, 2.0, 5.0, 10.0, 5.0, 8.0, 10.0, 10.0, 2.0, 5.0, 3.0, 9.0, 5.0, 6.0, 9.0, 7.0, 1.0, 10.0, 1.0, 5.0, 8.0, 6.0, 8.0, 6.0, 4.0, 6.0, 1.0, 7.0, 10.0, 3.0, 3.0, 7.0, 6.0, 3.0, 4.0, 1.0, 2.0, 7.0, 5.0, 2.0, 2.0, 1.0, 2.0, 9.0, 4.0, 2.0, 10.0, 6.0, 7.0, 4.0, 1.0, 1.0, 4.0, 9.0, 9.0, 6.0, 3.0, 10.0, 2.0, 10.0, 6.0, 2.0, 7.0, 2.0, 2.0, 5.0, 6.0, 1.0, 8.0, 4.0, 9.0, 9.0, 7.0, 8.0, 2.0, 5.0, 7.0, 3.0, 4.0, 8.0, 1.0, 10.0, 5.0, 9.0, 1.0, 3.0, 3.0, 1.0, 4.0, 8.0, 8.0, 2.0, 1.0, 4.0, 7.0, 8.0, 6.0, 2.0, 8.0, 5.0, 10.0, 10.0, 2.0, 9.0, 7.0, 3.0, 9.0, 2.0, 10.0, 1.0, 7.0, 10.0, 4.0, 5.0, 3.0, 6.0, 6.0, 5.0, 4.0, 2.0, 6.0, 5.0, 10.0, 8.0, 1.0, 1.0, 2.0]
global b_y = 10
global p = [0.339, 0.623, 0.703, 0.466, 0.238, 0.622, 0.545, 0.429, 0.483, 0.942, 0.909, 0.648, 0.942, 0.26, 0.669, 0.672, 0.358, 0.856, 0.822, 0.738, 0.943, 0.296, 0.729, 0.926, 0.093, 0.638, 0.494, 0.323, 0.611, 0.034, 0.353, 0.467, 0.714, 0.71, 0.272, 0.097, 0.786, 0.429, 0.137, 0.101, 0.114, 0.27, 0.794, 0.027, 0.181, 0.107, 0.2, 0.813, 0.53, 0.809, 0.216, 0.894, 0.763, 0.108, 0.716, 0.095, 0.981, 0.177, 0.171, 0.917, 0.915, 0.133, 0.273, 0.239, 0.766, 0.538, 0.259, 0.128, 0.16, 0.357, 0.375, 0.35, 0.304, 0.56, 0.884, 0.762, 0.296, 0.031, 0.634, 0.116, 0.764, 0.449, 0.314, 0.326, 0.341, 0.993, 0.969, 0.258, 0.255, 0.695, 0.767, 0.144, 0.556, 0.515, 0.995, 0.425, 0.346, 0.705, 0.731, 0.894, 0.708, 0.232, 0.14, 0.636, 0.572, 0.226, 0.852, 0.996, 0.505, 0.979, 0.039, 0.277, 0.38, 0.881, 0.101, 0.009, 0.897, 0.239, 0.038, 0.248, 0.305, 0.747, 0.044, 0.027, 0.705, 0.597, 0.053, 0.818, 0.409, 0.277, 0.063, 0.655, 0.673, 0.323, 0.093, 0.799, 0.52, 0.451, 0.46, 0.734, 0.802, 0.718, 0.586, 0.435, 0.651, 0.414, 0.487, 0.321, 0.938, 0.165, 0.467, 0.905, 0.063, 0.983, 0.34, 0.894, 0.591, 0.27, 0.308, 0.858, 0.507, 0.48, 0.818, 0.21, 0.835, 0.996, 0.673, 0.381, 0.437, 0.442, 0.285, 0.992]
global q = [0.404, 0.79, 0.953, 0.891, 0.686, 0.751, 0.894, 0.778, 0.544, 0.979, 0.955, 0.977, 0.942, 0.54, 0.764, 0.828, 0.86, 0.915, 0.832, 0.739, 0.954, 0.876, 0.807, 0.963, 0.848, 0.872, 0.939, 0.329, 0.729, 0.183, 0.611, 0.467, 0.885, 0.915, 0.488, 0.741, 0.904, 0.753, 0.875, 0.488, 0.435, 0.293, 0.845, 0.227, 0.887, 0.94, 0.377, 0.841, 0.71, 0.969, 0.317, 0.938, 0.819, 0.398, 0.98, 0.374, 0.996, 0.776, 0.656, 0.994, 0.942, 0.407, 0.738, 0.95, 0.777, 0.993, 0.885, 0.518, 0.833, 0.726, 0.857, 0.363, 0.403, 0.72, 0.909, 0.947, 0.344, 0.059, 0.681, 0.843, 0.766, 0.819, 0.79, 0.814, 0.883, 0.995, 0.988, 0.758, 0.743, 0.948, 0.807, 0.385, 0.958, 0.869, 0.995, 0.472, 0.406, 0.903, 0.971, 0.918, 0.758, 0.43, 0.868, 0.966, 0.767, 0.467, 0.968, 0.997, 0.791, 0.981, 0.611, 0.59, 0.614, 0.987, 0.792, 0.055, 0.962, 0.812, 0.153, 0.382, 0.636, 0.752, 0.773, 0.068, 0.949, 0.606, 0.642, 0.827, 0.857, 0.368, 0.314, 0.725, 0.76, 0.456, 0.337, 0.921, 0.871, 0.531, 0.769, 0.848, 0.923, 0.953, 0.94, 0.948, 0.899, 0.96, 0.968, 0.805, 0.96, 0.302, 0.716, 0.93, 0.423, 0.984, 0.746, 0.987, 0.745, 0.488, 0.863, 0.926, 0.53, 0.556, 0.857, 0.761, 0.917, 0.999, 0.956, 0.476, 0.671, 0.61, 0.635, 0.993]
global origin = 1
global destination = 40