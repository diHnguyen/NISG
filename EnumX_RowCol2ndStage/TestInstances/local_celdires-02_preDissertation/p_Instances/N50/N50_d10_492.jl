global arcs = [1 4; 1 11; 1 30; 1 47; 2 15; 2 26; 2 33; 2 35; 2 39; 2 45; 3 8; 3 22; 3 42; 4 3; 4 6; 4 14; 4 18; 4 39; 4 43; 4 44; 5 11; 5 14; 5 15; 5 45; 6 5; 6 7; 6 23; 6 47; 7 4; 7 14; 7 35; 8 19; 8 21; 8 24; 8 28; 8 41; 9 4; 9 18; 9 27; 9 37; 9 38; 9 43; 9 44; 9 47; 10 11; 10 21; 10 47; 11 5; 11 14; 11 30; 11 42; 12 5; 12 13; 12 24; 12 25; 12 34; 12 39; 12 46; 13 4; 13 25; 13 28; 14 15; 14 19; 14 29; 14 47; 14 50; 15 3; 15 12; 15 34; 15 39; 15 42; 16 6; 16 13; 16 25; 16 29; 16 33; 16 35; 16 37; 16 38; 16 42; 16 48; 17 5; 17 10; 17 26; 17 28; 17 37; 17 41; 17 43; 18 9; 18 28; 18 43; 19 4; 19 5; 19 7; 19 18; 19 45; 20 12; 20 18; 20 36; 20 37; 20 50; 21 5; 21 6; 21 19; 21 40; 22 31; 22 32; 22 40; 22 50; 23 2; 23 26; 23 30; 23 48; 24 2; 24 9; 24 12; 24 22; 24 23; 24 26; 24 34; 24 35; 24 43; 24 48; 25 14; 25 22; 25 23; 25 28; 25 39; 25 41; 25 47; 26 8; 26 13; 26 17; 26 24; 27 4; 27 13; 27 16; 27 18; 27 19; 27 28; 27 34; 27 39; 28 7; 28 15; 28 24; 28 26; 28 29; 28 38; 29 5; 29 11; 29 25; 29 35; 29 43; 30 2; 30 4; 30 34; 30 44; 30 49; 31 4; 31 28; 31 36; 31 41; 31 43; 32 3; 32 16; 32 30; 32 38; 32 44; 33 15; 33 37; 34 33; 34 44; 35 5; 35 13; 35 20; 35 27; 35 37; 35 43; 35 44; 35 46; 35 48; 36 10; 36 11; 36 13; 36 14; 36 29; 36 33; 36 45; 37 2; 37 3; 37 8; 37 18; 37 29; 37 33; 38 29; 38 43; 38 47; 39 5; 39 6; 39 10; 39 11; 39 49; 40 6; 40 22; 40 34; 40 36; 40 37; 41 9; 41 21; 41 27; 41 38; 41 40; 42 9; 42 25; 42 26; 42 46; 43 3; 43 8; 43 11; 43 12; 43 25; 43 35; 43 50; 44 9; 44 25; 44 39; 45 6; 45 19; 45 44; 45 48; 46 4; 46 26; 46 49; 47 9; 47 13; 47 45; 47 48; 47 49; 48 5; 48 6; 48 12; 48 14; 48 38; 48 39; 48 42; 49 10; 49 17; 49 19; 49 20]
global d_x = [10.0, 10.0, 1.0, 3.0, 9.0, 6.0, 6.0, 3.0, 2.0, 8.0, 2.0, 1.0, 7.0, 8.0, 7.0, 5.0, 5.0, 10.0, 8.0, 5.0, 6.0, 10.0, 3.0, 6.0, 10.0, 9.0, 7.0, 5.0, 4.0, 4.0, 6.0, 5.0, 7.0, 7.0, 9.0, 5.0, 6.0, 2.0, 8.0, 5.0, 4.0, 6.0, 1.0, 7.0, 9.0, 9.0, 1.0, 1.0, 10.0, 3.0, 1.0, 1.0, 9.0, 8.0, 6.0, 9.0, 1.0, 9.0, 9.0, 7.0, 3.0, 6.0, 2.0, 1.0, 8.0, 10.0, 10.0, 1.0, 6.0, 3.0, 1.0, 9.0, 2.0, 10.0, 10.0, 5.0, 1.0, 2.0, 3.0, 10.0, 3.0, 10.0, 10.0, 1.0, 10.0, 1.0, 10.0, 8.0, 4.0, 2.0, 6.0, 1.0, 5.0, 3.0, 4.0, 6.0, 3.0, 2.0, 7.0, 2.0, 7.0, 1.0, 10.0, 2.0, 6.0, 3.0, 3.0, 1.0, 10.0, 8.0, 2.0, 10.0, 2.0, 8.0, 10.0, 2.0, 6.0, 8.0, 1.0, 2.0, 7.0, 5.0, 2.0, 3.0, 9.0, 1.0, 5.0, 2.0, 2.0, 7.0, 1.0, 10.0, 1.0, 3.0, 2.0, 10.0, 1.0, 1.0, 5.0, 3.0, 9.0, 8.0, 6.0, 10.0, 9.0, 3.0, 4.0, 10.0, 7.0, 5.0, 7.0, 6.0, 9.0, 9.0, 3.0, 8.0, 8.0, 10.0, 2.0, 10.0, 8.0, 2.0, 6.0, 10.0, 9.0, 5.0, 4.0, 10.0, 5.0, 7.0, 4.0, 4.0, 3.0, 9.0, 8.0, 9.0, 1.0, 1.0, 8.0, 2.0, 5.0, 9.0, 9.0, 7.0, 9.0, 4.0, 5.0, 6.0, 1.0, 9.0, 5.0, 9.0, 3.0, 9.0, 5.0, 6.0, 5.0, 9.0, 9.0, 1.0, 6.0, 7.0, 3.0, 9.0, 3.0, 4.0, 4.0, 4.0, 6.0, 5.0, 3.0, 10.0, 8.0, 1.0, 7.0, 1.0, 1.0, 8.0, 5.0, 3.0, 10.0, 2.0, 7.0, 4.0, 4.0, 7.0, 8.0, 10.0, 1.0, 8.0, 4.0, 9.0, 2.0, 7.0, 3.0, 3.0, 7.0, 9.0, 1.0, 7.0, 7.0, 9.0, 10.0, 8.0, 6.0, 7.0, 9.0, 4.0, 6.0]
global b_x = 5
global d_y = [1.0, 6.0, 6.0, 4.0, 5.0, 10.0, 1.0, 3.0, 7.0, 2.0, 10.0, 10.0, 1.0, 4.0, 5.0, 2.0, 7.0, 6.0, 3.0, 8.0, 2.0, 3.0, 10.0, 5.0, 8.0, 4.0, 1.0, 2.0, 4.0, 2.0, 9.0, 8.0, 7.0, 6.0, 1.0, 4.0, 6.0, 5.0, 6.0, 9.0, 3.0, 8.0, 10.0, 10.0, 9.0, 4.0, 7.0, 2.0, 10.0, 7.0, 2.0, 8.0, 4.0, 5.0, 9.0, 2.0, 3.0, 8.0, 6.0, 9.0, 9.0, 6.0, 8.0, 8.0, 5.0, 5.0, 4.0, 10.0, 2.0, 9.0, 1.0, 8.0, 9.0, 2.0, 8.0, 2.0, 8.0, 10.0, 10.0, 4.0, 6.0, 1.0, 4.0, 1.0, 2.0, 9.0, 10.0, 7.0, 2.0, 7.0, 5.0, 1.0, 1.0, 2.0, 9.0, 10.0, 4.0, 3.0, 8.0, 9.0, 1.0, 1.0, 6.0, 9.0, 8.0, 4.0, 2.0, 7.0, 10.0, 3.0, 4.0, 9.0, 4.0, 2.0, 5.0, 1.0, 7.0, 5.0, 1.0, 8.0, 8.0, 8.0, 8.0, 2.0, 8.0, 5.0, 10.0, 6.0, 1.0, 9.0, 7.0, 3.0, 3.0, 10.0, 1.0, 2.0, 4.0, 1.0, 9.0, 6.0, 1.0, 9.0, 5.0, 2.0, 7.0, 6.0, 1.0, 10.0, 7.0, 7.0, 5.0, 1.0, 7.0, 3.0, 5.0, 1.0, 5.0, 8.0, 6.0, 2.0, 2.0, 5.0, 6.0, 2.0, 6.0, 1.0, 2.0, 10.0, 3.0, 5.0, 1.0, 6.0, 2.0, 2.0, 9.0, 3.0, 7.0, 1.0, 6.0, 3.0, 4.0, 4.0, 7.0, 10.0, 5.0, 5.0, 2.0, 9.0, 7.0, 8.0, 1.0, 7.0, 2.0, 6.0, 2.0, 5.0, 8.0, 8.0, 9.0, 8.0, 9.0, 8.0, 3.0, 4.0, 5.0, 3.0, 6.0, 2.0, 1.0, 10.0, 7.0, 1.0, 5.0, 8.0, 2.0, 4.0, 4.0, 2.0, 10.0, 6.0, 10.0, 7.0, 5.0, 10.0, 10.0, 9.0, 9.0, 4.0, 1.0, 6.0, 3.0, 7.0, 8.0, 7.0, 2.0, 9.0, 7.0, 6.0, 8.0, 9.0, 2.0, 7.0, 9.0, 5.0, 5.0, 10.0, 1.0, 4.0, 1.0]
global b_y = 10
global p = [0.553, 0.361, 0.552, 0.177, 0.03, 0.894, 0.041, 0.569, 0.725, 0.237, 0.552, 0.575, 0.532, 0.631, 0.037, 0.887, 0.968, 0.954, 0.757, 0.714, 0.965, 0.176, 0.683, 0.63, 0.614, 0.844, 0.312, 0.237, 0.355, 0.388, 0.927, 0.378, 0.435, 0.66, 0.208, 0.498, 0.721, 0.041, 0.224, 0.507, 0.304, 0.674, 0.08, 0.225, 0.092, 0.583, 0.66, 0.297, 0.578, 0.963, 0.537, 0.061, 0.484, 0.433, 0.001, 0.582, 0.255, 0.687, 0.84, 0.829, 0.428, 0.661, 0.792, 0.699, 0.421, 0.08, 0.507, 0.359, 0.33, 0.171, 0.146, 0.445, 0.095, 0.339, 0.779, 0.645, 0.413, 0.259, 0.881, 0.722, 0.963, 0.621, 0.199, 0.858, 0.945, 0.23, 0.404, 0.939, 0.961, 0.264, 0.175, 0.545, 0.46, 0.89, 0.669, 0.76, 0.441, 0.377, 0.878, 0.17, 0.849, 0.644, 0.287, 0.473, 0.063, 0.497, 0.159, 0.545, 0.488, 0.277, 0.116, 0.588, 0.744, 0.464, 0.803, 0.373, 0.509, 0.863, 0.428, 0.395, 0.367, 0.696, 0.738, 0.852, 0.161, 0.262, 0.7, 0.716, 0.295, 0.196, 0.004, 0.13, 0.435, 0.199, 0.653, 0.679, 0.433, 0.73, 0.909, 0.672, 0.443, 0.509, 0.402, 0.189, 0.457, 0.064, 0.906, 0.632, 0.199, 0.524, 0.826, 0.095, 0.372, 0.371, 0.268, 0.839, 0.495, 0.658, 0.535, 0.25, 0.389, 0.187, 0.445, 0.967, 0.011, 0.087, 0.051, 0.285, 0.383, 0.773, 0.542, 0.073, 0.839, 0.431, 0.495, 0.006, 0.088, 0.193, 0.5, 0.986, 0.899, 0.043, 0.524, 0.628, 0.412, 0.648, 0.1, 0.658, 0.089, 0.435, 0.407, 0.609, 0.875, 0.127, 0.806, 0.818, 0.163, 0.884, 0.168, 0.334, 0.204, 0.638, 0.861, 0.908, 0.417, 0.303, 0.578, 0.875, 0.65, 0.237, 0.039, 0.328, 0.683, 0.234, 0.2, 0.181, 0.664, 0.43, 0.051, 0.146, 0.665, 0.654, 0.128, 0.412, 0.108, 0.479, 0.339, 0.104, 0.773, 0.241, 0.639, 0.791, 0.239, 0.283, 0.577, 0.261, 0.737, 0.876, 0.33, 0.991, 0.718, 0.044, 0.208, 0.555, 0.956, 0.516, 0.671, 0.096, 0.7]
global q = [0.839, 0.415, 0.804, 0.848, 0.807, 0.998, 0.96, 0.878, 0.999, 0.285, 0.742, 0.723, 0.981, 0.671, 0.889, 0.965, 0.973, 0.957, 0.778, 0.756, 0.968, 0.19, 0.924, 0.926, 0.846, 0.897, 0.71, 0.94, 0.512, 0.452, 0.956, 0.901, 0.754, 0.876, 0.552, 0.843, 0.73, 0.958, 0.564, 0.793, 0.496, 0.734, 0.121, 0.354, 0.834, 0.732, 0.857, 0.612, 0.958, 0.972, 0.747, 0.632, 0.487, 0.912, 0.57, 0.727, 0.392, 0.995, 0.94, 0.931, 0.436, 0.874, 0.811, 0.913, 0.422, 0.273, 0.874, 0.675, 0.673, 0.997, 0.317, 0.927, 0.452, 0.502, 0.942, 0.665, 0.898, 0.801, 0.983, 0.826, 0.989, 0.766, 0.671, 0.959, 0.956, 0.542, 0.602, 0.943, 0.997, 0.78, 0.655, 0.698, 0.497, 0.991, 0.808, 0.804, 0.486, 0.769, 0.973, 0.859, 0.972, 0.836, 0.313, 0.859, 0.378, 0.584, 0.795, 0.545, 0.593, 0.802, 0.267, 0.847, 0.925, 0.73, 0.921, 0.514, 0.915, 0.992, 0.604, 0.861, 0.888, 0.708, 0.897, 0.864, 0.539, 0.928, 0.953, 0.998, 0.683, 0.363, 0.29, 0.79, 0.819, 0.563, 0.99, 0.821, 0.613, 0.8, 0.929, 0.906, 0.585, 0.84, 0.949, 0.926, 0.905, 0.14, 0.979, 0.913, 0.597, 0.855, 0.954, 0.98, 0.998, 0.72, 0.782, 0.89, 0.52, 0.72, 0.693, 0.922, 0.554, 0.952, 0.597, 0.985, 0.625, 0.381, 0.973, 0.462, 0.713, 0.991, 0.578, 0.859, 0.849, 0.856, 0.868, 0.181, 0.791, 0.356, 0.912, 0.989, 0.959, 0.728, 0.661, 0.888, 0.546, 0.892, 0.231, 0.872, 0.335, 0.706, 0.865, 0.915, 0.914, 0.535, 0.864, 0.958, 0.351, 0.908, 0.441, 0.667, 0.496, 0.882, 0.871, 0.952, 0.965, 0.646, 0.962, 0.953, 0.921, 0.956, 0.487, 0.652, 0.736, 0.97, 0.249, 0.244, 0.729, 0.802, 0.996, 0.534, 0.67, 0.916, 0.745, 0.426, 0.757, 0.827, 0.445, 0.509, 0.956, 0.841, 0.854, 0.84, 0.695, 0.666, 0.894, 0.658, 0.758, 0.983, 0.828, 0.994, 0.783, 0.304, 0.454, 0.96, 0.971, 0.734, 0.764, 0.941, 0.978]
global origin = 1
global destination = 50