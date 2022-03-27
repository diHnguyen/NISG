global arcs = [1 28; 1 30; 1 36; 1 39; 2 5; 2 16; 2 20; 2 21; 2 36; 2 40; 2 48; 3 36; 3 39; 3 45; 4 5; 4 14; 4 33; 4 38; 4 43; 5 8; 5 23; 5 30; 5 33; 5 46; 6 4; 6 13; 6 38; 6 45; 6 47; 6 49; 7 8; 7 26; 7 39; 7 49; 8 6; 8 24; 8 30; 8 38; 8 43; 9 5; 9 15; 9 16; 9 18; 9 25; 9 30; 9 46; 10 16; 10 19; 10 36; 10 41; 10 42; 10 44; 11 2; 11 15; 11 26; 11 27; 11 29; 11 32; 11 39; 11 46; 12 19; 12 27; 12 34; 12 36; 12 48; 13 11; 13 15; 13 20; 13 26; 13 28; 13 30; 13 43; 14 18; 14 19; 14 25; 14 30; 14 36; 14 41; 14 43; 14 47; 15 14; 15 19; 15 20; 15 21; 15 23; 15 33; 16 14; 16 17; 16 21; 17 6; 17 12; 17 15; 17 35; 17 41; 17 46; 18 2; 18 5; 18 11; 18 27; 18 34; 18 46; 18 48; 19 4; 19 10; 19 16; 19 24; 19 35; 20 3; 20 38; 20 41; 21 7; 21 9; 21 10; 21 12; 21 22; 21 41; 21 42; 22 13; 22 14; 22 16; 22 39; 22 41; 23 8; 23 11; 23 21; 23 24; 23 25; 23 44; 24 6; 24 17; 24 38; 24 43; 25 3; 25 5; 25 16; 25 28; 25 35; 25 47; 25 49; 26 5; 26 7; 26 21; 26 27; 26 36; 26 50; 27 22; 27 24; 27 25; 27 28; 27 33; 27 34; 27 35; 27 37; 27 42; 27 44; 27 46; 28 3; 28 9; 28 31; 29 11; 29 36; 29 47; 30 2; 30 10; 30 11; 30 13; 30 14; 30 19; 30 50; 31 10; 31 27; 31 28; 31 35; 31 48; 32 10; 32 24; 32 25; 32 31; 32 37; 32 40; 33 13; 33 19; 33 21; 33 26; 33 38; 34 10; 34 12; 34 13; 34 16; 34 17; 34 20; 34 46; 35 7; 35 28; 35 30; 36 5; 36 11; 36 31; 36 40; 36 49; 37 25; 37 27; 38 12; 38 21; 38 40; 38 41; 38 44; 38 50; 39 7; 39 20; 39 24; 39 41; 39 42; 40 5; 40 7; 40 8; 40 26; 40 27; 40 36; 40 47; 41 13; 41 36; 41 47; 42 18; 42 26; 42 48; 43 3; 43 15; 44 19; 44 21; 44 25; 44 39; 45 5; 45 16; 45 17; 45 22; 45 26; 45 41; 45 50; 46 12; 46 23; 46 25; 46 29; 46 34; 46 35; 46 37; 47 9; 47 18; 47 20; 47 27; 47 30; 48 11; 48 49; 49 4; 49 8; 49 17; 49 27]
global d_x = [4.0, 5.0, 3.0, 5.0, 4.0, 5.0, 10.0, 2.0, 6.0, 5.0, 10.0, 10.0, 6.0, 1.0, 2.0, 1.0, 1.0, 10.0, 8.0, 1.0, 9.0, 9.0, 10.0, 3.0, 1.0, 5.0, 5.0, 1.0, 2.0, 8.0, 8.0, 3.0, 8.0, 10.0, 6.0, 7.0, 8.0, 2.0, 8.0, 5.0, 2.0, 4.0, 4.0, 1.0, 1.0, 3.0, 9.0, 6.0, 7.0, 10.0, 1.0, 5.0, 7.0, 1.0, 2.0, 7.0, 9.0, 6.0, 8.0, 3.0, 9.0, 6.0, 4.0, 6.0, 8.0, 1.0, 6.0, 9.0, 7.0, 7.0, 7.0, 8.0, 1.0, 3.0, 8.0, 4.0, 4.0, 5.0, 5.0, 8.0, 6.0, 5.0, 10.0, 7.0, 9.0, 8.0, 4.0, 8.0, 2.0, 10.0, 2.0, 1.0, 9.0, 4.0, 2.0, 9.0, 7.0, 6.0, 4.0, 9.0, 5.0, 3.0, 4.0, 7.0, 10.0, 4.0, 1.0, 4.0, 3.0, 4.0, 10.0, 4.0, 2.0, 4.0, 4.0, 8.0, 9.0, 2.0, 6.0, 5.0, 1.0, 10.0, 6.0, 3.0, 10.0, 10.0, 4.0, 1.0, 3.0, 3.0, 9.0, 8.0, 9.0, 9.0, 7.0, 3.0, 2.0, 7.0, 3.0, 5.0, 4.0, 7.0, 3.0, 1.0, 8.0, 1.0, 10.0, 5.0, 8.0, 4.0, 7.0, 4.0, 10.0, 10.0, 2.0, 10.0, 10.0, 1.0, 5.0, 7.0, 6.0, 10.0, 2.0, 1.0, 3.0, 9.0, 9.0, 9.0, 5.0, 9.0, 8.0, 10.0, 2.0, 1.0, 10.0, 9.0, 2.0, 6.0, 9.0, 8.0, 4.0, 4.0, 4.0, 10.0, 2.0, 5.0, 7.0, 3.0, 6.0, 3.0, 3.0, 5.0, 8.0, 5.0, 7.0, 10.0, 10.0, 7.0, 6.0, 3.0, 10.0, 3.0, 7.0, 4.0, 6.0, 4.0, 1.0, 6.0, 4.0, 5.0, 8.0, 6.0, 5.0, 4.0, 5.0, 8.0, 8.0, 3.0, 1.0, 9.0, 8.0, 9.0, 7.0, 1.0, 1.0, 1.0, 4.0, 6.0, 6.0, 7.0, 8.0, 9.0, 9.0, 1.0, 9.0, 8.0, 8.0, 9.0, 5.0, 4.0, 8.0, 5.0, 5.0, 7.0, 4.0, 9.0, 10.0, 10.0, 10.0, 10.0, 10.0, 2.0, 3.0, 8.0, 5.0, 1.0, 7.0]
global b_x = 5
global d_y = [6.0, 6.0, 10.0, 8.0, 3.0, 3.0, 8.0, 9.0, 2.0, 1.0, 2.0, 9.0, 8.0, 1.0, 8.0, 10.0, 5.0, 7.0, 10.0, 2.0, 2.0, 4.0, 10.0, 10.0, 10.0, 8.0, 10.0, 3.0, 8.0, 8.0, 5.0, 2.0, 2.0, 3.0, 10.0, 4.0, 2.0, 7.0, 9.0, 2.0, 7.0, 8.0, 7.0, 5.0, 5.0, 2.0, 7.0, 4.0, 10.0, 8.0, 10.0, 10.0, 2.0, 2.0, 9.0, 8.0, 8.0, 6.0, 5.0, 10.0, 8.0, 5.0, 6.0, 10.0, 5.0, 9.0, 6.0, 10.0, 6.0, 9.0, 10.0, 7.0, 3.0, 4.0, 7.0, 6.0, 7.0, 6.0, 6.0, 5.0, 1.0, 5.0, 7.0, 2.0, 8.0, 9.0, 2.0, 3.0, 5.0, 10.0, 2.0, 2.0, 7.0, 9.0, 2.0, 10.0, 7.0, 3.0, 2.0, 5.0, 10.0, 9.0, 4.0, 2.0, 7.0, 10.0, 2.0, 6.0, 10.0, 3.0, 5.0, 4.0, 2.0, 9.0, 5.0, 9.0, 9.0, 5.0, 3.0, 1.0, 5.0, 10.0, 7.0, 6.0, 8.0, 10.0, 10.0, 5.0, 7.0, 1.0, 4.0, 6.0, 2.0, 2.0, 1.0, 5.0, 6.0, 1.0, 4.0, 4.0, 10.0, 3.0, 1.0, 4.0, 4.0, 3.0, 6.0, 7.0, 6.0, 2.0, 10.0, 6.0, 1.0, 8.0, 6.0, 9.0, 5.0, 3.0, 8.0, 9.0, 10.0, 6.0, 10.0, 7.0, 6.0, 8.0, 6.0, 9.0, 7.0, 3.0, 6.0, 1.0, 7.0, 2.0, 8.0, 3.0, 8.0, 5.0, 5.0, 10.0, 6.0, 1.0, 7.0, 2.0, 4.0, 5.0, 9.0, 7.0, 10.0, 6.0, 7.0, 2.0, 9.0, 8.0, 7.0, 5.0, 4.0, 1.0, 10.0, 2.0, 3.0, 7.0, 6.0, 2.0, 7.0, 1.0, 3.0, 2.0, 8.0, 10.0, 6.0, 1.0, 6.0, 6.0, 5.0, 7.0, 4.0, 6.0, 6.0, 8.0, 2.0, 3.0, 5.0, 1.0, 6.0, 7.0, 5.0, 2.0, 8.0, 3.0, 1.0, 2.0, 10.0, 7.0, 6.0, 8.0, 2.0, 9.0, 4.0, 6.0, 7.0, 4.0, 5.0, 3.0, 3.0, 10.0, 6.0, 2.0, 9.0, 2.0, 10.0, 4.0, 10.0, 6.0, 4.0, 5.0, 10.0]
global b_y = 10
global p = [0.722, 0.668, 0.126, 0.222, 0.884, 0.205, 0.03, 0.27, 0.108, 0.377, 0.018, 0.905, 0.992, 0.608, 0.737, 0.254, 0.772, 0.689, 0.596, 0.81, 0.554, 0.777, 0.061, 0.277, 0.775, 0.555, 0.547, 0.052, 0.872, 0.015, 0.748, 0.342, 0.221, 0.345, 0.11, 0.767, 0.343, 0.37, 0.558, 0.703, 0.392, 0.404, 0.209, 0.537, 0.713, 0.862, 0.163, 0.48, 0.033, 0.842, 0.176, 0.326, 0.366, 0.879, 0.947, 0.796, 0.064, 0.264, 0.347, 0.413, 0.968, 0.831, 0.991, 0.392, 0.993, 0.839, 0.601, 0.438, 0.344, 0.632, 0.634, 0.22, 0.702, 0.959, 0.811, 0.999, 0.113, 0.971, 0.395, 0.284, 0.062, 0.952, 0.438, 0.463, 0.473, 0.757, 0.62, 0.542, 0.938, 0.766, 0.031, 0.147, 0.641, 0.622, 0.335, 0.565, 0.839, 0.111, 0.195, 0.929, 0.784, 0.021, 0.773, 0.897, 0.911, 0.553, 0.671, 0.131, 0.884, 0.098, 0.039, 0.404, 0.533, 0.074, 0.846, 0.461, 0.632, 0.784, 0.158, 0.52, 0.358, 0.803, 0.661, 0.3, 0.17, 0.567, 0.998, 0.937, 0.672, 0.581, 0.819, 0.791, 0.71, 0.841, 0.832, 0.103, 0.544, 0.234, 0.032, 0.517, 0.784, 0.648, 0.581, 0.55, 0.542, 0.794, 0.547, 0.528, 0.511, 0.116, 0.278, 0.823, 0.724, 0.27, 0.744, 0.97, 0.342, 0.163, 0.087, 0.683, 0.404, 0.781, 0.781, 0.36, 0.341, 0.245, 0.705, 0.565, 0.247, 0.621, 0.507, 0.608, 0.767, 0.611, 0.996, 0.92, 0.023, 0.055, 0.4, 0.087, 0.219, 0.216, 0.866, 0.392, 0.893, 0.461, 0.817, 0.361, 0.507, 0.465, 0.076, 0.663, 0.934, 0.351, 0.791, 0.645, 0.084, 0.337, 0.622, 0.684, 0.938, 0.44, 0.748, 0.92, 0.455, 0.377, 0.399, 0.711, 0.686, 0.502, 0.04, 0.091, 0.299, 0.982, 0.548, 0.635, 0.649, 0.272, 0.608, 0.083, 0.177, 0.65, 0.051, 0.212, 0.621, 0.655, 0.981, 0.692, 0.812, 0.941, 0.95, 0.956, 0.174, 0.955, 0.586, 0.783, 0.142, 0.911, 0.883, 0.665, 0.177, 0.477, 0.187, 0.449, 0.123, 0.445, 0.5, 0.434, 0.048, 0.191, 0.455, 0.327, 0.41, 0.703, 0.059, 0.289, 0.439]
global q = [0.819, 0.979, 0.317, 0.787, 0.968, 0.407, 0.415, 0.443, 0.588, 0.862, 0.544, 0.933, 0.997, 0.732, 0.887, 0.693, 0.967, 0.998, 0.703, 0.85, 0.872, 0.822, 0.781, 0.355, 0.814, 0.968, 0.737, 0.4, 0.946, 0.526, 0.976, 0.942, 0.537, 0.678, 0.767, 0.964, 0.841, 0.945, 0.768, 0.876, 0.63, 0.481, 0.214, 0.725, 0.899, 0.993, 0.723, 0.89, 0.593, 0.91, 0.865, 0.56, 0.744, 0.921, 0.983, 0.92, 0.689, 0.723, 0.831, 0.564, 0.971, 0.834, 0.992, 0.603, 0.993, 0.983, 0.779, 0.502, 0.543, 0.723, 0.996, 0.733, 0.996, 0.968, 0.939, 0.999, 0.192, 0.973, 0.822, 0.33, 0.603, 0.998, 0.995, 0.486, 0.581, 0.767, 0.761, 0.77, 0.943, 0.792, 0.91, 0.894, 0.838, 0.954, 0.884, 0.656, 0.9, 0.416, 0.594, 0.934, 0.812, 0.622, 0.918, 0.992, 0.922, 0.815, 0.748, 0.915, 0.927, 0.247, 0.465, 0.404, 0.885, 0.483, 0.998, 0.827, 0.837, 0.942, 0.695, 0.998, 0.753, 0.939, 0.78, 0.933, 0.346, 0.899, 0.999, 0.941, 0.781, 0.661, 0.851, 0.83, 0.847, 0.944, 0.929, 0.414, 0.98, 0.982, 0.803, 0.884, 0.988, 0.85, 0.836, 0.941, 0.672, 0.883, 0.846, 0.636, 0.53, 0.237, 0.678, 0.883, 0.887, 0.357, 0.919, 0.99, 0.952, 0.214, 0.489, 0.941, 0.932, 0.864, 0.91, 0.469, 0.531, 0.804, 0.729, 0.833, 0.383, 0.625, 0.522, 0.945, 0.817, 0.66, 0.999, 0.951, 0.153, 0.262, 0.591, 0.434, 0.433, 0.641, 0.898, 0.826, 0.994, 0.531, 0.996, 0.383, 0.713, 0.91, 0.454, 0.766, 0.962, 0.666, 0.896, 0.669, 0.653, 0.873, 0.927, 0.904, 0.95, 0.585, 0.849, 0.996, 0.584, 0.657, 0.728, 0.952, 0.795, 0.517, 0.982, 0.727, 0.455, 0.983, 0.843, 0.716, 0.781, 0.54, 0.942, 0.333, 0.338, 0.967, 0.703, 0.318, 0.843, 0.732, 0.989, 0.899, 0.923, 0.968, 0.956, 0.964, 0.69, 0.964, 0.768, 0.964, 0.937, 0.939, 0.897, 0.995, 0.838, 0.551, 0.252, 0.699, 0.675, 0.485, 0.996, 0.569, 0.674, 0.437, 0.825, 0.962, 0.86, 0.762, 0.344, 0.899, 0.991]
global origin = 1
global destination = 50