global arcs = [1 10; 1 17; 1 31; 1 42; 2 3; 2 16; 2 17; 2 19; 2 31; 2 40; 2 46; 2 53; 2 55; 3 12; 3 24; 3 26; 3 27; 3 55; 4 2; 4 8; 4 22; 4 30; 5 11; 5 20; 5 21; 5 31; 5 43; 5 45; 5 46; 5 60; 6 26; 6 29; 6 35; 6 41; 6 42; 6 44; 6 47; 6 49; 7 34; 7 56; 7 59; 7 60; 8 23; 8 43; 8 60; 9 16; 9 28; 9 50; 9 55; 10 9; 10 25; 10 29; 10 34; 10 48; 10 50; 10 55; 11 13; 11 23; 11 30; 11 31; 11 58; 12 3; 12 6; 12 30; 12 37; 12 39; 13 2; 13 10; 13 29; 13 31; 13 51; 13 53; 13 56; 13 57; 14 5; 14 12; 14 16; 14 25; 14 30; 14 38; 14 48; 14 50; 14 58; 15 20; 15 23; 15 26; 15 27; 15 29; 15 35; 15 41; 16 6; 16 22; 16 39; 16 57; 17 4; 17 15; 17 22; 17 37; 18 20; 18 34; 18 44; 18 57; 19 10; 19 12; 19 15; 19 23; 19 54; 19 60; 20 10; 20 11; 20 31; 20 34; 20 37; 21 5; 21 18; 21 38; 21 42; 21 47; 21 50; 21 57; 21 59; 21 60; 22 18; 22 19; 22 29; 23 16; 23 32; 23 41; 23 49; 24 13; 24 17; 24 26; 24 40; 24 47; 24 60; 25 5; 25 21; 25 22; 25 32; 25 43; 26 19; 26 21; 26 31; 26 33; 26 48; 26 51; 26 54; 26 56; 27 14; 27 36; 27 45; 27 56; 27 58; 28 17; 28 26; 28 38; 28 48; 28 57; 29 7; 29 8; 29 10; 29 14; 29 15; 29 17; 29 19; 29 22; 29 25; 29 33; 29 45; 29 51; 30 9; 30 46; 30 60; 31 13; 31 14; 31 24; 31 38; 31 40; 31 54; 31 56; 31 57; 32 8; 32 48; 32 51; 33 6; 33 25; 33 44; 33 46; 34 4; 34 19; 34 28; 34 35; 34 55; 34 58; 35 10; 35 21; 35 29; 35 45; 35 47; 35 59; 36 14; 36 18; 36 21; 36 29; 36 57; 36 58; 37 13; 37 17; 37 25; 37 52; 37 59; 38 46; 38 50; 38 60; 39 2; 39 17; 39 25; 39 28; 39 30; 39 40; 39 54; 40 15; 40 21; 40 36; 40 58; 41 9; 41 20; 41 49; 42 3; 42 17; 42 32; 42 34; 42 35; 42 39; 42 46; 42 54; 42 59; 43 10; 43 14; 43 22; 43 29; 43 33; 43 35; 43 46; 43 48; 43 53; 44 13; 44 25; 44 45; 45 15; 45 34; 45 35; 46 10; 46 16; 46 26; 46 30; 46 32; 46 41; 47 10; 47 14; 47 16; 47 19; 47 24; 47 25; 47 32; 47 53; 48 8; 48 9; 48 15; 48 39; 48 40; 48 44; 49 24; 49 37; 49 38; 49 45; 49 51; 50 3; 50 17; 50 35; 50 38; 50 53; 50 54; 51 20; 51 21; 51 27; 51 43; 52 26; 52 48; 53 2; 53 4; 53 9; 53 10; 53 14; 53 19; 53 23; 53 48; 53 54; 54 4; 54 7; 54 10; 54 15; 54 31; 54 59; 55 21; 55 37; 55 38; 55 42; 55 44; 55 45; 55 51; 55 53; 55 58; 55 59; 56 21; 56 29; 56 30; 56 41; 56 44; 56 45; 56 46; 57 2; 57 18; 57 21; 57 40; 57 53; 58 10; 58 18; 58 32; 58 35; 58 44; 58 45; 58 46; 58 48; 58 52; 58 53; 59 9; 59 11; 59 22; 59 25; 59 29; 59 39]
global d_x = [4.0, 5.0, 4.0, 5.0, 1.0, 7.0, 10.0, 8.0, 3.0, 10.0, 3.0, 7.0, 3.0, 4.0, 8.0, 2.0, 5.0, 2.0, 9.0, 5.0, 10.0, 6.0, 8.0, 5.0, 7.0, 5.0, 5.0, 10.0, 6.0, 8.0, 8.0, 4.0, 2.0, 8.0, 7.0, 8.0, 10.0, 7.0, 6.0, 10.0, 6.0, 7.0, 7.0, 1.0, 3.0, 4.0, 6.0, 8.0, 4.0, 2.0, 5.0, 5.0, 2.0, 2.0, 9.0, 9.0, 6.0, 6.0, 1.0, 2.0, 2.0, 6.0, 4.0, 3.0, 3.0, 4.0, 8.0, 1.0, 2.0, 2.0, 10.0, 1.0, 3.0, 3.0, 7.0, 5.0, 10.0, 6.0, 9.0, 7.0, 3.0, 1.0, 8.0, 5.0, 6.0, 7.0, 4.0, 5.0, 3.0, 3.0, 9.0, 5.0, 8.0, 3.0, 1.0, 9.0, 8.0, 6.0, 1.0, 3.0, 7.0, 8.0, 7.0, 6.0, 5.0, 8.0, 4.0, 10.0, 2.0, 7.0, 9.0, 3.0, 6.0, 8.0, 5.0, 5.0, 8.0, 3.0, 4.0, 5.0, 8.0, 1.0, 1.0, 4.0, 9.0, 3.0, 7.0, 3.0, 1.0, 6.0, 9.0, 4.0, 9.0, 4.0, 5.0, 7.0, 3.0, 10.0, 7.0, 10.0, 2.0, 4.0, 3.0, 8.0, 8.0, 9.0, 1.0, 3.0, 6.0, 5.0, 8.0, 3.0, 7.0, 7.0, 7.0, 1.0, 4.0, 3.0, 3.0, 8.0, 7.0, 1.0, 7.0, 8.0, 4.0, 9.0, 8.0, 5.0, 6.0, 6.0, 9.0, 7.0, 8.0, 9.0, 3.0, 9.0, 5.0, 6.0, 2.0, 1.0, 3.0, 3.0, 6.0, 8.0, 6.0, 1.0, 3.0, 9.0, 6.0, 4.0, 5.0, 9.0, 10.0, 1.0, 6.0, 10.0, 4.0, 4.0, 6.0, 7.0, 10.0, 4.0, 1.0, 6.0, 3.0, 7.0, 10.0, 6.0, 5.0, 10.0, 8.0, 5.0, 3.0, 5.0, 8.0, 3.0, 1.0, 8.0, 10.0, 3.0, 10.0, 5.0, 3.0, 5.0, 3.0, 10.0, 7.0, 4.0, 2.0, 5.0, 9.0, 1.0, 6.0, 4.0, 7.0, 3.0, 2.0, 3.0, 2.0, 2.0, 9.0, 6.0, 10.0, 9.0, 4.0, 8.0, 2.0, 8.0, 4.0, 8.0, 4.0, 9.0, 5.0, 3.0, 3.0, 8.0, 10.0, 2.0, 1.0, 9.0, 10.0, 6.0, 8.0, 2.0, 10.0, 3.0, 10.0, 4.0, 6.0, 5.0, 4.0, 4.0, 8.0, 7.0, 9.0, 3.0, 2.0, 7.0, 10.0, 2.0, 4.0, 6.0, 9.0, 2.0, 2.0, 10.0, 5.0, 5.0, 2.0, 4.0, 10.0, 1.0, 1.0, 3.0, 5.0, 10.0, 9.0, 5.0, 7.0, 9.0, 9.0, 9.0, 10.0, 7.0, 8.0, 5.0, 7.0, 3.0, 3.0, 9.0, 1.0, 4.0, 2.0, 10.0, 2.0, 1.0, 10.0, 2.0, 9.0, 9.0, 5.0, 2.0, 9.0, 2.0, 6.0, 5.0, 10.0, 3.0, 3.0, 7.0, 4.0, 7.0, 5.0, 5.0, 2.0, 9.0, 6.0, 10.0, 4.0, 7.0, 6.0, 6.0]
global b_x = 5
global d_y = [10.0, 5.0, 10.0, 5.0, 1.0, 3.0, 1.0, 4.0, 1.0, 3.0, 1.0, 9.0, 4.0, 6.0, 4.0, 8.0, 10.0, 4.0, 5.0, 4.0, 9.0, 4.0, 6.0, 8.0, 1.0, 1.0, 2.0, 6.0, 1.0, 7.0, 7.0, 9.0, 5.0, 10.0, 6.0, 2.0, 1.0, 10.0, 10.0, 10.0, 3.0, 4.0, 2.0, 5.0, 9.0, 8.0, 10.0, 5.0, 3.0, 1.0, 8.0, 5.0, 3.0, 10.0, 4.0, 8.0, 10.0, 3.0, 10.0, 7.0, 3.0, 7.0, 3.0, 8.0, 4.0, 10.0, 1.0, 5.0, 6.0, 8.0, 2.0, 10.0, 4.0, 9.0, 6.0, 10.0, 9.0, 6.0, 2.0, 6.0, 3.0, 6.0, 4.0, 1.0, 4.0, 6.0, 4.0, 2.0, 3.0, 10.0, 9.0, 7.0, 6.0, 9.0, 1.0, 7.0, 3.0, 6.0, 9.0, 9.0, 6.0, 8.0, 8.0, 9.0, 3.0, 1.0, 6.0, 3.0, 6.0, 9.0, 5.0, 1.0, 10.0, 6.0, 7.0, 7.0, 7.0, 10.0, 10.0, 10.0, 1.0, 2.0, 10.0, 10.0, 7.0, 7.0, 5.0, 5.0, 7.0, 7.0, 10.0, 4.0, 8.0, 2.0, 4.0, 3.0, 3.0, 3.0, 8.0, 8.0, 1.0, 3.0, 7.0, 1.0, 7.0, 4.0, 3.0, 8.0, 9.0, 6.0, 2.0, 3.0, 3.0, 10.0, 7.0, 4.0, 10.0, 3.0, 8.0, 8.0, 3.0, 8.0, 8.0, 1.0, 5.0, 10.0, 4.0, 6.0, 7.0, 10.0, 5.0, 2.0, 6.0, 10.0, 6.0, 6.0, 8.0, 9.0, 8.0, 6.0, 4.0, 4.0, 7.0, 2.0, 1.0, 7.0, 2.0, 1.0, 3.0, 1.0, 5.0, 3.0, 10.0, 3.0, 4.0, 10.0, 9.0, 3.0, 2.0, 6.0, 7.0, 10.0, 1.0, 5.0, 10.0, 4.0, 7.0, 5.0, 4.0, 4.0, 9.0, 3.0, 8.0, 7.0, 8.0, 10.0, 3.0, 10.0, 2.0, 1.0, 3.0, 3.0, 4.0, 10.0, 3.0, 5.0, 9.0, 7.0, 5.0, 7.0, 7.0, 9.0, 3.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 3.0, 5.0, 4.0, 8.0, 4.0, 2.0, 7.0, 6.0, 6.0, 6.0, 3.0, 6.0, 5.0, 3.0, 9.0, 7.0, 3.0, 4.0, 6.0, 7.0, 6.0, 1.0, 1.0, 6.0, 6.0, 1.0, 1.0, 4.0, 10.0, 10.0, 1.0, 3.0, 9.0, 3.0, 3.0, 2.0, 10.0, 10.0, 8.0, 1.0, 4.0, 3.0, 10.0, 7.0, 9.0, 5.0, 4.0, 9.0, 3.0, 7.0, 6.0, 8.0, 1.0, 2.0, 10.0, 10.0, 8.0, 7.0, 8.0, 7.0, 2.0, 7.0, 5.0, 8.0, 6.0, 3.0, 4.0, 2.0, 10.0, 3.0, 6.0, 2.0, 10.0, 10.0, 7.0, 8.0, 6.0, 8.0, 6.0, 3.0, 4.0, 6.0, 3.0, 10.0, 7.0, 9.0, 10.0, 9.0, 1.0, 1.0, 7.0, 1.0, 10.0, 6.0, 3.0, 5.0, 6.0, 10.0, 6.0, 8.0, 1.0, 8.0, 5.0]
global b_y = 10
global p = [0.157, 0.825, 0.279, 0.479, 0.736, 0.069, 0.931, 0.896, 0.107, 0.633, 0.29, 0.527, 0.81, 0.756, 0.934, 0.272, 0.644, 0.328, 0.648, 0.814, 0.927, 0.691, 0.787, 0.665, 0.492, 0.787, 0.749, 0.9, 0.509, 0.597, 0.213, 0.663, 0.717, 0.239, 0.665, 0.03, 0.082, 0.763, 0.212, 0.784, 0.126, 0.879, 0.805, 0.151, 0.889, 0.236, 0.213, 0.973, 0.73, 0.798, 0.262, 0.079, 0.765, 0.499, 0.652, 0.222, 0.591, 0.269, 0.436, 0.558, 0.677, 0.13, 0.901, 0.93, 0.037, 0.241, 0.745, 0.431, 0.001, 0.195, 0.221, 0.776, 0.474, 0.165, 0.678, 0.444, 0.717, 0.787, 0.781, 0.951, 0.48, 0.719, 0.125, 0.874, 0.005, 0.451, 0.54, 0.745, 0.493, 0.449, 0.635, 0.662, 0.055, 0.988, 0.517, 0.393, 0.465, 0.789, 0.295, 0.219, 0.08, 0.102, 0.932, 0.254, 0.002, 0.581, 0.525, 0.273, 0.046, 0.049, 0.15, 0.215, 0.174, 0.653, 0.714, 0.209, 0.628, 0.31, 0.951, 0.076, 0.215, 0.511, 0.283, 0.405, 0.389, 0.711, 0.761, 0.454, 0.14, 0.296, 0.667, 0.84, 0.492, 0.068, 0.554, 0.591, 0.608, 0.055, 0.574, 0.047, 0.011, 0.692, 0.839, 0.295, 0.949, 0.939, 0.9, 0.645, 0.052, 0.683, 0.115, 0.155, 0.738, 0.026, 0.246, 0.381, 0.792, 0.263, 0.353, 0.404, 0.72, 0.559, 0.897, 0.054, 0.323, 0.706, 0.618, 0.625, 0.746, 0.452, 0.23, 0.325, 0.734, 0.442, 0.743, 0.828, 0.753, 0.82, 0.417, 0.837, 0.713, 0.475, 0.481, 0.026, 0.854, 0.634, 0.226, 0.996, 0.51, 0.16, 0.889, 0.68, 0.961, 0.667, 0.61, 0.122, 0.627, 0.691, 0.128, 0.225, 0.406, 0.795, 0.449, 0.404, 0.072, 0.398, 0.851, 0.955, 0.644, 0.043, 0.561, 0.229, 0.965, 0.582, 0.528, 0.245, 0.948, 0.189, 0.747, 0.065, 0.184, 0.518, 0.151, 0.968, 0.155, 0.791, 0.848, 0.713, 0.205, 0.747, 0.632, 0.99, 0.534, 0.583, 0.859, 0.192, 0.315, 0.57, 0.899, 0.7, 0.117, 0.193, 0.63, 0.26, 0.16, 0.04, 0.804, 0.755, 0.442, 0.924, 0.189, 0.527, 0.807, 0.541, 0.713, 0.901, 0.393, 0.951, 0.182, 0.649, 0.969, 0.446, 0.198, 0.202, 0.028, 0.298, 0.674, 0.659, 0.987, 0.766, 0.416, 0.843, 0.439, 0.819, 0.435, 0.867, 0.174, 0.054, 0.167, 0.137, 0.11, 0.005, 0.964, 0.664, 0.809, 0.978, 0.5, 0.835, 0.652, 0.858, 0.505, 0.242, 0.28, 0.838, 0.722, 0.898, 0.66, 0.689, 0.869, 0.158, 0.142, 0.223, 0.547, 0.862, 0.77, 0.509, 0.002, 0.795, 0.244, 0.531, 0.555, 0.261, 0.496, 0.811, 0.601, 0.316, 0.685, 0.503, 0.438, 0.479, 0.812, 0.449, 0.561, 0.121, 0.43, 0.334, 0.86, 0.92, 0.798, 0.1, 0.861, 0.884, 0.037, 0.779, 0.345, 0.273, 0.989, 0.922, 0.83, 0.269, 0.125, 0.352]
global q = [0.377, 0.904, 0.777, 0.72, 0.819, 0.563, 0.965, 0.99, 0.292, 0.897, 0.445, 0.846, 0.948, 0.898, 0.968, 0.595, 0.925, 0.473, 0.924, 0.912, 0.995, 0.794, 0.864, 0.766, 0.959, 0.995, 0.819, 0.96, 0.6, 0.599, 0.997, 0.68, 0.718, 0.646, 0.678, 0.33, 0.672, 0.988, 0.822, 0.81, 0.673, 0.94, 0.854, 0.428, 0.964, 0.675, 0.512, 0.989, 0.801, 0.912, 0.899, 0.782, 0.776, 0.6, 0.963, 0.55, 0.909, 0.439, 0.596, 0.96, 0.821, 0.239, 0.995, 0.992, 0.149, 0.976, 0.774, 0.613, 0.741, 0.447, 0.426, 0.901, 0.897, 0.676, 0.934, 0.849, 0.853, 0.931, 0.785, 0.98, 0.728, 0.734, 0.945, 0.996, 0.181, 0.952, 0.988, 0.957, 0.811, 0.814, 0.92, 0.974, 0.748, 0.989, 0.67, 0.713, 0.708, 0.791, 0.613, 0.349, 0.659, 0.933, 0.939, 0.63, 0.88, 0.991, 0.969, 0.773, 0.23, 0.373, 0.602, 0.501, 0.86, 0.876, 0.791, 0.969, 0.79, 0.918, 0.991, 0.299, 0.561, 0.971, 0.922, 0.774, 0.892, 0.732, 0.764, 0.887, 0.397, 0.319, 0.769, 0.932, 0.755, 0.981, 0.682, 0.708, 0.71, 0.345, 0.588, 0.865, 0.044, 0.902, 0.897, 0.584, 0.955, 0.959, 0.981, 0.816, 0.198, 0.919, 0.676, 0.625, 0.888, 0.097, 0.83, 0.962, 0.803, 0.577, 0.5, 0.589, 0.816, 0.628, 0.989, 0.959, 0.471, 0.977, 0.717, 0.975, 0.752, 0.889, 0.833, 0.952, 0.827, 0.454, 0.952, 0.912, 0.879, 0.942, 0.863, 0.986, 0.989, 0.533, 0.635, 0.3, 0.948, 0.734, 0.831, 0.998, 0.896, 0.206, 0.899, 0.812, 0.998, 0.724, 0.713, 0.489, 0.775, 0.851, 0.341, 0.427, 0.772, 0.893, 0.457, 0.61, 0.468, 0.829, 0.965, 0.973, 0.963, 0.449, 0.721, 0.853, 0.991, 0.602, 0.981, 0.482, 0.977, 0.316, 0.752, 0.883, 0.937, 0.566, 0.739, 0.976, 0.302, 0.865, 0.924, 0.773, 0.45, 0.78, 0.854, 0.992, 0.62, 0.803, 0.942, 0.299, 0.316, 0.689, 0.99, 0.994, 0.852, 0.432, 0.809, 0.351, 0.982, 0.26, 0.848, 0.992, 0.981, 0.949, 0.225, 0.741, 0.949, 0.663, 0.989, 0.957, 0.577, 0.958, 0.576, 0.89, 0.971, 0.472, 0.3, 0.632, 0.652, 0.58, 0.91, 0.672, 0.997, 0.977, 0.54, 0.99, 0.5, 0.94, 0.9, 0.966, 0.444, 0.351, 0.206, 0.389, 0.549, 0.585, 0.992, 0.839, 0.854, 0.978, 0.69, 0.946, 0.782, 0.958, 0.705, 0.858, 0.772, 0.916, 0.727, 0.962, 0.963, 0.97, 0.889, 0.329, 0.889, 0.374, 0.718, 0.957, 0.859, 0.99, 0.016, 0.941, 0.917, 0.935, 0.923, 0.77, 0.996, 0.982, 0.982, 0.909, 0.951, 0.818, 0.704, 0.584, 0.937, 0.498, 0.592, 0.808, 0.457, 0.641, 0.871, 0.96, 0.832, 0.796, 0.954, 0.903, 0.209, 0.833, 0.873, 0.844, 0.996, 0.948, 0.873, 0.618, 0.174, 0.722]
global origin = 1
global destination = 60