global arcs = [1 35; 1 38; 1 46; 1 49; 1 53; 2 3; 2 8; 2 11; 2 22; 2 24; 2 47; 3 23; 3 52; 3 58; 4 15; 4 24; 4 29; 4 39; 5 9; 5 10; 5 17; 5 37; 5 39; 5 51; 5 53; 6 8; 6 11; 6 24; 6 40; 6 50; 7 9; 8 3; 8 5; 8 12; 8 29; 8 40; 8 59; 9 4; 9 6; 9 25; 9 28; 9 48; 9 49; 9 58; 10 22; 10 24; 10 31; 10 54; 11 15; 11 24; 11 36; 12 11; 12 20; 12 27; 12 28; 12 33; 12 35; 12 42; 12 56; 12 59; 13 19; 13 27; 13 33; 13 34; 13 54; 14 12; 14 16; 14 38; 14 51; 14 52; 15 10; 15 23; 15 39; 15 59; 16 7; 16 20; 16 25; 16 26; 16 27; 16 38; 16 55; 16 60; 17 18; 17 38; 18 5; 18 7; 18 25; 18 31; 18 36; 18 47; 19 16; 19 23; 19 38; 19 42; 19 44; 19 49; 19 51; 19 58; 19 59; 20 4; 20 5; 20 11; 20 17; 20 23; 20 34; 20 42; 20 46; 21 22; 21 23; 21 31; 21 37; 21 41; 21 48; 22 9; 22 15; 22 20; 22 51; 22 55; 23 27; 23 52; 24 6; 24 10; 24 11; 24 22; 24 33; 24 57; 25 4; 25 9; 25 16; 25 23; 25 30; 25 34; 25 39; 25 44; 26 27; 26 28; 26 34; 26 38; 27 4; 27 9; 27 13; 27 20; 27 29; 27 55; 27 57; 27 58; 27 59; 28 4; 28 9; 28 12; 28 13; 28 18; 28 23; 28 37; 28 56; 28 57; 29 58; 29 60; 30 10; 30 13; 30 39; 31 4; 31 5; 31 10; 31 14; 31 28; 31 41; 31 56; 31 58; 32 3; 32 9; 32 19; 32 29; 32 31; 32 33; 32 40; 32 47; 32 57; 32 59; 33 38; 33 43; 34 12; 34 15; 34 23; 34 33; 34 43; 34 59; 35 7; 35 23; 35 28; 36 10; 36 27; 36 34; 36 55; 36 58; 37 27; 37 29; 38 16; 38 21; 38 22; 38 26; 38 33; 39 2; 39 11; 39 13; 39 21; 39 24; 39 26; 39 28; 39 60; 40 21; 40 31; 40 48; 41 15; 41 35; 41 45; 41 60; 42 6; 42 12; 42 13; 42 30; 42 34; 42 40; 42 41; 42 58; 42 59; 43 13; 43 16; 43 29; 43 47; 44 9; 44 10; 44 11; 44 12; 44 21; 44 25; 44 36; 44 39; 44 60; 45 47; 45 49; 45 55; 46 3; 46 12; 46 35; 46 41; 46 45; 47 12; 47 14; 47 32; 47 33; 47 38; 47 48; 47 50; 47 57; 48 2; 48 23; 48 24; 48 49; 49 25; 49 40; 49 51; 49 57; 50 15; 50 18; 50 31; 50 33; 50 37; 50 54; 51 4; 51 5; 51 8; 51 24; 51 27; 51 37; 52 5; 52 7; 52 18; 52 33; 52 47; 53 8; 53 14; 53 15; 53 22; 53 38; 53 45; 53 55; 53 56; 54 7; 54 14; 54 21; 54 22; 54 32; 54 42; 55 9; 55 16; 55 19; 55 22; 55 29; 55 43; 56 33; 56 45; 56 48; 56 53; 56 58; 57 4; 57 6; 57 8; 57 19; 57 23; 57 38; 57 44; 57 54; 58 16; 58 31; 58 33; 58 41; 59 11; 59 38; 59 44; 59 49; 59 54; 59 57]
global d_x = [5.0, 7.0, 2.0, 4.0, 1.0, 7.0, 1.0, 1.0, 9.0, 5.0, 2.0, 9.0, 10.0, 7.0, 10.0, 4.0, 4.0, 8.0, 3.0, 10.0, 3.0, 6.0, 2.0, 10.0, 5.0, 9.0, 8.0, 8.0, 10.0, 1.0, 4.0, 1.0, 3.0, 4.0, 8.0, 5.0, 4.0, 5.0, 8.0, 9.0, 5.0, 5.0, 7.0, 1.0, 10.0, 5.0, 3.0, 6.0, 10.0, 7.0, 2.0, 7.0, 6.0, 5.0, 6.0, 6.0, 5.0, 6.0, 6.0, 8.0, 2.0, 4.0, 3.0, 1.0, 7.0, 7.0, 8.0, 7.0, 8.0, 6.0, 6.0, 6.0, 9.0, 6.0, 1.0, 7.0, 4.0, 8.0, 4.0, 2.0, 9.0, 5.0, 3.0, 3.0, 5.0, 4.0, 9.0, 1.0, 1.0, 7.0, 7.0, 1.0, 6.0, 7.0, 10.0, 2.0, 1.0, 8.0, 1.0, 1.0, 2.0, 3.0, 3.0, 1.0, 3.0, 2.0, 9.0, 7.0, 7.0, 4.0, 10.0, 7.0, 3.0, 6.0, 3.0, 8.0, 9.0, 9.0, 9.0, 9.0, 8.0, 5.0, 10.0, 5.0, 4.0, 7.0, 5.0, 1.0, 3.0, 4.0, 2.0, 6.0, 3.0, 7.0, 5.0, 3.0, 8.0, 2.0, 7.0, 9.0, 3.0, 10.0, 9.0, 3.0, 1.0, 1.0, 5.0, 9.0, 9.0, 2.0, 3.0, 10.0, 10.0, 9.0, 9.0, 9.0, 5.0, 7.0, 1.0, 5.0, 9.0, 2.0, 7.0, 4.0, 8.0, 8.0, 4.0, 4.0, 9.0, 6.0, 3.0, 7.0, 3.0, 7.0, 2.0, 7.0, 5.0, 4.0, 3.0, 6.0, 4.0, 1.0, 2.0, 9.0, 7.0, 8.0, 7.0, 10.0, 1.0, 2.0, 7.0, 7.0, 7.0, 7.0, 1.0, 7.0, 7.0, 5.0, 1.0, 9.0, 5.0, 4.0, 9.0, 6.0, 4.0, 10.0, 1.0, 4.0, 8.0, 2.0, 3.0, 7.0, 3.0, 6.0, 10.0, 1.0, 4.0, 10.0, 3.0, 9.0, 3.0, 4.0, 4.0, 10.0, 10.0, 1.0, 2.0, 7.0, 9.0, 5.0, 1.0, 7.0, 9.0, 4.0, 7.0, 6.0, 4.0, 6.0, 10.0, 4.0, 7.0, 9.0, 3.0, 9.0, 1.0, 10.0, 6.0, 1.0, 9.0, 2.0, 1.0, 3.0, 10.0, 4.0, 6.0, 10.0, 6.0, 2.0, 10.0, 5.0, 2.0, 1.0, 6.0, 10.0, 6.0, 4.0, 10.0, 1.0, 8.0, 5.0, 3.0, 2.0, 4.0, 5.0, 10.0, 9.0, 4.0, 2.0, 2.0, 10.0, 6.0, 9.0, 4.0, 4.0, 4.0, 10.0, 9.0, 9.0, 10.0, 4.0, 4.0, 8.0, 3.0, 4.0, 4.0, 1.0, 6.0, 8.0, 3.0, 9.0, 5.0, 5.0, 6.0, 2.0, 8.0, 7.0, 8.0, 3.0, 10.0, 10.0, 10.0, 5.0, 2.0, 8.0, 5.0, 8.0, 4.0, 8.0, 2.0, 1.0, 1.0, 2.0, 8.0]
global b_x = 5
global d_y = [1.0, 10.0, 6.0, 10.0, 9.0, 2.0, 7.0, 9.0, 4.0, 2.0, 7.0, 3.0, 9.0, 8.0, 8.0, 5.0, 5.0, 2.0, 5.0, 8.0, 8.0, 2.0, 5.0, 2.0, 1.0, 1.0, 2.0, 3.0, 6.0, 3.0, 7.0, 8.0, 6.0, 2.0, 5.0, 6.0, 7.0, 6.0, 5.0, 4.0, 8.0, 10.0, 2.0, 7.0, 6.0, 3.0, 6.0, 2.0, 4.0, 10.0, 1.0, 1.0, 9.0, 9.0, 10.0, 9.0, 4.0, 10.0, 9.0, 7.0, 5.0, 7.0, 10.0, 8.0, 5.0, 6.0, 3.0, 8.0, 6.0, 7.0, 8.0, 8.0, 5.0, 7.0, 2.0, 3.0, 10.0, 5.0, 1.0, 9.0, 5.0, 8.0, 7.0, 6.0, 1.0, 1.0, 5.0, 6.0, 6.0, 8.0, 1.0, 1.0, 10.0, 9.0, 8.0, 2.0, 7.0, 2.0, 7.0, 8.0, 4.0, 6.0, 5.0, 3.0, 4.0, 4.0, 8.0, 3.0, 6.0, 2.0, 10.0, 8.0, 1.0, 7.0, 1.0, 8.0, 4.0, 9.0, 10.0, 3.0, 5.0, 1.0, 6.0, 2.0, 4.0, 7.0, 2.0, 2.0, 10.0, 2.0, 10.0, 1.0, 8.0, 9.0, 2.0, 8.0, 7.0, 1.0, 2.0, 6.0, 6.0, 4.0, 3.0, 5.0, 7.0, 9.0, 9.0, 3.0, 3.0, 5.0, 9.0, 4.0, 5.0, 4.0, 2.0, 8.0, 10.0, 6.0, 8.0, 10.0, 10.0, 8.0, 2.0, 10.0, 6.0, 1.0, 7.0, 5.0, 5.0, 5.0, 10.0, 1.0, 3.0, 9.0, 3.0, 3.0, 8.0, 3.0, 10.0, 10.0, 1.0, 3.0, 5.0, 7.0, 6.0, 7.0, 7.0, 10.0, 8.0, 6.0, 4.0, 9.0, 4.0, 9.0, 4.0, 7.0, 10.0, 2.0, 7.0, 1.0, 9.0, 7.0, 1.0, 10.0, 10.0, 7.0, 4.0, 9.0, 9.0, 6.0, 7.0, 2.0, 4.0, 6.0, 10.0, 1.0, 2.0, 2.0, 10.0, 8.0, 6.0, 6.0, 3.0, 1.0, 6.0, 7.0, 8.0, 10.0, 5.0, 7.0, 5.0, 2.0, 9.0, 6.0, 2.0, 8.0, 4.0, 3.0, 7.0, 7.0, 7.0, 2.0, 8.0, 1.0, 6.0, 10.0, 4.0, 2.0, 3.0, 8.0, 6.0, 2.0, 6.0, 3.0, 7.0, 10.0, 10.0, 2.0, 7.0, 1.0, 5.0, 9.0, 2.0, 7.0, 3.0, 8.0, 7.0, 3.0, 10.0, 1.0, 9.0, 6.0, 4.0, 6.0, 4.0, 3.0, 10.0, 5.0, 9.0, 6.0, 4.0, 4.0, 8.0, 6.0, 6.0, 1.0, 1.0, 6.0, 2.0, 10.0, 3.0, 3.0, 3.0, 5.0, 1.0, 5.0, 8.0, 6.0, 1.0, 2.0, 9.0, 3.0, 8.0, 9.0, 4.0, 9.0, 5.0, 7.0, 6.0, 2.0, 1.0, 7.0, 1.0, 1.0, 8.0, 6.0, 4.0, 4.0, 2.0, 10.0, 8.0, 1.0, 5.0]
global b_y = 10
global p = [0.98, 0.681, 0.496, 0.811, 0.72, 0.797, 0.366, 0.461, 0.154, 0.466, 0.953, 0.858, 0.298, 0.874, 0.1, 0.771, 0.214, 0.071, 0.16, 0.602, 0.187, 0.429, 0.901, 0.826, 0.761, 0.825, 0.726, 0.65, 0.291, 0.239, 0.08, 0.98, 0.072, 0.799, 0.993, 0.74, 0.772, 0.637, 0.497, 0.578, 0.64, 0.085, 0.892, 0.366, 0.391, 0.008, 0.15, 0.232, 0.787, 0.158, 0.852, 0.857, 0.138, 0.065, 0.655, 0.72, 0.725, 0.235, 0.049, 0.996, 0.652, 0.791, 0.727, 0.842, 0.1, 0.991, 0.03, 0.091, 0.458, 0.294, 0.985, 0.794, 0.912, 0.347, 0.445, 0.465, 0.1, 0.709, 0.228, 0.914, 0.505, 0.514, 0.173, 0.549, 0.258, 0.09, 0.052, 0.977, 0.778, 0.39, 0.561, 0.623, 0.786, 0.035, 0.675, 0.933, 0.365, 0.464, 0.175, 0.13, 0.368, 0.196, 0.646, 0.097, 0.637, 0.947, 0.301, 0.433, 0.086, 0.302, 0.61, 0.794, 0.388, 0.352, 0.935, 0.662, 0.407, 0.207, 0.536, 0.534, 0.199, 0.139, 0.301, 0.402, 0.344, 0.647, 0.871, 0.221, 0.228, 0.903, 0.836, 0.781, 0.278, 0.627, 0.578, 0.959, 0.159, 0.689, 0.378, 0.623, 0.328, 0.926, 0.478, 0.195, 0.917, 0.334, 0.088, 0.574, 0.731, 0.749, 0.11, 0.251, 0.928, 0.843, 0.29, 0.71, 0.321, 0.014, 0.085, 0.841, 0.447, 0.174, 0.905, 0.164, 0.694, 0.819, 0.734, 0.032, 0.273, 0.787, 0.925, 0.626, 0.569, 0.468, 0.823, 0.651, 0.829, 0.985, 0.663, 0.434, 0.504, 0.465, 0.655, 0.926, 0.083, 0.468, 0.406, 0.793, 0.332, 0.438, 0.623, 0.185, 0.178, 0.026, 0.578, 0.287, 0.815, 0.792, 0.134, 0.103, 0.762, 0.674, 0.054, 0.92, 0.954, 0.161, 0.776, 0.833, 0.81, 0.13, 0.017, 0.248, 0.155, 0.125, 0.498, 0.097, 0.548, 0.654, 0.676, 0.999, 0.743, 0.274, 0.993, 0.243, 0.034, 0.064, 0.882, 0.373, 0.081, 0.763, 0.409, 0.488, 0.21, 0.829, 0.747, 0.101, 0.801, 0.376, 0.446, 0.175, 0.315, 0.251, 0.382, 0.834, 0.612, 0.98, 0.642, 0.603, 0.258, 0.489, 0.111, 0.184, 0.218, 0.292, 0.742, 0.81, 0.148, 0.218, 0.334, 0.623, 0.08, 0.097, 0.594, 0.522, 0.061, 0.927, 0.206, 0.199, 0.528, 0.806, 0.887, 0.122, 0.202, 0.863, 0.888, 0.394, 0.172, 0.793, 0.439, 0.214, 0.39, 0.867, 0.18, 0.369, 0.899, 0.232, 0.151, 0.375, 0.871, 0.493, 0.221, 0.956, 0.625, 0.644, 0.557, 0.021, 0.959, 0.014, 0.21, 0.851, 0.574, 0.9, 0.354, 0.865, 0.21, 0.82, 0.604, 0.594, 0.009, 0.382, 0.006, 0.965, 0.646, 0.334, 0.299, 0.028, 0.326, 0.146, 0.724, 0.192, 0.724, 0.101, 0.148]
global q = [0.981, 0.8, 0.527, 0.817, 0.823, 0.905, 0.377, 0.63, 0.9, 0.675, 0.959, 0.871, 0.823, 0.914, 0.111, 0.961, 0.649, 0.672, 0.838, 0.782, 0.712, 0.568, 0.989, 0.909, 0.965, 0.832, 0.767, 0.914, 0.636, 0.831, 0.867, 0.994, 0.23, 0.981, 0.997, 0.864, 0.858, 0.677, 0.876, 0.831, 0.691, 0.912, 0.912, 0.707, 0.813, 0.921, 0.214, 0.512, 0.813, 0.718, 0.993, 0.886, 0.584, 0.621, 0.9, 0.989, 0.868, 0.553, 0.431, 0.998, 0.755, 0.978, 0.733, 0.91, 0.383, 0.998, 0.382, 0.192, 0.803, 0.485, 0.996, 0.889, 0.963, 0.869, 0.976, 0.841, 0.847, 0.875, 0.904, 0.995, 0.529, 0.617, 0.963, 0.837, 0.376, 0.719, 0.591, 0.981, 0.996, 0.973, 0.88, 0.948, 0.98, 0.973, 0.969, 0.962, 0.74, 0.473, 0.732, 0.259, 0.722, 0.733, 0.807, 0.315, 0.969, 0.993, 0.793, 0.898, 0.997, 0.693, 0.94, 0.951, 0.6, 0.862, 0.956, 0.861, 0.471, 0.286, 0.836, 0.772, 0.934, 0.545, 0.661, 0.621, 0.989, 0.872, 0.946, 0.595, 0.48, 0.98, 0.853, 0.794, 0.605, 0.809, 0.934, 0.962, 0.898, 0.899, 0.77, 0.82, 0.416, 0.988, 0.847, 0.253, 0.993, 0.389, 0.698, 0.932, 0.943, 0.959, 0.441, 0.783, 0.981, 0.868, 0.558, 0.754, 0.755, 0.193, 0.321, 0.919, 0.992, 0.782, 0.955, 0.265, 0.896, 0.919, 0.854, 0.179, 0.357, 0.807, 0.992, 0.928, 0.91, 0.51, 0.89, 0.678, 0.848, 0.994, 0.72, 0.829, 0.709, 0.471, 0.84, 0.982, 0.391, 0.739, 0.993, 0.926, 0.903, 0.871, 0.669, 0.527, 0.566, 0.294, 0.77, 0.678, 0.993, 0.81, 0.604, 0.698, 0.779, 0.93, 0.242, 0.954, 0.964, 0.992, 0.94, 0.858, 0.976, 0.609, 0.035, 0.776, 0.592, 0.372, 0.665, 0.311, 0.71, 0.983, 0.938, 0.999, 0.79, 0.827, 0.994, 0.288, 0.991, 0.37, 0.938, 0.783, 0.61, 0.915, 0.789, 0.773, 0.948, 0.955, 0.969, 0.408, 0.856, 0.469, 0.861, 0.46, 0.383, 0.405, 0.913, 0.889, 0.92, 0.997, 0.683, 0.998, 0.585, 0.54, 0.814, 0.505, 0.974, 0.73, 0.781, 0.973, 0.53, 0.815, 0.928, 0.726, 0.933, 0.893, 0.644, 0.874, 0.491, 0.965, 0.547, 0.671, 0.841, 0.858, 0.945, 0.402, 0.577, 0.935, 0.892, 0.832, 0.4, 0.945, 0.734, 0.352, 0.749, 0.928, 0.26, 0.69, 0.905, 0.66, 0.528, 0.965, 0.915, 0.832, 0.339, 0.963, 0.786, 0.884, 0.66, 0.856, 0.959, 0.416, 0.494, 0.889, 0.939, 0.907, 0.808, 0.932, 0.78, 0.944, 0.776, 0.84, 0.084, 0.884, 0.358, 0.977, 0.996, 0.99, 0.476, 0.686, 0.374, 0.309, 0.785, 0.642, 0.965, 0.914, 0.225]
global origin = 1
global destination = 60