global arcs = [1 2; 1 24; 1 25; 1 32; 1 46; 1 58; 2 5; 2 11; 2 13; 2 49; 2 56; 3 2; 3 17; 3 26; 3 47; 4 3; 4 16; 4 19; 4 26; 4 30; 4 46; 4 48; 4 49; 4 52; 5 7; 5 14; 5 15; 5 21; 5 36; 5 52; 5 53; 6 13; 6 16; 6 24; 6 33; 6 43; 6 45; 6 57; 6 60; 7 23; 7 30; 7 35; 7 36; 7 50; 7 54; 8 10; 8 15; 8 19; 8 35; 8 38; 8 43; 8 47; 9 31; 9 42; 9 59; 10 24; 10 31; 10 40; 10 46; 10 57; 11 12; 11 16; 11 36; 11 49; 12 3; 12 10; 12 17; 13 14; 13 17; 13 20; 13 30; 13 43; 13 47; 13 48; 14 39; 14 43; 14 45; 15 5; 15 14; 15 23; 15 34; 15 39; 15 46; 16 9; 16 10; 16 20; 17 3; 17 48; 17 57; 18 10; 18 12; 18 19; 18 28; 18 33; 18 41; 18 47; 18 50; 18 52; 18 59; 19 9; 19 17; 19 31; 19 46; 19 50; 19 56; 20 3; 20 7; 20 17; 20 18; 20 32; 20 47; 20 59; 21 5; 21 17; 21 20; 21 25; 21 60; 22 6; 22 32; 22 41; 22 55; 23 3; 23 4; 23 48; 24 4; 24 14; 24 22; 24 25; 24 32; 24 36; 24 47; 24 59; 25 32; 25 47; 25 52; 25 56; 26 3; 26 9; 26 13; 26 14; 26 15; 26 19; 26 21; 26 36; 26 39; 26 53; 26 58; 27 6; 27 44; 27 45; 27 57; 28 9; 28 22; 28 39; 28 45; 28 54; 28 60; 29 3; 29 26; 29 36; 29 46; 29 56; 30 3; 30 4; 30 10; 30 18; 30 42; 30 50; 30 60; 31 10; 31 14; 31 21; 31 36; 31 41; 31 42; 31 45; 31 46; 31 49; 32 18; 32 23; 32 26; 32 34; 32 40; 32 43; 33 15; 33 29; 33 40; 33 43; 33 47; 33 50; 34 12; 34 50; 34 58; 35 14; 35 30; 35 34; 35 40; 35 50; 35 53; 35 58; 36 4; 36 12; 36 45; 36 60; 37 2; 37 5; 37 34; 37 53; 38 17; 38 20; 38 30; 38 33; 38 47; 38 52; 39 3; 39 5; 39 18; 39 21; 39 24; 39 28; 39 43; 39 52; 40 9; 40 25; 40 28; 40 38; 40 41; 40 45; 40 46; 40 51; 41 9; 41 11; 41 26; 41 28; 42 19; 42 21; 42 35; 42 37; 42 51; 42 52; 42 57; 42 58; 43 8; 43 20; 43 26; 43 34; 43 48; 44 43; 44 46; 45 13; 45 14; 45 15; 45 20; 45 21; 45 25; 45 26; 46 4; 46 20; 46 23; 46 31; 46 43; 46 52; 46 58; 46 60; 47 16; 47 19; 47 51; 47 52; 47 53; 47 54; 48 6; 48 9; 48 57; 48 60; 49 16; 49 22; 49 34; 49 53; 50 6; 50 8; 50 12; 50 13; 50 20; 50 29; 50 39; 50 41; 50 45; 50 48; 50 58; 51 6; 51 20; 51 43; 51 48; 51 53; 51 55; 51 58; 52 32; 52 40; 52 56; 53 4; 53 22; 53 32; 53 60; 54 3; 54 15; 54 17; 54 25; 55 5; 55 22; 55 37; 55 51; 56 3; 56 4; 56 9; 56 12; 56 14; 56 21; 56 23; 56 59; 57 7; 57 9; 57 22; 57 27; 57 32; 57 35; 57 47; 57 60; 58 5; 58 22; 58 23; 58 55; 59 23; 59 27; 59 48]
global d_x = [5.0, 5.0, 7.0, 5.0, 2.0, 9.0, 4.0, 3.0, 7.0, 4.0, 5.0, 6.0, 8.0, 5.0, 2.0, 2.0, 8.0, 7.0, 10.0, 8.0, 3.0, 2.0, 1.0, 1.0, 9.0, 6.0, 1.0, 6.0, 8.0, 2.0, 1.0, 3.0, 10.0, 10.0, 2.0, 1.0, 2.0, 10.0, 1.0, 5.0, 7.0, 4.0, 1.0, 8.0, 10.0, 8.0, 3.0, 3.0, 10.0, 9.0, 6.0, 5.0, 5.0, 2.0, 4.0, 6.0, 10.0, 4.0, 8.0, 6.0, 8.0, 9.0, 5.0, 7.0, 7.0, 6.0, 5.0, 6.0, 9.0, 3.0, 4.0, 2.0, 7.0, 3.0, 8.0, 10.0, 2.0, 3.0, 6.0, 1.0, 9.0, 7.0, 4.0, 10.0, 3.0, 9.0, 2.0, 8.0, 3.0, 7.0, 6.0, 10.0, 4.0, 5.0, 1.0, 2.0, 9.0, 9.0, 2.0, 9.0, 4.0, 10.0, 4.0, 5.0, 8.0, 9.0, 8.0, 5.0, 9.0, 3.0, 6.0, 1.0, 10.0, 7.0, 10.0, 5.0, 3.0, 4.0, 6.0, 8.0, 8.0, 3.0, 8.0, 5.0, 5.0, 8.0, 4.0, 2.0, 10.0, 1.0, 7.0, 9.0, 9.0, 9.0, 3.0, 4.0, 10.0, 5.0, 4.0, 6.0, 3.0, 8.0, 4.0, 2.0, 10.0, 6.0, 2.0, 2.0, 4.0, 6.0, 5.0, 3.0, 9.0, 9.0, 8.0, 4.0, 4.0, 4.0, 10.0, 9.0, 9.0, 1.0, 5.0, 8.0, 10.0, 9.0, 1.0, 9.0, 10.0, 7.0, 3.0, 5.0, 9.0, 4.0, 5.0, 1.0, 6.0, 10.0, 3.0, 8.0, 5.0, 7.0, 10.0, 2.0, 4.0, 8.0, 2.0, 4.0, 10.0, 4.0, 8.0, 2.0, 8.0, 6.0, 8.0, 3.0, 5.0, 7.0, 6.0, 8.0, 3.0, 8.0, 1.0, 2.0, 2.0, 8.0, 2.0, 5.0, 2.0, 6.0, 9.0, 9.0, 8.0, 8.0, 5.0, 10.0, 8.0, 10.0, 4.0, 5.0, 9.0, 1.0, 6.0, 10.0, 6.0, 7.0, 3.0, 5.0, 2.0, 9.0, 2.0, 2.0, 9.0, 8.0, 9.0, 1.0, 6.0, 4.0, 3.0, 8.0, 4.0, 6.0, 2.0, 6.0, 3.0, 1.0, 3.0, 8.0, 4.0, 3.0, 5.0, 3.0, 10.0, 1.0, 3.0, 3.0, 6.0, 7.0, 7.0, 1.0, 9.0, 2.0, 8.0, 2.0, 6.0, 8.0, 5.0, 5.0, 9.0, 7.0, 6.0, 2.0, 5.0, 5.0, 5.0, 4.0, 9.0, 4.0, 10.0, 7.0, 3.0, 9.0, 1.0, 7.0, 8.0, 4.0, 6.0, 1.0, 10.0, 2.0, 7.0, 9.0, 9.0, 9.0, 2.0, 1.0, 9.0, 2.0, 5.0, 4.0, 4.0, 10.0, 3.0, 4.0, 6.0, 10.0, 2.0, 6.0, 3.0, 2.0, 4.0, 5.0, 5.0, 8.0, 6.0, 1.0, 1.0, 2.0, 3.0, 6.0, 2.0, 6.0, 6.0, 8.0, 1.0, 4.0, 4.0, 7.0, 9.0, 10.0, 7.0, 9.0, 6.0, 7.0]
global b_x = 5
global d_y = [6.0, 7.0, 5.0, 10.0, 7.0, 6.0, 5.0, 2.0, 3.0, 9.0, 9.0, 6.0, 6.0, 6.0, 9.0, 10.0, 10.0, 2.0, 4.0, 4.0, 4.0, 8.0, 7.0, 10.0, 5.0, 6.0, 7.0, 1.0, 9.0, 4.0, 3.0, 1.0, 6.0, 2.0, 2.0, 7.0, 6.0, 2.0, 4.0, 6.0, 1.0, 2.0, 10.0, 10.0, 8.0, 6.0, 10.0, 2.0, 5.0, 7.0, 10.0, 1.0, 9.0, 4.0, 1.0, 2.0, 10.0, 4.0, 1.0, 3.0, 9.0, 5.0, 3.0, 4.0, 9.0, 10.0, 6.0, 5.0, 10.0, 8.0, 5.0, 2.0, 10.0, 1.0, 8.0, 1.0, 4.0, 2.0, 7.0, 4.0, 2.0, 5.0, 10.0, 9.0, 5.0, 9.0, 8.0, 1.0, 6.0, 4.0, 5.0, 1.0, 1.0, 2.0, 2.0, 6.0, 9.0, 3.0, 2.0, 3.0, 1.0, 1.0, 7.0, 1.0, 10.0, 5.0, 3.0, 7.0, 4.0, 7.0, 5.0, 5.0, 6.0, 4.0, 7.0, 7.0, 3.0, 2.0, 10.0, 5.0, 3.0, 9.0, 6.0, 2.0, 4.0, 7.0, 2.0, 8.0, 7.0, 8.0, 6.0, 4.0, 4.0, 2.0, 10.0, 8.0, 8.0, 1.0, 7.0, 7.0, 1.0, 5.0, 4.0, 7.0, 10.0, 8.0, 5.0, 3.0, 6.0, 8.0, 8.0, 1.0, 8.0, 5.0, 1.0, 5.0, 9.0, 4.0, 4.0, 10.0, 8.0, 9.0, 9.0, 1.0, 9.0, 2.0, 5.0, 5.0, 8.0, 3.0, 8.0, 10.0, 4.0, 1.0, 7.0, 8.0, 7.0, 1.0, 4.0, 7.0, 9.0, 4.0, 2.0, 6.0, 5.0, 7.0, 7.0, 6.0, 6.0, 2.0, 9.0, 8.0, 2.0, 4.0, 6.0, 2.0, 8.0, 4.0, 4.0, 6.0, 2.0, 7.0, 8.0, 6.0, 1.0, 8.0, 3.0, 8.0, 10.0, 1.0, 10.0, 7.0, 7.0, 4.0, 3.0, 5.0, 7.0, 9.0, 5.0, 10.0, 9.0, 7.0, 3.0, 5.0, 1.0, 1.0, 7.0, 8.0, 1.0, 1.0, 4.0, 8.0, 2.0, 5.0, 1.0, 2.0, 9.0, 3.0, 4.0, 5.0, 2.0, 8.0, 8.0, 4.0, 1.0, 6.0, 5.0, 1.0, 9.0, 7.0, 5.0, 2.0, 6.0, 6.0, 7.0, 4.0, 5.0, 3.0, 7.0, 4.0, 5.0, 10.0, 6.0, 8.0, 3.0, 9.0, 1.0, 1.0, 7.0, 7.0, 5.0, 1.0, 3.0, 10.0, 6.0, 5.0, 1.0, 10.0, 8.0, 4.0, 9.0, 10.0, 9.0, 10.0, 8.0, 2.0, 3.0, 2.0, 1.0, 2.0, 1.0, 3.0, 1.0, 8.0, 9.0, 9.0, 9.0, 4.0, 7.0, 3.0, 3.0, 8.0, 10.0, 9.0, 7.0, 5.0, 9.0, 10.0, 1.0, 7.0, 2.0, 2.0, 4.0, 9.0, 6.0, 8.0, 1.0, 4.0, 2.0, 1.0, 5.0, 4.0, 7.0, 3.0, 1.0, 6.0, 10.0, 1.0, 1.0, 6.0, 4.0, 9.0, 2.0, 7.0]
global b_y = 10
global p = [0.954, 0.561, 0.938, 0.348, 0.716, 0.049, 0.915, 0.458, 0.334, 0.866, 0.984, 0.844, 0.638, 0.183, 0.522, 0.986, 0.095, 0.696, 0.436, 0.062, 0.764, 0.603, 0.807, 0.071, 0.494, 0.658, 0.297, 0.67, 0.644, 0.099, 0.806, 0.123, 0.058, 0.916, 0.643, 0.376, 0.517, 0.924, 0.178, 0.117, 0.623, 0.091, 0.788, 0.496, 0.548, 0.026, 0.76, 0.182, 0.486, 0.922, 0.148, 0.59, 0.789, 0.125, 0.359, 0.81, 0.171, 0.932, 0.278, 0.567, 0.735, 0.565, 0.468, 0.625, 0.854, 0.149, 0.962, 0.134, 0.253, 0.353, 0.179, 0.331, 0.082, 0.478, 0.948, 0.84, 0.348, 0.253, 0.236, 0.473, 0.713, 0.114, 0.191, 0.43, 0.773, 0.374, 0.235, 0.433, 0.021, 0.789, 0.629, 0.97, 0.586, 0.353, 0.898, 0.835, 0.08, 0.136, 0.133, 0.033, 0.701, 0.288, 0.587, 0.587, 0.384, 0.155, 0.741, 0.766, 0.349, 0.417, 0.853, 0.881, 0.788, 0.156, 0.056, 0.332, 0.653, 0.509, 0.599, 0.425, 0.148, 0.917, 0.87, 0.232, 0.819, 0.143, 0.563, 0.73, 0.331, 0.902, 0.236, 0.893, 0.634, 0.796, 0.161, 0.277, 0.274, 0.609, 0.004, 0.622, 0.479, 0.726, 0.504, 0.828, 0.5, 0.042, 0.761, 0.311, 0.893, 0.245, 0.803, 0.813, 0.199, 0.185, 0.437, 0.978, 0.146, 0.557, 0.16, 0.694, 0.913, 0.801, 0.65, 0.409, 0.697, 0.933, 0.34, 0.946, 0.026, 0.395, 0.011, 0.86, 0.544, 0.772, 0.706, 0.614, 0.922, 0.863, 0.484, 0.969, 0.08, 0.59, 0.979, 0.005, 0.304, 0.71, 0.89, 0.616, 0.045, 0.273, 0.099, 0.675, 0.36, 0.404, 0.987, 0.332, 0.491, 0.069, 0.491, 0.351, 0.078, 0.599, 0.034, 0.439, 0.827, 0.638, 0.669, 0.893, 0.42, 0.454, 0.63, 0.673, 0.471, 0.413, 0.925, 0.149, 0.471, 0.744, 0.181, 0.309, 0.716, 0.315, 0.176, 0.879, 0.807, 0.386, 0.72, 0.748, 0.069, 0.935, 0.424, 0.548, 0.455, 0.814, 0.932, 0.137, 0.257, 0.777, 0.563, 0.621, 0.267, 0.982, 0.389, 0.777, 0.435, 0.433, 0.984, 0.675, 0.009, 0.949, 0.539, 0.176, 0.225, 0.099, 0.138, 0.809, 0.002, 0.349, 0.509, 0.277, 0.604, 0.89, 0.881, 0.903, 0.154, 0.954, 0.989, 0.165, 0.512, 0.132, 0.054, 0.855, 0.713, 0.316, 0.615, 0.733, 0.021, 0.231, 0.79, 0.75, 0.487, 0.114, 0.073, 0.301, 0.537, 0.401, 0.664, 0.933, 0.347, 0.757, 0.66, 0.192, 0.237, 0.226, 0.195, 0.044, 0.694, 0.322, 0.407, 0.424, 0.55, 0.098, 0.861, 0.213, 0.011, 0.423, 0.257, 0.17, 0.331, 0.098, 0.871, 0.596, 0.863, 0.321, 0.754, 0.509, 0.642, 0.404, 0.347, 0.873, 0.346, 0.639, 0.048, 0.366, 0.135, 0.073, 0.266, 0.159, 0.76, 0.033, 0.002, 0.563, 0.838, 0.313]
global q = [0.965, 0.589, 0.947, 0.648, 0.783, 0.669, 0.932, 0.796, 0.977, 0.975, 0.989, 0.962, 0.809, 0.736, 0.879, 0.999, 0.212, 0.701, 0.83, 0.602, 0.775, 0.777, 0.918, 0.575, 0.86, 0.881, 0.789, 0.896, 0.996, 0.119, 0.91, 0.508, 0.441, 0.988, 0.833, 0.45, 0.858, 0.939, 0.328, 0.476, 0.848, 0.684, 0.869, 0.508, 0.608, 0.814, 0.898, 0.798, 0.916, 0.974, 0.797, 0.711, 0.858, 0.204, 0.798, 0.827, 0.288, 0.988, 0.346, 0.963, 0.758, 0.631, 0.945, 0.907, 0.992, 0.167, 0.979, 0.583, 0.821, 0.866, 0.744, 0.987, 0.501, 0.68, 0.962, 0.982, 0.617, 0.629, 0.846, 0.801, 0.965, 0.888, 0.845, 0.736, 0.775, 0.733, 0.885, 0.527, 0.347, 0.87, 0.953, 0.979, 0.892, 0.743, 0.99, 0.964, 0.945, 0.639, 0.253, 0.211, 0.852, 0.522, 0.979, 0.892, 0.915, 0.924, 0.922, 0.889, 0.769, 0.642, 0.879, 0.972, 0.824, 0.757, 0.748, 0.391, 0.797, 0.619, 0.963, 0.715, 0.611, 0.968, 0.887, 0.26, 0.999, 0.479, 0.795, 0.877, 0.454, 0.97, 0.958, 0.976, 0.976, 0.916, 0.379, 0.356, 0.605, 0.668, 0.312, 0.952, 0.873, 0.821, 0.89, 0.835, 0.801, 0.953, 0.95, 0.648, 0.902, 0.395, 0.957, 0.863, 0.421, 0.81, 0.471, 0.992, 0.411, 0.877, 0.895, 0.967, 0.95, 0.93, 0.929, 0.582, 0.715, 0.987, 0.57, 0.958, 0.518, 0.515, 0.151, 0.961, 0.841, 0.811, 0.906, 0.99, 0.973, 0.97, 0.497, 0.97, 0.701, 0.652, 0.993, 0.757, 0.947, 0.854, 0.896, 0.821, 0.422, 0.899, 0.805, 0.946, 0.603, 0.854, 0.99, 0.699, 0.921, 0.783, 0.902, 0.851, 0.538, 0.729, 0.195, 0.569, 0.831, 0.645, 0.947, 0.982, 0.8, 0.498, 0.665, 0.692, 0.72, 0.796, 0.952, 0.281, 0.952, 0.823, 0.632, 0.815, 0.909, 0.826, 0.56, 0.967, 0.919, 0.587, 0.738, 0.95, 0.714, 0.962, 0.886, 0.889, 0.814, 0.967, 0.98, 0.931, 0.519, 0.991, 0.657, 0.73, 0.638, 0.996, 0.827, 0.927, 0.992, 0.51, 0.998, 0.931, 0.521, 0.982, 0.554, 0.332, 0.645, 0.868, 0.65, 0.862, 0.215, 0.884, 0.713, 0.362, 0.997, 0.939, 0.947, 0.97, 0.285, 0.959, 0.996, 0.956, 0.996, 0.598, 0.089, 0.997, 0.85, 0.767, 0.895, 0.939, 0.368, 0.368, 0.941, 0.878, 0.784, 0.664, 0.463, 0.438, 0.867, 0.458, 0.705, 0.971, 0.591, 0.905, 0.901, 0.54, 0.289, 0.742, 0.223, 0.923, 0.934, 0.72, 0.486, 0.695, 0.745, 0.297, 0.971, 0.517, 0.336, 0.602, 0.804, 0.658, 0.361, 0.501, 0.944, 0.65, 0.923, 0.361, 0.824, 0.919, 0.783, 0.949, 0.703, 0.976, 0.385, 0.87, 0.129, 0.997, 0.156, 0.213, 0.443, 0.472, 0.798, 0.499, 0.189, 0.653, 0.997, 0.719]
global origin = 1
global destination = 60