global arcs = [1 31; 1 35; 1 37; 1 40; 1 46; 2 11; 2 13; 2 32; 2 35; 2 52; 3 2; 3 5; 3 25; 3 31; 3 47; 3 49; 3 52; 3 55; 3 56; 4 10; 4 19; 4 24; 4 30; 4 37; 4 40; 4 46; 5 16; 5 18; 5 20; 5 28; 5 37; 5 43; 5 50; 5 53; 5 57; 6 12; 6 14; 6 25; 6 29; 6 40; 6 55; 6 58; 7 15; 7 25; 7 32; 7 42; 7 47; 8 45; 8 57; 9 3; 9 6; 9 19; 9 27; 9 42; 10 12; 10 23; 10 28; 10 29; 10 42; 10 46; 11 19; 11 34; 11 39; 12 6; 12 46; 12 56; 13 2; 13 36; 13 45; 13 50; 14 4; 14 8; 14 36; 14 42; 15 8; 15 55; 15 57; 15 60; 16 27; 16 31; 16 34; 16 35; 16 47; 17 3; 17 9; 17 11; 17 14; 17 15; 17 27; 17 49; 18 21; 18 32; 18 57; 19 7; 19 14; 19 16; 19 24; 19 48; 19 53; 19 59; 20 23; 20 26; 20 36; 20 51; 20 60; 21 2; 21 4; 21 6; 21 12; 21 20; 21 35; 21 38; 21 59; 22 2; 22 45; 22 54; 22 58; 22 59; 23 22; 23 33; 23 38; 23 51; 23 59; 24 14; 24 20; 24 25; 24 34; 24 43; 24 46; 24 54; 25 3; 25 11; 25 21; 25 26; 25 34; 25 40; 25 52; 25 55; 26 4; 26 5; 26 16; 26 24; 26 34; 27 5; 27 15; 27 21; 27 32; 27 34; 27 41; 27 48; 28 46; 29 26; 29 44; 30 7; 30 10; 30 17; 30 19; 30 49; 30 54; 30 58; 30 60; 31 3; 31 8; 31 20; 31 33; 31 39; 31 41; 31 56; 32 14; 32 39; 33 48; 33 53; 34 5; 34 27; 34 36; 34 41; 35 16; 35 27; 35 39; 36 5; 36 8; 36 17; 36 21; 36 44; 36 51; 36 60; 37 18; 37 31; 37 34; 37 55; 37 56; 38 3; 38 21; 38 43; 38 55; 39 11; 39 24; 39 27; 39 30; 39 32; 39 43; 39 47; 39 48; 39 58; 40 4; 40 17; 40 18; 40 27; 40 45; 40 48; 40 51; 40 52; 40 54; 40 55; 40 57; 41 4; 41 8; 41 9; 41 48; 41 59; 41 60; 42 4; 42 15; 42 22; 42 38; 43 15; 43 17; 43 21; 43 35; 43 44; 43 49; 43 54; 44 2; 44 22; 44 23; 44 39; 44 41; 44 47; 44 50; 44 54; 45 2; 45 22; 45 25; 45 30; 45 41; 46 21; 46 24; 46 30; 46 53; 47 14; 47 20; 47 24; 47 42; 47 58; 48 6; 48 13; 48 17; 48 26; 48 49; 49 6; 49 17; 49 31; 49 41; 49 44; 49 54; 50 5; 50 36; 50 39; 50 49; 50 57; 51 12; 51 16; 51 17; 51 18; 51 24; 51 38; 51 59; 52 10; 52 14; 52 30; 52 38; 52 47; 52 50; 52 54; 52 58; 53 4; 53 5; 53 25; 53 29; 53 34; 53 39; 53 41; 54 23; 54 28; 54 35; 54 45; 54 50; 55 2; 55 6; 55 16; 55 17; 55 24; 55 37; 56 18; 56 19; 56 31; 56 35; 56 46; 56 48; 56 50; 56 51; 56 59; 57 7; 57 20; 57 34; 57 36; 57 52; 57 58; 58 14; 58 30; 58 56; 59 12; 59 17; 59 31; 59 44; 59 47; 59 60]
global d_x = [5.0, 9.0, 1.0, 6.0, 2.0, 9.0, 6.0, 10.0, 4.0, 1.0, 10.0, 10.0, 5.0, 9.0, 10.0, 6.0, 10.0, 8.0, 1.0, 10.0, 9.0, 3.0, 2.0, 4.0, 5.0, 4.0, 6.0, 3.0, 9.0, 4.0, 6.0, 6.0, 6.0, 10.0, 10.0, 10.0, 6.0, 6.0, 7.0, 1.0, 9.0, 10.0, 10.0, 1.0, 10.0, 5.0, 5.0, 3.0, 4.0, 9.0, 9.0, 7.0, 2.0, 9.0, 9.0, 8.0, 3.0, 9.0, 6.0, 1.0, 5.0, 4.0, 8.0, 1.0, 7.0, 3.0, 1.0, 2.0, 8.0, 2.0, 8.0, 7.0, 3.0, 8.0, 10.0, 4.0, 6.0, 7.0, 4.0, 9.0, 5.0, 1.0, 8.0, 7.0, 4.0, 1.0, 2.0, 7.0, 6.0, 8.0, 3.0, 6.0, 3.0, 8.0, 2.0, 10.0, 9.0, 2.0, 7.0, 8.0, 3.0, 1.0, 5.0, 3.0, 10.0, 5.0, 8.0, 6.0, 9.0, 1.0, 9.0, 4.0, 5.0, 4.0, 2.0, 8.0, 1.0, 9.0, 2.0, 1.0, 2.0, 2.0, 9.0, 9.0, 1.0, 4.0, 4.0, 10.0, 10.0, 2.0, 8.0, 6.0, 9.0, 5.0, 6.0, 10.0, 10.0, 8.0, 2.0, 1.0, 8.0, 5.0, 2.0, 1.0, 2.0, 6.0, 9.0, 4.0, 10.0, 10.0, 5.0, 9.0, 6.0, 3.0, 3.0, 3.0, 4.0, 2.0, 5.0, 8.0, 5.0, 10.0, 1.0, 6.0, 4.0, 8.0, 8.0, 8.0, 8.0, 4.0, 6.0, 6.0, 9.0, 5.0, 8.0, 6.0, 3.0, 10.0, 4.0, 8.0, 5.0, 1.0, 9.0, 5.0, 8.0, 1.0, 5.0, 7.0, 4.0, 10.0, 2.0, 4.0, 8.0, 4.0, 9.0, 6.0, 9.0, 2.0, 3.0, 4.0, 10.0, 1.0, 3.0, 4.0, 2.0, 5.0, 9.0, 8.0, 7.0, 8.0, 4.0, 10.0, 5.0, 9.0, 10.0, 10.0, 7.0, 1.0, 9.0, 3.0, 8.0, 4.0, 6.0, 8.0, 9.0, 7.0, 3.0, 4.0, 4.0, 2.0, 2.0, 1.0, 4.0, 9.0, 9.0, 1.0, 4.0, 1.0, 2.0, 5.0, 6.0, 4.0, 2.0, 7.0, 7.0, 7.0, 3.0, 4.0, 5.0, 5.0, 2.0, 8.0, 8.0, 4.0, 6.0, 7.0, 8.0, 10.0, 5.0, 5.0, 3.0, 3.0, 7.0, 6.0, 3.0, 1.0, 5.0, 7.0, 1.0, 5.0, 9.0, 9.0, 4.0, 7.0, 4.0, 5.0, 9.0, 1.0, 3.0, 7.0, 1.0, 9.0, 1.0, 7.0, 3.0, 9.0, 1.0, 10.0, 2.0, 5.0, 1.0, 1.0, 5.0, 2.0, 9.0, 4.0, 4.0, 3.0, 3.0, 7.0, 10.0, 10.0, 8.0, 1.0, 5.0, 5.0, 8.0, 1.0, 2.0, 3.0, 5.0, 1.0, 9.0, 3.0, 6.0, 9.0, 7.0, 9.0, 2.0, 9.0, 10.0, 9.0, 4.0, 2.0, 3.0, 3.0, 8.0]
global b_x = 5
global d_y = [10.0, 5.0, 2.0, 7.0, 9.0, 4.0, 6.0, 10.0, 4.0, 10.0, 10.0, 1.0, 5.0, 8.0, 2.0, 2.0, 1.0, 4.0, 2.0, 8.0, 10.0, 10.0, 5.0, 10.0, 7.0, 3.0, 5.0, 10.0, 5.0, 9.0, 1.0, 5.0, 6.0, 7.0, 3.0, 9.0, 7.0, 8.0, 6.0, 2.0, 3.0, 7.0, 1.0, 4.0, 10.0, 8.0, 3.0, 2.0, 10.0, 10.0, 10.0, 7.0, 3.0, 2.0, 5.0, 1.0, 2.0, 5.0, 4.0, 8.0, 1.0, 5.0, 1.0, 7.0, 9.0, 10.0, 6.0, 5.0, 3.0, 9.0, 5.0, 3.0, 9.0, 2.0, 3.0, 1.0, 9.0, 2.0, 7.0, 8.0, 4.0, 7.0, 9.0, 8.0, 7.0, 3.0, 2.0, 7.0, 10.0, 8.0, 10.0, 5.0, 8.0, 3.0, 9.0, 2.0, 9.0, 7.0, 3.0, 5.0, 8.0, 6.0, 8.0, 4.0, 10.0, 8.0, 10.0, 1.0, 4.0, 7.0, 4.0, 6.0, 7.0, 9.0, 4.0, 2.0, 4.0, 2.0, 1.0, 9.0, 2.0, 8.0, 6.0, 10.0, 7.0, 5.0, 8.0, 7.0, 1.0, 5.0, 7.0, 10.0, 8.0, 6.0, 6.0, 4.0, 2.0, 10.0, 8.0, 10.0, 2.0, 5.0, 6.0, 4.0, 8.0, 8.0, 1.0, 9.0, 7.0, 9.0, 5.0, 7.0, 2.0, 8.0, 8.0, 9.0, 2.0, 8.0, 5.0, 4.0, 9.0, 5.0, 1.0, 8.0, 4.0, 2.0, 9.0, 3.0, 2.0, 7.0, 3.0, 3.0, 9.0, 2.0, 7.0, 10.0, 1.0, 1.0, 7.0, 9.0, 6.0, 1.0, 8.0, 4.0, 7.0, 10.0, 10.0, 6.0, 6.0, 9.0, 5.0, 7.0, 8.0, 9.0, 3.0, 9.0, 7.0, 6.0, 6.0, 6.0, 10.0, 7.0, 1.0, 7.0, 1.0, 8.0, 4.0, 2.0, 3.0, 3.0, 9.0, 3.0, 3.0, 10.0, 1.0, 2.0, 8.0, 2.0, 7.0, 6.0, 10.0, 4.0, 9.0, 9.0, 10.0, 2.0, 7.0, 6.0, 1.0, 8.0, 9.0, 6.0, 2.0, 6.0, 3.0, 1.0, 7.0, 4.0, 1.0, 4.0, 1.0, 4.0, 6.0, 5.0, 5.0, 4.0, 4.0, 10.0, 1.0, 1.0, 6.0, 4.0, 7.0, 7.0, 4.0, 4.0, 1.0, 7.0, 10.0, 3.0, 5.0, 1.0, 7.0, 2.0, 4.0, 10.0, 9.0, 8.0, 5.0, 8.0, 7.0, 2.0, 8.0, 7.0, 9.0, 10.0, 2.0, 3.0, 8.0, 8.0, 10.0, 2.0, 9.0, 2.0, 10.0, 5.0, 9.0, 2.0, 9.0, 1.0, 3.0, 8.0, 6.0, 9.0, 1.0, 4.0, 7.0, 3.0, 10.0, 7.0, 6.0, 8.0, 5.0, 2.0, 1.0, 3.0, 8.0, 4.0, 7.0, 8.0, 7.0, 5.0, 4.0, 2.0, 5.0, 6.0, 2.0, 4.0, 4.0, 7.0, 5.0, 2.0, 9.0, 6.0, 6.0, 10.0, 10.0]
global b_y = 10
global p = [0.269, 0.697, 0.258, 0.201, 0.211, 0.393, 0.282, 0.662, 0.269, 0.232, 0.34, 0.11, 0.36, 0.568, 0.311, 0.656, 0.945, 0.674, 0.073, 0.633, 0.755, 0.992, 0.277, 0.245, 0.318, 0.882, 0.81, 0.621, 0.241, 0.387, 0.947, 0.913, 0.66, 0.024, 0.144, 0.22, 0.325, 0.849, 0.194, 0.797, 0.397, 0.031, 0.307, 0.572, 0.719, 0.041, 0.209, 0.042, 0.237, 0.198, 0.052, 0.125, 0.338, 0.671, 0.651, 0.851, 0.837, 0.314, 0.581, 0.985, 0.936, 0.459, 0.784, 0.192, 0.512, 0.927, 0.297, 0.194, 0.009, 0.28, 0.428, 0.518, 0.461, 0.812, 0.529, 0.124, 0.684, 0.038, 0.07, 0.857, 0.696, 0.91, 0.836, 0.355, 0.26, 0.869, 0.246, 0.101, 0.685, 0.669, 0.486, 0.348, 0.461, 0.175, 0.029, 0.825, 0.638, 0.869, 0.122, 0.759, 0.383, 0.407, 0.403, 0.496, 0.474, 0.156, 0.912, 0.064, 0.713, 0.244, 0.8, 0.513, 0.323, 0.114, 0.615, 0.32, 0.028, 0.385, 0.285, 0.8, 0.103, 0.325, 0.504, 0.219, 0.299, 0.922, 0.325, 0.446, 0.896, 0.588, 0.322, 0.019, 0.371, 0.851, 0.783, 0.405, 0.943, 0.741, 0.004, 0.782, 0.867, 0.404, 0.023, 0.705, 0.483, 0.864, 0.715, 0.877, 0.449, 0.252, 0.072, 0.436, 0.431, 0.527, 0.266, 0.836, 0.929, 0.65, 0.532, 0.938, 0.357, 0.75, 0.867, 0.538, 0.857, 0.618, 0.014, 0.266, 0.253, 0.346, 0.961, 0.779, 0.892, 0.735, 0.523, 0.353, 0.86, 0.807, 0.326, 0.466, 0.429, 0.125, 0.258, 0.849, 0.287, 0.645, 0.869, 0.064, 0.364, 0.397, 0.691, 0.416, 0.464, 0.824, 0.152, 0.844, 0.293, 0.938, 0.821, 0.715, 0.273, 0.572, 0.639, 0.603, 0.344, 0.781, 0.124, 0.561, 0.597, 0.904, 0.761, 0.12, 0.82, 0.344, 0.043, 0.894, 0.739, 0.595, 0.352, 0.708, 0.276, 0.514, 0.008, 0.761, 0.545, 0.563, 0.672, 0.518, 0.429, 0.721, 0.298, 0.625, 0.393, 0.727, 0.139, 0.933, 0.696, 0.727, 0.871, 0.635, 0.538, 0.775, 0.008, 0.563, 0.927, 0.179, 0.065, 0.421, 0.039, 0.934, 0.772, 0.005, 0.281, 0.945, 0.419, 0.436, 0.501, 0.68, 0.613, 0.901, 0.942, 0.028, 0.35, 0.233, 0.311, 0.558, 0.052, 0.524, 0.912, 0.019, 0.797, 0.281, 0.861, 0.774, 0.996, 0.98, 0.869, 0.721, 0.163, 0.737, 0.226, 0.336, 0.632, 0.733, 0.857, 0.679, 0.112, 0.46, 0.142, 0.781, 0.205, 0.676, 0.465, 0.076, 0.054, 0.514, 0.545, 0.126, 0.555, 0.144, 0.163, 0.083, 0.996, 0.641, 0.503, 0.714, 0.305, 0.888, 0.527, 0.029, 0.862, 0.384, 0.962, 0.182, 0.698, 0.946, 0.847, 0.396, 0.105, 0.423, 0.486, 0.283, 0.953, 0.269, 0.231, 0.531, 0.622]
global q = [0.734, 0.896, 0.57, 0.399, 0.752, 0.835, 0.322, 0.684, 0.892, 0.328, 0.769, 0.697, 0.996, 0.599, 0.49, 0.852, 0.945, 0.689, 0.755, 0.886, 0.912, 0.995, 0.324, 0.35, 0.815, 0.974, 0.947, 0.664, 0.427, 0.62, 0.959, 0.921, 0.826, 0.283, 0.91, 0.666, 0.519, 0.881, 0.754, 0.849, 0.96, 0.195, 0.679, 0.687, 0.932, 0.152, 0.827, 0.299, 0.519, 0.948, 0.803, 0.945, 0.863, 0.684, 0.996, 0.954, 0.968, 0.509, 0.79, 0.99, 0.971, 0.963, 0.95, 0.652, 0.688, 0.932, 0.94, 0.554, 0.766, 0.298, 0.84, 0.812, 0.877, 0.977, 0.674, 0.8, 0.941, 0.653, 0.246, 0.916, 0.984, 0.965, 0.892, 0.881, 0.75, 0.978, 0.647, 0.481, 0.945, 0.68, 0.786, 0.707, 0.756, 0.966, 0.047, 0.843, 0.955, 0.926, 0.226, 0.774, 0.916, 0.912, 0.856, 0.609, 0.926, 0.772, 0.971, 0.462, 0.882, 0.388, 0.997, 0.842, 0.725, 0.263, 0.787, 0.329, 0.167, 0.886, 0.929, 0.909, 0.43, 0.884, 0.73, 0.355, 0.526, 0.976, 0.368, 0.633, 0.922, 0.921, 0.77, 0.409, 0.489, 0.86, 0.939, 0.997, 0.961, 0.934, 0.725, 0.954, 0.945, 0.61, 0.643, 0.896, 0.606, 0.908, 0.99, 0.996, 0.559, 0.589, 0.512, 0.757, 0.498, 0.967, 0.922, 0.978, 0.948, 0.845, 0.835, 0.958, 0.358, 0.907, 0.94, 0.6, 0.924, 0.681, 0.251, 0.848, 0.457, 0.787, 0.961, 0.893, 0.914, 0.741, 0.739, 0.789, 0.992, 0.888, 0.809, 0.982, 0.72, 0.711, 0.545, 0.995, 0.716, 0.788, 0.913, 0.767, 0.811, 0.901, 0.918, 0.857, 0.788, 0.951, 0.277, 0.91, 0.754, 0.969, 0.893, 0.954, 0.775, 0.979, 0.741, 0.788, 0.52, 0.996, 0.946, 0.943, 0.881, 0.92, 0.868, 0.206, 0.874, 0.886, 0.831, 0.966, 0.861, 0.704, 0.42, 0.796, 0.558, 0.934, 0.046, 0.946, 0.744, 0.714, 0.793, 0.994, 0.889, 0.842, 0.36, 0.977, 0.97, 0.882, 0.989, 0.956, 0.853, 0.919, 0.876, 0.985, 0.829, 0.864, 0.384, 0.999, 0.973, 0.817, 0.701, 0.609, 0.164, 0.994, 0.99, 0.216, 0.586, 0.95, 0.719, 0.792, 0.69, 0.701, 0.622, 0.931, 0.962, 0.362, 0.48, 0.896, 0.951, 0.927, 0.096, 0.793, 0.99, 0.115, 0.895, 0.834, 0.878, 0.882, 0.998, 0.997, 0.913, 0.942, 0.692, 0.906, 0.898, 0.71, 0.826, 0.961, 0.88, 0.81, 0.583, 0.608, 0.581, 0.878, 0.833, 0.688, 0.763, 0.706, 0.855, 0.832, 0.555, 0.89, 0.9, 0.284, 0.372, 0.334, 0.997, 0.935, 0.742, 0.922, 0.985, 0.896, 0.752, 0.3, 0.952, 0.886, 0.974, 0.859, 0.991, 0.985, 0.976, 0.517, 0.458, 0.945, 0.983, 0.792, 0.993, 0.996, 0.598, 0.702, 0.967]
global origin = 1
global destination = 60