global arcs = [1 11; 1 22; 1 24; 1 34; 1 52; 1 57; 2 17; 2 33; 2 46; 3 2; 3 6; 3 8; 3 14; 3 20; 3 33; 3 42; 3 52; 4 10; 4 25; 4 27; 4 29; 4 41; 4 42; 4 51; 4 56; 5 3; 5 29; 5 37; 5 40; 5 46; 5 55; 6 27; 6 39; 6 41; 6 46; 6 53; 7 5; 7 35; 7 40; 7 54; 7 59; 8 11; 8 21; 8 27; 8 31; 8 39; 9 8; 9 23; 9 29; 9 30; 9 36; 9 48; 9 59; 10 47; 11 14; 11 33; 11 47; 11 48; 11 54; 12 10; 12 15; 12 16; 12 18; 12 23; 12 25; 12 32; 12 41; 13 17; 13 29; 13 30; 13 32; 13 43; 13 46; 13 51; 14 6; 14 7; 14 18; 14 28; 14 33; 14 47; 14 50; 14 54; 15 21; 15 39; 15 50; 15 56; 16 20; 16 32; 16 35; 16 40; 16 43; 16 54; 17 2; 17 3; 17 13; 17 15; 17 50; 17 51; 18 6; 18 11; 18 47; 18 51; 18 60; 19 3; 19 10; 19 11; 19 16; 19 39; 19 52; 20 49; 20 59; 21 9; 21 17; 21 36; 21 41; 21 42; 21 48; 21 51; 21 58; 22 6; 22 9; 23 12; 23 26; 23 32; 23 48; 23 50; 24 4; 24 7; 24 17; 24 25; 24 27; 24 32; 24 37; 24 44; 24 46; 25 19; 25 24; 25 26; 25 27; 25 32; 26 11; 26 14; 26 20; 26 50; 26 57; 27 12; 27 17; 27 19; 27 23; 27 43; 27 45; 27 54; 28 16; 28 18; 28 31; 29 21; 29 48; 29 58; 29 60; 30 4; 30 33; 30 35; 30 42; 30 44; 30 50; 31 5; 31 21; 31 29; 31 39; 31 58; 31 59; 32 10; 32 12; 32 21; 32 35; 32 44; 33 7; 33 19; 33 29; 33 32; 33 52; 33 54; 33 58; 34 2; 34 8; 34 9; 34 20; 34 31; 34 33; 34 36; 34 41; 34 43; 34 48; 35 12; 35 13; 35 21; 35 34; 35 40; 35 46; 36 15; 36 24; 36 33; 36 39; 36 49; 36 53; 37 30; 37 32; 37 45; 37 47; 38 2; 38 4; 38 35; 38 39; 38 45; 38 47; 38 48; 38 51; 38 54; 39 18; 39 22; 39 23; 39 29; 39 35; 39 41; 39 43; 39 48; 39 53; 40 19; 40 22; 40 25; 40 34; 41 4; 41 7; 41 28; 41 44; 41 54; 42 7; 42 27; 42 35; 43 41; 44 4; 44 8; 44 26; 44 28; 44 36; 44 46; 45 2; 45 6; 45 18; 45 33; 45 41; 45 58; 46 12; 46 36; 46 38; 47 6; 47 41; 47 45; 47 56; 48 3; 48 24; 48 31; 48 37; 48 40; 48 44; 48 45; 48 52; 49 14; 49 39; 49 42; 49 48; 49 50; 49 51; 50 5; 50 6; 50 18; 50 19; 50 31; 50 35; 50 43; 51 9; 51 19; 51 25; 51 46; 51 53; 51 58; 52 28; 52 55; 52 56; 53 7; 53 10; 53 38; 53 40; 53 45; 53 46; 54 7; 54 21; 54 22; 54 23; 54 27; 54 45; 55 30; 55 31; 56 19; 56 20; 56 44; 56 59; 57 10; 57 16; 57 28; 57 35; 57 49; 57 60; 58 53; 59 15; 59 16; 59 51; 59 56]
global d_x = [8.0, 7.0, 8.0, 2.0, 9.0, 7.0, 10.0, 4.0, 6.0, 10.0, 6.0, 6.0, 5.0, 6.0, 7.0, 6.0, 8.0, 4.0, 7.0, 1.0, 10.0, 3.0, 9.0, 1.0, 9.0, 7.0, 7.0, 2.0, 10.0, 3.0, 8.0, 6.0, 6.0, 5.0, 10.0, 5.0, 3.0, 10.0, 1.0, 7.0, 2.0, 5.0, 9.0, 2.0, 7.0, 10.0, 10.0, 8.0, 8.0, 8.0, 4.0, 7.0, 5.0, 5.0, 3.0, 9.0, 7.0, 1.0, 4.0, 9.0, 1.0, 3.0, 10.0, 3.0, 5.0, 10.0, 1.0, 10.0, 2.0, 8.0, 5.0, 1.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 8.0, 4.0, 7.0, 10.0, 3.0, 8.0, 4.0, 9.0, 5.0, 4.0, 4.0, 9.0, 8.0, 7.0, 1.0, 5.0, 1.0, 9.0, 2.0, 4.0, 5.0, 6.0, 2.0, 6.0, 10.0, 1.0, 9.0, 1.0, 8.0, 2.0, 10.0, 6.0, 7.0, 3.0, 5.0, 2.0, 9.0, 8.0, 3.0, 5.0, 6.0, 9.0, 10.0, 9.0, 2.0, 6.0, 1.0, 7.0, 4.0, 9.0, 2.0, 9.0, 4.0, 6.0, 3.0, 4.0, 6.0, 9.0, 3.0, 6.0, 10.0, 2.0, 7.0, 9.0, 6.0, 7.0, 7.0, 7.0, 10.0, 6.0, 10.0, 5.0, 2.0, 8.0, 10.0, 6.0, 6.0, 9.0, 4.0, 10.0, 7.0, 1.0, 3.0, 3.0, 10.0, 5.0, 8.0, 6.0, 8.0, 1.0, 1.0, 9.0, 10.0, 9.0, 6.0, 1.0, 10.0, 7.0, 1.0, 5.0, 7.0, 4.0, 3.0, 10.0, 3.0, 4.0, 3.0, 2.0, 2.0, 2.0, 10.0, 1.0, 4.0, 3.0, 10.0, 4.0, 9.0, 1.0, 5.0, 9.0, 9.0, 6.0, 6.0, 8.0, 6.0, 2.0, 2.0, 10.0, 3.0, 9.0, 5.0, 2.0, 4.0, 3.0, 2.0, 1.0, 9.0, 6.0, 5.0, 10.0, 9.0, 1.0, 6.0, 2.0, 5.0, 4.0, 9.0, 5.0, 8.0, 6.0, 6.0, 9.0, 7.0, 10.0, 10.0, 5.0, 6.0, 5.0, 9.0, 8.0, 4.0, 4.0, 5.0, 7.0, 9.0, 6.0, 4.0, 7.0, 8.0, 3.0, 4.0, 2.0, 2.0, 1.0, 8.0, 4.0, 6.0, 6.0, 1.0, 7.0, 6.0, 5.0, 5.0, 3.0, 5.0, 1.0, 8.0, 6.0, 3.0, 2.0, 6.0, 8.0, 7.0, 8.0, 3.0, 1.0, 5.0, 7.0, 1.0, 7.0, 10.0, 3.0, 10.0, 10.0, 6.0, 7.0, 8.0, 2.0, 2.0, 1.0, 8.0, 2.0, 9.0, 2.0, 5.0, 6.0, 3.0, 9.0, 7.0, 8.0, 4.0, 3.0, 10.0, 8.0, 6.0, 7.0, 7.0, 10.0, 9.0, 2.0, 3.0, 8.0, 3.0, 2.0, 2.0, 3.0, 8.0, 7.0, 4.0, 3.0]
global b_x = 5
global d_y = [10.0, 2.0, 2.0, 6.0, 2.0, 8.0, 5.0, 8.0, 2.0, 3.0, 3.0, 5.0, 4.0, 9.0, 8.0, 2.0, 10.0, 6.0, 2.0, 10.0, 5.0, 1.0, 7.0, 3.0, 2.0, 6.0, 10.0, 4.0, 9.0, 6.0, 3.0, 3.0, 5.0, 10.0, 2.0, 2.0, 7.0, 10.0, 10.0, 8.0, 10.0, 9.0, 8.0, 1.0, 2.0, 6.0, 1.0, 10.0, 3.0, 2.0, 7.0, 5.0, 1.0, 7.0, 1.0, 5.0, 10.0, 2.0, 2.0, 2.0, 9.0, 4.0, 5.0, 5.0, 9.0, 1.0, 4.0, 2.0, 3.0, 2.0, 8.0, 5.0, 8.0, 4.0, 7.0, 10.0, 1.0, 4.0, 7.0, 10.0, 5.0, 5.0, 2.0, 4.0, 3.0, 10.0, 3.0, 9.0, 10.0, 10.0, 1.0, 8.0, 8.0, 2.0, 6.0, 10.0, 6.0, 7.0, 6.0, 8.0, 6.0, 2.0, 9.0, 9.0, 4.0, 10.0, 2.0, 9.0, 9.0, 2.0, 3.0, 5.0, 4.0, 4.0, 1.0, 8.0, 10.0, 10.0, 8.0, 6.0, 9.0, 5.0, 8.0, 10.0, 4.0, 1.0, 4.0, 3.0, 9.0, 6.0, 7.0, 1.0, 2.0, 10.0, 10.0, 2.0, 3.0, 2.0, 2.0, 5.0, 5.0, 4.0, 1.0, 3.0, 1.0, 1.0, 4.0, 5.0, 2.0, 6.0, 7.0, 2.0, 7.0, 2.0, 3.0, 7.0, 5.0, 8.0, 1.0, 1.0, 10.0, 3.0, 9.0, 9.0, 5.0, 5.0, 1.0, 8.0, 6.0, 10.0, 9.0, 5.0, 4.0, 9.0, 5.0, 4.0, 8.0, 10.0, 8.0, 6.0, 8.0, 1.0, 7.0, 5.0, 7.0, 6.0, 6.0, 9.0, 1.0, 7.0, 6.0, 10.0, 2.0, 7.0, 7.0, 4.0, 9.0, 8.0, 4.0, 8.0, 9.0, 2.0, 10.0, 10.0, 9.0, 5.0, 9.0, 10.0, 6.0, 1.0, 8.0, 10.0, 1.0, 3.0, 3.0, 3.0, 10.0, 9.0, 4.0, 3.0, 8.0, 1.0, 5.0, 3.0, 1.0, 5.0, 6.0, 1.0, 9.0, 5.0, 3.0, 9.0, 6.0, 7.0, 6.0, 7.0, 8.0, 7.0, 9.0, 1.0, 2.0, 10.0, 4.0, 5.0, 6.0, 2.0, 6.0, 5.0, 10.0, 9.0, 8.0, 2.0, 5.0, 7.0, 1.0, 10.0, 8.0, 6.0, 7.0, 8.0, 3.0, 5.0, 1.0, 4.0, 3.0, 2.0, 4.0, 2.0, 4.0, 4.0, 7.0, 3.0, 2.0, 6.0, 10.0, 1.0, 10.0, 10.0, 5.0, 6.0, 5.0, 2.0, 2.0, 7.0, 2.0, 4.0, 5.0, 4.0, 9.0, 4.0, 4.0, 8.0, 4.0, 7.0, 8.0, 8.0, 10.0, 2.0, 6.0, 4.0, 7.0, 3.0, 7.0, 9.0, 4.0, 7.0, 1.0, 5.0, 7.0, 6.0, 8.0, 2.0, 3.0, 2.0, 1.0, 1.0, 1.0, 3.0]
global b_y = 10
global p = [0.369, 0.234, 0.915, 0.322, 0.562, 0.298, 0.687, 0.624, 0.164, 0.018, 0.604, 0.81, 0.421, 0.42, 0.06, 0.061, 0.806, 0.594, 0.053, 0.445, 0.529, 0.047, 0.491, 0.239, 0.649, 0.842, 0.231, 0.345, 0.235, 0.621, 0.947, 0.623, 0.992, 0.32, 0.18, 0.998, 0.49, 0.943, 0.524, 0.878, 0.645, 0.433, 0.424, 0.376, 0.485, 0.679, 0.178, 0.955, 0.182, 0.901, 0.011, 0.66, 0.289, 0.127, 0.349, 0.933, 0.106, 0.339, 0.176, 0.196, 0.09, 0.132, 0.529, 0.512, 0.741, 0.858, 0.783, 0.84, 0.014, 0.079, 0.592, 0.845, 0.198, 0.225, 0.12, 0.151, 0.397, 0.357, 0.991, 0.283, 0.504, 0.641, 0.327, 0.241, 0.864, 0.063, 0.383, 0.689, 0.948, 0.937, 0.961, 0.596, 0.146, 0.029, 0.05, 0.13, 0.881, 0.912, 0.907, 0.997, 0.576, 0.799, 0.23, 0.526, 0.619, 0.033, 0.581, 0.611, 0.837, 0.009, 0.955, 0.926, 0.377, 0.916, 0.489, 0.393, 0.169, 0.757, 0.767, 0.562, 0.26, 0.783, 0.945, 0.943, 0.346, 0.926, 0.302, 0.286, 0.724, 0.991, 0.526, 0.351, 0.886, 0.835, 0.194, 0.84, 0.437, 0.711, 0.183, 0.957, 0.669, 0.721, 0.891, 0.004, 0.032, 0.871, 0.539, 0.078, 0.85, 0.535, 0.412, 0.614, 0.136, 0.545, 0.458, 0.389, 0.884, 0.236, 0.061, 0.42, 0.665, 0.957, 0.071, 0.357, 0.502, 0.234, 0.25, 0.925, 0.721, 0.267, 0.369, 0.857, 0.231, 0.879, 0.282, 0.084, 0.605, 0.242, 0.565, 0.895, 0.773, 0.679, 0.424, 0.917, 0.345, 0.101, 0.864, 0.573, 0.094, 0.874, 0.414, 0.926, 0.231, 0.011, 0.728, 0.553, 0.483, 0.548, 0.907, 0.574, 0.241, 0.792, 0.995, 0.117, 0.532, 0.748, 0.195, 0.039, 0.371, 0.666, 0.062, 0.784, 0.012, 0.377, 0.353, 0.148, 0.188, 0.366, 0.462, 0.435, 0.575, 0.709, 0.144, 0.582, 0.25, 0.724, 0.399, 0.778, 0.749, 0.139, 0.503, 0.807, 0.241, 0.428, 0.043, 0.644, 0.148, 0.812, 0.246, 0.826, 0.724, 0.991, 0.055, 0.877, 0.105, 0.943, 0.432, 0.436, 0.93, 0.155, 0.394, 0.293, 0.764, 0.598, 0.515, 0.183, 0.985, 0.043, 0.952, 0.742, 0.709, 0.306, 0.936, 0.332, 0.533, 0.381, 0.239, 0.082, 0.523, 0.107, 0.731, 0.222, 0.061, 0.353, 0.756, 0.779, 0.096, 0.645, 0.599, 0.914, 0.309, 0.255, 0.462, 0.711, 0.272, 0.713, 0.754, 0.21, 0.881, 0.058, 0.038, 0.434, 0.618, 0.602, 0.734, 0.179, 0.357, 0.922, 0.761, 0.227, 0.692, 0.658, 0.012, 0.338, 0.593, 0.789, 0.164, 0.788, 0.602, 0.972, 0.755, 0.154, 0.463, 0.132, 0.246, 0.225, 0.714, 0.416]
global q = [0.567, 0.299, 0.923, 0.983, 0.726, 0.998, 0.758, 0.811, 0.601, 0.196, 0.626, 0.911, 0.739, 0.424, 0.524, 0.269, 0.9, 0.87, 0.203, 0.85, 0.845, 0.753, 0.499, 0.458, 0.799, 0.981, 0.802, 0.744, 0.703, 0.955, 0.978, 0.652, 0.999, 0.772, 0.25, 0.999, 0.907, 0.959, 0.839, 0.967, 0.819, 0.735, 0.657, 0.737, 0.64, 0.831, 0.52, 0.992, 0.689, 0.929, 0.487, 0.707, 0.836, 0.998, 0.76, 0.936, 0.111, 0.713, 0.54, 0.452, 0.637, 0.395, 0.641, 0.818, 0.874, 0.958, 0.931, 0.908, 0.207, 0.618, 0.92, 0.88, 0.579, 0.555, 0.266, 0.704, 0.514, 0.469, 0.997, 0.54, 0.958, 0.726, 0.346, 0.804, 0.974, 0.335, 0.574, 0.758, 0.988, 0.991, 0.971, 0.997, 0.573, 0.305, 0.134, 0.754, 0.94, 0.942, 0.928, 0.997, 0.711, 0.956, 0.575, 0.597, 0.956, 0.861, 0.625, 0.707, 0.983, 0.208, 0.996, 0.973, 0.523, 0.986, 0.761, 0.468, 0.753, 0.832, 0.925, 0.879, 0.293, 0.947, 0.957, 0.965, 0.449, 0.978, 0.974, 0.561, 0.773, 0.991, 0.674, 0.923, 0.901, 0.946, 0.438, 0.948, 0.492, 0.833, 0.783, 0.998, 0.98, 0.754, 0.955, 0.549, 0.812, 0.993, 0.873, 0.567, 0.922, 0.981, 0.639, 0.616, 0.672, 0.698, 0.508, 0.428, 0.97, 0.45, 0.886, 0.987, 0.779, 0.997, 0.74, 0.985, 0.977, 0.669, 0.884, 0.926, 0.997, 0.292, 0.953, 0.892, 0.238, 0.991, 0.843, 0.866, 0.944, 0.479, 0.941, 0.924, 0.895, 0.721, 0.533, 0.935, 0.395, 0.478, 0.892, 0.746, 0.977, 0.936, 0.791, 0.941, 0.901, 0.602, 0.802, 0.65, 0.794, 0.655, 0.97, 0.725, 0.462, 0.977, 0.996, 0.67, 0.811, 0.882, 0.274, 0.599, 0.495, 0.869, 0.284, 0.922, 0.641, 0.911, 0.514, 0.639, 0.902, 0.667, 0.485, 0.8, 0.825, 0.802, 0.414, 0.607, 0.562, 0.974, 0.938, 0.877, 0.853, 0.929, 0.631, 0.891, 0.574, 0.936, 0.336, 0.799, 0.637, 0.951, 0.815, 0.826, 0.783, 0.991, 0.83, 0.979, 0.139, 0.973, 0.465, 0.772, 0.948, 0.707, 0.488, 0.454, 0.841, 0.732, 0.745, 0.45, 0.99, 0.724, 0.983, 0.936, 0.719, 0.479, 0.946, 0.487, 0.709, 0.624, 0.499, 0.22, 0.885, 0.537, 0.958, 0.28, 0.51, 0.951, 0.963, 0.849, 0.491, 0.71, 0.934, 0.946, 0.64, 0.544, 0.582, 0.89, 0.504, 0.774, 0.996, 0.469, 0.914, 0.684, 0.254, 0.517, 0.894, 0.938, 0.774, 0.217, 0.974, 0.926, 0.805, 0.373, 0.928, 0.91, 0.391, 0.781, 0.898, 0.964, 0.705, 0.999, 0.663, 0.993, 0.991, 0.636, 0.799, 0.985, 0.918, 0.295, 0.99, 0.652]
global origin = 1
global destination = 60