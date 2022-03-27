global arcs = [1 2; 1 5; 1 8; 1 10; 1 11; 1 33; 1 34; 1 35; 1 44; 1 45; 2 16; 2 34; 2 47; 2 50; 3 7; 3 11; 3 16; 3 21; 3 26; 3 27; 3 32; 3 37; 3 48; 4 26; 4 30; 4 32; 4 34; 4 35; 4 43; 5 14; 5 30; 5 40; 5 41; 5 49; 6 5; 6 9; 6 17; 6 18; 6 30; 6 43; 7 16; 7 22; 7 26; 7 27; 7 42; 7 44; 8 34; 8 36; 9 3; 9 12; 9 25; 9 32; 9 33; 10 32; 10 34; 11 5; 11 41; 11 42; 11 44; 11 50; 12 2; 12 3; 12 9; 12 18; 12 24; 12 36; 12 41; 13 7; 13 30; 13 39; 14 22; 14 24; 14 48; 15 9; 15 12; 15 41; 15 46; 15 49; 16 34; 16 42; 16 45; 16 49; 17 15; 17 21; 17 37; 18 22; 18 29; 18 32; 18 37; 18 48; 19 2; 19 10; 19 11; 19 24; 19 40; 20 6; 20 14; 20 15; 20 44; 21 25; 21 29; 21 36; 21 37; 21 39; 21 43; 22 4; 22 47; 22 49; 23 24; 23 33; 23 37; 23 47; 24 16; 24 18; 24 27; 24 29; 24 31; 24 34; 24 41; 24 42; 24 47; 25 12; 25 28; 25 31; 25 34; 25 36; 25 44; 26 7; 26 45; 27 18; 27 30; 27 43; 28 10; 28 29; 28 35; 28 42; 29 5; 29 43; 30 2; 30 6; 30 8; 30 15; 30 18; 30 29; 30 37; 30 44; 30 47; 30 49; 31 2; 31 17; 31 26; 31 35; 31 38; 31 48; 32 4; 32 30; 32 31; 32 39; 32 47; 33 28; 33 46; 34 11; 34 13; 34 33; 34 49; 35 13; 35 16; 35 37; 36 9; 36 14; 36 16; 36 19; 36 24; 36 31; 37 13; 37 19; 37 22; 37 25; 37 35; 37 44; 38 20; 38 35; 38 42; 38 43; 38 48; 38 49; 39 8; 39 11; 39 15; 39 18; 39 26; 39 43; 40 4; 40 8; 40 12; 40 44; 40 47; 41 5; 41 8; 41 9; 41 14; 41 25; 41 27; 41 28; 41 38; 41 43; 42 18; 42 22; 42 30; 42 40; 43 8; 43 10; 43 12; 43 14; 43 36; 43 49; 44 10; 44 13; 44 18; 44 25; 44 33; 44 37; 44 41; 44 43; 44 49; 45 26; 45 36; 45 38; 45 50; 46 6; 46 15; 46 33; 46 38; 46 44; 46 50; 47 2; 47 6; 47 17; 47 23; 47 27; 47 46; 47 50; 48 13; 48 39; 48 42; 49 8; 49 15; 49 17; 49 23; 49 29; 49 42; 49 47; 49 48]
global d_x = [7.0, 6.0, 10.0, 10.0, 9.0, 5.0, 9.0, 1.0, 8.0, 8.0, 7.0, 6.0, 3.0, 9.0, 10.0, 1.0, 3.0, 4.0, 5.0, 1.0, 8.0, 6.0, 5.0, 10.0, 10.0, 1.0, 3.0, 4.0, 6.0, 7.0, 3.0, 7.0, 6.0, 1.0, 2.0, 7.0, 8.0, 4.0, 5.0, 5.0, 8.0, 7.0, 1.0, 6.0, 4.0, 5.0, 4.0, 4.0, 7.0, 5.0, 5.0, 10.0, 6.0, 1.0, 3.0, 10.0, 6.0, 5.0, 7.0, 3.0, 6.0, 3.0, 3.0, 8.0, 10.0, 1.0, 6.0, 4.0, 6.0, 7.0, 9.0, 1.0, 2.0, 9.0, 8.0, 4.0, 10.0, 9.0, 3.0, 9.0, 5.0, 6.0, 5.0, 3.0, 3.0, 9.0, 2.0, 9.0, 6.0, 1.0, 5.0, 7.0, 10.0, 2.0, 5.0, 3.0, 3.0, 1.0, 3.0, 4.0, 1.0, 6.0, 1.0, 2.0, 1.0, 8.0, 2.0, 10.0, 7.0, 8.0, 4.0, 8.0, 10.0, 9.0, 2.0, 9.0, 4.0, 9.0, 10.0, 1.0, 8.0, 4.0, 9.0, 4.0, 9.0, 9.0, 2.0, 5.0, 8.0, 9.0, 6.0, 5.0, 9.0, 3.0, 3.0, 9.0, 9.0, 7.0, 1.0, 8.0, 8.0, 10.0, 10.0, 8.0, 6.0, 5.0, 9.0, 5.0, 8.0, 5.0, 8.0, 5.0, 10.0, 2.0, 6.0, 7.0, 1.0, 7.0, 2.0, 5.0, 1.0, 7.0, 8.0, 10.0, 1.0, 2.0, 2.0, 2.0, 5.0, 8.0, 7.0, 1.0, 9.0, 2.0, 2.0, 4.0, 4.0, 2.0, 5.0, 8.0, 7.0, 6.0, 5.0, 10.0, 3.0, 4.0, 3.0, 2.0, 3.0, 1.0, 4.0, 9.0, 8.0, 10.0, 10.0, 1.0, 2.0, 10.0, 10.0, 7.0, 5.0, 1.0, 2.0, 1.0, 7.0, 3.0, 3.0, 3.0, 6.0, 7.0, 6.0, 9.0, 2.0, 1.0, 4.0, 6.0, 6.0, 2.0, 4.0, 5.0, 1.0, 9.0, 4.0, 7.0, 6.0, 3.0, 9.0, 8.0, 6.0, 4.0, 9.0, 6.0, 6.0, 2.0, 9.0, 1.0, 8.0, 1.0, 9.0, 7.0, 3.0, 6.0, 5.0, 1.0, 7.0, 9.0, 2.0, 2.0, 9.0, 2.0, 1.0, 4.0, 2.0]
global b_x = 5
global d_y = [8.0, 4.0, 1.0, 7.0, 1.0, 4.0, 4.0, 3.0, 7.0, 1.0, 5.0, 3.0, 7.0, 8.0, 2.0, 4.0, 6.0, 2.0, 7.0, 9.0, 2.0, 6.0, 7.0, 1.0, 7.0, 5.0, 3.0, 2.0, 10.0, 2.0, 5.0, 8.0, 8.0, 7.0, 7.0, 10.0, 1.0, 3.0, 6.0, 1.0, 10.0, 9.0, 1.0, 8.0, 3.0, 8.0, 6.0, 6.0, 8.0, 9.0, 3.0, 4.0, 8.0, 6.0, 10.0, 2.0, 4.0, 5.0, 7.0, 9.0, 7.0, 6.0, 2.0, 4.0, 10.0, 6.0, 2.0, 9.0, 9.0, 3.0, 3.0, 4.0, 1.0, 5.0, 6.0, 9.0, 5.0, 8.0, 2.0, 9.0, 8.0, 1.0, 6.0, 1.0, 1.0, 5.0, 6.0, 2.0, 9.0, 5.0, 8.0, 6.0, 2.0, 3.0, 4.0, 3.0, 7.0, 9.0, 2.0, 6.0, 9.0, 5.0, 2.0, 5.0, 9.0, 8.0, 10.0, 5.0, 6.0, 4.0, 5.0, 5.0, 2.0, 8.0, 9.0, 3.0, 4.0, 7.0, 3.0, 5.0, 9.0, 3.0, 3.0, 2.0, 9.0, 4.0, 8.0, 9.0, 9.0, 1.0, 5.0, 2.0, 3.0, 9.0, 1.0, 4.0, 7.0, 4.0, 8.0, 8.0, 3.0, 4.0, 5.0, 1.0, 8.0, 8.0, 2.0, 2.0, 1.0, 3.0, 10.0, 8.0, 5.0, 2.0, 7.0, 7.0, 5.0, 3.0, 5.0, 3.0, 10.0, 6.0, 3.0, 7.0, 3.0, 8.0, 4.0, 8.0, 3.0, 7.0, 6.0, 3.0, 3.0, 3.0, 6.0, 3.0, 2.0, 1.0, 9.0, 10.0, 3.0, 3.0, 5.0, 3.0, 9.0, 2.0, 1.0, 5.0, 3.0, 9.0, 6.0, 5.0, 7.0, 4.0, 7.0, 6.0, 4.0, 3.0, 10.0, 2.0, 6.0, 8.0, 8.0, 5.0, 1.0, 3.0, 6.0, 3.0, 7.0, 6.0, 10.0, 3.0, 8.0, 1.0, 2.0, 7.0, 9.0, 5.0, 10.0, 8.0, 1.0, 8.0, 7.0, 3.0, 10.0, 8.0, 5.0, 9.0, 8.0, 9.0, 3.0, 10.0, 5.0, 1.0, 4.0, 3.0, 8.0, 8.0, 4.0, 2.0, 9.0, 6.0, 6.0, 1.0, 4.0, 8.0, 10.0, 4.0, 3.0, 2.0, 2.0, 1.0, 1.0]
global b_y = 10
global p = [0.855, 0.05, 0.57, 0.788, 0.812, 0.546, 0.636, 0.692, 0.861, 0.586, 0.257, 0.252, 0.52, 0.086, 0.952, 0.084, 0.099, 0.525, 0.754, 0.12, 0.857, 0.213, 0.004, 0.85, 0.172, 0.858, 0.319, 0.193, 0.248, 0.678, 0.124, 0.193, 0.476, 0.382, 0.747, 0.355, 0.006, 0.392, 0.22, 0.745, 0.912, 0.324, 0.859, 0.316, 0.746, 0.377, 0.766, 0.573, 0.309, 0.852, 0.128, 0.002, 0.236, 0.213, 0.416, 0.268, 0.312, 0.835, 0.69, 0.316, 0.75, 0.319, 0.855, 0.573, 0.711, 0.186, 0.109, 0.003, 0.435, 0.654, 0.9, 0.962, 0.783, 0.467, 0.2, 0.284, 0.906, 0.405, 0.093, 0.669, 0.86, 0.759, 0.975, 0.144, 0.069, 0.905, 0.057, 0.211, 0.612, 0.298, 0.218, 0.454, 0.389, 0.17, 0.663, 0.67, 0.795, 0.954, 0.379, 0.257, 0.081, 0.263, 0.434, 0.657, 0.325, 0.626, 0.337, 0.506, 0.467, 0.365, 0.449, 0.208, 0.343, 0.741, 0.596, 0.289, 0.347, 0.46, 0.294, 0.3, 0.27, 0.665, 0.384, 0.782, 0.941, 0.865, 0.036, 0.344, 0.42, 0.503, 0.177, 0.289, 0.764, 0.646, 0.604, 0.01, 0.479, 0.194, 0.238, 0.582, 0.185, 0.802, 0.373, 0.927, 0.902, 0.258, 0.707, 0.986, 0.61, 0.313, 0.027, 0.527, 0.301, 0.941, 0.13, 0.961, 0.724, 0.613, 0.211, 0.498, 0.391, 0.764, 0.14, 0.028, 0.878, 0.922, 0.463, 0.679, 0.939, 0.167, 0.07, 0.177, 0.741, 0.846, 0.635, 0.829, 0.329, 0.014, 0.144, 0.952, 0.578, 0.275, 0.441, 0.199, 0.765, 0.096, 0.422, 0.24, 0.218, 0.667, 0.118, 0.663, 0.581, 0.795, 0.624, 0.11, 0.502, 0.653, 0.868, 0.283, 0.546, 0.867, 0.665, 0.65, 0.686, 0.007, 0.004, 0.045, 0.021, 0.901, 0.799, 0.326, 0.231, 0.163, 0.988, 0.95, 0.688, 0.637, 0.603, 0.674, 0.362, 0.883, 0.446, 0.927, 0.242, 0.741, 0.349, 0.762, 0.044, 0.726, 0.137, 0.136, 0.523, 0.553, 0.768, 0.926, 0.455, 0.287, 0.035, 0.285, 0.286, 0.893, 0.607, 0.025, 0.172, 0.285, 0.941, 0.678, 0.232, 0.398, 0.604, 0.404, 0.956]
global q = [0.938, 0.486, 0.958, 0.813, 0.889, 0.947, 0.812, 0.875, 0.897, 0.788, 0.423, 0.339, 0.859, 0.903, 0.966, 0.52, 0.48, 0.944, 0.849, 0.624, 0.95, 0.507, 0.779, 0.937, 0.439, 0.956, 0.753, 0.992, 0.819, 0.952, 0.168, 0.889, 0.896, 0.816, 0.916, 0.93, 0.349, 0.72, 0.984, 0.987, 0.924, 0.915, 0.877, 0.493, 0.789, 0.94, 0.792, 0.908, 0.561, 0.947, 0.345, 0.432, 0.36, 0.325, 0.828, 0.786, 0.823, 0.934, 0.833, 0.425, 0.788, 0.384, 0.876, 0.866, 0.896, 0.545, 0.723, 0.347, 0.586, 0.696, 0.963, 0.978, 0.798, 0.522, 0.692, 0.666, 0.958, 0.573, 0.81, 0.958, 0.976, 0.799, 0.988, 0.299, 0.85, 0.993, 0.071, 0.514, 0.652, 0.471, 0.517, 0.884, 0.768, 0.962, 0.97, 0.961, 0.974, 0.993, 0.743, 0.518, 0.191, 0.639, 0.948, 0.907, 0.455, 0.944, 0.63, 0.988, 0.498, 0.405, 0.865, 0.946, 0.636, 0.749, 0.929, 0.523, 0.643, 0.575, 0.511, 0.534, 0.971, 0.751, 0.826, 0.883, 0.948, 0.899, 0.669, 0.864, 0.487, 0.764, 0.995, 0.505, 0.781, 0.721, 0.878, 0.098, 0.768, 0.277, 0.307, 0.789, 0.273, 0.912, 0.848, 0.933, 0.927, 0.617, 0.932, 0.986, 0.788, 0.97, 0.202, 0.944, 0.95, 0.982, 0.367, 0.976, 0.903, 0.799, 0.693, 0.777, 0.79, 0.808, 0.39, 0.094, 0.904, 0.988, 0.492, 0.696, 0.953, 0.952, 0.563, 0.824, 0.872, 0.931, 0.855, 0.834, 0.998, 0.422, 0.735, 0.973, 0.957, 0.681, 0.715, 0.449, 0.867, 0.313, 0.616, 0.502, 0.878, 0.713, 0.12, 0.925, 0.866, 0.858, 0.875, 0.391, 0.75, 0.957, 0.994, 0.463, 0.65, 0.897, 0.824, 0.804, 0.766, 0.723, 0.446, 0.603, 0.585, 0.97, 0.924, 0.391, 0.74, 0.564, 0.995, 0.987, 0.786, 0.708, 0.947, 0.961, 0.653, 0.917, 0.77, 0.941, 0.674, 0.809, 0.531, 0.832, 0.302, 0.811, 0.913, 0.943, 0.933, 0.853, 0.784, 0.967, 0.5, 0.423, 0.476, 0.722, 0.715, 0.971, 0.631, 0.108, 0.795, 0.445, 0.98, 0.788, 0.71, 0.94, 0.821, 0.957, 0.971]
global origin = 1
global destination = 50