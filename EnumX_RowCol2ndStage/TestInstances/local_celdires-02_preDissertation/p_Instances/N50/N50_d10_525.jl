global arcs = [1 21; 1 29; 1 34; 2 14; 2 20; 2 23; 2 29; 2 37; 3 9; 3 13; 3 18; 3 21; 3 26; 3 28; 3 32; 3 34; 3 39; 4 47; 4 50; 5 6; 6 4; 6 17; 6 25; 6 35; 7 8; 7 12; 7 31; 8 33; 8 35; 8 37; 8 43; 8 44; 9 2; 9 23; 9 28; 9 30; 9 47; 10 17; 10 22; 10 30; 10 35; 10 50; 11 3; 11 6; 11 10; 11 14; 11 17; 11 23; 11 37; 11 45; 11 47; 12 31; 12 32; 12 36; 12 47; 13 20; 13 42; 13 44; 13 45; 14 6; 14 29; 14 34; 14 43; 14 46; 15 13; 15 30; 15 38; 15 39; 15 48; 16 40; 17 27; 17 34; 18 8; 18 32; 18 40; 18 42; 19 4; 19 18; 19 23; 19 39; 19 49; 20 3; 20 6; 20 12; 20 25; 20 32; 20 38; 20 43; 21 8; 21 15; 21 43; 22 2; 22 7; 22 19; 22 31; 22 34; 22 40; 22 41; 22 42; 23 13; 23 21; 24 6; 24 7; 24 8; 24 19; 24 32; 24 47; 25 20; 25 21; 25 26; 25 35; 25 38; 25 49; 26 4; 26 11; 26 31; 26 37; 26 38; 27 33; 27 43; 28 15; 28 38; 28 42; 28 47; 29 7; 29 11; 29 32; 29 33; 29 40; 29 43; 30 14; 30 42; 31 12; 31 21; 31 29; 31 32; 31 40; 31 44; 31 46; 32 11; 32 25; 32 38; 33 3; 33 20; 33 21; 33 48; 34 2; 34 28; 34 35; 35 16; 35 17; 35 27; 35 42; 36 3; 36 6; 36 16; 36 45; 36 49; 37 2; 37 31; 38 4; 38 15; 38 26; 38 34; 38 49; 39 7; 39 8; 39 15; 39 19; 39 22; 39 27; 39 28; 39 31; 39 33; 39 38; 39 42; 39 50; 40 10; 40 14; 40 27; 40 34; 40 48; 41 6; 41 18; 41 42; 41 45; 41 50; 42 15; 42 20; 42 23; 42 43; 43 5; 43 22; 43 24; 43 26; 43 30; 43 50; 44 9; 44 16; 44 25; 45 10; 45 17; 45 43; 46 6; 46 8; 46 31; 46 35; 46 49; 47 2; 47 16; 47 17; 47 21; 47 23; 47 38; 47 39; 47 41; 48 7; 48 32; 48 36; 48 44; 49 13; 49 19; 49 22; 49 33]
global d_x = [9.0, 10.0, 6.0, 6.0, 1.0, 6.0, 9.0, 8.0, 10.0, 10.0, 7.0, 4.0, 1.0, 7.0, 3.0, 5.0, 8.0, 10.0, 8.0, 5.0, 6.0, 1.0, 10.0, 4.0, 5.0, 7.0, 5.0, 1.0, 5.0, 7.0, 3.0, 6.0, 5.0, 3.0, 8.0, 9.0, 9.0, 6.0, 4.0, 2.0, 7.0, 7.0, 7.0, 10.0, 9.0, 4.0, 4.0, 4.0, 8.0, 1.0, 6.0, 4.0, 10.0, 3.0, 7.0, 5.0, 7.0, 7.0, 3.0, 3.0, 5.0, 6.0, 5.0, 4.0, 10.0, 8.0, 4.0, 9.0, 9.0, 1.0, 10.0, 3.0, 4.0, 3.0, 1.0, 8.0, 7.0, 5.0, 2.0, 3.0, 6.0, 2.0, 5.0, 9.0, 1.0, 10.0, 10.0, 7.0, 5.0, 7.0, 4.0, 8.0, 5.0, 10.0, 7.0, 4.0, 1.0, 5.0, 3.0, 2.0, 6.0, 9.0, 9.0, 3.0, 6.0, 7.0, 3.0, 3.0, 8.0, 2.0, 1.0, 10.0, 7.0, 7.0, 5.0, 9.0, 5.0, 1.0, 3.0, 5.0, 6.0, 7.0, 2.0, 8.0, 1.0, 10.0, 9.0, 2.0, 3.0, 2.0, 5.0, 6.0, 6.0, 8.0, 6.0, 5.0, 1.0, 2.0, 6.0, 2.0, 6.0, 3.0, 5.0, 4.0, 8.0, 7.0, 6.0, 10.0, 4.0, 1.0, 9.0, 7.0, 2.0, 10.0, 1.0, 10.0, 6.0, 10.0, 7.0, 2.0, 7.0, 6.0, 4.0, 9.0, 5.0, 10.0, 6.0, 2.0, 10.0, 9.0, 8.0, 5.0, 4.0, 5.0, 8.0, 1.0, 10.0, 2.0, 1.0, 2.0, 6.0, 3.0, 8.0, 8.0, 1.0, 1.0, 10.0, 7.0, 10.0, 1.0, 8.0, 10.0, 4.0, 5.0, 5.0, 2.0, 2.0, 1.0, 1.0, 9.0, 5.0, 4.0, 10.0, 6.0, 4.0, 10.0, 10.0, 3.0, 7.0, 7.0, 5.0, 4.0, 1.0, 8.0, 2.0, 9.0, 8.0, 6.0, 3.0, 1.0, 7.0, 2.0, 9.0, 7.0]
global b_x = 5
global d_y = [5.0, 7.0, 3.0, 7.0, 7.0, 7.0, 10.0, 7.0, 2.0, 1.0, 1.0, 5.0, 7.0, 5.0, 6.0, 3.0, 9.0, 6.0, 1.0, 8.0, 6.0, 4.0, 7.0, 5.0, 4.0, 6.0, 7.0, 6.0, 5.0, 3.0, 8.0, 10.0, 2.0, 5.0, 2.0, 7.0, 1.0, 5.0, 5.0, 8.0, 3.0, 7.0, 8.0, 4.0, 2.0, 4.0, 2.0, 2.0, 5.0, 7.0, 8.0, 8.0, 10.0, 2.0, 9.0, 9.0, 5.0, 8.0, 8.0, 6.0, 10.0, 2.0, 3.0, 5.0, 5.0, 3.0, 6.0, 5.0, 4.0, 5.0, 6.0, 3.0, 4.0, 5.0, 7.0, 2.0, 3.0, 10.0, 1.0, 9.0, 7.0, 5.0, 8.0, 1.0, 4.0, 4.0, 3.0, 2.0, 9.0, 9.0, 9.0, 5.0, 8.0, 7.0, 2.0, 6.0, 10.0, 5.0, 10.0, 2.0, 7.0, 1.0, 7.0, 8.0, 3.0, 9.0, 2.0, 3.0, 7.0, 4.0, 5.0, 7.0, 6.0, 7.0, 2.0, 10.0, 6.0, 4.0, 6.0, 1.0, 6.0, 8.0, 8.0, 9.0, 3.0, 2.0, 2.0, 9.0, 6.0, 6.0, 8.0, 4.0, 1.0, 7.0, 8.0, 4.0, 10.0, 4.0, 1.0, 3.0, 3.0, 8.0, 5.0, 4.0, 5.0, 6.0, 1.0, 2.0, 9.0, 10.0, 9.0, 5.0, 5.0, 3.0, 4.0, 9.0, 10.0, 10.0, 1.0, 10.0, 1.0, 8.0, 8.0, 6.0, 5.0, 2.0, 3.0, 4.0, 10.0, 5.0, 6.0, 5.0, 9.0, 7.0, 1.0, 7.0, 3.0, 9.0, 4.0, 3.0, 7.0, 4.0, 10.0, 2.0, 10.0, 1.0, 3.0, 5.0, 7.0, 3.0, 2.0, 3.0, 8.0, 2.0, 10.0, 5.0, 7.0, 6.0, 7.0, 10.0, 7.0, 5.0, 7.0, 1.0, 8.0, 4.0, 8.0, 9.0, 8.0, 9.0, 5.0, 6.0, 4.0, 9.0, 2.0, 6.0, 3.0, 9.0, 2.0, 2.0, 9.0, 3.0, 1.0, 4.0]
global b_y = 10
global p = [0.239, 0.196, 0.66, 0.138, 0.211, 0.427, 0.661, 0.361, 0.712, 0.453, 0.205, 0.05, 0.217, 0.322, 0.752, 0.034, 0.277, 0.705, 0.418, 0.534, 0.216, 0.261, 0.103, 0.54, 0.306, 0.658, 0.592, 0.764, 0.653, 0.86, 0.84, 0.068, 0.329, 0.877, 0.837, 0.97, 0.344, 0.979, 0.728, 0.572, 0.212, 0.591, 0.945, 0.842, 0.75, 0.457, 0.345, 0.211, 0.206, 0.194, 0.615, 0.211, 0.409, 0.609, 0.067, 0.814, 0.122, 0.985, 0.072, 0.912, 0.103, 0.666, 0.433, 0.672, 0.149, 0.113, 0.551, 0.467, 0.968, 0.736, 0.448, 0.77, 0.283, 0.196, 0.191, 0.406, 0.833, 0.218, 0.234, 0.317, 0.053, 0.292, 0.059, 0.905, 0.338, 0.505, 0.716, 0.803, 0.103, 0.615, 0.926, 0.281, 0.624, 0.308, 0.029, 0.35, 0.253, 0.658, 0.525, 0.443, 0.931, 0.633, 0.912, 0.631, 0.622, 0.376, 0.741, 0.037, 0.279, 0.427, 0.221, 0.528, 0.244, 0.714, 0.953, 0.202, 0.953, 0.463, 0.247, 0.768, 0.926, 0.423, 0.822, 0.8, 0.382, 0.091, 0.146, 0.35, 0.686, 0.894, 0.608, 0.107, 0.746, 0.509, 0.388, 0.984, 0.132, 0.7, 0.177, 0.88, 0.479, 0.034, 0.349, 0.218, 0.815, 0.829, 0.726, 0.753, 0.295, 0.242, 0.175, 0.814, 0.237, 0.63, 0.414, 0.177, 0.75, 0.998, 0.339, 0.511, 0.089, 0.747, 0.576, 0.138, 0.46, 0.537, 0.346, 0.528, 0.195, 0.752, 0.271, 0.979, 0.129, 0.15, 0.86, 0.109, 0.691, 0.616, 0.921, 0.734, 0.564, 0.939, 0.194, 0.601, 0.815, 0.656, 0.952, 0.529, 0.096, 0.012, 0.059, 0.6, 0.193, 0.686, 0.512, 0.241, 0.663, 0.319, 0.299, 0.78, 0.01, 0.622, 0.706, 0.621, 0.73, 0.452, 0.486, 0.283, 0.664, 0.021, 0.597, 0.287, 0.575, 0.996, 0.634, 0.077, 0.021, 0.079, 0.708, 0.88, 0.597, 0.375, 0.823, 0.884]
global q = [0.998, 0.342, 0.948, 0.184, 0.26, 0.766, 0.95, 0.812, 0.972, 0.719, 0.873, 0.1, 0.957, 0.622, 0.928, 0.584, 0.454, 0.987, 0.976, 0.909, 0.526, 0.573, 0.268, 0.653, 0.757, 0.964, 0.607, 0.992, 0.882, 0.956, 0.962, 0.831, 0.643, 0.943, 0.959, 0.998, 0.714, 0.979, 0.808, 0.735, 0.503, 0.718, 0.965, 0.941, 0.798, 0.96, 0.527, 0.488, 0.242, 0.667, 0.982, 0.901, 0.743, 0.838, 0.661, 0.897, 0.226, 0.998, 0.436, 0.941, 0.54, 0.798, 0.52, 0.943, 0.16, 0.301, 0.701, 0.517, 0.977, 0.76, 0.471, 0.919, 0.959, 0.403, 0.389, 0.624, 0.917, 0.756, 0.37, 0.374, 0.817, 0.961, 0.783, 0.917, 0.363, 0.529, 0.85, 0.858, 0.83, 0.912, 0.974, 0.802, 0.749, 0.366, 0.621, 0.748, 0.969, 0.667, 0.689, 0.835, 0.989, 0.757, 0.957, 0.862, 0.731, 0.597, 0.847, 0.344, 0.928, 0.565, 0.718, 0.586, 0.428, 0.825, 0.997, 0.656, 0.969, 0.978, 0.91, 0.782, 0.99, 0.878, 0.944, 0.886, 0.713, 0.463, 0.748, 0.956, 0.9, 0.925, 0.767, 0.552, 0.945, 0.617, 0.918, 0.993, 0.298, 0.745, 0.731, 0.947, 0.741, 0.907, 0.824, 0.624, 0.848, 0.991, 0.839, 0.823, 0.969, 0.954, 0.385, 0.908, 0.521, 0.899, 0.446, 0.532, 0.775, 0.998, 0.456, 0.743, 0.82, 0.981, 0.692, 0.189, 0.512, 0.963, 0.871, 0.701, 0.583, 0.802, 0.962, 0.992, 0.237, 0.655, 0.979, 0.189, 0.817, 0.643, 0.999, 0.979, 0.74, 0.991, 0.573, 0.65, 0.901, 0.664, 0.981, 0.915, 0.561, 0.99, 0.373, 0.729, 0.904, 0.959, 0.81, 0.703, 0.783, 0.761, 0.416, 0.86, 0.434, 0.757, 0.911, 0.755, 0.775, 0.762, 0.565, 0.306, 0.848, 0.893, 0.895, 0.644, 0.73, 0.999, 0.917, 0.696, 0.615, 0.303, 0.769, 0.904, 0.7, 0.379, 0.895, 0.927]
global origin = 1
global destination = 50