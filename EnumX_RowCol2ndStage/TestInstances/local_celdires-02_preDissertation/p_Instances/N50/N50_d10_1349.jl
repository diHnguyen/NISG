global arcs = [1 6; 1 39; 1 42; 1 47; 2 17; 2 36; 2 37; 3 2; 3 4; 3 17; 3 30; 3 32; 3 49; 4 12; 4 22; 4 35; 4 47; 5 7; 5 10; 5 28; 5 38; 5 43; 5 49; 6 14; 6 18; 6 23; 6 33; 7 4; 7 12; 7 34; 7 42; 8 3; 8 35; 8 44; 8 45; 8 48; 9 7; 9 16; 9 25; 9 33; 9 36; 9 40; 9 46; 10 8; 10 13; 10 14; 10 26; 10 41; 10 48; 10 49; 10 50; 11 27; 11 34; 12 8; 12 36; 12 39; 12 44; 13 2; 13 10; 14 3; 14 7; 14 16; 15 20; 15 24; 15 31; 15 32; 15 46; 15 48; 16 4; 16 8; 16 19; 16 24; 16 48; 17 15; 17 16; 17 43; 17 50; 18 16; 18 19; 18 22; 18 29; 19 12; 19 16; 19 26; 19 42; 19 43; 20 3; 20 5; 20 9; 20 29; 20 31; 20 33; 20 40; 21 6; 21 9; 21 20; 21 25; 21 30; 21 31; 21 48; 22 45; 22 49; 23 2; 23 8; 23 18; 23 34; 23 49; 24 2; 24 3; 24 32; 24 40; 24 43; 25 20; 25 23; 25 37; 25 48; 25 50; 26 5; 26 6; 26 9; 26 11; 26 18; 26 23; 26 35; 26 40; 27 10; 27 21; 27 26; 27 34; 28 2; 28 23; 28 24; 28 41; 28 45; 28 49; 29 19; 29 35; 29 40; 30 11; 30 32; 30 37; 30 49; 31 6; 31 12; 31 13; 31 16; 31 29; 31 32; 31 35; 31 42; 32 9; 32 16; 32 34; 32 36; 32 41; 32 43; 32 45; 33 36; 33 47; 34 2; 34 18; 34 25; 34 48; 35 23; 35 38; 35 41; 36 5; 36 9; 36 16; 36 21; 36 44; 37 8; 37 25; 37 44; 37 46; 38 15; 38 24; 38 39; 39 7; 39 9; 39 16; 39 24; 39 25; 39 27; 39 34; 39 35; 39 44; 39 46; 39 49; 40 8; 40 15; 40 50; 41 20; 41 34; 41 40; 41 46; 41 48; 42 5; 42 6; 42 8; 42 13; 43 20; 43 24; 43 25; 44 9; 44 13; 44 32; 44 48; 45 4; 45 22; 45 37; 45 47; 46 2; 46 11; 46 12; 46 13; 46 19; 46 23; 46 32; 46 37; 46 42; 47 5; 47 20; 47 23; 47 38; 48 5; 48 6; 48 16; 48 18; 48 49; 49 3; 49 9; 49 13; 49 38; 49 47; 49 48]
global d_x = [3.0, 3.0, 10.0, 7.0, 9.0, 2.0, 6.0, 5.0, 1.0, 5.0, 8.0, 2.0, 10.0, 9.0, 1.0, 5.0, 2.0, 8.0, 4.0, 4.0, 6.0, 3.0, 9.0, 3.0, 10.0, 8.0, 5.0, 10.0, 3.0, 5.0, 4.0, 3.0, 2.0, 4.0, 6.0, 5.0, 7.0, 8.0, 10.0, 4.0, 8.0, 2.0, 7.0, 5.0, 4.0, 2.0, 9.0, 3.0, 4.0, 6.0, 3.0, 6.0, 10.0, 3.0, 9.0, 1.0, 2.0, 2.0, 9.0, 2.0, 4.0, 6.0, 7.0, 6.0, 4.0, 9.0, 1.0, 3.0, 5.0, 6.0, 8.0, 5.0, 10.0, 4.0, 6.0, 2.0, 4.0, 9.0, 4.0, 3.0, 1.0, 2.0, 5.0, 10.0, 8.0, 7.0, 6.0, 8.0, 2.0, 7.0, 5.0, 4.0, 6.0, 10.0, 6.0, 1.0, 10.0, 5.0, 4.0, 5.0, 9.0, 3.0, 3.0, 9.0, 10.0, 7.0, 2.0, 5.0, 5.0, 8.0, 2.0, 10.0, 9.0, 7.0, 4.0, 8.0, 9.0, 7.0, 5.0, 2.0, 5.0, 6.0, 9.0, 6.0, 7.0, 4.0, 10.0, 6.0, 7.0, 10.0, 8.0, 4.0, 6.0, 2.0, 6.0, 8.0, 10.0, 8.0, 8.0, 4.0, 4.0, 2.0, 2.0, 4.0, 5.0, 10.0, 2.0, 6.0, 10.0, 3.0, 3.0, 7.0, 7.0, 9.0, 7.0, 8.0, 1.0, 8.0, 3.0, 6.0, 1.0, 10.0, 10.0, 7.0, 8.0, 2.0, 2.0, 6.0, 10.0, 2.0, 8.0, 7.0, 3.0, 10.0, 8.0, 1.0, 9.0, 10.0, 7.0, 3.0, 2.0, 4.0, 9.0, 10.0, 6.0, 7.0, 2.0, 6.0, 6.0, 6.0, 6.0, 5.0, 4.0, 3.0, 10.0, 5.0, 1.0, 3.0, 1.0, 6.0, 5.0, 10.0, 10.0, 7.0, 1.0, 10.0, 7.0, 4.0, 2.0, 1.0, 8.0, 8.0, 5.0, 7.0, 2.0, 1.0, 10.0, 2.0, 2.0, 3.0, 4.0, 5.0, 4.0, 10.0, 1.0, 1.0, 5.0, 2.0, 7.0, 5.0, 7.0, 9.0, 2.0, 6.0, 10.0, 2.0]
global b_x = 5
global d_y = [4.0, 8.0, 3.0, 7.0, 10.0, 9.0, 4.0, 8.0, 1.0, 6.0, 7.0, 2.0, 8.0, 2.0, 9.0, 1.0, 1.0, 5.0, 9.0, 9.0, 4.0, 9.0, 9.0, 10.0, 2.0, 1.0, 1.0, 4.0, 2.0, 7.0, 3.0, 4.0, 6.0, 5.0, 8.0, 7.0, 8.0, 1.0, 1.0, 1.0, 7.0, 10.0, 6.0, 9.0, 4.0, 10.0, 4.0, 6.0, 10.0, 3.0, 1.0, 2.0, 7.0, 10.0, 3.0, 6.0, 5.0, 3.0, 1.0, 8.0, 6.0, 10.0, 10.0, 1.0, 7.0, 5.0, 2.0, 6.0, 5.0, 8.0, 9.0, 3.0, 2.0, 6.0, 9.0, 3.0, 6.0, 6.0, 1.0, 3.0, 2.0, 2.0, 4.0, 6.0, 6.0, 6.0, 4.0, 1.0, 7.0, 6.0, 6.0, 10.0, 6.0, 2.0, 6.0, 1.0, 8.0, 1.0, 7.0, 9.0, 8.0, 8.0, 8.0, 8.0, 7.0, 5.0, 9.0, 6.0, 4.0, 1.0, 6.0, 3.0, 1.0, 7.0, 7.0, 3.0, 3.0, 5.0, 5.0, 10.0, 5.0, 6.0, 5.0, 3.0, 9.0, 10.0, 4.0, 1.0, 1.0, 1.0, 7.0, 9.0, 7.0, 6.0, 3.0, 4.0, 9.0, 1.0, 7.0, 10.0, 10.0, 10.0, 8.0, 2.0, 4.0, 9.0, 1.0, 2.0, 8.0, 6.0, 7.0, 6.0, 10.0, 4.0, 1.0, 10.0, 2.0, 3.0, 6.0, 9.0, 4.0, 1.0, 8.0, 9.0, 6.0, 2.0, 5.0, 7.0, 7.0, 2.0, 9.0, 1.0, 8.0, 3.0, 8.0, 3.0, 2.0, 6.0, 2.0, 3.0, 2.0, 10.0, 6.0, 5.0, 4.0, 3.0, 8.0, 9.0, 4.0, 5.0, 10.0, 2.0, 3.0, 8.0, 6.0, 3.0, 9.0, 1.0, 8.0, 1.0, 5.0, 9.0, 9.0, 7.0, 4.0, 10.0, 2.0, 9.0, 9.0, 10.0, 6.0, 10.0, 3.0, 7.0, 6.0, 8.0, 10.0, 8.0, 8.0, 10.0, 8.0, 4.0, 3.0, 10.0, 9.0, 8.0, 4.0, 5.0, 10.0, 4.0, 4.0, 2.0, 6.0, 8.0, 9.0, 3.0]
global b_y = 10
global p = [0.905, 0.936, 0.559, 0.2, 0.199, 0.301, 0.806, 0.657, 0.04, 0.56, 0.827, 0.108, 0.363, 0.208, 0.648, 0.224, 0.514, 0.874, 0.547, 0.955, 0.901, 0.551, 0.534, 0.728, 0.757, 0.129, 0.692, 0.328, 0.987, 0.976, 0.142, 0.171, 0.929, 0.709, 0.128, 0.328, 0.49, 0.97, 0.799, 0.048, 0.173, 0.798, 0.284, 0.553, 0.155, 0.258, 0.354, 0.286, 0.198, 0.589, 0.659, 0.468, 0.907, 0.726, 0.091, 0.293, 0.55, 0.188, 0.407, 0.466, 0.544, 0.775, 0.011, 0.971, 0.889, 0.801, 0.445, 0.482, 0.595, 0.168, 0.368, 0.346, 0.92, 0.375, 0.194, 0.096, 0.27, 0.359, 0.319, 0.179, 0.848, 0.092, 0.314, 0.148, 0.16, 0.807, 0.325, 0.894, 0.649, 0.339, 0.663, 0.305, 0.899, 0.287, 0.462, 0.852, 0.208, 0.743, 0.657, 0.066, 0.639, 0.635, 0.711, 0.425, 0.123, 0.836, 0.121, 0.437, 0.003, 0.428, 0.948, 0.764, 0.995, 0.624, 0.341, 0.96, 0.9, 0.961, 0.511, 0.162, 0.03, 0.084, 0.973, 0.352, 0.892, 0.887, 0.764, 0.13, 0.018, 0.253, 0.282, 0.688, 0.927, 0.246, 0.809, 0.975, 0.5, 0.52, 0.82, 0.978, 0.122, 0.433, 0.831, 0.151, 0.284, 0.771, 0.256, 0.203, 0.573, 0.269, 0.717, 0.251, 0.922, 0.935, 0.775, 0.106, 0.091, 0.738, 0.034, 0.801, 0.834, 0.839, 0.958, 0.142, 0.701, 0.044, 0.646, 0.071, 0.035, 0.877, 0.472, 0.859, 0.511, 0.07, 0.205, 0.229, 0.338, 0.748, 0.537, 0.863, 0.096, 0.148, 0.375, 0.938, 0.304, 0.95, 0.136, 0.64, 0.096, 0.745, 0.927, 0.367, 0.46, 0.934, 0.555, 0.689, 0.711, 0.485, 0.3, 0.761, 0.265, 0.242, 0.28, 0.139, 0.58, 0.354, 0.246, 0.98, 0.317, 0.612, 0.934, 0.003, 0.078, 0.128, 0.863, 0.506, 0.09, 0.661, 0.352, 0.581, 0.747, 0.563, 0.991, 0.618, 0.546, 0.646, 0.757, 0.532, 0.554, 0.831, 0.899, 0.902, 0.873, 0.173, 0.916, 0.404]
global q = [0.966, 0.977, 0.756, 0.912, 0.813, 0.791, 0.817, 0.76, 0.375, 0.882, 0.864, 0.688, 0.684, 0.678, 0.909, 0.815, 0.935, 0.879, 0.949, 0.992, 0.931, 0.784, 0.653, 0.986, 0.805, 0.2, 0.822, 0.961, 0.995, 0.997, 0.718, 0.979, 0.956, 0.809, 0.408, 0.748, 0.532, 0.981, 0.886, 0.162, 0.337, 0.846, 0.978, 0.969, 0.508, 0.514, 0.365, 0.393, 0.222, 0.616, 0.842, 0.559, 0.945, 0.759, 0.898, 0.53, 0.896, 0.465, 0.632, 0.716, 0.976, 0.778, 0.076, 0.997, 0.969, 0.897, 0.696, 0.973, 0.926, 0.899, 0.814, 0.657, 0.937, 0.979, 0.367, 0.115, 0.482, 0.739, 0.508, 0.355, 0.942, 0.563, 0.378, 0.209, 0.258, 0.9, 0.369, 0.988, 0.805, 0.544, 0.706, 0.406, 0.964, 0.784, 0.604, 0.969, 0.804, 0.831, 0.85, 0.673, 0.778, 0.687, 0.976, 0.783, 0.393, 0.863, 0.171, 0.653, 0.109, 0.46, 0.962, 0.898, 0.995, 0.99, 0.412, 0.968, 0.914, 0.966, 0.585, 0.27, 0.09, 0.329, 0.997, 0.912, 0.905, 0.968, 0.942, 0.196, 0.481, 0.966, 0.723, 0.905, 0.944, 0.864, 0.989, 0.977, 0.703, 0.544, 0.842, 0.988, 0.361, 0.716, 0.896, 0.357, 0.626, 0.964, 0.328, 0.994, 0.702, 0.3, 0.901, 0.731, 0.926, 0.995, 0.829, 0.827, 0.17, 0.778, 0.134, 0.849, 0.873, 0.843, 0.97, 0.645, 0.837, 0.274, 0.772, 0.752, 0.373, 0.926, 0.672, 0.904, 0.858, 0.096, 0.205, 0.643, 0.592, 0.859, 0.811, 0.972, 0.764, 0.972, 0.586, 0.948, 0.71, 0.975, 0.957, 0.732, 0.152, 0.965, 0.947, 0.498, 0.869, 0.946, 0.725, 0.946, 0.728, 0.764, 0.49, 0.774, 0.789, 0.334, 0.352, 0.312, 0.81, 0.696, 0.277, 0.998, 0.513, 0.822, 0.955, 0.686, 0.402, 0.582, 0.979, 0.784, 0.247, 0.821, 0.7, 0.616, 0.947, 0.939, 0.997, 0.922, 0.986, 0.715, 0.855, 0.879, 0.787, 0.947, 0.936, 0.945, 0.877, 0.264, 0.937, 0.782]
global origin = 1
global destination = 50