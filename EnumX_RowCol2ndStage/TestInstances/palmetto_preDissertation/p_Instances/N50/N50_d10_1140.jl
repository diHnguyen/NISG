global arcs = [1 4; 1 5; 1 18; 1 34; 2 9; 2 24; 2 26; 2 29; 2 31; 2 44; 2 45; 3 5; 3 12; 3 15; 3 33; 4 5; 4 12; 4 26; 4 36; 4 41; 5 12; 5 16; 5 23; 5 29; 5 31; 6 24; 6 42; 7 12; 7 20; 7 47; 8 26; 8 34; 8 36; 9 2; 9 13; 9 29; 9 37; 10 6; 10 40; 10 42; 10 44; 10 50; 11 7; 11 15; 11 22; 11 25; 11 33; 11 48; 12 27; 12 32; 13 14; 13 18; 13 19; 13 20; 13 21; 13 22; 13 38; 14 3; 14 4; 14 15; 14 30; 14 37; 15 16; 15 23; 15 33; 15 35; 15 41; 15 47; 16 31; 17 18; 17 27; 17 28; 17 30; 17 33; 17 40; 17 45; 18 9; 18 15; 18 37; 18 47; 19 10; 19 29; 19 30; 19 39; 19 43; 20 13; 20 32; 20 48; 21 23; 21 24; 21 34; 21 45; 21 46; 21 48; 22 13; 22 14; 22 18; 22 33; 22 39; 22 44; 22 46; 23 5; 23 11; 23 16; 23 31; 23 47; 24 8; 24 12; 24 14; 24 16; 24 43; 24 49; 25 3; 25 11; 25 14; 25 29; 26 11; 26 29; 26 31; 26 47; 26 50; 27 2; 27 16; 27 19; 27 35; 27 42; 27 45; 28 3; 28 12; 28 23; 28 25; 28 31; 28 44; 29 3; 29 6; 29 30; 29 33; 29 35; 30 17; 30 27; 30 40; 30 45; 31 8; 31 15; 31 22; 31 25; 31 28; 31 38; 32 4; 32 35; 32 37; 32 45; 33 6; 34 4; 34 14; 34 19; 34 23; 35 8; 35 17; 35 22; 36 25; 36 30; 36 31; 37 5; 37 10; 37 26; 37 43; 37 44; 37 48; 38 3; 38 22; 38 47; 39 3; 39 5; 39 7; 39 24; 39 30; 39 44; 39 46; 39 50; 40 6; 40 9; 40 17; 40 25; 40 31; 40 37; 40 39; 40 50; 41 2; 41 6; 41 9; 41 50; 42 11; 42 12; 42 13; 42 20; 42 21; 42 32; 42 37; 42 45; 42 48; 42 49; 43 12; 43 13; 43 21; 43 50; 44 5; 44 21; 44 28; 44 29; 44 30; 44 42; 45 4; 45 11; 45 13; 45 20; 45 36; 45 40; 45 44; 45 46; 46 7; 46 38; 46 42; 47 2; 47 12; 47 28; 47 31; 47 36; 47 40; 47 48; 48 3; 48 14; 48 23; 48 41; 48 50; 49 39; 49 43]
global d_x = [1.0, 6.0, 10.0, 10.0, 6.0, 4.0, 5.0, 9.0, 4.0, 4.0, 9.0, 10.0, 8.0, 2.0, 7.0, 4.0, 5.0, 9.0, 2.0, 3.0, 7.0, 7.0, 8.0, 10.0, 2.0, 9.0, 6.0, 9.0, 2.0, 7.0, 7.0, 6.0, 6.0, 3.0, 1.0, 6.0, 9.0, 1.0, 2.0, 8.0, 1.0, 8.0, 6.0, 3.0, 10.0, 6.0, 7.0, 4.0, 4.0, 1.0, 10.0, 4.0, 9.0, 9.0, 7.0, 1.0, 2.0, 7.0, 4.0, 3.0, 9.0, 8.0, 1.0, 9.0, 9.0, 1.0, 10.0, 10.0, 2.0, 3.0, 10.0, 8.0, 1.0, 10.0, 6.0, 1.0, 7.0, 5.0, 9.0, 10.0, 7.0, 1.0, 3.0, 3.0, 4.0, 1.0, 8.0, 8.0, 4.0, 3.0, 7.0, 6.0, 8.0, 3.0, 8.0, 7.0, 7.0, 7.0, 4.0, 3.0, 4.0, 9.0, 2.0, 1.0, 9.0, 5.0, 1.0, 5.0, 2.0, 8.0, 6.0, 10.0, 5.0, 1.0, 2.0, 10.0, 6.0, 4.0, 6.0, 7.0, 9.0, 9.0, 4.0, 2.0, 3.0, 1.0, 8.0, 7.0, 5.0, 2.0, 10.0, 1.0, 3.0, 3.0, 6.0, 5.0, 6.0, 5.0, 4.0, 3.0, 1.0, 2.0, 8.0, 6.0, 3.0, 4.0, 10.0, 1.0, 2.0, 5.0, 7.0, 7.0, 4.0, 10.0, 8.0, 3.0, 10.0, 9.0, 10.0, 5.0, 5.0, 10.0, 4.0, 8.0, 8.0, 5.0, 5.0, 10.0, 8.0, 8.0, 9.0, 4.0, 3.0, 4.0, 10.0, 1.0, 2.0, 5.0, 7.0, 1.0, 8.0, 6.0, 2.0, 4.0, 9.0, 8.0, 4.0, 3.0, 3.0, 10.0, 10.0, 5.0, 1.0, 9.0, 8.0, 10.0, 5.0, 9.0, 10.0, 9.0, 5.0, 6.0, 1.0, 7.0, 6.0, 3.0, 5.0, 2.0, 2.0, 10.0, 3.0, 9.0, 10.0, 1.0, 5.0, 8.0, 6.0, 10.0, 7.0, 6.0, 9.0, 6.0, 7.0, 10.0, 10.0, 8.0, 4.0, 9.0, 4.0, 4.0, 5.0, 6.0, 4.0, 6.0, 2.0, 1.0, 7.0]
global b_x = 5
global d_y = [1.0, 5.0, 1.0, 3.0, 7.0, 4.0, 1.0, 6.0, 10.0, 7.0, 10.0, 4.0, 10.0, 1.0, 5.0, 1.0, 3.0, 9.0, 9.0, 6.0, 5.0, 1.0, 1.0, 2.0, 6.0, 1.0, 10.0, 1.0, 2.0, 7.0, 6.0, 3.0, 9.0, 9.0, 7.0, 7.0, 5.0, 5.0, 4.0, 3.0, 10.0, 8.0, 3.0, 5.0, 3.0, 3.0, 8.0, 3.0, 4.0, 3.0, 2.0, 3.0, 6.0, 10.0, 10.0, 8.0, 10.0, 3.0, 1.0, 5.0, 3.0, 7.0, 3.0, 2.0, 8.0, 4.0, 8.0, 7.0, 7.0, 10.0, 9.0, 3.0, 8.0, 2.0, 7.0, 3.0, 7.0, 4.0, 7.0, 6.0, 6.0, 3.0, 3.0, 10.0, 5.0, 7.0, 2.0, 5.0, 8.0, 1.0, 6.0, 1.0, 4.0, 6.0, 4.0, 6.0, 4.0, 7.0, 10.0, 4.0, 8.0, 3.0, 10.0, 9.0, 7.0, 2.0, 3.0, 5.0, 4.0, 5.0, 2.0, 2.0, 6.0, 1.0, 6.0, 3.0, 1.0, 3.0, 6.0, 10.0, 3.0, 10.0, 1.0, 4.0, 4.0, 2.0, 3.0, 1.0, 10.0, 7.0, 1.0, 8.0, 2.0, 3.0, 6.0, 5.0, 4.0, 9.0, 7.0, 7.0, 9.0, 8.0, 8.0, 8.0, 7.0, 10.0, 9.0, 5.0, 1.0, 9.0, 8.0, 10.0, 6.0, 7.0, 9.0, 3.0, 6.0, 4.0, 2.0, 10.0, 9.0, 5.0, 1.0, 1.0, 2.0, 4.0, 2.0, 5.0, 8.0, 3.0, 1.0, 1.0, 10.0, 8.0, 7.0, 3.0, 8.0, 4.0, 1.0, 6.0, 9.0, 5.0, 5.0, 4.0, 7.0, 9.0, 10.0, 8.0, 10.0, 2.0, 1.0, 2.0, 8.0, 3.0, 2.0, 1.0, 7.0, 2.0, 3.0, 7.0, 7.0, 1.0, 4.0, 2.0, 9.0, 4.0, 6.0, 4.0, 3.0, 1.0, 5.0, 5.0, 9.0, 9.0, 7.0, 1.0, 7.0, 1.0, 4.0, 2.0, 6.0, 10.0, 10.0, 2.0, 4.0, 4.0, 10.0, 5.0, 1.0, 10.0, 2.0, 6.0, 5.0, 5.0, 2.0, 3.0, 4.0]
global b_y = 10
global p = [0.188, 0.434, 0.797, 0.808, 0.082, 0.249, 0.353, 0.598, 0.702, 0.235, 0.093, 0.176, 0.457, 0.186, 0.304, 0.228, 0.491, 0.822, 0.536, 0.805, 0.666, 0.551, 0.297, 0.49, 0.152, 0.361, 0.787, 0.946, 0.76, 0.595, 0.683, 0.181, 0.378, 0.199, 0.859, 0.99, 0.593, 0.451, 0.789, 0.122, 0.278, 0.738, 0.995, 0.166, 0.065, 0.123, 0.426, 0.983, 0.305, 0.546, 0.324, 0.522, 0.602, 0.554, 0.423, 0.746, 0.751, 0.34, 0.526, 0.7, 0.184, 0.584, 0.096, 0.658, 0.316, 0.418, 0.698, 0.474, 0.766, 0.215, 0.759, 0.202, 0.208, 0.222, 0.367, 0.617, 0.14, 0.765, 0.405, 0.216, 0.008, 0.919, 0.797, 0.613, 0.323, 0.683, 0.215, 0.291, 0.271, 0.716, 0.089, 0.13, 0.833, 0.927, 0.389, 0.815, 0.043, 0.931, 0.053, 0.705, 0.534, 0.97, 0.116, 0.932, 0.037, 0.198, 0.128, 0.513, 0.982, 0.2, 0.82, 0.723, 0.061, 0.642, 0.916, 0.702, 0.665, 0.902, 0.141, 0.903, 0.174, 0.17, 0.783, 0.731, 0.765, 0.144, 0.673, 0.715, 0.375, 0.432, 0.074, 0.588, 0.219, 0.488, 0.088, 0.101, 0.852, 0.518, 0.394, 0.147, 0.581, 0.6, 0.255, 0.48, 0.764, 0.797, 0.107, 0.442, 0.896, 0.283, 0.677, 0.253, 0.381, 0.071, 0.486, 0.783, 0.666, 0.645, 0.373, 0.843, 0.226, 0.823, 0.891, 0.771, 0.193, 0.18, 0.642, 0.193, 0.476, 0.53, 0.195, 0.16, 0.267, 0.956, 0.657, 0.686, 0.824, 0.548, 0.459, 0.543, 0.249, 0.51, 0.076, 0.77, 0.284, 0.748, 0.977, 0.893, 0.513, 0.749, 0.563, 0.649, 0.317, 0.222, 0.69, 0.347, 0.522, 0.06, 0.109, 0.386, 0.605, 0.813, 0.691, 0.563, 0.928, 0.4, 0.313, 0.043, 0.871, 0.183, 0.78, 0.894, 0.292, 0.253, 0.369, 0.553, 0.8, 0.993, 0.333, 0.103, 0.011, 0.91, 0.099, 0.106, 0.718, 0.627, 0.112, 0.369, 0.918, 0.009, 0.574, 0.29, 0.057, 0.29, 0.633, 0.879, 0.823]
global q = [0.331, 0.796, 0.829, 0.892, 0.761, 0.922, 0.422, 0.948, 0.891, 0.779, 0.345, 0.879, 0.644, 0.305, 0.654, 0.931, 0.548, 0.893, 0.906, 0.88, 0.803, 0.734, 0.893, 0.883, 0.974, 0.542, 0.903, 0.979, 0.76, 0.71, 0.85, 0.266, 0.986, 0.651, 0.909, 0.994, 0.752, 0.908, 0.852, 0.211, 0.44, 0.973, 0.997, 0.393, 0.261, 0.252, 0.446, 0.999, 0.4, 0.769, 0.792, 0.899, 0.943, 0.669, 0.494, 0.843, 0.836, 0.461, 0.567, 0.725, 0.827, 0.639, 0.688, 0.662, 0.836, 0.929, 0.845, 0.509, 0.954, 0.856, 0.896, 0.222, 0.79, 0.684, 0.725, 0.99, 0.772, 0.855, 0.422, 0.995, 0.686, 0.932, 0.909, 0.674, 0.488, 0.807, 0.541, 0.997, 0.634, 0.765, 0.615, 0.497, 0.969, 0.988, 0.434, 0.947, 0.526, 0.944, 0.715, 0.942, 0.695, 0.981, 0.959, 0.965, 0.934, 0.265, 0.813, 0.617, 0.992, 0.805, 0.987, 0.762, 0.963, 0.718, 0.994, 0.923, 0.833, 0.914, 0.267, 0.934, 0.759, 0.531, 0.856, 0.74, 0.889, 0.751, 0.929, 0.789, 0.4, 0.569, 0.115, 0.655, 0.564, 0.533, 0.249, 0.706, 0.937, 0.767, 0.454, 0.271, 0.844, 0.8, 0.449, 0.662, 0.901, 0.877, 0.152, 0.739, 0.928, 0.577, 0.698, 0.486, 0.863, 0.493, 0.773, 0.922, 0.93, 0.715, 0.541, 0.843, 0.629, 0.905, 0.909, 0.941, 0.448, 0.916, 0.745, 0.535, 0.806, 0.62, 0.737, 0.86, 0.648, 0.991, 0.804, 0.816, 0.995, 0.851, 0.665, 0.702, 0.636, 0.608, 0.626, 0.828, 0.289, 0.928, 0.988, 0.921, 0.659, 0.896, 0.849, 0.795, 0.876, 0.981, 0.75, 0.396, 0.622, 0.437, 0.909, 0.648, 0.718, 0.925, 0.844, 0.754, 0.948, 0.877, 0.908, 0.778, 0.974, 0.65, 0.819, 0.963, 0.309, 0.906, 0.843, 0.956, 0.891, 0.995, 0.924, 0.668, 0.449, 0.954, 0.696, 0.52, 0.865, 0.715, 0.674, 0.571, 0.922, 0.926, 0.87, 0.602, 0.208, 0.908, 0.978, 0.927, 0.969]
global origin = 1
global destination = 50