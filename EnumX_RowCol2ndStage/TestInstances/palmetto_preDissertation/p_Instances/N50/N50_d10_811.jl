global arcs = [1 11; 1 13; 1 22; 1 23; 1 41; 1 48; 2 7; 2 19; 2 37; 3 22; 3 27; 3 36; 3 42; 3 44; 3 47; 4 3; 4 8; 4 28; 4 34; 4 40; 4 44; 4 47; 5 36; 5 37; 5 39; 6 3; 6 30; 6 36; 6 38; 6 40; 6 44; 6 47; 7 9; 7 19; 7 27; 7 32; 7 36; 7 44; 7 45; 8 6; 8 39; 8 47; 9 12; 9 14; 10 17; 10 44; 10 45; 11 9; 11 16; 11 20; 11 24; 11 34; 11 41; 12 5; 12 9; 12 24; 12 49; 13 6; 13 7; 13 32; 13 41; 14 4; 14 11; 14 45; 15 8; 15 36; 15 38; 15 40; 16 9; 16 24; 16 27; 16 32; 16 36; 17 3; 17 9; 17 27; 18 3; 18 24; 18 36; 19 4; 19 15; 19 16; 19 25; 19 36; 19 44; 19 45; 20 15; 20 31; 20 33; 20 39; 20 40; 20 45; 21 19; 21 20; 21 23; 21 25; 21 33; 21 45; 22 2; 22 15; 22 41; 22 43; 22 47; 23 3; 23 4; 23 6; 23 7; 23 9; 23 18; 23 28; 23 41; 23 44; 24 17; 24 32; 24 40; 24 42; 24 49; 25 11; 25 19; 25 20; 25 50; 26 22; 26 23; 26 25; 26 31; 27 13; 27 20; 27 33; 27 39; 27 40; 27 48; 28 4; 28 13; 28 34; 28 40; 28 41; 28 43; 28 48; 29 14; 29 21; 29 22; 29 24; 29 43; 30 5; 30 7; 30 13; 30 21; 30 25; 30 43; 30 47; 31 19; 31 43; 31 46; 32 2; 32 6; 32 7; 32 9; 32 13; 32 16; 32 45; 33 17; 33 21; 33 25; 33 29; 33 35; 33 40; 33 43; 34 5; 34 10; 34 16; 34 26; 34 39; 34 46; 34 49; 35 2; 35 15; 35 28; 35 29; 35 30; 36 4; 36 12; 36 31; 36 44; 37 11; 37 16; 37 27; 37 31; 38 16; 38 29; 38 37; 38 43; 39 14; 39 27; 39 31; 39 50; 40 5; 40 6; 40 9; 40 11; 40 14; 40 22; 40 39; 40 46; 40 48; 41 6; 41 21; 41 38; 41 39; 42 8; 42 41; 42 45; 42 46; 42 50; 43 28; 43 48; 44 21; 44 48; 45 8; 45 16; 45 37; 46 21; 46 28; 47 29; 47 35; 47 45; 48 8; 48 9; 48 24; 48 29; 48 36; 49 10; 49 13; 49 25; 49 36]
global d_x = [3.0, 5.0, 7.0, 5.0, 3.0, 1.0, 2.0, 1.0, 10.0, 7.0, 8.0, 7.0, 2.0, 10.0, 10.0, 3.0, 4.0, 8.0, 7.0, 8.0, 8.0, 1.0, 3.0, 7.0, 2.0, 9.0, 3.0, 5.0, 6.0, 5.0, 5.0, 5.0, 9.0, 8.0, 8.0, 5.0, 8.0, 1.0, 4.0, 8.0, 9.0, 4.0, 5.0, 1.0, 7.0, 8.0, 5.0, 8.0, 4.0, 2.0, 1.0, 5.0, 1.0, 3.0, 7.0, 3.0, 9.0, 1.0, 1.0, 10.0, 3.0, 9.0, 5.0, 1.0, 7.0, 10.0, 7.0, 7.0, 1.0, 9.0, 9.0, 2.0, 2.0, 4.0, 8.0, 5.0, 10.0, 6.0, 8.0, 9.0, 9.0, 5.0, 7.0, 1.0, 9.0, 10.0, 5.0, 6.0, 10.0, 6.0, 3.0, 5.0, 6.0, 1.0, 8.0, 2.0, 8.0, 1.0, 7.0, 10.0, 8.0, 4.0, 1.0, 8.0, 8.0, 3.0, 2.0, 10.0, 3.0, 8.0, 6.0, 1.0, 4.0, 7.0, 1.0, 7.0, 10.0, 2.0, 5.0, 1.0, 2.0, 2.0, 4.0, 8.0, 2.0, 8.0, 1.0, 3.0, 1.0, 3.0, 4.0, 5.0, 3.0, 4.0, 3.0, 2.0, 2.0, 8.0, 1.0, 9.0, 1.0, 6.0, 4.0, 7.0, 10.0, 10.0, 9.0, 8.0, 10.0, 3.0, 4.0, 8.0, 2.0, 2.0, 5.0, 5.0, 3.0, 8.0, 4.0, 6.0, 6.0, 5.0, 9.0, 5.0, 2.0, 3.0, 5.0, 1.0, 8.0, 7.0, 7.0, 7.0, 6.0, 7.0, 1.0, 8.0, 5.0, 4.0, 6.0, 2.0, 10.0, 3.0, 8.0, 4.0, 9.0, 5.0, 5.0, 8.0, 3.0, 5.0, 10.0, 6.0, 10.0, 5.0, 2.0, 5.0, 10.0, 6.0, 4.0, 4.0, 5.0, 6.0, 3.0, 1.0, 10.0, 3.0, 6.0, 5.0, 8.0, 9.0, 7.0, 8.0, 5.0, 9.0, 5.0, 4.0, 2.0, 7.0, 3.0, 3.0, 1.0, 1.0, 10.0, 8.0, 7.0, 5.0, 8.0, 10.0, 4.0, 2.0, 5.0, 10.0, 1.0, 7.0]
global b_x = 5
global d_y = [3.0, 2.0, 7.0, 3.0, 8.0, 6.0, 1.0, 4.0, 1.0, 6.0, 6.0, 8.0, 2.0, 8.0, 10.0, 7.0, 9.0, 7.0, 2.0, 9.0, 1.0, 6.0, 7.0, 10.0, 8.0, 3.0, 6.0, 4.0, 4.0, 5.0, 2.0, 2.0, 9.0, 6.0, 1.0, 7.0, 4.0, 8.0, 7.0, 10.0, 7.0, 8.0, 9.0, 2.0, 4.0, 6.0, 5.0, 10.0, 3.0, 5.0, 2.0, 7.0, 10.0, 6.0, 5.0, 1.0, 5.0, 1.0, 9.0, 9.0, 7.0, 9.0, 2.0, 9.0, 6.0, 5.0, 10.0, 6.0, 9.0, 4.0, 8.0, 2.0, 4.0, 3.0, 9.0, 2.0, 1.0, 10.0, 10.0, 6.0, 10.0, 9.0, 3.0, 8.0, 7.0, 8.0, 8.0, 3.0, 3.0, 5.0, 6.0, 9.0, 6.0, 9.0, 7.0, 7.0, 4.0, 6.0, 7.0, 6.0, 8.0, 10.0, 9.0, 8.0, 4.0, 7.0, 5.0, 4.0, 8.0, 4.0, 8.0, 3.0, 4.0, 2.0, 2.0, 7.0, 3.0, 8.0, 7.0, 9.0, 10.0, 1.0, 3.0, 10.0, 3.0, 5.0, 4.0, 9.0, 4.0, 3.0, 6.0, 10.0, 1.0, 10.0, 10.0, 7.0, 2.0, 5.0, 6.0, 1.0, 6.0, 7.0, 8.0, 5.0, 7.0, 10.0, 1.0, 2.0, 5.0, 8.0, 6.0, 1.0, 5.0, 10.0, 6.0, 7.0, 10.0, 5.0, 9.0, 6.0, 7.0, 5.0, 9.0, 8.0, 9.0, 6.0, 8.0, 1.0, 2.0, 3.0, 8.0, 4.0, 2.0, 4.0, 1.0, 2.0, 7.0, 1.0, 3.0, 8.0, 2.0, 3.0, 4.0, 8.0, 4.0, 9.0, 9.0, 1.0, 7.0, 8.0, 8.0, 4.0, 1.0, 7.0, 8.0, 5.0, 2.0, 10.0, 3.0, 6.0, 7.0, 9.0, 9.0, 6.0, 5.0, 8.0, 8.0, 8.0, 10.0, 6.0, 9.0, 2.0, 8.0, 2.0, 10.0, 10.0, 4.0, 2.0, 10.0, 2.0, 9.0, 10.0, 9.0, 9.0, 6.0, 9.0, 5.0, 4.0, 2.0, 10.0, 3.0, 7.0, 3.0, 2.0]
global b_y = 10
global p = [0.079, 0.695, 0.864, 0.429, 0.115, 0.54, 0.021, 0.251, 0.671, 0.822, 0.019, 0.701, 0.591, 0.457, 0.064, 0.084, 0.576, 0.612, 0.192, 0.27, 0.813, 0.29, 0.101, 0.22, 0.692, 0.346, 0.784, 0.864, 0.145, 0.154, 0.349, 0.365, 0.599, 0.349, 0.474, 0.028, 0.644, 0.131, 0.941, 0.013, 0.184, 0.715, 0.815, 0.619, 0.201, 0.714, 0.416, 0.366, 0.9, 0.583, 0.057, 0.88, 0.314, 0.78, 0.621, 0.687, 0.748, 0.735, 0.95, 0.844, 0.984, 0.331, 0.456, 0.151, 0.161, 0.876, 0.143, 0.051, 0.545, 0.643, 0.871, 0.65, 0.492, 0.765, 0.43, 0.668, 0.43, 0.194, 0.078, 0.989, 0.166, 0.735, 0.228, 0.693, 0.036, 0.648, 0.34, 0.263, 0.314, 0.953, 0.793, 0.619, 0.285, 0.099, 0.806, 0.436, 0.546, 0.566, 0.513, 0.061, 0.139, 0.664, 0.493, 0.232, 0.893, 0.505, 0.718, 0.766, 0.608, 0.714, 0.746, 0.389, 0.744, 0.789, 0.056, 0.212, 0.605, 0.446, 0.565, 0.553, 0.83, 0.233, 0.6, 0.403, 0.049, 0.915, 0.156, 0.468, 0.229, 0.277, 0.038, 0.395, 0.861, 0.602, 0.313, 0.627, 0.601, 0.354, 0.017, 0.762, 0.27, 0.112, 0.171, 0.291, 0.211, 0.649, 0.735, 0.627, 0.382, 0.185, 0.634, 0.127, 0.627, 0.148, 0.377, 0.723, 0.906, 0.74, 0.107, 0.015, 0.943, 0.623, 0.138, 0.059, 0.072, 0.884, 0.236, 0.844, 0.019, 0.173, 0.853, 0.189, 0.62, 0.798, 0.912, 0.807, 0.489, 0.885, 0.117, 0.307, 0.837, 0.785, 0.941, 0.806, 0.157, 0.305, 0.144, 0.255, 0.968, 0.978, 0.08, 0.982, 0.31, 0.417, 0.125, 0.709, 0.033, 0.207, 0.881, 0.109, 0.754, 0.645, 0.608, 0.513, 0.65, 0.454, 0.582, 0.649, 0.573, 0.743, 0.542, 0.964, 0.701, 0.851, 0.854, 0.252, 0.89, 0.103, 0.18, 0.455, 0.277, 0.032, 0.708, 0.29, 0.754, 0.89, 0.767, 0.501, 0.749, 0.548, 0.92, 0.48, 0.386, 0.705]
global q = [0.96, 0.784, 0.92, 0.848, 0.417, 0.559, 0.907, 0.418, 0.859, 0.911, 0.03, 0.734, 0.928, 0.642, 0.944, 0.81, 0.819, 0.911, 0.35, 0.433, 0.95, 0.841, 0.475, 0.57, 0.975, 0.938, 0.877, 0.912, 0.771, 0.978, 0.924, 0.606, 0.838, 0.624, 0.637, 0.534, 0.874, 0.678, 0.958, 0.453, 0.921, 0.787, 0.849, 0.668, 0.691, 0.969, 0.79, 0.781, 0.978, 0.974, 0.491, 0.894, 0.604, 0.784, 0.865, 0.694, 0.943, 0.934, 0.985, 0.87, 0.994, 0.382, 0.936, 0.314, 0.39, 0.893, 0.183, 0.788, 0.734, 0.699, 0.874, 0.97, 0.843, 0.979, 0.499, 0.71, 0.905, 0.803, 0.546, 0.994, 0.56, 0.971, 0.729, 0.713, 0.757, 0.838, 0.44, 0.285, 0.798, 0.982, 0.955, 0.714, 0.723, 0.953, 0.9, 0.788, 0.577, 0.658, 0.679, 0.259, 0.985, 0.686, 0.718, 0.565, 0.944, 0.674, 0.861, 0.913, 0.965, 0.794, 0.811, 0.986, 0.778, 0.864, 0.507, 0.519, 0.694, 0.775, 0.574, 0.797, 0.902, 0.901, 0.604, 0.92, 0.405, 0.948, 0.508, 0.566, 0.804, 0.59, 0.414, 0.813, 0.921, 0.887, 0.467, 0.669, 0.852, 0.925, 0.474, 0.787, 0.705, 0.712, 0.679, 0.308, 0.691, 0.7, 0.942, 0.79, 0.424, 0.22, 0.947, 0.459, 0.794, 0.313, 0.616, 0.842, 0.986, 0.917, 0.313, 0.025, 0.957, 0.815, 0.768, 0.718, 0.472, 0.97, 0.451, 0.846, 0.337, 0.825, 0.86, 0.982, 0.995, 0.925, 0.973, 0.989, 0.63, 0.921, 0.899, 0.451, 0.999, 0.982, 0.986, 0.811, 0.977, 0.555, 0.388, 0.957, 0.983, 0.985, 0.549, 0.996, 0.983, 0.896, 0.847, 0.78, 0.316, 0.6, 0.997, 0.922, 0.844, 0.727, 0.806, 0.718, 0.712, 0.647, 0.938, 0.731, 0.879, 0.985, 0.618, 0.976, 0.979, 0.957, 0.939, 0.416, 0.897, 0.578, 0.365, 0.611, 0.688, 0.691, 0.878, 0.467, 0.863, 0.951, 0.806, 0.699, 0.813, 0.641, 0.965, 0.829, 0.581, 0.996]
global origin = 1
global destination = 50