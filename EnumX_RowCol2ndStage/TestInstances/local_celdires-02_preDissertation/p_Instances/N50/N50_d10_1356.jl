global arcs = [1 4; 1 16; 1 17; 1 24; 1 27; 1 31; 2 5; 2 10; 2 25; 2 27; 2 32; 2 36; 3 12; 3 16; 3 29; 4 15; 4 17; 4 32; 4 34; 4 37; 4 48; 5 6; 5 9; 5 19; 5 40; 5 41; 5 42; 5 44; 6 24; 6 47; 7 5; 7 14; 7 18; 7 20; 7 31; 8 5; 8 10; 8 16; 8 19; 8 29; 8 31; 9 5; 9 10; 9 27; 9 33; 9 45; 10 5; 10 9; 10 15; 10 16; 10 24; 10 26; 10 32; 10 38; 10 40; 10 41; 10 47; 11 18; 11 20; 11 43; 11 44; 12 5; 12 7; 12 10; 12 31; 12 47; 13 9; 13 16; 13 22; 14 17; 14 23; 14 29; 14 37; 14 40; 15 4; 15 6; 15 19; 15 42; 15 46; 16 11; 16 18; 16 33; 17 3; 17 4; 17 7; 17 13; 17 35; 17 46; 18 32; 18 48; 19 5; 19 10; 19 14; 19 17; 20 6; 20 33; 21 23; 21 26; 21 31; 21 36; 22 6; 22 23; 22 26; 22 33; 22 35; 22 38; 23 8; 23 15; 23 20; 24 3; 24 27; 24 33; 24 34; 24 37; 24 39; 24 44; 25 5; 25 7; 25 19; 25 24; 25 27; 25 29; 25 33; 25 39; 25 42; 26 12; 26 13; 26 23; 26 49; 26 50; 27 12; 27 13; 27 23; 27 28; 27 29; 27 37; 27 44; 28 27; 28 29; 28 46; 29 4; 29 18; 29 31; 29 33; 29 37; 29 40; 29 43; 29 50; 30 5; 30 14; 30 18; 30 27; 30 43; 30 50; 31 7; 31 16; 31 35; 31 47; 31 49; 32 5; 32 10; 32 43; 33 13; 33 26; 34 3; 34 4; 34 27; 34 38; 34 43; 34 50; 35 4; 35 33; 36 2; 36 5; 36 14; 36 17; 36 21; 36 29; 36 30; 36 44; 37 13; 37 15; 37 29; 38 8; 38 17; 38 21; 39 7; 39 10; 39 30; 39 37; 39 48; 40 4; 40 12; 40 25; 40 34; 41 4; 41 10; 41 25; 41 42; 41 49; 42 29; 42 36; 42 41; 42 46; 42 50; 43 18; 43 44; 43 47; 43 48; 44 2; 44 15; 44 17; 44 23; 45 7; 45 14; 45 39; 45 50; 46 2; 46 9; 46 16; 46 18; 46 23; 46 24; 46 48; 47 4; 47 23; 47 33; 47 34; 47 39; 47 40; 47 50; 48 11; 48 18; 48 26; 48 28; 48 30; 48 32; 49 5; 49 15; 49 22; 49 38; 49 50]
global d_x = [8.0, 8.0, 1.0, 2.0, 10.0, 5.0, 4.0, 8.0, 7.0, 2.0, 6.0, 5.0, 4.0, 7.0, 6.0, 1.0, 7.0, 6.0, 8.0, 10.0, 10.0, 6.0, 9.0, 6.0, 10.0, 8.0, 8.0, 1.0, 9.0, 7.0, 3.0, 3.0, 2.0, 6.0, 4.0, 3.0, 5.0, 4.0, 9.0, 4.0, 1.0, 2.0, 1.0, 3.0, 7.0, 5.0, 1.0, 7.0, 6.0, 8.0, 4.0, 2.0, 9.0, 9.0, 10.0, 3.0, 7.0, 8.0, 7.0, 10.0, 1.0, 5.0, 9.0, 3.0, 3.0, 2.0, 4.0, 2.0, 9.0, 6.0, 10.0, 1.0, 6.0, 7.0, 8.0, 2.0, 4.0, 5.0, 2.0, 8.0, 7.0, 10.0, 7.0, 4.0, 3.0, 8.0, 6.0, 4.0, 5.0, 7.0, 7.0, 2.0, 4.0, 5.0, 6.0, 1.0, 1.0, 1.0, 5.0, 10.0, 4.0, 9.0, 2.0, 7.0, 3.0, 6.0, 9.0, 2.0, 3.0, 8.0, 1.0, 7.0, 9.0, 6.0, 4.0, 1.0, 6.0, 9.0, 10.0, 7.0, 5.0, 8.0, 2.0, 2.0, 2.0, 9.0, 6.0, 5.0, 10.0, 5.0, 3.0, 8.0, 1.0, 3.0, 5.0, 6.0, 7.0, 8.0, 7.0, 10.0, 2.0, 5.0, 3.0, 1.0, 6.0, 1.0, 8.0, 3.0, 9.0, 5.0, 2.0, 3.0, 4.0, 5.0, 2.0, 9.0, 2.0, 3.0, 5.0, 8.0, 4.0, 4.0, 9.0, 6.0, 3.0, 7.0, 6.0, 7.0, 10.0, 2.0, 9.0, 2.0, 1.0, 9.0, 5.0, 7.0, 4.0, 2.0, 6.0, 4.0, 6.0, 10.0, 8.0, 7.0, 1.0, 1.0, 7.0, 1.0, 7.0, 2.0, 6.0, 8.0, 7.0, 1.0, 7.0, 10.0, 1.0, 8.0, 9.0, 4.0, 10.0, 1.0, 9.0, 3.0, 9.0, 5.0, 5.0, 4.0, 10.0, 1.0, 2.0, 9.0, 3.0, 2.0, 4.0, 5.0, 9.0, 9.0, 1.0, 8.0, 9.0, 5.0, 2.0, 5.0, 5.0, 2.0, 9.0, 1.0, 10.0, 7.0, 3.0, 7.0, 7.0, 1.0, 2.0, 6.0, 9.0, 8.0, 8.0, 10.0, 7.0, 8.0]
global b_x = 5
global d_y = [4.0, 2.0, 2.0, 7.0, 6.0, 9.0, 6.0, 9.0, 6.0, 4.0, 6.0, 3.0, 4.0, 7.0, 5.0, 5.0, 7.0, 6.0, 2.0, 2.0, 2.0, 3.0, 5.0, 1.0, 4.0, 9.0, 1.0, 6.0, 8.0, 9.0, 7.0, 8.0, 3.0, 4.0, 3.0, 2.0, 5.0, 5.0, 5.0, 5.0, 1.0, 8.0, 4.0, 8.0, 2.0, 3.0, 6.0, 10.0, 1.0, 8.0, 9.0, 10.0, 8.0, 3.0, 8.0, 10.0, 9.0, 6.0, 9.0, 2.0, 3.0, 3.0, 8.0, 5.0, 9.0, 5.0, 9.0, 1.0, 2.0, 8.0, 3.0, 9.0, 10.0, 5.0, 9.0, 8.0, 1.0, 9.0, 9.0, 4.0, 6.0, 10.0, 8.0, 4.0, 2.0, 3.0, 5.0, 6.0, 4.0, 5.0, 1.0, 2.0, 7.0, 8.0, 10.0, 1.0, 3.0, 7.0, 8.0, 5.0, 5.0, 9.0, 2.0, 8.0, 3.0, 2.0, 8.0, 5.0, 1.0, 5.0, 2.0, 2.0, 5.0, 4.0, 5.0, 9.0, 6.0, 5.0, 6.0, 3.0, 1.0, 1.0, 6.0, 2.0, 7.0, 9.0, 9.0, 5.0, 7.0, 4.0, 8.0, 4.0, 9.0, 2.0, 10.0, 7.0, 9.0, 4.0, 2.0, 4.0, 10.0, 8.0, 9.0, 8.0, 4.0, 2.0, 2.0, 10.0, 4.0, 7.0, 2.0, 9.0, 5.0, 7.0, 5.0, 9.0, 10.0, 10.0, 6.0, 9.0, 5.0, 1.0, 1.0, 1.0, 8.0, 5.0, 2.0, 5.0, 6.0, 2.0, 9.0, 7.0, 10.0, 7.0, 1.0, 10.0, 2.0, 8.0, 9.0, 9.0, 4.0, 1.0, 3.0, 6.0, 1.0, 8.0, 2.0, 5.0, 6.0, 8.0, 7.0, 3.0, 9.0, 2.0, 2.0, 8.0, 1.0, 7.0, 10.0, 8.0, 9.0, 2.0, 8.0, 5.0, 9.0, 2.0, 6.0, 1.0, 9.0, 4.0, 4.0, 10.0, 1.0, 1.0, 3.0, 1.0, 4.0, 9.0, 5.0, 1.0, 9.0, 6.0, 10.0, 3.0, 6.0, 8.0, 3.0, 3.0, 3.0, 4.0, 2.0, 3.0, 10.0, 4.0, 1.0, 9.0, 10.0, 4.0, 10.0, 2.0, 4.0, 8.0]
global b_y = 10
global p = [0.612, 0.772, 0.731, 0.378, 0.648, 0.93, 0.307, 0.931, 0.528, 0.474, 0.765, 0.628, 0.519, 0.243, 0.133, 0.072, 0.764, 0.598, 0.565, 0.94, 0.177, 0.351, 0.262, 0.534, 0.789, 0.435, 0.382, 0.519, 0.149, 0.755, 0.748, 0.312, 0.459, 0.853, 0.464, 0.896, 0.738, 0.377, 0.512, 0.232, 0.544, 0.792, 0.752, 0.03, 0.926, 0.628, 0.755, 0.669, 0.092, 0.591, 0.838, 0.037, 0.545, 0.745, 0.428, 0.325, 0.669, 0.039, 0.315, 0.678, 0.051, 0.108, 0.136, 0.951, 0.523, 0.038, 0.033, 0.258, 0.691, 0.484, 0.185, 0.146, 0.433, 0.149, 0.192, 0.893, 0.294, 0.09, 0.245, 0.055, 0.455, 0.655, 0.805, 0.89, 0.046, 0.531, 0.177, 0.238, 0.895, 0.64, 0.81, 0.275, 0.942, 0.961, 0.586, 0.031, 0.709, 0.348, 0.82, 0.168, 0.284, 0.512, 0.618, 0.917, 0.462, 0.498, 0.331, 0.101, 0.269, 0.005, 0.783, 0.683, 0.953, 0.283, 0.168, 0.136, 0.034, 0.165, 0.799, 0.42, 0.406, 0.975, 0.38, 0.608, 0.565, 0.613, 0.292, 0.402, 0.599, 0.646, 0.193, 0.12, 0.421, 0.738, 0.507, 0.508, 0.707, 0.113, 0.276, 0.06, 0.941, 0.612, 0.809, 0.118, 0.699, 0.873, 0.46, 0.12, 0.295, 0.937, 0.1, 0.438, 0.061, 0.191, 0.34, 0.158, 0.558, 0.336, 0.561, 0.903, 0.436, 0.207, 0.706, 0.905, 0.521, 0.218, 0.724, 0.197, 0.917, 0.793, 0.89, 0.93, 0.239, 0.138, 0.492, 0.533, 0.46, 0.887, 0.764, 0.818, 0.445, 0.237, 0.635, 0.93, 0.557, 0.064, 0.768, 0.418, 0.067, 0.298, 0.895, 0.128, 0.078, 0.45, 0.379, 0.208, 0.924, 0.298, 0.02, 0.478, 0.04, 0.095, 0.372, 0.491, 0.199, 0.431, 0.919, 0.579, 0.818, 0.122, 0.694, 0.802, 0.104, 0.079, 0.474, 0.72, 0.837, 0.116, 0.287, 0.192, 0.132, 0.444, 0.083, 0.971, 0.363, 0.706, 0.42, 0.988, 0.423, 0.933, 0.599, 0.049, 0.98, 0.613, 0.8, 0.29, 0.523, 0.729, 0.1, 0.189, 0.201, 0.906]
global q = [0.621, 0.947, 0.954, 0.608, 0.744, 0.953, 0.405, 0.937, 0.716, 0.698, 0.855, 0.716, 0.795, 0.476, 0.17, 0.549, 0.799, 0.819, 0.667, 0.986, 0.221, 0.999, 0.609, 0.924, 0.923, 0.455, 0.653, 0.911, 0.534, 0.783, 0.79, 0.506, 0.816, 0.932, 0.656, 0.977, 0.913, 0.867, 0.979, 0.963, 0.715, 0.88, 0.88, 0.985, 0.974, 0.96, 0.804, 0.885, 0.639, 0.996, 0.93, 0.346, 0.745, 0.94, 0.894, 0.837, 0.832, 0.812, 0.907, 0.807, 0.881, 0.582, 0.801, 0.957, 0.721, 0.36, 0.742, 0.512, 0.959, 0.847, 0.205, 0.264, 0.584, 0.151, 0.246, 0.936, 0.5, 0.68, 0.709, 0.144, 0.842, 0.766, 0.902, 0.935, 0.737, 0.944, 0.961, 0.948, 0.941, 0.69, 0.983, 0.793, 0.957, 0.982, 0.804, 0.696, 0.897, 0.63, 0.847, 0.953, 0.348, 0.678, 0.872, 0.992, 0.882, 0.934, 0.334, 0.197, 0.485, 0.677, 0.956, 0.971, 0.963, 0.341, 0.995, 0.349, 0.842, 0.322, 0.999, 0.462, 0.481, 0.994, 0.42, 0.757, 0.818, 0.695, 0.988, 0.948, 0.885, 0.854, 0.433, 0.724, 0.531, 0.81, 0.676, 0.647, 0.746, 0.425, 0.334, 0.448, 0.964, 0.689, 0.873, 0.695, 0.88, 0.902, 0.731, 0.209, 0.79, 0.96, 0.326, 0.53, 0.177, 0.467, 0.617, 0.403, 0.711, 0.915, 0.729, 0.945, 0.99, 0.568, 0.821, 0.923, 0.598, 0.54, 0.734, 0.782, 0.963, 0.846, 0.978, 0.962, 0.751, 0.9, 0.621, 0.656, 0.961, 0.921, 0.915, 0.898, 0.992, 0.45, 0.764, 0.93, 0.879, 0.314, 0.995, 0.978, 0.809, 0.962, 0.926, 0.18, 0.464, 0.817, 0.557, 0.607, 0.977, 0.711, 0.728, 0.731, 0.638, 0.324, 0.689, 0.976, 0.65, 0.589, 0.931, 0.809, 0.845, 0.827, 0.862, 0.834, 0.199, 0.953, 0.991, 0.946, 0.938, 0.646, 0.314, 0.712, 0.516, 0.625, 0.801, 0.979, 0.747, 0.857, 0.898, 0.993, 0.481, 0.936, 0.895, 0.821, 0.998, 0.861, 0.973, 0.658, 0.82, 0.747, 0.926, 0.82, 0.539, 0.938]
global origin = 1
global destination = 50